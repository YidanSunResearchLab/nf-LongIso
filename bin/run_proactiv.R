#!/usr/bin/env Rscript
# NextLongIso: promoter activity (proActiv) from long-read junctions.
#
# Usage:
#   run_proactiv.R --gtf=cleaned_annotations.gtf --design=design.csv --species=Homo_sapiens \
#              --reference=control --altprom-args-file=altprom_args.txt --cpus=8
#
# BAM files (<sample>.sorted.bam) are read from the working directory. Junctions are
# converted to STAR SJ.out.tab format, which proActiv accepts as input.
# --altprom-args-file contains extra R arguments for getAlternativePromoters(), e.g.
#   minAbs = 0.25, minRel = 0.05, maxPval = 0.05, promoterFC = 2.0, geneFC = 1.5
#
# Outputs: promoter_activity.tsv, promoter_counts.tsv, alternative_promoters.tsv (if any)

suppressPackageStartupMessages({
    library(proActiv)
    library(GenomicFeatures)
    library(GenomicAlignments)
    library(SummarizedExperiment)
})

# Shared helpers (parse_args, design handling) from bin/nli_design.R
.script <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "nli_design.R"))
opt <- parse_args(commandArgs(trailingOnly = TRUE))
for (k in c("gtf", "design", "species", "reference", "altprom-args-file", "cpus")) {
    if (is.null(opt[[k]])) stop(sprintf("Missing required argument --%s", k))
}

# ---- 1. Promoter annotation ---------------------------------------------------------------------
message("Preparing promoter annotation (species: ", opt$species, ")...")
txdb           <- makeTxDbFromGFF(opt$gtf)
promoter_annot <- preparePromoterAnnotation(txdb = txdb, species = opt$species)

# ---- 2. Junctions from BAM files -----------------------------------------------------------------
bam_files <- sort(list.files(pattern = "[.]sorted[.]bam$"))
message(sprintf("Found %d BAM files", length(bam_files)))
junc_files   <- character(0)
sample_names <- character(0)
for (bam_file in bam_files) {
    sname <- sub("[.]sorted[.]bam$", "", bam_file)
    gal   <- readGAlignments(bam_file, param = ScanBamParam(
                 flag = scanBamFlag(isSecondaryAlignment = FALSE, isUnmappedQuery = FALSE)))
    juncs <- summarizeJunctions(gal)
    if (length(juncs) == 0) {
        warning(sprintf("No junctions found in %s - sample skipped", bam_file))
        next
    }
    df <- as.data.frame(juncs)
    star_df <- data.frame(
        chrom     = as.character(df$seqnames),
        start     = df$start,
        end       = df$end,
        strand    = ifelse(as.character(df$strand) == "+", 1L, ifelse(as.character(df$strand) == "-", 2L, 0L)),
        motif     = 1L,
        annotated = 1L,
        unique    = pmax(as.integer(df$score), 1L),
        multi     = 0L,
        overhang  = 10L
    )
    out_name <- paste0(sname, ".SJ.out.tab")
    write.table(star_df, out_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    junc_files   <- c(junc_files, out_name)
    sample_names <- c(sample_names, sname)
    message(sprintf("  %s: %d junctions", sname, nrow(star_df)))
}
if (length(junc_files) == 0) stop("No junction files created - aborting.")

# ---- 3. Sample -> condition (exact sample names from the validated samplesheet) ---------------
design <- read.csv(opt$design, stringsAsFactors = FALSE)
unknown <- setdiff(sample_names, design$sample)
if (length(unknown) > 0) stop(sprintf("BAM files without a samplesheet entry: %s", paste(unknown, collapse = ", ")))
conditions <- design$condition[match(sample_names, design$sample)]
print(data.frame(Sample = sample_names, Condition = conditions))

# ---- 4. proActiv -------------------------------------------------------------------------------
result <- proActiv(files = junc_files, promoterAnnotation = promoter_annot,
                   condition = conditions, ncores = as.integer(opt$cpus))

pos_info <- as.data.frame(rowData(result))
for (col in names(which(sapply(pos_info, is.list)))) {
    pos_info[[col]] <- sapply(pos_info[[col]], function(x) {
        if (is.null(x) || length(x) == 0) NA_character_ else paste(x, collapse = ",")
    })
}
assay_list <- assays(result)
write.table(cbind(pos_info, as.data.frame(assay_list[["absolutePromoterActivity"]])),
            "promoter_activity.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cbind(pos_info, as.data.frame(assay_list[["promoterCounts"]])),
            "promoter_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 5. Alternative promoters --------------------------------------------------------------------
if (length(unique(conditions)) < 2) {
    message("Only one condition - skipping alternative promoter detection")
    quit(save = "no", status = 0)
}
reference <- opt$reference
if (!reference %in% conditions) stop(sprintf("Reference condition '%s' not found among samples", reference))

arg_text <- trimws(paste(readLines(opt[["altprom-args-file"]], warn = FALSE), collapse = " "))
extra    <- if (nzchar(arg_text)) eval(parse(text = sprintf("list(%s)", arg_text))) else list()
dup      <- intersect(names(extra), c("result", "referenceCondition"))
if (length(dup) > 0) stop(sprintf("--proactiv_altprom_args must not set %s (set by the pipeline)", paste(dup, collapse = ", ")))
message(sprintf("Running getAlternativePromoters (reference: %s; %s)", reference, arg_text))

# As in v0.1.0, an internal failure of getAlternativePromoters (e.g. too little data) is
# reported but does not fail the pipeline; invalid user arguments fail above.
alt <- tryCatch(
    do.call(getAlternativePromoters, c(list(result = result, referenceCondition = reference), extra)),
    error = function(e) { message("getAlternativePromoters failed: ", conditionMessage(e)); NULL })
if (is.null(alt)) quit(save = "no", status = 0)
up   <- alt$upReg
down <- alt$downReg
message(sprintf("Alternative promoters up in %s: %d; down in %s: %d", reference, nrow(up), reference, nrow(down)))
if (nrow(up) + nrow(down) > 0) {
    if (nrow(up) > 0)   up$direction   <- "up_in_reference"
    if (nrow(down) > 0) down$direction <- "down_in_reference"
    combined <- rbind(up, down)
    rdata    <- as.data.frame(rowData(result))
    combined <- merge(combined,
                      rdata[, c("promoterId", "seqnames", "start", "strand", "internalPromoter", "promoterPosition")],
                      by = "promoterId", all.x = TRUE)
    combined <- combined[order(combined$padjAbs), ]
    write.table(combined, "alternative_promoters.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
} else {
    message("No alternative promoters found with the current thresholds (--proactiv_altprom_args)")
}
