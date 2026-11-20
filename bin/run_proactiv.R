#!/usr/bin/env Rscript
# NextLongIso: promoter activity (proActiv) from long-read junctions.
#
# Usage:
#   run_proactiv.R --gtf=cleaned_annotations.gtf --design=design.csv --species=Homo_sapiens \
#              --reference=control --altprom-args-file=altprom_args.txt --cpus=8
#
# Junction tables (<sample>.junctions.tsv from EXTRACT_JUNCTIONS: chr, start (0-based), end, strand,
# reads, motif, annotated) are read from the working directory and converted to STAR SJ.out.tab files,
# which proActiv accepts as input. The strand comes from the minimap2 ts tag (+ -> 1, - -> 2);
# junctions of unknown strand ('.') are not used: proActiv counts a junction for the intron with the
# same coordinates on the strand given in the file, and STAR's code 0 would be read as '+'.
# --altprom-args-file contains extra R arguments for getAlternativePromoters(), e.g.
#   minAbs = 0.25, minRel = 0.05, maxPval = 0.05, promoterFC = 2.0, geneFC = 1.5
#
# Outputs: promoter_activity.tsv, promoter_counts.tsv, alternative_promoters.tsv (if any)

suppressPackageStartupMessages({
    library(proActiv)
    library(GenomicFeatures)
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

# ---- 2. Junction tables -> STAR SJ.out.tab ----------------------------------------------------------
# STAR intron motif codes (genome orientation): 0 non-canonical, 1 GT/AG, 2 CT/AC, 3 GC/AG, 4 CT/GC,
# 5 AT/AC, 6 GT/AT. The junction tables give the motif in transcript orientation, so a minus-strand
# GT-AG junction is code 2, etc.
motif_code <- function(motif, strand) {
    base <- c("GT-AG" = 1L, "GC-AG" = 3L, "AT-AC" = 5L)[motif]
    ifelse(is.na(base), 0L, base + as.integer(strand == "-"))
}
tsv_files <- sort(list.files(pattern = "[.]junctions[.]tsv$"))
message(sprintf("Found %d junction tables", length(tsv_files)))
junc_files   <- character(0)
sample_names <- character(0)
for (tsv in tsv_files) {
    sname <- sub("[.]junctions[.]tsv$", "", tsv)
    j <- read.delim(tsv, colClasses = c(chr = "character", strand = "character", motif = "character",
                                        annotated = "character"), na.strings = character(0))
    unstranded <- j$strand == "."
    j <- j[!unstranded, , drop = FALSE]
    if (nrow(j) == 0) {
        warning(sprintf("No stranded junctions in %s - sample skipped", tsv))
        next
    }
    star_df <- data.frame(
        chrom     = j$chr,
        start     = j$start + 1L,                       # first intron base, 1-based
        end       = j$end,                              # last intron base, 1-based
        strand    = ifelse(j$strand == "+", 1L, 2L),
        motif     = motif_code(j$motif, j$strand),
        annotated = as.integer(j$annotated == "yes"),
        unique    = as.integer(j$reads),
        multi     = 0L,
        overhang  = 10L                                 # not used by proActiv
    )
    out_name <- paste0(sname, ".SJ.out.tab")
    write.table(star_df, out_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    junc_files   <- c(junc_files, out_name)
    sample_names <- c(sample_names, sname)
    message(sprintf("  %s: %d junctions (+: %d, -: %d; %d of unknown strand not used)", sname, nrow(star_df),
                    sum(star_df$strand == 1L), sum(star_df$strand == 2L), sum(unstranded)))
}
if (length(junc_files) == 0) stop("No junction files created - aborting.")

# ---- 3. Sample -> condition (exact sample names from the validated samplesheet) ---------------
design <- read_design(opt$design)
unknown <- setdiff(sample_names, design$sample)
if (length(unknown) > 0) stop(sprintf("Junction tables without a samplesheet entry: %s", paste(unknown, collapse = ", ")))
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
