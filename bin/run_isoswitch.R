#!/usr/bin/env Rscript
# NextLongIso: isoform switching (IsoformSwitchAnalyzeR) and alternative TSS/TTS usage from Bambu
# isoforms. (Poly(A)-site usage is analysed by the APA process: apa_sites.R, usage_dexseq.R.)
#
# Usage:
#   run_isoswitch.R --counts=counts_transcript.txt --abundance=CPM_transcript.txt \
#               --gtf=extended_annotations.gtf --design=design.csv --design-formula="~ condition" \
#               --reference=control \
#               --alpha=0.05 --dif-cutoff=0.1
#
# Every non-reference condition is compared with the reference (condition_1 = reference,
# condition_2 = other condition; dIF and dPAU are condition_2 - condition_1).
#
# Outputs
#   isoform_switch_results.csv        isoformFeatures for all comparisons (required)
#   significant_switches.csv          q < alpha and |dIF| > dif-cutoff (required)
#   all_AS_events.csv, atss_results.csv, atts_results.csv
#   switchlist.rds
# With fewer than 2 conditions the two required files contain a single 'Message' column.

.libPaths(unique(c(.libPaths(), "/usr/local/lib/R/site-library", "/usr/lib/R/library")))
suppressPackageStartupMessages({
    library(IsoformSwitchAnalyzeR)
    library(GenomicRanges)
})

# Shared helpers (parse_args, design handling) from bin/nli_design.R
.script <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "nli_design.R"))
opt <- parse_args(commandArgs(trailingOnly = TRUE))
for (k in c("counts", "abundance", "gtf", "design", "design-formula", "reference", "alpha", "dif-cutoff")) {
    if (is.null(opt[[k]])) stop(sprintf("Missing required argument --%s", k))
}
alpha      <- as.numeric(opt$alpha)
dif_cutoff <- as.numeric(opt[["dif-cutoff"]])

strip_suffix  <- function(x) sub("[.]sorted([.]bam)?$", "", x)
write_message <- function(file, msg) write.csv(data.frame(Message = msg), file, row.names = FALSE)

# ---- 1. Load ---------------------------------------------------------------------------------
counts_data <- read.delim(opt$counts,    check.names = FALSE, stringsAsFactors = FALSE)
abund_data  <- read.delim(opt$abundance, check.names = FALSE, stringsAsFactors = FALSE)
colnames(counts_data)[1] <- "isoform_id"
colnames(abund_data)[1]  <- "isoform_id"
colnames(counts_data) <- strip_suffix(colnames(counts_data))
colnames(abund_data)  <- strip_suffix(colnames(abund_data))

design  <- read.csv(opt$design, stringsAsFactors = FALSE)
missing <- setdiff(design$sample, intersect(colnames(counts_data), colnames(abund_data)))
if (length(missing) > 0) stop(sprintf("Samples missing from the count tables: %s", paste(missing, collapse = ", ")))

conditions <- unique(design$condition)
if (length(conditions) < 2) {
    msg <- sprintf("Isoform switch analysis skipped: the samplesheet has %d condition(s); at least 2 are required.",
                   length(conditions))
    message(msg)
    write_message("isoform_switch_results.csv", msg)
    write_message("significant_switches.csv", msg)
    quit(save = "no", status = 0)
}
reference <- opt$reference
if (!reference %in% conditions) stop(sprintf("Reference condition '%s' not found in the design", reference))

# The DEXSeq-based switch test needs biological replicates: only conditions with >= 2 samples
# are compared, and the reference itself must have >= 2 samples.
n_rep  <- table(design$condition)
single <- names(n_rep)[n_rep < 2]
others <- sort(setdiff(names(n_rep)[n_rep >= 2], reference))
if (n_rep[[reference]] < 2 || length(others) == 0) {
    msg <- sprintf(paste0("Isoform switch analysis skipped: it needs >= 2 replicates in the reference ",
                          "condition ('%s', %d sample(s)) and in at least one other condition. ",
                          "Conditions with < 2 replicates: %s."),
                   reference, n_rep[[reference]], if (length(single)) paste(single, collapse = ", ") else "none")
    message(msg)
    write_message("isoform_switch_results.csv", msg)
    write_message("significant_switches.csv", msg)
    quit(save = "no", status = 0)
}
if (length(single) > 0) {
    message("Conditions with < 2 replicates are excluded from isoform switching: ", paste(single, collapse = ", "))
}
design <- design[design$condition %in% c(reference, others), ]
covariates <- design_covariates(opt[["design-formula"]])
if (length(covariates) > 0) {
    # Covariates are passed as extra designMatrix columns (used by the DEXSeq-based switch test)
    cd <- design_coldata(design, opt[["design-formula"]], reference)
    check_full_rank(design_formula(opt[["design-formula"]], cd), cd)
    message("Covariates: ", paste(covariates, collapse = ", "))
}
comparisons <- data.frame(condition_1 = reference, condition_2 = others, stringsAsFactors = FALSE)
message("Comparisons: ", paste(comparisons$condition_2, "vs", comparisons$condition_1, collapse = "; "))

counts_data <- counts_data[, c("isoform_id", design$sample)]
abund_data  <- abund_data[,  c("isoform_id", design$sample)]
counts_data[, design$sample] <- round(counts_data[, design$sample])

# ---- 2. Isoform switching and alternative splicing ---------------------------------------------
switchList <- importRdata(
    isoformCountMatrix   = counts_data,
    isoformRepExpression = abund_data,
    designMatrix         = data.frame(sampleID = design$sample, condition = design$condition,
                                      design[, covariates, drop = FALSE], stringsAsFactors = FALSE),
    isoformExonAnnoation = opt$gtf,
    comparisonsToMake    = comparisons,
    showProgress         = FALSE,
    fixStringTieAnnotationProblem = TRUE
)
switchList <- preFilter(switchList)
switchList <- isoformSwitchTestDEXSeq(switchList, reduceToSwitchingGenes = FALSE)
switchList <- analyzeAlternativeSplicing(switchList, onlySwitchingGenes = FALSE)

as_analysis <- switchList$AlternativeSplicingAnalysis
if (!is.null(as_analysis) && nrow(as_analysis) > 0) {
    message("Total AS events: ", nrow(as_analysis))
    write.csv(as_analysis, "all_AS_events.csv", row.names = FALSE)
    atss <- as_analysis[!is.na(as_analysis$ATSS) & as_analysis$ATSS > 0, ]
    atts <- as_analysis[!is.na(as_analysis$ATTS) & as_analysis$ATTS > 0, ]
    message("ATSS events: ", nrow(atss), "; ATTS events: ", nrow(atts))
    if (nrow(atss) > 0) write.csv(atss, "atss_results.csv", row.names = FALSE)
    if (nrow(atts) > 0) write.csv(atts, "atts_results.csv", row.names = FALSE)
} else {
    message("No alternative splicing results")
}

feat <- switchList$isoformFeatures
clean_id <- function(x) trimws(gsub('"', '', sub(".*:", "", as.character(x))))
feat$iso_clean <- clean_id(feat$isoform_id)

# ---- 3. Switch tables --------------------------------------------------------------------------
saveRDS(switchList, "switchlist.rds")
write.csv(feat, "isoform_switch_results.csv", row.names = FALSE)
sig_switch <- feat[which(feat$isoform_switch_q_value < alpha & abs(feat$dIF) > dif_cutoff), ]
write.csv(sig_switch, "significant_switches.csv", row.names = FALSE)
message(sprintf("Significant switches (q < %g, |dIF| > %g): %d", alpha, dif_cutoff, nrow(sig_switch)))
