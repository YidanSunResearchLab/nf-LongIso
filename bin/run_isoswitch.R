#!/usr/bin/env Rscript
# NextLongIso: isoform switching (IsoformSwitchAnalyzeR), alternative TSS/TTS usage and
# poly(A)-site (APA) usage from Bambu isoforms.
#
# Usage:
#   run_isoswitch.R --counts=counts_transcript.txt --abundance=CPM_transcript.txt \
#               --gtf=extended_annotations.gtf --design=design.csv --design-formula="~ condition" \
#               --reference=control \
#               --alpha=0.05 --dif-cutoff=0.1 --seed=42
#
# Every non-reference condition is compared with the reference (condition_1 = reference,
# condition_2 = other condition; dIF and dPAU are condition_2 - condition_1).
#
# Outputs
#   isoform_switch_results.csv        isoformFeatures for all comparisons (required)
#   significant_switches.csv          q < alpha and |dIF| > dif-cutoff (required)
#   all_AS_events.csv, atss_results.csv, atts_results.csv
#   clustered_base_level_PAS_usage.csv  PAS-cluster usage per comparison
#   APA_gene_summary.csv              one row per (comparison, gene) with significant APA
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
for (k in c("counts", "abundance", "gtf", "design", "design-formula", "reference", "alpha", "dif-cutoff", "seed")) {
    if (is.null(opt[[k]])) stop(sprintf("Missing required argument --%s", k))
}
alpha      <- as.numeric(opt$alpha)
dif_cutoff <- as.numeric(opt[["dif-cutoff"]])
seed       <- as.integer(opt$seed)

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

# ---- 3. APA: cluster 3' ends and compare PAS-cluster usage per comparison ---------------------
exons_df <- as.data.frame(switchList$exons)
pas_final <- NULL
if (nrow(exons_df) > 0) {
    message("--- Clustering 3' ends ---")
    pas_raw <- do.call(rbind, lapply(split(exons_df, exons_df$isoform_id), function(x) {
        data.frame(
            isoform_id = as.character(x$isoform_id[1]),
            seqnames   = as.character(x$seqnames[1]),
            strand     = as.character(x$strand[1]),
            pas_coord  = if (x$strand[1] == "+") max(x$end) else min(x$start),
            stringsAsFactors = FALSE
        )
    }))
    pas_gr       <- makeGRangesFromDataFrame(pas_raw, start.field = "pas_coord", end.field = "pas_coord")
    pas_clusters <- resize(reduce(resize(pas_gr, width = 25, fix = "center")), width = 1, fix = "center")
    hits_df <- as.data.frame(findOverlaps(pas_gr, pas_clusters, maxgap = 25))
    hits_df <- hits_df[!duplicated(hits_df$queryHits), ]
    pas_raw$row_idx <- seq_len(nrow(pas_raw))
    pas_mapped <- merge(pas_raw, data.frame(
        row_idx        = hits_df$queryHits,
        pas_cluster_id = paste0("PAS_", hits_df$subjectHits),
        cluster_coord  = start(pas_clusters[hits_df$subjectHits]),
        stringsAsFactors = FALSE
    ), by = "row_idx")
    pas_mapped$iso_clean <- clean_id(pas_mapped$isoform_id)

    set.seed(seed)
    per_comparison <- list()
    for (k in seq_len(nrow(comparisons))) {
        c1 <- comparisons$condition_1[k]
        c2 <- comparisons$condition_2[k]
        # isoformFeatures has one row per isoform and comparison: use only this comparison's rows
        fc  <- feat[feat$condition_1 == c1 & feat$condition_2 == c2, ]
        idx <- match(pas_mapped$iso_clean, fc$iso_clean)
        pm  <- cbind(pas_mapped[!is.na(idx), ],
                     fc[idx[!is.na(idx)], c("gene_id", "IF1", "IF2", "gene_value_1", "gene_value_2")])
        if (nrow(pm) == 0) {
            message(sprintf("APA %s vs %s: no isoforms with PAS clusters", c2, c1))
            next
        }
        pm$gene_id <- ifelse(is.na(pm$gene_id), pm$iso_clean, as.character(pm$gene_id))
        pm$c1 <- pm$IF1 * pm$gene_value_1   # expression of the isoform in condition_1
        pm$c2 <- pm$IF2 * pm$gene_value_2   # expression of the isoform in condition_2

        agg <- aggregate(cbind(IF1, IF2, c1, c2) ~ gene_id + seqnames + cluster_coord + strand + pas_cluster_id,
                         data = pm, FUN = sum, na.rm = TRUE)
        if (nrow(agg) == 0) next
        agg$dPAU        <- agg$IF2 - agg$IF1
        agg$pas_p_value <- NA_real_
        agg$n_sites     <- NA_integer_

        for (g in unique(agg$gene_id)) {
            rows <- which(agg$gene_id == g)
            agg$n_sites[rows] <- length(rows)
            if (length(rows) < 2) next
            mat <- as.matrix(agg[rows, c("c1", "c2")])
            if (all(mat >= 0) && sum(mat) > 10) {
                agg$pas_p_value[rows] <- tryCatch({
                    mat_int <- round(mat)
                    mat_int <- mat_int[rowSums(mat_int) > 0, , drop = FALSE]
                    if (nrow(mat_int) > 1) fisher.test(mat_int, simulate.p.value = TRUE, B = 10000)$p.value else NA_real_
                }, error = function(e) NA_real_)
            }
        }
        agg$pas_q_value <- NA_real_
        ok <- !is.na(agg$pas_p_value)
        if (any(ok)) agg$pas_q_value[ok] <- p.adjust(agg$pas_p_value[ok], method = "BH")

        agg <- cbind(condition_1 = c1, condition_2 = c2, agg, stringsAsFactors = FALSE)
        message(sprintf("APA %s vs %s: %d genes tested, %d PAS clusters with q < %g",
                        c2, c1, length(unique(agg$gene_id[agg$n_sites > 1])),
                        sum(agg$pas_q_value < alpha, na.rm = TRUE), alpha))
        per_comparison[[k]] <- agg
    }
    per_comparison <- Filter(Negate(is.null), per_comparison)
    if (length(per_comparison) > 0) pas_final <- do.call(rbind, per_comparison)
}

if (!is.null(pas_final) && nrow(pas_final) > 0) {
    write.csv(pas_final[order(pas_final$pas_q_value, na.last = TRUE), ],
              "clustered_base_level_PAS_usage.csv", row.names = FALSE)

    # Per comparison and gene: the significant PAS cluster with the largest |dPAU|, so that
    # dPAU, q-value and cluster coordinates all come from the same site.
    sig <- pas_final[!is.na(pas_final$pas_q_value) & pas_final$pas_q_value < alpha, ]
    if (nrow(sig) > 0) {
        # Ties in |dPAU| (common for two-site genes, whose dPAU values are +x and -x) are broken
        # deterministically by q-value and then cluster coordinate
        sig <- sig[order(sig$condition_2, sig$gene_id, -abs(sig$dPAU), sig$pas_q_value, sig$cluster_coord), ]
        top <- sig[!duplicated(sig[, c("condition_1", "condition_2", "gene_id")]),
                   c("condition_1", "condition_2", "gene_id", "seqnames", "strand", "pas_cluster_id",
                     "cluster_coord", "dPAU", "pas_p_value", "pas_q_value", "n_sites")]
        write.csv(top, "APA_gene_summary.csv", row.names = FALSE)
        message("APA_gene_summary.csv: ", nrow(top), " gene(s)")
    } else {
        message("No significant APA (q < ", alpha, "); APA_gene_summary.csv not written")
    }
} else {
    message("No PAS usage could be computed; APA outputs not written")
}

# ---- 4. Switch tables --------------------------------------------------------------------------
saveRDS(switchList, "switchlist.rds")
write.csv(feat, "isoform_switch_results.csv", row.names = FALSE)
sig_switch <- feat[which(feat$isoform_switch_q_value < alpha & abs(feat$dIF) > dif_cutoff), ]
write.csv(sig_switch, "significant_switches.csv", row.names = FALSE)
message(sprintf("Significant switches (q < %g, |dIF| > %g): %d", alpha, dif_cutoff, nrow(sig_switch)))
