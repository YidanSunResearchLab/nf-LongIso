#!/usr/bin/env Rscript
# NextLongIso: differential TE expression with DESeq2.
#
# TE counts are modelled together with the transcript counts, so that size factors are
# estimated from the whole transcriptome (same model as v0.1.0).
#
# Usage:
#   te_deseq2.R --te-counts TE_family_counts.tsv --tx-counts counts_transcript.txt \
#               --design design.csv --type Family --reference control
#
# Outputs (<type> = Family or Instance)
#   DE_TE_<type>_results.csv                  (exactly 2 conditions)
#   DE_TE_<type>_results_<cond>_vs_<ref>.csv  (one per non-reference condition when > 2)
#   norm_counts_<type>.csv, pca_<type>.pdf, volcano_<type>.pdf (one page per comparison)

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
})

parse_args <- function(args) {
    if (length(args) %% 2 != 0) stop("Arguments must be given as --key value pairs")
    keys <- sub("^--", "", args[c(TRUE, FALSE)])
    setNames(as.list(args[c(FALSE, TRUE)]), keys)
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
for (k in c("te-counts", "tx-counts", "design", "type", "reference")) {
    if (is.null(opt[[k]])) stop(sprintf("Missing required argument --%s", k))
}
type <- opt$type
if (!type %in% c("Family", "Instance")) stop("--type must be 'Family' or 'Instance'")

strip_suffix <- function(x) sub("[.]sorted([.]bam)?$", "", x)

# ---- 1. Load ---------------------------------------------------------------------------------
te_data <- read.delim(opt[["te-counts"]], check.names = FALSE, stringsAsFactors = FALSE)
tx_data <- read.delim(opt[["tx-counts"]], check.names = FALSE, stringsAsFactors = FALSE)
design  <- read.csv(opt$design, stringsAsFactors = FALSE)

if (type == "Family") {
    rownames(te_data) <- make.unique(as.character(te_data[, 1]))
    te_meta <- te_data[, 1, drop = FALSE]
    te_mat  <- te_data[, -1, drop = FALSE]
} else {
    rownames(te_data) <- make.unique(paste(te_data[, 1], te_data[, 2], te_data[, 3], te_data[, 4], sep = "_"))
    te_meta <- te_data[, 1:4]
    te_mat  <- te_data[, -(1:4), drop = FALSE]
}
te_ids <- rownames(te_mat)
te_meta$ID_JOIN <- te_ids

rownames(tx_data) <- tx_data$TXNAME
tx_mat <- tx_data[, -(1:2), drop = FALSE]
colnames(te_mat) <- strip_suffix(colnames(te_mat))
colnames(tx_mat) <- strip_suffix(colnames(tx_mat))

missing <- setdiff(design$sample, intersect(colnames(te_mat), colnames(tx_mat)))
if (length(missing) > 0) stop(sprintf("Samples missing from the count tables: %s", paste(missing, collapse = ", ")))
te_mat <- te_mat[, design$sample, drop = FALSE]
tx_mat <- tx_mat[, design$sample, drop = FALSE]

combined   <- rbind(as.matrix(tx_mat), as.matrix(te_mat))
coldata    <- data.frame(row.names = design$sample, condition = factor(design$condition))
conditions <- levels(coldata$condition)
reference  <- opt$reference
if (length(conditions) > 1) {
    if (!reference %in% conditions) stop(sprintf("Reference condition '%s' not found in the design", reference))
    coldata$condition <- relevel(coldata$condition, ref = reference)
}

# ---- 2. DESeq2 -------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(combined), colData = coldata,
                              design = if (length(conditions) > 1) ~ condition else ~ 1)

if (length(conditions) > 1) {
    dds <- DESeq(dds)
    pdf(sprintf("volcano_%s.pdf", type))
    for (cond in setdiff(conditions, reference)) {
        res <- as.data.frame(results(dds, contrast = c("condition", cond, reference)))[te_ids, ]
        res$ID_JOIN <- rownames(res)
        fname <- if (length(conditions) == 2) sprintf("DE_TE_%s_results.csv", type) else
            sprintf("DE_TE_%s_results_%s_vs_%s.csv", type, cond, reference)
        write.csv(merge(te_meta, res, by = "ID_JOIN"), fname, row.names = FALSE)
        print(ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = padj < 0.05)) +
              geom_point(alpha = 0.4) + theme_minimal() + ggtitle(sprintf("TE %s: %s vs %s", type, cond, reference)))
        message(sprintf("TE %s, %s vs %s: %d with padj < 0.05", type, cond, reference, sum(res$padj < 0.05, na.rm = TRUE)))
    }
    invisible(dev.off())
} else {
    dds <- estimateSizeFactors(dds)
}

# ---- 3. Normalised counts and PCA ------------------------------------------------------------
te_nc <- as.data.frame(counts(dds, normalized = TRUE)[te_ids, , drop = FALSE])
te_nc$ID_JOIN <- rownames(te_nc)
write.csv(merge(te_meta, te_nc, by = "ID_JOIN"), sprintf("norm_counts_%s.csv", type), row.names = FALSE)

vsd <- tryCatch(vst(dds, blind = TRUE),
                error = function(e) varianceStabilizingTransformation(dds, blind = TRUE))[te_ids, ]
pdf(sprintf("pca_%s.pdf", type))
if (ncol(dds) > 1) print(plotPCA(vsd, intgroup = "condition")) else plot.new()
invisible(dev.off())
