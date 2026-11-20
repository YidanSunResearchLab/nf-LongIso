#!/usr/bin/env Rscript
# NextLongIso: differential TE expression with DESeq2.
#
# TE counts are modelled together with the transcript counts, so that size factors are
# estimated from the whole transcriptome (same model as v0.1.0).
#
# Usage:
#   run_te_deseq2.R --te-counts=TE_family_counts.tsv --tx-counts=counts_transcript.txt \
#               --design=design.csv --design-formula="~ condition" --type=Family --reference=control
#
# Outputs (<type> = Family or Instance)
#   DE_TE_<type>_results.csv                  (exactly 2 conditions)
#   DE_TE_<type>_results_<cond>_vs_<ref>.csv  (one per non-reference condition when > 2)
#   norm_counts_<type>.csv, pca_<type>.pdf, volcano_<type>.pdf (one page per comparison)

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
})

# Shared helpers (parse_args, design handling) from bin/nli_design.R
.script <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "nli_design.R"))
opt <- parse_args(commandArgs(trailingOnly = TRUE))
for (k in c("te-counts", "tx-counts", "design", "design-formula", "type", "reference")) {
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
reference  <- opt$reference
coldata    <- design_coldata(design, opt[["design-formula"]], reference)
conditions <- levels(coldata$condition)
formula    <- design_formula(opt[["design-formula"]], coldata)
check_full_rank(formula, coldata)

# ---- 2. DESeq2 -------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(combined), colData = coldata, design = formula)

if (length(conditions) > 1) {
    dds <- DESeq(dds)
    pdf(sprintf("volcano_%s.pdf", type))
    for (cond in setdiff(conditions, reference)) {
        res <- as.data.frame(results(dds, contrast = c("condition", cond, reference)))[te_ids, ]
        res$ID_JOIN <- rownames(res)
        fname <- if (length(conditions) == 2) sprintf("DE_TE_%s_results.csv", type) else
            sprintf("DE_TE_%s_results_%s_vs_%s.csv", type, cond, reference)
        write.csv(merge(te_meta, res, by = "ID_JOIN"), fname, row.names = FALSE)
        p <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = padj < 0.05)) +
             geom_point(alpha = 0.4) + theme_minimal()
        # One page per comparison; with two conditions the plot is drawn exactly as in v0.1.0 (no title)
        if (length(conditions) > 2) p <- p + ggtitle(sprintf("TE %s: %s vs %s", type, cond, reference))
        print(p)
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
