#!/usr/bin/env Rscript
# NextLongIso: differential expression with DESeq2 at gene or transcript level.
#
# Usage:
#   deseq2.R --counts counts_gene.txt --design design.csv --level gene --reference control
#
# Inputs
#   --counts     Bambu count table: counts_gene.txt (GENEID, samples...) for --level gene,
#                counts_transcript.txt (TXNAME, GENEID, samples...) for --level transcript
#   --design     CSV with columns sample,condition
#   --level      gene | transcript
#   --reference  reference condition (denominator of every fold change)
#
# Outputs
#   DE_results.csv                    (exactly 2 conditions)
#   DE_results_<cond>_vs_<ref>.csv    (one per non-reference condition when > 2 conditions)
#   normalized_counts.csv, pca_plot.pdf
#   deseq2_<level>_mqc.tsv            (MultiQC custom-content table of significant counts)

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
for (k in c("counts", "design", "level", "reference")) {
    if (is.null(opt[[k]])) stop(sprintf("Missing required argument --%s", k))
}
if (!opt$level %in% c("gene", "transcript")) stop("--level must be 'gene' or 'transcript'")

# Sample names are validated by the samplesheet schema; only the exact BAM suffix is removed
strip_suffix <- function(x) sub("[.]sorted([.]bam)?$", "", x)

# ---- 1. Load data -------------------------------------------------------------------------
tab <- read.delim(opt$counts, check.names = FALSE, stringsAsFactors = FALSE)
if (opt$level == "transcript") {
    ids      <- tab$TXNAME
    gene_map <- setNames(tab$GENEID, tab$TXNAME)
    mat      <- tab[, -(1:2), drop = FALSE]
} else {
    ids      <- tab$GENEID
    gene_map <- NULL
    mat      <- tab[, -1, drop = FALSE]
}
colnames(mat) <- strip_suffix(colnames(mat))
rownames(mat) <- ids

design <- read.csv(opt$design, stringsAsFactors = FALSE)
missing <- setdiff(design$sample, colnames(mat))
if (length(missing) > 0) {
    stop(sprintf("Samples missing from the count table: %s (columns: %s)",
                 paste(missing, collapse = ", "), paste(colnames(mat), collapse = ", ")))
}
mat <- round(as.matrix(mat[, design$sample, drop = FALSE]))

coldata    <- data.frame(row.names = design$sample, condition = factor(design$condition))
conditions <- levels(coldata$condition)
reference  <- opt$reference
if (length(conditions) > 1) {
    if (!reference %in% conditions) stop(sprintf("Reference condition '%s' not found in the design", reference))
    coldata$condition <- relevel(coldata$condition, ref = reference)
}
message(sprintf("Level: %s; conditions: %s; reference: %s",
                opt$level, paste(conditions, collapse = ", "), reference))

with_ids <- function(df) {
    out <- data.frame(id = rownames(df), stringsAsFactors = FALSE)
    if (!is.null(gene_map)) out$gene_id <- unname(gene_map[out$id])
    cbind(out, df, row.names = NULL)
}

# ---- 2. Model and PCA -----------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = mat, colData = coldata,
                              design = if (length(conditions) > 1) ~ condition else ~ 1)

vsd <- tryCatch(vst(dds, blind = TRUE),
                error = function(e) varianceStabilizingTransformation(dds, blind = TRUE))
pdf("pca_plot.pdf", width = 7, height = 5)
if (ncol(dds) > 1) {
    print(plotPCA(vsd, intgroup = "condition") + theme_minimal())
} else {
    plot.new(); text(0.5, 0.5, "PCA requires > 1 sample")
}
invisible(dev.off())

# ---- 3. Differential expression --------------------------------------------------------------
level_label <- if (opt$level == "gene") "genes" else "transcripts"
mqc_rows <- list()
if (length(conditions) > 1) {
    dds <- DESeq(dds)
    others <- setdiff(conditions, reference)
    for (cond in others) {
        res <- as.data.frame(results(dds, contrast = c("condition", cond, reference)))
        out <- with_ids(res)
        out <- out[order(out$padj, out$id, na.last = TRUE), ]
        fname <- if (length(conditions) == 2) "DE_results.csv" else sprintf("DE_results_%s_vs_%s.csv", cond, reference)
        write.csv(out, fname, row.names = FALSE)
        sig <- !is.na(res$padj) & res$padj < 0.05
        mqc_rows[[cond]] <- data.frame(
            Comparison = sprintf("%s_vs_%s", cond, reference),
            Tested     = sum(!is.na(res$padj)),
            Significant = sum(sig),
            Up         = sum(sig & res$log2FoldChange > 0),
            Down       = sum(sig & res$log2FoldChange < 0),
            check.names = FALSE
        )
        message(sprintf("%s vs %s: %d %s with padj < 0.05", cond, reference, sum(sig), level_label))
    }
} else {
    dds <- estimateSizeFactors(dds)
    write.csv(data.frame(Message = "Single condition - no DE"), "DE_results.csv", row.names = FALSE)
}

# ---- 4. Normalised counts --------------------------------------------------------------------
write.csv(with_ids(as.data.frame(counts(dds, normalized = TRUE))), "normalized_counts.csv", row.names = FALSE)

# ---- 5. MultiQC custom content ----------------------------------------------------------------
mqc_file <- sprintf("deseq2_%s_mqc.tsv", opt$level)
writeLines(c(
    sprintf("# id: 'deseq2_%s'", opt$level),
    sprintf("# section_name: 'DESeq2: %s level'", opt$level),
    sprintf("# description: 'Differentially expressed %s (padj < 0.05) per comparison; reference condition: %s.'", level_label, reference),
    "# plot_type: 'table'"
), mqc_file)
mqc <- if (length(mqc_rows) > 0) do.call(rbind, mqc_rows) else
    data.frame(Comparison = "none (single condition)", Tested = NA, Significant = NA, Up = NA, Down = NA)
suppressWarnings(write.table(mqc, mqc_file, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE))
