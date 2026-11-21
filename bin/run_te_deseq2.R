#!/usr/bin/env Rscript
# NextLongIso: differential TE expression with DESeq2.
#
# Size factors are estimated from the gene-level counts (Bambu counts_gene.txt) and assigned to the
# TE DESeq2 object; TE counts are modelled on their own. (v0.1.0 estimated size factors from
# transcript and TE counts combined, although TE counts are themselves sums of transcript counts.)
#
# Usage:
#   run_te_deseq2.R --te-counts=TE_family_counts.tsv --gene-counts=counts_gene.txt \
#               --design=design.csv --design-formula="~ condition" --type=Family --reference=control
#               [--numeric-covariates=age,rin]   (covariates of the design that are numeric)
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
for (k in c("te-counts", "gene-counts", "design", "design-formula", "type", "reference")) {
    if (is.null(opt[[k]])) stop(sprintf("Missing required argument --%s", k))
}
type <- opt$type
if (!type %in% c("Family", "Instance")) stop("--type must be 'Family' or 'Instance'")

strip_suffix <- function(x) sub("[.]sorted([.]bam)?$", "", x)

# ---- 1. Load ---------------------------------------------------------------------------------
te_data   <- read.delim(opt[["te-counts"]], check.names = FALSE, stringsAsFactors = FALSE)
gene_data <- read.delim(opt[["gene-counts"]], check.names = FALSE, stringsAsFactors = FALSE)
design    <- read_design(opt$design)

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

rownames(gene_data) <- gene_data$GENEID
gene_mat <- gene_data[, -1, drop = FALSE]
colnames(te_mat)   <- strip_suffix(colnames(te_mat))
colnames(gene_mat) <- strip_suffix(colnames(gene_mat))

missing <- setdiff(design$sample, intersect(colnames(te_mat), colnames(gene_mat)))
if (length(missing) > 0) stop(sprintf("Samples missing from the count tables: %s", paste(missing, collapse = ", ")))
te_mat   <- round(as.matrix(te_mat[, design$sample, drop = FALSE]))
gene_mat <- round(as.matrix(gene_mat[, design$sample, drop = FALSE]))

reference  <- opt$reference
coldata    <- design_coldata(design, opt[["design-formula"]], reference, numeric_covariates(opt[["numeric-covariates"]]))
conditions <- levels(coldata$condition)
formula    <- design_formula(opt[["design-formula"]], coldata)
check_full_rank(formula, coldata)

# Why the test cannot be run (written to the results file instead of failing)
no_de <- NULL
if (nrow(te_mat) == 0 || sum(te_mat) == 0) {
    no_de <- sprintf(paste0("No TE %s has any read: check that the chromosome names of --te_gtf match --genome ",
                            "and --gtf (e.g. 'chr1' vs '1')"), tolower(type))
} else if (length(conditions) < 2) {
    no_de <- "Single condition - no DE"
} else if (residual_df(formula, coldata) < 1) {
    no_de <- sprintf("No DE: the design %s leaves no residual degrees of freedom (no replicates)", deparse(formula))
}
if (!is.null(no_de)) message(no_de)
if (nrow(te_mat) == 0 || sum(te_mat) == 0) {
    # No TE read at all: placeholder outputs (DESeq2 cannot build an object from an empty or all-zero
    # table); the normalised counts are all 0
    write.csv(data.frame(Message = no_de), sprintf("DE_TE_%s_results.csv", type), row.names = FALSE)
    write.csv(merge(te_meta, data.frame(ID_JOIN = te_ids, te_mat, check.names = FALSE), by = "ID_JOIN"),
              sprintf("norm_counts_%s.csv", type), row.names = FALSE)
    pdf(sprintf("pca_%s.pdf", type)); plot.new(); text(0.5, 0.5, no_de, cex = 0.7); invisible(dev.off())
    quit(save = "no", status = 0)
}

# ---- 2. Size factors from gene-level counts -----------------------------------------------------
dds_gene <- estimateSizeFactors(DESeqDataSetFromMatrix(countData = gene_mat, colData = coldata, design = formula))
sf <- sizeFactors(dds_gene)
message("Size factors (gene level): ", paste(sprintf("%s=%.3f", names(sf), sf), collapse = ", "))

dds <- DESeqDataSetFromMatrix(countData = te_mat, colData = coldata, design = formula)
sizeFactors(dds) <- sf

# ---- 3. DESeq2 -------------------------------------------------------------------------------
if (!is.null(no_de)) {
    write.csv(data.frame(Message = no_de), sprintf("DE_TE_%s_results.csv", type), row.names = FALSE)
} else {
    dds <- run_deseq(dds)   # keeps the pre-assigned size factors
    pdf(sprintf("volcano_%s.pdf", type))
    for (cond in setdiff(conditions, reference)) {
        res <- as.data.frame(results(dds, contrast = c("condition", cond, reference), alpha = 0.05))[te_ids, ]
        res$ID_JOIN <- rownames(res)
        fname <- if (length(conditions) == 2) sprintf("DE_TE_%s_results.csv", type) else
            sprintf("DE_TE_%s_results_%s_vs_%s.csv", type, cond, reference)
        write.csv(merge(te_meta, res, by = "ID_JOIN"), fname, row.names = FALSE)
        p <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = padj < 0.05)) +
             geom_point(alpha = 0.4) + theme_minimal()
        # One page per comparison; titled when there is more than one comparison
        if (length(conditions) > 2) p <- p + ggtitle(sprintf("TE %s: %s vs %s", type, cond, reference))
        print(p)
        message(sprintf("TE %s, %s vs %s: %d with padj < 0.05", type, cond, reference, sum(res$padj < 0.05, na.rm = TRUE)))
    }
    invisible(dev.off())
}

# ---- 4. Normalised counts and PCA ------------------------------------------------------------
te_nc <- as.data.frame(counts(dds, normalized = TRUE)[te_ids, , drop = FALSE])
te_nc$ID_JOIN <- rownames(te_nc)
write.csv(merge(te_meta, te_nc, by = "ID_JOIN"), sprintf("norm_counts_%s.csv", type), row.names = FALSE)

vsd <- tryCatch(vst(dds, blind = TRUE),
                error = function(e) tryCatch(varianceStabilizingTransformation(dds, blind = TRUE),
                                             error = function(e2) NULL))
pdf(sprintf("pca_%s.pdf", type))
if (ncol(dds) > 1 && !is.null(vsd)) print(plotPCA(vsd, intgroup = "condition")) else plot.new()
invisible(dev.off())
