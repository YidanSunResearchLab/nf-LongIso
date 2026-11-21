#!/usr/bin/env Rscript
# NextLongIso: differential feature usage within genes with DEXSeq (used for PAS clusters and TSS
# clusters).
#
# Usage:
#   usage_dexseq.R --counts=pas_counts.tsv --design=design.csv --design-formula="~ condition" \
#                  --reference=control --min-count=10 --prefix=apa --cpus=4
#
# Inputs
#   --counts     TSV: feature_id, gene_id, one column per sample (raw counts)
#   --design     design.csv (sample, condition, covariates)
#   --design-formula  e.g. "~ batch + condition"; the DEXSeq models are
#                full    ~ sample + exon + <covariate>:exon + condition:exon
#                reduced ~ sample + exon + <covariate>:exon
#   --numeric-covariates  optional: covariates of the design that are numeric (comma-separated)
#   --min-count  features with fewer reads over all samples are removed; then genes with < 2
#                features are removed
#   --usage-label  name of the usage measure in the output columns (default 'usage'; 'PAU' for APA)
#
# Output: <prefix>_usage.tsv, one row per (comparison, feature):
#   condition_1 (reference), condition_2, feature_id, gene_id, pvalue, padj (feature level),
#   gene_qvalue (perGeneQValue), log2fold (DEXSeq), mean_<label>_1, mean_<label>_2,
#   d<label> = mean_<label>_2 - mean_<label>_1, <label>_<sample> (feature share of the gene's reads)
# Each comparison (condition vs reference) is checked separately: it is skipped, with a warning, when
# either condition has < 2 replicates, its design is not of full rank or leaves no residual degrees of
# freedom, or no gene has >= 2 features with reads. When no comparison can be tested (or there are
# < 2 conditions) the file has a single 'Message' column giving the reasons.

.script <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "nli_design.R"))
opt <- parse_args(commandArgs(trailingOnly = TRUE))
for (k in c("counts", "design", "design-formula", "reference", "min-count", "prefix")) {
    if (is.null(opt[[k]])) stop(sprintf("Missing required argument --%s", k))
}
cpus      <- as.integer(if (is.null(opt$cpus)) 1 else opt$cpus)
label     <- if (is.null(opt[["usage-label"]])) "usage" else opt[["usage-label"]]
min_count <- as.numeric(opt[["min-count"]])
out_file  <- sprintf("%s_usage.tsv", opt$prefix)
reference <- opt$reference

placeholder <- function(msg) {
    message(msg)
    write.table(data.frame(Message = msg), out_file, sep = "\t", quote = FALSE, row.names = FALSE)
    quit(save = "no", status = 0)
}

# ---- 1. Inputs and replicate checks ---------------------------------------------------------------
tab    <- read.delim(opt$counts, check.names = FALSE, stringsAsFactors = FALSE)
design <- read_design(opt$design)
missing <- setdiff(design$sample, colnames(tab))
if (length(missing) > 0) stop(sprintf("Samples missing from %s: %s", opt$counts, paste(missing, collapse = ", ")))

n_rep <- table(design$condition)
if (length(n_rep) < 2) placeholder(sprintf("Differential usage not tested: %d condition(s); at least 2 are required.", length(n_rep)))
if (!reference %in% names(n_rep)) stop(sprintf("Reference condition '%s' not found", reference))

# ---- 2. Filters ----------------------------------------------------------------------------------
mat <- as.matrix(tab[, design$sample, drop = FALSE])
mat[is.na(mat)] <- 0
mat <- round(mat)
keep <- rowSums(mat) >= min_count
tab <- tab[keep, , drop = FALSE]; mat <- mat[keep, , drop = FALSE]
n_feat <- table(tab$gene_id)
keep <- tab$gene_id %in% names(n_feat)[n_feat >= 2]
tab <- tab[keep, , drop = FALSE]; mat <- mat[keep, , drop = FALSE]
message(sprintf("Features after filters (>= %g reads, genes with >= 2 features): %d in %d genes",
                min_count, nrow(tab), length(unique(tab$gene_id))))
if (nrow(tab) == 0) placeholder(sprintf("Differential usage not tested: no gene has >= 2 features with >= %g reads.", min_count))

# Per-sample usage (feature share of the gene's reads)
gene_tot <- rowsum(mat, tab$gene_id)[tab$gene_id, , drop = FALSE]
usage <- mat / gene_tot
usage[!is.finite(usage)] <- NA

# ---- 3. DEXSeq per comparison ---------------------------------------------------------------------
suppressPackageStartupMessages(library(DEXSeq))
bp <- if (cpus > 1) BiocParallel::MulticoreParam(cpus) else BiocParallel::SerialParam()
covariates <- design_covariates(opt[["design-formula"]])
others <- sort(setdiff(names(n_rep), reference))

fit_dexseq <- function(dxd, full, reduced) {
    dxd <- estimateSizeFactors(dxd)
    fit <- NULL
    for (ft in c("parametric", "local", "mean")) {
        fit <- tryCatch(estimateDispersions(dxd, formula = full, fitType = ft, BPPARAM = bp), error = function(e) {
            message(sprintf("Dispersion fit '%s' failed: %s", ft, conditionMessage(e))); NULL })
        if (!is.null(fit)) { if (ft != "parametric") warning(sprintf("DEXSeq dispersion fit used fitType = '%s'", ft)); break }
    }
    if (is.null(fit)) return(NULL)
    fit <- testForDEU(fit, fullModel = full, reducedModel = reduced, BPPARAM = bp)
    estimateExonFoldChanges(fit, fitExpToVar = "condition", BPPARAM = bp)
}

results <- list()
skipped <- character(0)   # comparisons that could not be tested, with the reason
skip <- function(cond, why) {
    msg <- sprintf("%s vs %s: %s", cond, reference, why)
    message("Skipped: ", msg); warning(msg)
    skipped <<- c(skipped, msg)
}
for (cond in others) {
    # Each comparison is checked on its own: a problem in one comparison does not stop the others
    if (n_rep[[reference]] < 2 || n_rep[[cond]] < 2) {
        skip(cond, sprintf("needs >= 2 replicates per condition (%s: %d, %s: %d)", reference, n_rep[[reference]],
                           cond, n_rep[[cond]]))
        next
    }
    smp <- design$sample[design$condition %in% c(reference, cond)]
    sd  <- design[match(smp, design$sample), , drop = FALSE]
    cd  <- design_coldata(sd, opt[["design-formula"]], reference, numeric_covariates(opt[["numeric-covariates"]]))
    covs <- covariates[vapply(covariates, function(v) n_distinct(cd[[v]]) > 1, logical(1))]
    if (length(covs) < length(covariates)) {
        warning(sprintf("%s vs %s: covariate(s) %s constant within this comparison and dropped",
                        cond, reference, paste(setdiff(covariates, covs), collapse = ", ")))
    }
    cd <- droplevels(cd[, c(covs, "condition"), drop = FALSE])
    cmp_formula <- as.formula(paste("~", paste(c(covs, "condition"), collapse = " + ")))
    problem <- rank_problem(cmp_formula, cd)
    if (!is.null(problem)) { skip(cond, problem); next }
    if (residual_df(cmp_formula, cd) < 1) {
        skip(cond, sprintf("the design %s leaves no residual degrees of freedom in this comparison", deparse(cmp_formula)))
        next
    }
    full    <- as.formula(paste("~ sample + exon", paste0(c(sprintf(" + %s:exon", covs), " + condition:exon"), collapse = "")))
    reduced <- as.formula(paste("~ sample + exon", paste0(sprintf(" + %s:exon", covs), collapse = "")))
    message(sprintf("%s vs %s: full %s; reduced %s", cond, reference, deparse(full), deparse(reduced)))

    sub_mat <- mat[, smp, drop = FALSE]
    # Features or genes without reads in this comparison cannot be tested
    ok <- rowSums(sub_mat) > 0
    ok <- ok & tab$gene_id %in% names(which(tapply(ok, tab$gene_id, sum) >= 2))
    if (sum(ok) == 0) { skip(cond, "no gene with >= 2 features with reads in this comparison"); next }
    dxd <- DEXSeqDataSet(countData = sub_mat[ok, , drop = FALSE], sampleData = cd, design = full,
                         featureID = tab$feature_id[ok], groupID = tab$gene_id[ok])
    dxd <- fit_dexseq(dxd, full, reduced)
    if (is.null(dxd)) { skip(cond, "dispersions could not be estimated"); next }
    res <- DEXSeqResults(dxd)
    q   <- perGeneQValue(res)
    lfc_col <- grep("^log2fold_", colnames(res), value = TRUE)[1]
    out <- data.frame(condition_1 = reference, condition_2 = cond,
                      feature_id = res$featureID, gene_id = res$groupID,
                      pvalue = res$pvalue, padj = res$padj,
                      gene_qvalue = unname(q[res$groupID]),
                      log2fold = if (is.na(lfc_col)) NA_real_ else res[[lfc_col]],
                      stringsAsFactors = FALSE)
    idx <- match(paste(out$gene_id, out$feature_id), paste(tab$gene_id, tab$feature_id))
    u <- usage[idx, smp, drop = FALSE]
    colnames(u) <- paste0(label, "_", smp)
    m1 <- rowMeans(usage[idx, sd$sample[sd$condition == reference], drop = FALSE], na.rm = TRUE)
    m2 <- rowMeans(usage[idx, sd$sample[sd$condition == cond], drop = FALSE], na.rm = TRUE)
    out[[paste0("mean_", label, "_1")]] <- m1
    out[[paste0("mean_", label, "_2")]] <- m2
    out[[paste0("d", label)]] <- m2 - m1
    results[[cond]] <- cbind(out, u)
    message(sprintf("%s vs %s: %d features in %d genes tested; %d genes with q < 0.05",
                    cond, reference, nrow(out), length(unique(out$gene_id)), sum(q < 0.05, na.rm = TRUE)))
}
if (length(results) == 0) placeholder(paste0("Differential usage not tested. ", paste(skipped, collapse = "; ")))
if (length(skipped) > 0) message("Comparisons not tested: ", paste(skipped, collapse = "; "))

# Comparisons may involve different samples: align the usage_<sample> columns
all_cols <- unique(unlist(lapply(results, colnames)))
results <- lapply(results, function(d) { for (c in setdiff(all_cols, colnames(d))) d[[c]] <- NA; d[, all_cols] })
final <- do.call(rbind, results)
final <- final[order(final$condition_2, final$gene_qvalue, final$padj, final$gene_id, final$feature_id, na.last = TRUE), ]
write.table(final, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
