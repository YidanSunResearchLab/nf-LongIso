# NextLongIso: shared helpers for experimental designs (sourced by the run_*.R scripts with
#   source(file.path(dirname(normalizePath(<this script>)), "nli_design.R"))).
#
# The design is a formula over samplesheet columns whose last term is `condition`,
# e.g. "~ condition" or "~ batch + sex + condition". Covariates are categorical (their values follow
# the sample-name rules, enforced by the pipeline before any R step runs).

# Arguments: --key=value (base R, no optparse)
parse_args <- function(args) {
    out <- list()
    for (a in args) {
        if (!grepl("^--[A-Za-z0-9_-]+=", a)) stop(sprintf("Arguments must be given as --key=value (got '%s')", a))
        out[[sub("^--([^=]+)=.*$", "\\1", a)]] <- sub("^--[^=]+=", "", a)
    }
    out
}

# Covariates of a design string, in formula order, without `condition`
design_covariates <- function(design) {
    terms <- trimws(strsplit(sub("^\\s*~", "", design), "+", fixed = TRUE)[[1]])
    if (length(terms) == 0 || terms[length(terms)] != "condition") {
        stop(sprintf("Design '%s': the last term must be 'condition'", design))
    }
    setdiff(terms, "condition")
}

# Sample table (data.frame, row names = samples) with factors: condition releveled to `reference`,
# covariates as factors. Only the columns used by the design are kept.
design_coldata <- function(design_table, design, reference) {
    covs <- design_covariates(design)
    missing <- setdiff(c("sample", "condition", covs), colnames(design_table))
    if (length(missing) > 0) stop(sprintf("Design columns missing from design.csv: %s", paste(missing, collapse = ", ")))
    cd <- data.frame(row.names = design_table$sample, condition = factor(design_table$condition))
    if (nlevels(cd$condition) > 1) {
        if (!reference %in% levels(cd$condition)) stop(sprintf("Reference condition '%s' not found", reference))
        cd$condition <- relevel(cd$condition, ref = reference)
    }
    for (v in covs) cd[[v]] <- factor(design_table[[v]])
    cd[, c(covs, "condition"), drop = FALSE]
}

# Formula used for the model: with a single condition, `condition` is dropped (intercept-only
# if there are no covariates either)
design_formula <- function(design, coldata) {
    covs <- design_covariates(design)
    terms <- if (nlevels(coldata$condition) > 1) c(covs, "condition") else covs
    as.formula(paste("~", if (length(terms)) paste(terms, collapse = " + ") else "1"))
}

# Stop early, naming the columns involved, if the design matrix is not of full rank
check_full_rank <- function(formula, coldata) {
    for (v in all.vars(formula)) {
        if (nlevels(coldata[[v]]) < 2 && v != "condition") {
            stop(sprintf("Design column '%s' has a single value (%s); remove it from --design", v, levels(coldata[[v]])))
        }
    }
    X <- model.matrix(formula, coldata)
    q <- qr(X)
    if (q$rank < ncol(X)) {
        dependent <- colnames(X)[q$pivot[(q$rank + 1):ncol(X)]]
        # Columns involved: every term whose removal makes the remaining design full rank
        labs <- labels(terms(formula))
        involved <- labs[vapply(labs, function(t) {
            keep <- attr(X, "assign") != match(t, labs)
            Xr <- X[, keep, drop = FALSE]
            qr(Xr)$rank == ncol(Xr)
        }, logical(1))]
        if (length(involved) == 0) involved <- labs
        stop(sprintf(paste0("The design %s is not of full rank (model column(s) %s are linear combinations of ",
                            "other columns). Samplesheet column(s) involved: %s - a covariate is confounded with ",
                            "condition or with another covariate (e.g. all samples of a condition in the same ",
                            "batch). Remove or recode one of them in --design."),
                     deparse(formula), paste(dependent, collapse = ", "), paste(involved, collapse = ", ")))
    }
    invisible(TRUE)
}


# DESeq() with the fallback recommended by DESeq2 when the dispersion trend cannot be fitted
# (all gene-wise estimates within two orders of magnitude, typical for very small tables or designs
# with few residual degrees of freedom): use the gene-wise dispersion estimates, then the Wald test.
run_deseq <- function(dds) {
    tryCatch(DESeq2::DESeq(dds), error = function(e) {
        if (!grepl("dispersion estimates are within 2 orders of magnitude", conditionMessage(e))) stop(e)
        warning("DESeq2 could not fit a dispersion trend; using gene-wise dispersion estimates")
        dds <- DESeq2::estimateDispersionsGeneEst(dds)
        DESeq2::dispersions(dds) <- S4Vectors::mcols(dds)$dispGeneEst
        DESeq2::nbinomWaldTest(dds)
    })
}
