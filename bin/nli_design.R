# NextLongIso: shared helpers for experimental designs (sourced by the run_*.R scripts with
#   source(file.path(dirname(normalizePath(<this script>)), "nli_design.R"))).
#
# The design is a formula over samplesheet columns whose last term is `condition`,
# e.g. "~ condition" or "~ batch + age + condition". Covariates are categorical (factors) unless they are
# listed in --numeric_covariates (numbers, centred and scaled to unit standard deviation). Values are
# checked by the pipeline before any R step runs.

# Arguments: --key=value (base R, no optparse)
parse_args <- function(args) {
    out <- list()
    for (a in args) {
        if (!grepl("^--[A-Za-z0-9_-]+=", a)) stop(sprintf("Arguments must be given as --key=value (got '%s')", a))
        out[[sub("^--([^=]+)=.*$", "\\1", a)]] <- sub("^--[^=]+=", "", a)
    }
    out
}

# Design table (design.csv: sample, condition, covariates). All values are read as text: names such as
# T, F, TRUE or NA are valid sample/condition/covariate values and must not become logical or missing.
read_design <- function(path) {
    read.csv(path, colClasses = "character", na.strings = character(0), stringsAsFactors = FALSE)
}

# Covariates of a design string, in formula order, without `condition`
design_covariates <- function(design) {
    terms <- trimws(strsplit(sub("^\\s*~", "", design), "+", fixed = TRUE)[[1]])
    if (length(terms) == 0 || terms[length(terms)] != "condition") {
        stop(sprintf("Design '%s': the last term must be 'condition'", design))
    }
    setdiff(terms, "condition")
}

# Numeric covariates from the --numeric-covariates argument (comma-separated; NULL or "" = none)
numeric_covariates <- function(x) {
    if (is.null(x) || !nzchar(trimws(x))) return(character(0))
    trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

# Number of distinct values of a design column (factor or numeric)
n_distinct <- function(x) length(unique(as.character(x)))

# Sample table (data.frame, row names = samples): condition as a factor releveled to `reference`;
# covariates as factors, except those in `numeric`, which are centred and scaled (as recommended by
# DESeq2 for GLM convergence; the condition effects do not change). Only the design columns are kept.
design_coldata <- function(design_table, design, reference, numeric = character(0)) {
    covs <- design_covariates(design)
    unknown <- setdiff(numeric, covs)
    if (length(unknown) > 0) stop(sprintf("Numeric covariate(s) %s are not in the design %s", paste(unknown, collapse = ", "), design))
    missing <- setdiff(c("sample", "condition", covs), colnames(design_table))
    if (length(missing) > 0) stop(sprintf("Design columns missing from design.csv: %s", paste(missing, collapse = ", ")))
    cd <- data.frame(row.names = design_table$sample, condition = factor(design_table$condition))
    if (nlevels(cd$condition) > 1) {
        if (!reference %in% levels(cd$condition)) stop(sprintf("Reference condition '%s' not found", reference))
        cd$condition <- relevel(cd$condition, ref = reference)
    }
    for (v in covs) {
        if (v %in% numeric) {
            x <- suppressWarnings(as.numeric(design_table[[v]]))
            if (anyNA(x)) stop(sprintf("Numeric covariate '%s' has non-numeric value(s): %s", v,
                                       paste(unique(design_table[[v]][is.na(x)]), collapse = ", ")))
            # A constant column (possible within a subset of samples) is left constant: rank_problem()
            # reports it, and the DEXSeq usage test drops it for that comparison
            cd[[v]] <- if (n_distinct(x) > 1) as.numeric(scale(x)) else rep(0, length(x))
        } else {
            cd[[v]] <- factor(design_table[[v]])
        }
    }
    cd[, c(covs, "condition"), drop = FALSE]
}

# Formula used for the model: with a single condition, `condition` is dropped (intercept-only
# if there are no covariates either)
design_formula <- function(design, coldata) {
    covs <- design_covariates(design)
    terms <- if (nlevels(coldata$condition) > 1) c(covs, "condition") else covs
    as.formula(paste("~", if (length(terms)) paste(terms, collapse = " + ") else "1"))
}

# NULL if the design matrix is of full rank, otherwise a message naming the samplesheet columns involved
rank_problem <- function(formula, coldata) {
    for (v in all.vars(formula)) {
        if (n_distinct(coldata[[v]]) < 2 && v != "condition") {
            return(sprintf("Design column '%s' has a single value (%s); remove it from --design", v, unique(as.character(coldata[[v]]))))
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
        return(sprintf(paste0("The design %s is not of full rank (model column(s) %s are linear combinations of ",
                            "other columns). Samplesheet column(s) involved: %s - a covariate is confounded with ",
                            "condition or with another covariate (e.g. all samples of a condition in the same ",
                            "batch). Remove or recode one of them in --design."),
                     deparse(formula), paste(dependent, collapse = ", "), paste(involved, collapse = ", ")))
    }
    NULL
}

# Stop early, naming the columns involved, if the design matrix is not of full rank
check_full_rank <- function(formula, coldata) {
    msg <- rank_problem(formula, coldata)
    if (!is.null(msg)) stop(msg)
    invisible(TRUE)
}

# Residual degrees of freedom of the design (samples - model coefficients). With 0, dispersions cannot
# be estimated (e.g. one sample per condition, or a covariate using up the replicates), so the
# differential tests are skipped with a message instead of failing.
residual_df <- function(formula, coldata) {
    X <- model.matrix(formula, coldata)
    nrow(X) - qr(X)$rank
}


# DESeq() with the fallback recommended by DESeq2 when the dispersion trend cannot be fitted
# (all gene-wise estimates within two orders of magnitude, typical for very small tables or designs
# with few residual degrees of freedom): use the gene-wise dispersion estimates, then the Wald test.
run_deseq <- function(dds) {
    tryCatch(DESeq2::DESeq(dds), error = function(e) {
        if (!grepl("dispersion estimates are within 2 orders of magnitude", conditionMessage(e))) stop(e)
        warning("DESeq2 could not fit a dispersion trend; using gene-wise dispersion estimates")
        # DESeq() failed after its own size-factor step, which is not kept in `dds`: estimate them here
        # unless they were set by the caller (TE analysis: size factors from the gene counts)
        if (is.null(DESeq2::sizeFactors(dds))) dds <- DESeq2::estimateSizeFactors(dds)
        dds <- DESeq2::estimateDispersionsGeneEst(dds)
        DESeq2::dispersions(dds) <- S4Vectors::mcols(dds)$dispGeneEst
        DESeq2::nbinomWaldTest(dds)
    })
}
