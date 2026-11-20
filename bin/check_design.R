#!/usr/bin/env Rscript
# NextLongIso: validate the experimental design before any analysis runs.
#
# Usage:
#   check_design.R --design=design.csv --design-formula="~ batch + condition" --reference=control
#
# Stops with a message naming the samplesheet columns involved when the design matrix is not of
# full rank (e.g. a covariate confounded with condition) or a covariate has a single value.

.script <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "nli_design.R"))

opt <- parse_args(commandArgs(trailingOnly = TRUE))
for (k in c("design", "design-formula", "reference")) if (is.null(opt[[k]])) stop(sprintf("Missing required argument --%s", k))

design  <- read_design(opt$design)
coldata <- design_coldata(design, opt[["design-formula"]], opt$reference)
formula <- design_formula(opt[["design-formula"]], coldata)
check_full_rank(formula, coldata)
message(sprintf("Design %s is of full rank for %d samples (%d conditions; reference: %s)",
                deparse(formula), nrow(coldata), nlevels(coldata$condition), opt$reference))
