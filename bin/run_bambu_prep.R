#!/usr/bin/env Rscript
# NextLongIso: prepare Bambu annotations from the reference GTF.
#
# Usage:
#   run_bambu_prep.R --gtf=annotation.gtf --out=bambu_annotations.rds

parse_args <- function(args) {
    # Base-R parser for --key=value arguments (no optparse dependency)
    out <- list()
    for (a in args) {
        if (!grepl("^--[A-Za-z0-9_-]+=", a)) stop(sprintf("Arguments must be given as --key=value (got '%s')", a))
        out[[sub("^--([^=]+)=.*$", "\\1", a)]] <- sub("^--[^=]+=", "", a)
    }
    out
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
for (k in c("gtf", "out")) if (is.null(opt[[k]])) stop(sprintf("Missing required argument --%s", k))

library(bambu)
annot <- prepareAnnotations(opt$gtf)
saveRDS(annot, opt$out)
