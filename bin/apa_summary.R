#!/usr/bin/env Rscript
# NextLongIso: gene-level APA summary from the DEXSeq PAS usage results.
#
# Usage:
#   apa_summary.R --usage=apa_usage.tsv --annotation=pas_annotation.tsv --alpha=0.05
#
# Output APA_gene_summary.csv, one row per (comparison, gene):
#   top_pas_id / top_pas_coord / top_dPAU / top_padj: the PAS with the largest |dPAU|
#   gene_qvalue: DEXSeq perGeneQValue; significant: gene_qvalue < alpha
#   pas_shift_bp: change of the PAU-weighted mean PAS position in transcription direction
#     (condition_2 - condition_1); > 0 = shift to distal sites (longer 3' ends / 3' UTRs),
#     < 0 = shift to proximal sites (shorter)
#   direction: distal / proximal / none

.script <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "nli_design.R"))
opt <- parse_args(commandArgs(trailingOnly = TRUE))
for (k in c("usage", "annotation", "alpha")) if (is.null(opt[[k]])) stop(sprintf("Missing required argument --%s", k))
alpha <- as.numeric(opt$alpha)

u <- read.delim(opt$usage, check.names = FALSE, stringsAsFactors = FALSE)
if ("Message" %in% colnames(u)) {
    write.csv(data.frame(Message = u$Message[1]), "APA_gene_summary.csv", row.names = FALSE)
    quit(save = "no", status = 0)
}
a <- read.delim(opt$annotation, stringsAsFactors = FALSE)
idx <- match(paste(u$gene_id, u$feature_id), paste(a$gene_id, a$pas_id))
u$pas_coord <- a$pas_coord[idx]; u$strand <- a$strand[idx]; u$chr <- a$chr[idx]

rows <- lapply(split(u, paste(u$condition_2, u$gene_id)), function(x) {
    x <- x[order(-abs(x$dPAU), x$padj, x$pas_coord), ]
    sgn <- ifelse(x$strand[1] == "+", 1, -1)   # position along the transcription direction
    w1 <- x$mean_PAU_1; w2 <- x$mean_PAU_2
    shift <- if (sum(w1, na.rm = TRUE) > 0 && sum(w2, na.rm = TRUE) > 0)
        sgn * (sum(w2 * x$pas_coord, na.rm = TRUE) / sum(w2, na.rm = TRUE) -
               sum(w1 * x$pas_coord, na.rm = TRUE) / sum(w1, na.rm = TRUE)) else NA_real_
    data.frame(condition_1 = x$condition_1[1], condition_2 = x$condition_2[1], gene_id = x$gene_id[1],
               chr = x$chr[1], strand = x$strand[1], n_sites = nrow(x), gene_qvalue = x$gene_qvalue[1],
               significant = !is.na(x$gene_qvalue[1]) && x$gene_qvalue[1] < alpha,
               top_pas_id = x$feature_id[1], top_pas_coord = x$pas_coord[1], top_dPAU = x$dPAU[1],
               top_padj = x$padj[1], pas_shift_bp = round(shift, 1),
               direction = if (is.na(shift) || abs(shift) < 1) "none" else if (shift > 0) "distal" else "proximal",
               stringsAsFactors = FALSE)
})
s <- do.call(rbind, rows)
s <- s[order(s$condition_2, s$gene_qvalue, -abs(s$top_dPAU), na.last = TRUE), ]
write.csv(s, "APA_gene_summary.csv", row.names = FALSE)
message(sprintf("APA genes tested: %d; significant (q < %g): %d (distal %d, proximal %d)",
                nrow(s), alpha, sum(s$significant), sum(s$significant & s$direction == "distal"),
                sum(s$significant & s$direction == "proximal")))
