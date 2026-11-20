#!/usr/bin/env Rscript
# NextLongIso: promoters (TSS clusters) from transcript 5' ends, promoter counts and activity.
#
# Usage:
#   tss_sites.R --exons=transcript_exons.tsv --counts=counts_transcript.txt --gene-counts=counts_gene.txt \
#               --design=design.csv --design-formula="~ condition" --reference=control \
#               --window=50 --max-width=100 --min-activity=0.25
#
# Inputs
#   --exons    TSV without header: transcript_id, gene_id, chr, start, end, strand (Bambu extended
#              annotation exons, 1-based)
#   --counts   Bambu counts_transcript.txt; --gene-counts Bambu counts_gene.txt (size factors)
#
# Method
#   * TSS of a transcript = its 5' end (strand-aware); single-exon transcripts are included.
#   * TSSs are clustered per gene and strand (window --window, maximal width --max-width).
#   * Promoter counts = sum of the Bambu read counts of the transcripts starting in the cluster.
#   * Absolute activity = log2(normalised count + 1), normalised with DESeq2 size factors estimated
#     from the gene-level counts (as proActiv, but with gene-level size factors).
#   * Relative activity = promoter count / gene count (sum of the gene's promoter counts).
#   * Class per condition (proActiv rules): mean absolute activity < --min-activity -> Inactive; the
#     most active promoter of the gene -> Major; other active promoters -> Minor.
#
# Outputs
#   promoter_annotation.tsv  promoter_id, gene_id, chr, strand, start, end, tss_coord, transcripts,
#                            n_transcripts, n_single_exon, annotated, internal_promoter,
#                            dist_to_annotated_tss
#   promoter_counts.tsv      feature_id, gene_id, one column per sample (input of usage_dexseq.R)
#   promoter_activity.tsv    annotation columns + absolute/relative activity per sample +
#                            mean activity and class per condition

suppressPackageStartupMessages(library(DESeq2))
.script <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "nli_design.R"))
source(file.path(dirname(normalizePath(.script)), "nli_clusters.R"))
opt <- parse_args(commandArgs(trailingOnly = TRUE))
for (k in c("exons", "counts", "gene-counts", "design", "design-formula", "reference", "window", "max-width", "min-activity")) {
    if (is.null(opt[[k]])) stop(sprintf("Missing required argument --%s", k))
}
window <- as.numeric(opt$window); max_width <- as.numeric(opt[["max-width"]]); min_act <- as.numeric(opt[["min-activity"]])
strip_suffix <- function(x) sub("[.]sorted([.]bam)?$", "", x)

# ---- 1. Transcript 5' ends -----------------------------------------------------------------------
ex <- read.delim(opt$exons, header = FALSE, stringsAsFactors = FALSE,
                 col.names = c("transcript_id", "gene_id", "chr", "start", "end", "strand"))
ex <- ex[ex$strand %in% c("+", "-"), ]
tx <- ex[!duplicated(ex$transcript_id), c("transcript_id", "gene_id", "chr", "strand")]
min_start <- tapply(ex$start, ex$transcript_id, min); max_end <- tapply(ex$end, ex$transcript_id, max)
tx$tss <- ifelse(tx$strand == "+", min_start[tx$transcript_id], max_end[tx$transcript_id])
tx$tes <- ifelse(tx$strand == "+", max_end[tx$transcript_id], min_start[tx$transcript_id])
tx$n_exons <- as.integer(table(ex$transcript_id)[tx$transcript_id])
tx$annotated <- !startsWith(tx$transcript_id, "BambuTx")

design <- read.csv(opt$design, stringsAsFactors = FALSE)
counts <- read.delim(opt$counts, check.names = FALSE, stringsAsFactors = FALSE)
colnames(counts) <- strip_suffix(colnames(counts))
cmat <- as.matrix(counts[, design$sample, drop = FALSE]); rownames(cmat) <- counts$TXNAME
cmat[is.na(cmat)] <- 0
tx <- tx[tx$transcript_id %in% rownames(cmat), ]
tx$reads <- rowSums(cmat)[tx$transcript_id]

# ---- 2. TSS clusters and annotation ----------------------------------------------------------------
tx$cluster <- cluster_ends(tx$gene_id, tx$strand, tx$tss, window, max_width)
tx$promoter_id <- paste0("P", tx$cluster)
key <- paste(tx$gene_id, tx$strand, tx$promoter_id, sep = "\r")
tx$tss_coord <- cluster_coord(key, tx$tss, tx$reads, tx$strand)
first <- !duplicated(key)
ann <- data.frame(
    promoter_id = tx$promoter_id[first], gene_id = tx$gene_id[first], chr = tx$chr[first], strand = tx$strand[first],
    start = as.numeric(tapply(tx$tss, key, min)[key[first]]), end = as.numeric(tapply(tx$tss, key, max)[key[first]]),
    tss_coord = tx$tss_coord[first],
    transcripts = tapply(tx$transcript_id, key, function(t) paste(sort(t), collapse = ","))[key[first]],
    n_transcripts = as.integer(table(key)[key[first]]),
    n_single_exon = as.integer(tapply(tx$n_exons == 1, key, sum)[key[first]]),
    annotated = ifelse(tapply(tx$annotated, key, any)[key[first]], "yes", "no"),
    stringsAsFactors = FALSE)
ann$key <- key[first]

# Internal promoter: the TSS lies inside the body of another transcript of the same gene that
# starts further upstream (in transcription direction)
s <- ifelse(tx$strand == "+", 1, -1)
body_from <- s * tx$tss; body_to <- s * tx$tes        # transcription-oriented coordinates
ann$internal_promoter <- "no"
by_gene <- split(seq_len(nrow(tx)), tx$gene_id)
gene_idx <- match(ann$gene_id, names(by_gene))      # positional lookup (named-list lookup is linear)
for (i in seq_len(nrow(ann))) {
    j <- by_gene[[gene_idx[i]]]
    j <- j[tx$strand[j] == ann$strand[i] & tx$promoter_id[j] != ann$promoter_id[i]]
    if (length(j) == 0) next
    p <- ifelse(ann$strand[i] == "+", 1, -1) * ann$tss_coord[i]
    if (any(body_from[j] < p & body_to[j] > p)) ann$internal_promoter[i] <- "yes"
}

# Distance to the nearest annotated TSS (same chromosome and strand)
a_tss <- tx[tx$annotated, c("chr", "strand", "tss")]
ann$dist_to_annotated_tss <- NA_real_
if (nrow(a_tss) > 0) {
    grp <- split(a_tss$tss, paste(a_tss$chr, a_tss$strand))
    k2 <- paste(ann$chr, ann$strand)
    for (g in intersect(unique(k2), names(grp))) {
        ref <- sort(grp[[g]]); i <- which(k2 == g)
        pos <- findInterval(ann$tss_coord[i], ref)
        lo <- ref[pmax(pos, 1)]; hi <- ref[pmin(pos + 1, length(ref))]
        ann$dist_to_annotated_tss[i] <- pmin(abs(ann$tss_coord[i] - lo), abs(ann$tss_coord[i] - hi))
    }
}
ann <- ann[order(ann$chr, ann$start, ann$gene_id, ann$promoter_id), ]

# ---- 3. Counts ---------------------------------------------------------------------------------------
pm <- rowsum(cmat[tx$transcript_id, , drop = FALSE], key)[ann$key, , drop = FALSE]
pcounts <- data.frame(feature_id = ann$promoter_id, gene_id = ann$gene_id, pm, check.names = FALSE, stringsAsFactors = FALSE)
write.table(pcounts, "promoter_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 4. Activity -------------------------------------------------------------------------------------
gene <- read.delim(opt[["gene-counts"]], check.names = FALSE, stringsAsFactors = FALSE)
colnames(gene) <- strip_suffix(colnames(gene))
gmat <- round(as.matrix(gene[, design$sample, drop = FALSE])); gmat[is.na(gmat)] <- 0
sf <- sizeFactors(estimateSizeFactors(DESeqDataSetFromMatrix(gmat, data.frame(row.names = design$sample, x = rep(1, nrow(design))), ~ 1)))
abs_act <- log2(sweep(pm, 2, sf[design$sample], "/") + 1)
gene_tot <- rowsum(pm, ann$gene_id)[ann$gene_id, , drop = FALSE]
rel_act <- pm / gene_tot; rel_act[!is.finite(rel_act)] <- NA
act <- ann[, setdiff(colnames(ann), "key")]
for (smp in design$sample) {
    act[[paste0("absolute_", smp)]] <- abs_act[, smp]
    act[[paste0("relative_", smp)]] <- rel_act[, smp]
}
for (cond in unique(design$condition)) {
    m <- rowMeans(abs_act[, design$sample[design$condition == cond], drop = FALSE])
    cls <- ifelse(m < min_act, "Inactive", "Minor")
    top <- tapply(seq_along(m), ann$gene_id, function(i) i[which.max(m[i])])
    major <- unlist(top); major <- major[m[major] >= min_act]
    cls[major] <- "Major"
    act[[paste0(cond, ".mean")]] <- m
    act[[paste0(cond, ".class")]] <- cls
}
write.table(act, "promoter_activity.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(act[, c("promoter_id", "gene_id", "chr", "strand", "start", "end", "tss_coord", "transcripts", "n_transcripts",
                    "n_single_exon", "annotated", "internal_promoter", "dist_to_annotated_tss")],
            "promoter_annotation.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
message(sprintf("Promoters: %d in %d genes (%d genes with >= 2 promoters; %d promoters only from single-exon transcripts; %d internal)",
                nrow(ann), length(unique(ann$gene_id)), sum(table(ann$gene_id) >= 2),
                sum(ann$n_single_exon == ann$n_transcripts), sum(ann$internal_promoter == "yes")))
