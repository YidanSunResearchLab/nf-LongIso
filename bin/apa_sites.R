#!/usr/bin/env Rscript
# NextLongIso: poly(A) sites (PAS) from transcript 3' ends, internal-priming check and PAS counts.
#
# Usage:
#   apa_sites.R --exons=transcript_exons.tsv --counts=counts_transcript.txt --genome=genome.fa \
#               --design=design.csv --window=25 --max-width=50 --ip-a-frac=0.7 --ip-mode=filter
#
# Inputs
#   --exons    TSV without header: transcript_id, gene_id, chr, start, end, strand (exons of the
#              Bambu extended annotation, 1-based)
#   --counts   Bambu counts_transcript.txt (TXNAME, GENEID, samples...)
#   --genome   genome FASTA with a .fai index next to it
#
# Method
#   * PAS of a transcript = its 3' end (strand-aware). Transcripts whose ID does not start with
#     'BambuTx' come from the reference annotation ('annotated' ends).
#   * Internal priming (novel transcripts only): the 20 nt downstream of the 3' end contain
#     >= --ip-a-frac A or a run of >= 6 A, unless an annotated 3' end lies within 50 bp (same strand)
#     or a poly(A) signal hexamer occurs 10-40 nt upstream. --ip-mode filter removes such ends before
#     clustering; --ip-mode flag keeps them and marks the cluster.
#   * Ends are clustered per gene and strand (window --window, maximal width --max-width).
#   * PAS counts = sum of the Bambu read counts of the transcripts in each cluster, per sample.
#
# Outputs
#   pas_transcripts.tsv  one row per transcript 3' end (cluster, annotation status, priming check)
#   pas_annotation.tsv   one row per PAS cluster
#   pas_counts.tsv       feature_id, gene_id, one column per sample (input of usage_dexseq.R)

suppressPackageStartupMessages({
    library(Rsamtools)
    library(Biostrings)
    library(GenomicRanges)
})
.script <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "nli_design.R"))
source(file.path(dirname(normalizePath(.script)), "nli_clusters.R"))
opt <- parse_args(commandArgs(trailingOnly = TRUE))
for (k in c("exons", "counts", "genome", "design", "window", "max-width", "ip-a-frac", "ip-mode")) {
    if (is.null(opt[[k]])) stop(sprintf("Missing required argument --%s", k))
}
window    <- as.numeric(opt$window)
max_width <- as.numeric(opt[["max-width"]])
a_frac    <- as.numeric(opt[["ip-a-frac"]])
ip_mode   <- opt[["ip-mode"]]
if (!ip_mode %in% c("filter", "flag")) stop("--ip-mode must be 'filter' or 'flag'")
PAS_SIGNALS <- c("AATAAA", "ATTAAA", "AGTAAA", "TATAAA", "CATAAA", "GATAAA",
                 "AATATA", "AATACA", "AATAGA", "AAAAAG", "ACTAAA", "AAGAAA")
strip_suffix <- function(x) sub("[.]sorted([.]bam)?$", "", x)

# ---- 1. Transcript 3' ends ---------------------------------------------------------------------
ex <- read.delim(opt$exons, header = FALSE, stringsAsFactors = FALSE,
                 col.names = c("transcript_id", "gene_id", "chr", "start", "end", "strand"))
ex <- ex[ex$strand %in% c("+", "-"), ]
ends <- ex[!duplicated(ex$transcript_id), c("transcript_id", "gene_id", "chr", "strand")]
max_end <- tapply(ex$end, ex$transcript_id, max)
min_start <- tapply(ex$start, ex$transcript_id, min)
ends$pos <- ifelse(ends$strand == "+", max_end[ends$transcript_id], min_start[ends$transcript_id])
ends$annotated <- !startsWith(ends$transcript_id, "BambuTx")

counts <- read.delim(opt$counts, check.names = FALSE, stringsAsFactors = FALSE)
colnames(counts) <- strip_suffix(colnames(counts))
design <- read.csv(opt$design, stringsAsFactors = FALSE)
cmat <- as.matrix(counts[, design$sample, drop = FALSE]); rownames(cmat) <- counts$TXNAME
ends$reads <- rowSums(cmat)[ends$transcript_id]; ends$reads[is.na(ends$reads)] <- 0

# ---- 2. Internal-priming check (novel 3' ends) ----------------------------------------------------
fa <- FaFile(opt$genome)
chr_len <- seqlengths(scanFaIndex(fa)); chr_len <- setNames(as.numeric(chr_len), names(chr_len))
get_seq <- function(chr, from, to, strand) {
    from <- pmax(1, from); to <- pmin(chr_len[chr], to)
    out <- rep("", length(chr)); ok <- !is.na(to) & to >= from
    if (any(ok)) {
        s <- scanFa(fa, GRanges(chr[ok], IRanges(from[ok], to[ok])))
        s[strand[ok] == "-"] <- reverseComplement(s[strand[ok] == "-"])
        out[ok] <- toupper(as.character(s))
    }
    out
}
plus <- ends$strand == "+"
down <- get_seq(ends$chr, ifelse(plus, ends$pos + 1, ends$pos - 20), ifelse(plus, ends$pos + 20, ends$pos - 1), ends$strand)
up   <- get_seq(ends$chr, ifelse(plus, ends$pos - 40, ends$pos + 10), ifelse(plus, ends$pos - 10, ends$pos + 40), ends$strand)
ends$downstream_A_frac <- ifelse(nchar(down) > 0, nchar(gsub("[^A]", "", down)) / nchar(down), NA)
ends$downstream_A6     <- grepl("AAAAAA", down)
ends$pas_signal <- vapply(up, function(s) { h <- PAS_SIGNALS[vapply(PAS_SIGNALS, grepl, logical(1), x = s, fixed = TRUE)]; if (length(h)) h[1] else NA_character_ }, character(1))

ann <- ends[ends$annotated, c("chr", "strand", "pos")]
near_annotated <- rep(FALSE, nrow(ends))
if (nrow(ann) > 0) {
    q <- GRanges(ends$chr, IRanges(ends$pos, width = 1), strand = ends$strand)
    s <- GRanges(ann$chr, IRanges(ann$pos, width = 1), strand = ann$strand)
    near_annotated <- overlapsAny(q, s, maxgap = 50)
}
a_rich <- (!is.na(ends$downstream_A_frac) & ends$downstream_A_frac >= a_frac) | ends$downstream_A6
ends$internal_priming <- !ends$annotated & a_rich & !near_annotated & is.na(ends$pas_signal)
message(sprintf("3' ends: %d (%d novel); novel ends with A-rich downstream sequence: %d; flagged as internal priming: %d",
                nrow(ends), sum(!ends$annotated), sum(!ends$annotated & a_rich), sum(ends$internal_priming)))

# ---- 3. Clusters ------------------------------------------------------------------------------------
use <- if (ip_mode == "filter") !ends$internal_priming else rep(TRUE, nrow(ends))
ends$cluster <- NA_integer_
ends$cluster[use] <- cluster_ends(ends$gene_id[use], ends$strand[use], ends$pos[use], window, max_width)
ends$pas_id <- ifelse(is.na(ends$cluster), NA, paste0("PAS", ends$cluster))
ends$used <- use

u <- ends[use, ]
key <- paste(u$gene_id, u$strand, u$pas_id, sep = "\r")
u$pas_coord <- cluster_coord(key, u$pos, u$reads, u$strand)
first <- !duplicated(key)
rep_idx <- which(u$pos == u$pas_coord)
rep_idx <- rep_idx[!duplicated(key[rep_idx])]
ann_tab <- data.frame(
    gene_id = u$gene_id[first], pas_id = u$pas_id[first], chr = u$chr[first], strand = u$strand[first],
    start = as.numeric(tapply(u$pos, key, min)[key[first]]), end = as.numeric(tapply(u$pos, key, max)[key[first]]),
    pas_coord = u$pas_coord[first], n_transcripts = as.integer(table(key)[key[first]]),
    transcripts = tapply(u$transcript_id, key, function(t) paste(sort(t), collapse = ","))[key[first]],
    annotated = ifelse(tapply(u$annotated, key, any)[key[first]], "yes", "no"),
    internal_priming = ifelse(tapply(u$internal_priming, key, any)[key[first]], "yes", "no"),
    pas_signal = u$pas_signal[rep_idx][match(key[first], key[rep_idx])],
    reads = as.numeric(tapply(u$reads, key, sum)[key[first]]), stringsAsFactors = FALSE)
ann_tab <- ann_tab[order(ann_tab$chr, ann_tab$start, ann_tab$gene_id, ann_tab$pas_id), ]
write.table(ann_tab, "pas_annotation.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

out_ends <- ends[order(ends$gene_id, ends$transcript_id), c("transcript_id", "gene_id", "chr", "strand", "pos",
    "annotated", "reads", "downstream_A_frac", "downstream_A6", "pas_signal", "internal_priming", "used", "pas_id")]
write.table(out_ends, "pas_transcripts.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# ---- 4. PAS x sample counts -------------------------------------------------------------------------
pm <- rowsum(cmat[u$transcript_id, , drop = FALSE], key)
idx <- match(rownames(pm), paste(ann_tab$gene_id, ann_tab$strand, ann_tab$pas_id, sep = "\r"))
pas_counts <- data.frame(feature_id = ann_tab$pas_id[idx], gene_id = ann_tab$gene_id[idx], pm, check.names = FALSE,
                         stringsAsFactors = FALSE)
pas_counts <- pas_counts[order(pas_counts$gene_id, pas_counts$feature_id), ]
write.table(pas_counts, "pas_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("PAS clusters: %d in %d genes (%d genes with >= 2 clusters); internal-priming mode: %s",
                nrow(ann_tab), length(unique(ann_tab$gene_id)), sum(table(ann_tab$gene_id) >= 2), ip_mode))
