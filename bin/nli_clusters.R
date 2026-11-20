# NextLongIso: clustering of transcript ends (sourced by apa_sites.R and tss_sites.R).
#
# cluster_ends(gene, strand, pos, window, max_width)
#   Clusters positions within each gene and strand. Positions are visited in genomic order; a
#   position joins the current cluster when it is <= `window` bp from the previous position AND
#   <= `max_width` bp from the first position of the cluster, otherwise it starts a new cluster.
#   The width limit stops single-linkage chaining from joining distant sites.
#   Returns an integer cluster index per input row, numbered in transcription order (5' -> 3')
#   within each gene and strand. One ordered pass over all rows (no per-gene loop).
cluster_ends <- function(gene, strand, pos, window, max_width) {
    key <- paste(gene, strand, sep = "\r")
    o <- order(key, pos)
    k <- key[o]; p <- pos[o]
    n <- length(p)
    lab <- integer(n)
    id <- 0L; first <- prev <- NA_real_
    for (j in seq_len(n)) {
        if (j == 1 || k[j] != k[j - 1]) {
            id <- 1L; first <- p[j]
        } else if (p[j] - prev > window || p[j] - first > max_width) {
            id <- id + 1L; first <- p[j]
        }
        lab[j] <- id; prev <- p[j]
    }
    # transcription order: on the minus strand the right-most cluster is the first one
    mx <- ave(lab, k, FUN = max)
    lab <- ifelse(strand[o] == "-", mx + 1L - lab, lab)
    out <- integer(n); out[o] <- lab
    out
}

# Representative coordinate per cluster: the position with the most reads; ties are broken by the
# number of transcripts ending there, then by the 5'-most position (in transcription direction).
cluster_coord <- function(cluster_key, pos, weight, strand) {
    agg <- aggregate(data.frame(w = weight, n = 1), by = list(key = cluster_key, pos = pos, strand = strand), FUN = sum)
    tx_pos <- ifelse(agg$strand == "-", -agg$pos, agg$pos)
    agg <- agg[order(agg$key, -agg$w, -agg$n, tx_pos), ]
    best <- agg[!duplicated(agg$key), ]
    best$pos[match(cluster_key, best$key)]
}
