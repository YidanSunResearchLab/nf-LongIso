# NextLongIso: introns of a GTF annotation. Input: exon lines reduced to
#   transcript_id, chr, start (1-based), end, strand   and sorted by transcript, then start.
# Output: chr, intron start0, intron end, strand (BED interval of the intron).  POSIX awk.
BEGIN { FS = OFS = "\t" }
{
    if ($1 == tx && $3 > prev_end) print $2, prev_end, $3 - 1, $5
    tx = $1; prev_end = $4
}
