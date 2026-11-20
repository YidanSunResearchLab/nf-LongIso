# NextLongIso: introns of a GTF annotation. Input: exon lines reduced to
#   transcript_id, chr, start (1-based), end, strand   and sorted by transcript, chr, then start.
# Output: chr, intron start0, intron end, strand (BED interval of the intron).  POSIX awk.
# Introns are only formed between exons of the same transcript on the same chromosome (a transcript
# ID reused on another chromosome, e.g. PAR genes in some annotations, is kept apart), and abutting
# exons (no intronic base between them) give no intron.
BEGIN { FS = OFS = "\t" }
{
    if ($1 == tx && $2 == chr && $3 > prev_end + 1) print $2, prev_end, $3 - 1, $5
    tx = $1; chr = $2; prev_end = $4
}
