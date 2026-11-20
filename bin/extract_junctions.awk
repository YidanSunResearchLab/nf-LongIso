# NextLongIso: count splice junctions in SAM records (samtools view output of primary, mapped reads).
# Output (unsorted, one line per junction and strand): chr, intron start0, intron end, strand, reads
#
# Transcript strand of a read:
#   ts:A tag (minimap2): '+' = transcript on the read's alignment strand, '-' = opposite strand
#   otherwise XS:A (other aligners), otherwise '.' (unknown); the read orientation alone is not used.
# Bit tests use int(flag / 2^k) % 2 (POSIX awk, no gawk and()).
BEGIN { OFS = "\t" }
{
    flag = $2; cigar = $6
    ts = ""; xs = ""
    for (i = 12; i <= NF; i++) {
        if ($i ~ /^ts:A:/) ts = substr($i, 6, 1)
        else if ($i ~ /^XS:A:/) xs = substr($i, 6, 1)
    }
    read_strand = (int(flag / 16) % 2 == 1) ? "-" : "+"
    if (ts == "+") strand = read_strand
    else if (ts == "-") strand = (read_strand == "+") ? "-" : "+"
    else if (xs == "+" || xs == "-") strand = xs
    else strand = "."

    ref_pos = $4
    while (match(cigar, /[0-9]+[MIDNSHPX=]/)) {
        len = substr(cigar, RSTART, RLENGTH - 1) + 0
        op  = substr(cigar, RSTART + RLENGTH - 1, 1)
        cigar = substr(cigar, RSTART + RLENGTH)
        if (op == "M" || op == "=" || op == "X" || op == "D") {
            ref_pos += len
        } else if (op == "N") {
            n[$3 OFS (ref_pos - 1) OFS (ref_pos + len - 1) OFS strand]++
            ref_pos += len
        }
    }
}
END { for (k in n) print k, n[k] }
