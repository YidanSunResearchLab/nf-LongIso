# NextLongIso: splice junctions from SAM records (samtools view -F 4 output), as BED6 lines
# (chr, intron start0, intron end, name, 1000, strand).
# Strand: XS:A tag if present, otherwise the read strand (FLAG bit 16).
# Bit tests use int(flag / 2^k) % 2 instead of gawk's and().  POSIX awk.
BEGIN { OFS = "\t" }
{
    chr = $3; pos = $4; flag = $2; cigar = $6
    strand = "."
    for (i = 12; i <= NF; i++) {
        if ($i ~ /^XS:A:/) { strand = substr($i, 6, 1); break }
    }
    if (strand == ".") strand = (int(flag / 16) % 2 == 1) ? "-" : "+"

    ref_pos = pos
    while (match(cigar, /[0-9]+[MIDNSHPX=]/)) {
        len = substr(cigar, RSTART, RLENGTH - 1) + 0
        op  = substr(cigar, RSTART + RLENGTH - 1, 1)
        cigar = substr(cigar, RSTART + RLENGTH)
        if (op == "M" || op == "=" || op == "X" || op == "D") {
            ref_pos += len
        } else if (op == "N") {
            junc_start = ref_pos - 1
            junc_end = ref_pos + len - 1
            print chr, junc_start, junc_end, "junc_" chr "_" junc_start "_" junc_end "_" strand, "1000", strand
            ref_pos += len
        }
    }
}
