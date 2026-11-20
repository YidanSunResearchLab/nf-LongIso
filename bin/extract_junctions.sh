#!/bin/sh
# NextLongIso: splice junctions of one BAM file.
#
# Usage: extract_junctions.sh <bam> <genome.fa (indexed)> <gtf_introns.tsv> <min_mapq> <prefix>
#
# Outputs
#   <prefix>.junctions.bed  BED6: chr, start, end, name, score = min(reads, 1000), strand
#   <prefix>.junctions.tsv  chr, start, end, strand, reads, motif, annotated
#     motif: GT-AG, GC-AG, AT-AC or non-canonical, in transcript orientation (for strand '.', the
#            orientation that gives a canonical motif, if any); read from the genome sequence;
#            NA when the contig is not in the genome index or the intron extends past its end
#     annotated: yes if the intron (chr, start, end) is an intron of the reference annotation
# Junction coordinates are BED intervals of the intron (0-based start, end exclusive).
# Requires samtools, a POSIX awk and <genome.fa>.fai. Runs under sh with pipefail when available.
set -eu
# A failing samtools or awk inside a pipe must fail the task (not give truncated junction tables)
if (set -o pipefail) 2>/dev/null; then set -o pipefail; fi
BAM=$1; GENOME=$2; INTRONS=$3; MAPQ=$4; PREFIX=$5
AWK_DIR=$(dirname "$(command -v extract_junctions.awk || echo "$0")")

# Primary, mapped reads (-F 2308: unmapped, secondary, supplementary) with MAPQ >= MAPQ
samtools view -F 2308 -q "$MAPQ" "$BAM" | awk -f "$AWK_DIR/extract_junctions.awk" \
    | sort -k1,1 -k2,2n -k3,3n -k4,4 > junctions.raw.tsv

# Splice-site dinucleotides: donor = first two intron bases, acceptor = last two (genome strand).
# Only junctions on contigs of the genome index and within the contig are looked up; each sequence is
# kept with its region name, so a missing region cannot shift the motifs of other junctions.
awk -F '\t' 'fileno == 1 { len[$1] = $2; next }
     ($1 in len) && $3 <= len[$1] { printf "%s:%d-%d\n%s:%d-%d\n", $1, $2 + 1, $2 + 2, $1, $3 - 1, $3 }' \
    fileno=1 "$GENOME.fai" fileno=2 junctions.raw.tsv > sites.txt
if [ -s sites.txt ]; then
    samtools faidx "$GENOME" -r sites.txt \
        | awk '/^>/ { if (name != "") print name "\t" seq; name = substr($0, 2); seq = ""; next }
               { seq = seq toupper($0) }
               END { if (name != "") print name "\t" seq }' > dinucleotides.txt
else
    : > dinucleotides.txt
fi

awk -F '\t' -v OFS='\t' -v prefix="$PREFIX" '
    function rc(s,    i, c, o) {
        o = ""
        for (i = length(s); i >= 1; i--) {
            c = substr(s, i, 1)
            o = o (c == "A" ? "T" : c == "T" ? "A" : c == "G" ? "C" : c == "C" ? "G" : "N")
        }
        return o
    }
    function canonical(m) { return m == "GT-AG" || m == "GC-AG" || m == "AT-AC" }
    fileno == 1 { dn[$1] = $2; next }
    fileno == 2 { annot[$1 SUBSEP $2 SUBSEP $3] = 1; next }
    fileno == 3 {
        donor    = $1 ":" ($2 + 1) "-" ($2 + 2)
        acceptor = $1 ":" ($3 - 1) "-" $3
        if (!(donor in dn) || !(acceptor in dn) || length(dn[donor]) != 2 || length(dn[acceptor]) != 2) {
            motif = "NA"; n_na++
        } else {
            fwd = dn[donor] "-" dn[acceptor]                 # motif read on the + strand
            rev = rc(dn[acceptor]) "-" rc(dn[donor])         # motif read on the - strand
            if ($4 == "+")      motif = canonical(fwd) ? fwd : "non-canonical"
            else if ($4 == "-") motif = canonical(rev) ? rev : "non-canonical"
            else                motif = canonical(fwd) ? fwd : (canonical(rev) ? rev : "non-canonical")
        }
        ann = (($1 SUBSEP $2 SUBSEP $3) in annot) ? "yes" : "no"
        print $1, $2, $3, $4, $5, motif, ann > (prefix ".junctions.body.tsv")
        print $1, $2, $3, "junc_" $1 "_" $2 "_" $3 "_" $4, ($5 < 1000 ? $5 : 1000), $4 > (prefix ".junctions.bed")
    }
    END {
        if (n_na > 0) printf "WARNING: %s: %d junction(s) on contigs missing from the genome index or past a contig end; motif NA\n", prefix, n_na > "/dev/stderr"
    }
' fileno=1 dinucleotides.txt fileno=2 "$INTRONS" fileno=3 junctions.raw.tsv

touch "$PREFIX.junctions.bed" "$PREFIX.junctions.body.tsv"
{ printf 'chr\tstart\tend\tstrand\treads\tmotif\tannotated\n'; cat "$PREFIX.junctions.body.tsv"; } > "$PREFIX.junctions.tsv"
rm -f "$PREFIX.junctions.body.tsv" junctions.raw.tsv sites.txt dinucleotides.txt

echo "$PREFIX: $(wc -l < "$PREFIX.junctions.bed") junctions; by strand: $(awk -F '\t' '{ n[$6]++ } END { for (s in n) printf "%s=%d ", s, n[s] }' "$PREFIX.junctions.bed")"
