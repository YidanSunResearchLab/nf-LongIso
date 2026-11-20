#!/bin/sh
# Unit test for bin/te_analysis.sh on a synthetic locus (chrT); expected outputs were derived by hand
# (exon mode) or produced by the v0.1.0 code (span mode).
#   R1  reference transcript (2 exons); a MIR overlaps the first exon            -> TSS
#   N1  4 exons; two novel intronic exons fully covered by AluY / L1M, both ends supported by
#       junctions -> two exonization rows in exon mode (v0.1.0/span: none, the ID occurs 4 times)
#   N2  single-exon transcript, 75 % covered by L2                               -> autonomous
#   N3  novel exon covered by MER5A, only spanned by an exon-skipping junction   -> not validated in
#       exon mode; validated by the v0.1.0 window test (span mode)
# Files are compared on the header line plus the sorted remaining rows (span mode writes counts in
# awk array order, as v0.1.0 did, which differs between awk implementations).
# Run inside a container that provides bedtools, e.g. quay.io/biocontainers/bedtools:2.31.1--h13024bc_3:
#   sh tests/te_unit/run_te_unit_test.sh
set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
SCRIPT="$HERE/../../bin/te_analysis.sh"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
status=0

normalise() { head -n 1 "$1"; tail -n +2 "$1" | LC_ALL=C sort; }

run() {  # run <dir> <mode> <strand> [junction files...]
    d=$1; m=$2; st=$3; shift 3
    mkdir -p "$WORK/$d"
    (cd "$WORK/$d" && sh "$SCRIPT" "$m" 20 "$st" 0.5 "$HERE/iso.gtf" "$HERE/te.gtf" \
        "$HERE/counts_transcript.txt" "$HERE/ref.gtf" "$@" > log.txt 2>&1) || { echo "FAIL  $d: script error"; cat "$WORK/$d/log.txt"; status=1; }
}

for mode in exon span; do
    run "$mode" "$mode" sense "$HERE/s1.junctions.bed" "$HERE/s2.junctions.bed"
    for f in TE_family_counts.tsv TE_instance_counts.tsv isoform_TE_annotation.tsv isoform_TE_exonization.tsv; do
        if [ -f "$WORK/$mode/$f" ] && [ "$(normalise "$HERE/expected/$mode/$f")" = "$(normalise "$WORK/$mode/$f")" ]; then
            echo "PASS  $mode/$f"
        else
            echo "FAIL  $mode/$f"; diff "$HERE/expected/$mode/$f" "$WORK/$mode/$f" || true; status=1
        fi
    done
done

# strand=both adds the antisense AluSx overlapping the last exon of R1, N1 and N3 (category PAS)
run both exon both "$HERE/s1.junctions.bed" "$HERE/s2.junctions.bed"
n=$(awk -F '\t' '$3 == "AluSx_dup5" && $13 == "antisense" && $14 == "PAS" && $15 == "last" && $16 == "yes"' "$WORK/both/isoform_TE_annotation.tsv" | wc -l)
if [ "$n" -eq 3 ]; then echo "PASS  both/antisense PAS rows (3)"; else echo "FAIL  both/antisense PAS rows: $n (expected 3)"; status=1; fi

# No junction files (empty junction_ends.txt): same annotation pairs and categories, no exonization,
# splice sites unsupported ('no') except single-exon N2 ('NA')
run nojunc exon sense
if [ "$(cut -f1-15 "$WORK/nojunc/isoform_TE_annotation.tsv" | normalise /dev/stdin)" = "$(cut -f1-15 "$HERE/expected/exon/isoform_TE_annotation.tsv" | normalise /dev/stdin)" ] \
   && [ "$(tail -n +2 "$WORK/nojunc/isoform_TE_exonization.tsv" | wc -l)" -eq 0 ] \
   && [ "$(awk -F '\t' 'NR > 1 { n[$16]++ } END { printf "NA=%d no=%d yes=%d", n["NA"], n["no"], n["yes"] }' "$WORK/nojunc/isoform_TE_annotation.tsv")" = "NA=1 no=6 yes=0" ]; then
    echo "PASS  nojunc (empty junction input)"
else
    echo "FAIL  nojunc (empty junction input)"; cat "$WORK/nojunc/isoform_TE_annotation.tsv" "$WORK/nojunc/isoform_TE_exonization.tsv"; status=1
fi
exit $status
