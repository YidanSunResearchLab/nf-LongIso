#!/bin/sh
# Unit test for bin/te_analysis.sh on a synthetic locus (chrT); expected outputs were derived by hand.
#   R1  reference transcript (2 exons); a MIR overlaps the first exon            -> TSS
#   N1  4 exons; two novel intronic exons fully covered by AluY / L1M, both ends supported by
#       junctions -> two exonization rows for one transcript (v0.1.0 kept only one)
#   N2  single-exon transcript, 75 % covered by L2                               -> autonomous
#   N3  novel exon covered by MER5A, only spanned by an exon-skipping junction   -> not validated
# Run from any directory inside a container that provides bedtools (e.g. biocontainers bedtools):
#   sh tests/te_unit/run_te_unit_test.sh
set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
SCRIPT="$HERE/../../bin/te_analysis.sh"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
status=0
for mode in exon span; do
    mkdir -p "$WORK/$mode"
    (cd "$WORK/$mode" && sh "$SCRIPT" "$mode" 20 sense 0.5 "$HERE/iso.gtf" "$HERE/te.gtf" \
        "$HERE/counts_transcript.txt" "$HERE/ref.gtf" "$HERE/s1.junctions.bed" "$HERE/s2.junctions.bed" > log.txt)
    for f in TE_family_counts.tsv TE_instance_counts.tsv isoform_TE_annotation.tsv isoform_TE_exonization.tsv; do
        if cmp -s "$HERE/expected/$mode/$f" "$WORK/$mode/$f"; then
            echo "PASS  $mode/$f"
        else
            echo "FAIL  $mode/$f"; diff "$HERE/expected/$mode/$f" "$WORK/$mode/$f" || true; status=1
        fi
    done
done
# strand=both adds the antisense AluSx overlapping the last exon of R1, N1 and N3 (category PAS)
mkdir -p "$WORK/both"
(cd "$WORK/both" && sh "$SCRIPT" exon 20 both 0.5 "$HERE/iso.gtf" "$HERE/te.gtf" \
    "$HERE/counts_transcript.txt" "$HERE/ref.gtf" "$HERE/s1.junctions.bed" "$HERE/s2.junctions.bed" > log.txt)
n=$(awk -F '\t' '$3 == "AluSx_dup5" && $13 == "antisense" && $14 == "PAS" && $15 == "yes"' "$WORK/both/isoform_TE_annotation.tsv" | wc -l)
if [ "$n" -eq 3 ]; then echo "PASS  both/antisense PAS rows (3)"; else echo "FAIL  both/antisense PAS rows: $n (expected 3)"; status=1; fi
exit $status
