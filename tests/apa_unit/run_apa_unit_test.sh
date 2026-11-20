#!/bin/sh
# Unit test for bin/apa_sites.R (internal-priming check and PAS clustering) on a synthetic 6-kb
# contig without any A in the background sequence (expected outputs checked by hand). Annotated
# transcripts are those listed in ref_ids.txt (ENST1, ENST2).
#   BambuTx1  +, novel 3' end at 1000, 20 A downstream, no PAS signal, no annotated end nearby -> flagged
#   BambuTx2  +, novel 3' end at 2000, 20 A downstream, AATAAA 25 nt upstream                  -> rescued
#   BambuTx3  +, novel 3' end at 3000, 20 A downstream, annotated end (ENST1) 30 bp away       -> rescued
#   ENST2     +, annotated 3' end at 4000 with 20 A downstream                                -> never flagged
#   BambuTx4  +, novel 3' end at 4500, 20 A downstream, only the weak variant AAGAAA upstream  -> flagged
#   BambuTx5  -, novel 3' end at 5000, 20 A downstream (genomic T), no signal                 -> flagged
#   BambuTx6  -, novel 3' end at 5500, 20 A downstream, AATAAA 15-20 nt upstream              -> rescued
# 'filter' removes BambuTx1 before clustering, 'flag' keeps it and marks its cluster.
# Run inside the IsoformSwitchAnalyzeR container (R with Rsamtools/Biostrings):
#   sh tests/apa_unit/run_apa_unit_test.sh
set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
status=0
for mode in filter flag; do
    mkdir -p "$WORK/$mode"
    (cd "$WORK/$mode" && Rscript "$HERE/../../bin/apa_sites.R" --exons="$HERE/ex.tsv" --counts="$HERE/counts.tsv" \
        --genome="$HERE/g.fa" --ref-ids="$HERE/ref_ids.txt" --design="$HERE/design.csv" --window=25 --max-width=50 --ip-a-frac=0.7 \
        --ip-mode="$mode" > log.txt 2>&1) || { echo "FAIL  $mode: script error"; cat "$WORK/$mode/log.txt"; status=1; continue; }
    for f in pas_transcripts.tsv pas_annotation.tsv pas_counts.tsv; do
        if cmp -s "$HERE/expected/$mode/$f" "$WORK/$mode/$f"; then echo "PASS  $mode/$f"
        else echo "FAIL  $mode/$f"; diff "$HERE/expected/$mode/$f" "$WORK/$mode/$f" || true; status=1; fi
    done
done
exit $status
