#!/bin/sh
# Download the test dataset of -profile test from Zenodo into a local directory.
#
# Usage: sh tests/fetch_test_data.sh <outdir>
#        nextflow run main.nf -profile test,docker -c <outdir>/local.config --outdir results_test
#
# -profile test reads the files directly from Zenodo. Zenodo limits the request rate per IP address
# (shared CI runners and some institutional networks are affected); a refused request then shows up
# as "the file or directory 'https://zenodo.org/...' does not exist" in the parameter validation.
# This script downloads the files with retries, checks them against tests/test_data.md5 and writes
#   <outdir>/local.config                  genome, gtf, te_gtf and samplesheet pointing to <outdir>
#   <outdir>/samplesheet_test*.csv          copies of tests/samplesheet_test*.csv with local paths
# so that any of the test samplesheets can be used with --samplesheet <outdir>/<name>.csv.

set -eu

if [ $# -ne 1 ]; then
    echo "Usage: $0 <outdir>" >&2
    exit 1
fi
REPO=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$1"
OUT=$(cd "$1" && pwd)

# Zenodo record used by conf/test.config
BASE=$(grep -o "https://zenodo.org/records/[0-9]*/files/" "$REPO/conf/test.config" | head -n 1)
if [ -z "$BASE" ]; then
    echo "No Zenodo URL found in conf/test.config" >&2
    exit 1
fi

while read -r md5 name; do
    if [ -s "$OUT/$name" ] && echo "$md5  $OUT/$name" | md5sum -c --status -; then
        echo "$name: present"
        continue
    fi
    echo "$name: downloading $BASE$name"
    curl -fsSL --retry 6 --retry-delay 20 --retry-all-errors -o "$OUT/$name" "$BASE$name"
    if ! echo "$md5  $OUT/$name" | md5sum -c --status -; then
        echo "$name: md5 mismatch (expected $md5)" >&2
        exit 1
    fi
done < "$REPO/tests/test_data.md5"

for ss in "$REPO"/tests/samplesheet_test*.csv; do
    sed "s#$BASE#$OUT/#" "$ss" > "$OUT/$(basename "$ss")"
done

cat > "$OUT/local.config" <<EOF
// Written by tests/fetch_test_data.sh: local copy of the test dataset ($BASE)
params {
    samplesheet = '$OUT/samplesheet_test.csv'
    genome      = '$OUT/genome.fa'
    gtf         = '$OUT/genes.gtf'
    te_gtf      = '$OUT/te.gtf'
}
EOF
echo "Test data in $OUT; use: -profile test,<engine> -c $OUT/local.config"
