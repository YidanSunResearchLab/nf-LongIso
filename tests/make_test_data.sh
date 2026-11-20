#!/bin/sh
# Build the tiny test dataset in tests/data/ used by -profile test.
#
# The dataset is cut from the demo run of NextLongIso:
#   * genome: eight 1-Mb GRCh38 windows (six on chr1, two on chr2) that contain reads from
#     all four demo samples. Windows on the same chromosome are joined with 1,000-N spacers
#     so that contigs keep standard names (chr1, chr2), which proActiv requires.
#   * annotation / TE annotation: GENCODE transcripts and TE records lying entirely inside
#     a window, with coordinates shifted to the new contigs. The gene GTF keeps only gene,
#     transcript and exon records and the attributes the pipeline uses.
#   * reads: primary alignments in the windows, taken from the demo BAMs and written back
#     to FASTQ (2 conditions x 2 replicates).
#
# Requirements: samtools, gzip and a POSIX awk (tested with BusyBox awk), e.g.
#   apptainer exec -B /path quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1 \
#       sh tests/make_test_data.sh <demo_bam_dir> <genome.fa> <genes.gtf> <te.gtf> <outdir>
#
# <demo_bam_dir> must contain {GM12878,H1}_demo_rep{1,2}.sorted.bam (+ .bai), i.e. the
# results/minimap2/ directory of a demo run.
#
# The published copy of this dataset is on Zenodo (see conf/test.config); its files are genome.fa,
# genes.gtf, te.gtf and reads/*.fastq.gz from <outdir> (uploaded without the reads/ folder). The
# samplesheet written here uses local paths (tests/data/...); tests/samplesheet_test.csv uses the URLs.
# <outdir> = tests/data is git-ignored.

set -eu

if [ $# -ne 5 ]; then
    echo "Usage: $0 <demo_bam_dir> <genome.fa> <genes.gtf> <te.gtf> <outdir>" >&2
    exit 1
fi
BAM_DIR=$1; GENOME=$2; GTF=$3; TE_GTF=$4; OUT=$5
SPACER=1000

mkdir -p "$OUT/reads"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# Windows (1-based, inclusive), in output order
cat > "$TMP/windows.txt" <<'EOF'
chr1 8000001 9000000
chr1 23000001 24000000
chr1 26000001 27000000
chr1 32000001 33000000
chr1 154000001 155000000
chr1 156000001 157000000
chr2 27000001 28000000
chr2 197000001 198000000
EOF

# Offsets: new coordinate = old - start + new_start
awk -v sp="$SPACER" '
    $1 != chr { chr = $1; pos = 1 }
    { print $1, $2, $3, pos; pos += ($3 - $2 + 1) + sp }
' "$TMP/windows.txt" > "$TMP/offsets.txt"

# 1. Genome
: > "$TMP/genome.raw"
for chr in $(awk '{ print $1 }' "$TMP/windows.txt" | uniq); do
    printf '>%s\n' "$chr" >> "$TMP/genome.raw"
    first=1
    awk -v c="$chr" '$1 == c { print $1 ":" $2 "-" $3 }' "$TMP/windows.txt" | while read -r region; do
        if [ $first -eq 0 ]; then
            awk -v n="$SPACER" 'BEGIN { s = ""; for (i = 0; i < n; i++) s = s "N"; print s }' >> "$TMP/genome.raw"
        fi
        samtools faidx "$GENOME" "$region" | awk 'NR > 1' >> "$TMP/genome.raw"
        first=0
    done
done
# Re-wrap sequence lines to 60 bp
awk '
    /^>/ { if (seq != "") flush(); print; next }
    { seq = seq $0 }
    function flush() { while (length(seq) > 60) { print substr(seq, 1, 60); seq = substr(seq, 61) } if (seq != "") print seq; seq = "" }
    END { flush() }
' "$TMP/genome.raw" > "$OUT/genome.fa"
samtools faidx "$OUT/genome.fa"

# 2. Gene annotation: keep transcripts entirely inside a window (and their exons/CDS/UTRs),
#    plus gene records entirely inside a window that keep at least one transcript.
awk '$1 == "chr1" || $1 == "chr2"' "$GTF" > "$TMP/genes.chr.gtf"
awk -F '\t' -v OFS='\t' '
    function attr(s, key,    r) {
        if (match(s, key " \"[^\"]+\"")) { r = substr(s, RSTART + length(key) + 2, RLENGTH - length(key) - 3); return r }
        return ""
    }
    function win(c, s, e,    i) {
        for (i = 1; i <= nw; i++) if (wc[i] == c && s >= ws[i] && e <= we[i]) return i
        return 0
    }
    BEGIN { nkeys = split("gene_id transcript_id gene_type gene_name transcript_type transcript_name", keys, " ") }
    fileno == 1 { split($0, a, " "); nw++; wc[nw] = a[1]; ws[nw] = a[2]; we[nw] = a[3]; wn[nw] = a[4]; next }
    fileno == 2 {
        if ($3 == "transcript" && win($1, $4, $5)) { keep_tx[attr($9, "transcript_id")] = 1; keep_gene[attr($9, "gene_id")] = 1 }
        next
    }
    fileno == 3 {
        if ($3 != "gene" && $3 != "transcript" && $3 != "exon") next
        i = win($1, $4, $5)
        if (!i) next
        if ($3 == "gene") { if (!(attr($9, "gene_id") in keep_gene)) next }
        else if (!(attr($9, "transcript_id") in keep_tx)) next
        $4 = $4 - ws[i] + wn[i]; $5 = $5 - ws[i] + wn[i]
        # Keep only the attributes the pipeline uses, to keep the file small
        a = ""
        for (k = 1; k <= nkeys; k++) {
            v = attr($9, keys[k])
            if (v != "") a = a keys[k] " \"" v "\"; "
        }
        if (match($9, /exon_number [0-9]+/)) a = a substr($9, RSTART, RLENGTH) "; "
        sub(/ $/, "", a)
        $9 = a
        print
    }
' fileno=1 "$TMP/offsets.txt" fileno=2 "$TMP/genes.chr.gtf" fileno=3 "$TMP/genes.chr.gtf" > "$OUT/genes.gtf"

# 3. TE annotation: records entirely inside a window
awk -F '\t' -v OFS='\t' '
    fileno == 1 { split($0, a, " "); nw++; wc[nw] = a[1]; ws[nw] = a[2]; we[nw] = a[3]; wn[nw] = a[4]; next }
    {
        for (i = 1; i <= nw; i++) {
            if ($1 == wc[i] && $4 >= ws[i] && $5 <= we[i]) {
                $4 = $4 - ws[i] + wn[i]; $5 = $5 - ws[i] + wn[i]; print; next
            }
        }
    }
' fileno=1 "$TMP/offsets.txt" fileno=2 "$TE_GTF" > "$OUT/te.gtf"

# 4. Reads: primary alignments overlapping the windows, as FASTQ
REGIONS=$(awk '{ printf "%s:%s-%s ", $1, $2, $3 }' "$TMP/windows.txt")
printf 'sample,fastq,condition\n' > "$OUT/samplesheet.csv"
for s in GM12878_demo_rep1 GM12878_demo_rep2 H1_demo_rep1 H1_demo_rep2; do
    name=$(echo "$s" | sed 's/_demo//')
    cond=control; case "$s" in H1_*) cond=treatment ;; esac
    # shellcheck disable=SC2086
    samtools view -b -F 2308 "$BAM_DIR/$s.sorted.bam" $REGIONS \
        | samtools sort -n -o "$TMP/$s.byname.bam" -
    samtools fastq "$TMP/$s.byname.bam" 2> /dev/null | gzip -n -9 > "$OUT/reads/$name.fastq.gz"
    printf '%s,tests/data/reads/%s.fastq.gz,%s\n' "$name" "$name" "$cond" >> "$OUT/samplesheet.csv"
    echo "$name: $(samtools view -c "$TMP/$s.byname.bam") reads"
done

echo "Test data written to $OUT:"
ls -l "$OUT" "$OUT/reads"
