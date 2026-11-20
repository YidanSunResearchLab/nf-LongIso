#!/bin/sh
# NextLongIso: transposable-element (TE) analysis of Bambu isoforms.
#
# Usage:
#   te_analysis.sh MODE MIN_OVERLAP STRAND AUTONOMOUS_FRAC ISO_GTF TE_GTF COUNTS REF_GTF JUNCTION_BED...
#
#   MODE             exon (default in the pipeline) | span (v0.1.0 counting rule, for comparison)
#   MIN_OVERLAP      minimum exon-TE overlap in bp (exon mode and annotation)
#   STRAND           sense (TE on the transcript strand) | both (exon mode and annotation)
#   AUTONOMOUS_FRAC  TE fraction of a transcript's exonic length for category 'autonomous'
#   ISO_GTF          Bambu extended_annotations.gtf
#   TE_GTF           TE annotation in TEtranscripts format: one record per TE instance with
#                    gene_id = subfamily, transcript_id = instance, family_id, class_id
#   COUNTS           Bambu counts_transcript.txt (TXNAME, GENEID, samples...)
#   REF_GTF          reference annotation (defines annotated exons for exonization)
#   JUNCTION_BED     per-sample junction BEDs (intron intervals, from EXTRACT_JUNCTIONS)
#
# Outputs
#   TE_family_counts.tsv        TE_Family, samples...            (same layout as v0.1.0)
#   TE_instance_counts.tsv      TE_id, TE_chr, TE_start, TE_end, samples...  (same layout as v0.1.0)
#   isoform_TE_annotation.tsv   one row per (transcript, TE instance) with exonic overlap
#   isoform_TE_exonization.tsv  strict exonization: pure-intronic novel exons >= 90 % covered by a
#                               same-strand TE, with a junction endpoint at both exon ends
#
# MODE only changes the TE family/instance counts; the annotation and exonization tables are the
# same in both modes.
#
# Only POSIX awk is used (tested with BusyBox awk in the bedtools container).

set -eu

if [ $# -lt 8 ]; then
    sed -n '2,30p' "$0" >&2
    exit 1
fi
MODE=$1; MIN_OV=$2; STRAND=$3; AUTO_FRAC=$4; ISO_GTF=$5; TE_GTF=$6; COUNTS=$7; REF_GTF=$8
shift 8
TAB=$(printf '\t')
export LC_ALL=C

case "$MODE" in exon|span) ;; *) echo "ERROR: mode must be 'exon' or 'span' (got '$MODE')" >&2; exit 1 ;; esac
case "$STRAND" in sense|both) ;; *) echo "ERROR: strand must be 'sense' or 'both' (got '$STRAND')" >&2; exit 1 ;; esac
echo "TE analysis: mode=$MODE min_overlap=${MIN_OV}bp strand=$STRAND autonomous_frac=$AUTO_FRAC junction_files=$#"
if [ "$MODE" = span ]; then
    echo "Note: span mode uses the v0.1.0 counting rule (same-strand overlap with the whole transcript span);"
    echo "      min_overlap and strand only apply to isoform_TE_annotation.tsv."
fi

# Multi-file awk programs identify their inputs with command-line assignments (fileno=1 a fileno=2 b);
# unlike counting FNR == 1, this stays correct when an input file is empty.

# GTF attribute lookup by name (no greedy sed): attr($9, "gene_id")
# Junction end lookup: a junction supports a transcript when it is on the transcript strand or its
# strand is unknown ('.'): has_end(je, chr, "S"|"E", position, transcript_strand)
AWK_ATTR='
function has_end(je, chr, type, pos, strand) {
    return ((chr SUBSEP type SUBSEP pos SUBSEP strand) in je) || ((chr SUBSEP type SUBSEP pos SUBSEP ".") in je)
}
function attr(s, key,    p) {
    if (match(s, "(^|[; ])" key " \"[^\"]*\"")) {
        p = substr(s, RSTART, RLENGTH)
        sub(/^[; ]*[^ ]+ "/, "", p)
        sub(/"$/, "", p)
        return p
    }
    return ""
}'

############################################################
# 1. Inputs as BED
############################################################
# Count table with exact '.sorted(.bam)' suffixes removed from sample names
awk 'BEGIN { FS = OFS = "\t" } NR == 1 { for (i = 1; i <= NF; i++) sub(/[.]sorted([.]bam)?$/, "", $i) } { print }' \
    "$COUNTS" > counts.tsv
SAMPLES=$(head -n 1 counts.tsv | cut -f3-)
NS=$(head -n 1 counts.tsv | awk -F '\t' '{ print NF - 2 }')

# Chromosome names shared by the TE GTF and the isoform GTF: none usually means a naming mismatch
# ('chr1' vs '1'), which would silently give empty TE tables
SHARED_CHR=$(awk -F '\t' '!/^#/ && NF >= 8 { print $1 }' "$TE_GTF" | sort -u > te_chr.txt
             awk -F '\t' '!/^#/ && NF >= 8 { print $1 }' "$ISO_GTF" | sort -u > iso_chr.txt
             comm -12 te_chr.txt iso_chr.txt | wc -l)
echo "Chromosomes: $(wc -l < te_chr.txt) in the TE GTF, $(wc -l < iso_chr.txt) in the isoform GTF, $SHARED_CHR shared"
if [ "$SHARED_CHR" -eq 0 ]; then
    echo "WARNING: the TE GTF and the isoform GTF have no chromosome name in common" \
         "(TE GTF: $(head -n 3 te_chr.txt | tr '\n' ' ')...; isoform GTF: $(head -n 3 iso_chr.txt | tr '\n' ' ')...)." \
         "All TE counts will be 0; use a TE GTF with the same chromosome names as --genome and --gtf." >&2
fi

# Reference exons (annotated exonic space)
awk -F '\t' -v OFS='\t' '$3 == "exon" && NF >= 8 && $4 >= 1 && $5 > $4 && $7 != "" { print $1, $4 - 1, $5, ".", ".", $7 }' \
    "$REF_GTF" > ref_exons_raw.bed
if [ ! -s ref_exons_raw.bed ]; then
    echo "ERROR: no exon records in $REF_GTF" >&2; exit 1
fi
bedtools sort -i ref_exons_raw.bed > ref_exons.bed

# Isoform exons and transcript spans: chr, start0, end, transcript_id, ., strand
awk -F '\t' -v OFS='\t' "$AWK_ATTR"'
    $3 == "exon" { print $1, $4 - 1, $5, attr($9, "transcript_id"), ".", $7 > "iso_exons.bed" }
    $3 == "transcript" { print $1, $4 - 1, $5, attr($9, "transcript_id"), ".", $7 > "isoforms.bed" }
' "$ISO_GTF"

# TE instances: chr, start0, end, instance, ., strand, subfamily, family, class
awk -F '\t' -v OFS='\t' "$AWK_ATTR"'
    ($3 == "exon" || $3 == "transcript") && $4 != "" && $5 != "" && $7 != "" {
        start = $4 - 1; if (start < 0) start = 0
        if (start >= $5) next
        print $1, start, $5, attr($9, "transcript_id"), ".", $7, attr($9, "gene_id"), attr($9, "family_id"), attr($9, "class_id")
    }
' "$TE_GTF" | sort -u > te_full.bed
cut -f1-6 te_full.bed | sort -u > te.bed
if [ ! -s te.bed ]; then
    echo "ERROR: no TE records parsed from $TE_GTF" >&2; exit 1
fi

# Observed junction endpoints (intron start0 = exon end; intron end = next exon start0)
if [ $# -gt 0 ]; then
    # chr, S|E, position, strand ('.' when the BED has no strand column or the strand is unknown)
    cat "$@" | awk -F '\t' -v OFS='\t' 'NF >= 3 { st = (NF >= 6 && ($6 == "+" || $6 == "-")) ? $6 : "."; print $1, "S", $2, st; print $1, "E", $3, st }' \
        | sort -u > junction_ends.txt
else
    : > junction_ends.txt
fi
echo "Junction endpoints: $(wc -l < junction_ends.txt)"

# Exons ranked in transcription order: chr, start0, end, tx, ., strand, rank, n_exons, genomic_index
sort -t "$TAB" -k4,4 -k2,2n iso_exons.bed > iso_exons.by_tx.bed
awk -F '\t' -v OFS='\t' '
    FNR == NR { n[$4]++; next }
    { g[$4]++; rank = ($6 == "-") ? n[$4] - g[$4] + 1 : g[$4]; print $0, rank, n[$4], g[$4] }
' iso_exons.by_tx.bed iso_exons.by_tx.bed > exons_ranked.bed
awk -F '\t' -v OFS='\t' '{ len[$4] += $3 - $2 } END { for (t in len) print t, len[t] }' iso_exons.bed > tx_exonic_length.txt

############################################################
# 2. Strict exonization
############################################################
# Novel exons outside all annotated exons and >= 90 % covered by a same-strand TE ...
bedtools intersect -v -a iso_exons.bed -b ref_exons.bed > iso_pure_intronic_exons.bed
bedtools intersect -s -f 0.9 -wa -wb -a iso_pure_intronic_exons.bed -b te.bed > exon_overlap_strict.txt

# ... and whose boundaries each coincide (+/- 5 bp) with an observed junction endpoint (a junction
# ending at the exon start and one starting at the exon end); every (transcript, exon, TE)
# combination is reported. Same test in both counting modes.
{
    printf 'transcript_id\tgene_id\texon_chr\texon_start\texon_end\tTE_id\tTE_chr\tTE_start\tTE_end\t%s\n' "$SAMPLES"
    awk -F '\t' -v OFS='\t' "$AWK_ATTR"'
        fileno == 1 { je[$1 SUBSEP $2 SUBSEP $3 SUBSEP $4] = 1; next }
        fileno == 2 {
            # exon: $1-$6, TE: $7-$12. Left boundary needs a junction END, right boundary a junction START.
            okl = 0; okr = 0
            for (d = -5; d <= 5; d++) {
                if (has_end(je, $1, "E", $2 + d, $6)) okl = 1
                if (has_end(je, $1, "S", $3 + d, $6)) okr = 1
            }
            if (okl && okr) { key = $4 SUBSEP $1 SUBSEP $2 SUBSEP $3 SUBSEP $10 SUBSEP $7 SUBSEP $8 SUBSEP $9; hit[key] = 1; tx[$4] = 1 }
            next
        }
        fileno == 3 && ($1 in tx) {
            for (k in hit) {
                split(k, a, SUBSEP)
                if (a[1] != $1) continue
                line = $1 OFS $2 OFS a[2] OFS a[3] OFS a[4] OFS a[5] OFS a[6] OFS a[7] OFS a[8]
                for (i = 3; i <= NF; i++) line = line OFS $i
                print line
            }
        }
    ' fileno=1 junction_ends.txt fileno=2 exon_overlap_strict.txt fileno=3 counts.tsv | sort
} > isoform_TE_exonization.tsv
echo "Strict exonization rows: $(( $(wc -l < isoform_TE_exonization.tsv) - 1 ))"

############################################################
# 3. Exon-level TE annotation (both modes)
############################################################
if [ "$STRAND" = sense ]; then STRAND_OPT="-s"; else STRAND_OPT=""; fi
# shellcheck disable=SC2086
bedtools intersect $STRAND_OPT -wo -a exons_ranked.bed -b te_full.bed > exon_te_overlap.txt

{
    printf 'transcript_id\tgene_id\tTE_id\tsubfamily\tfamily\tclass\tTE_chr\tTE_start\tTE_end\tTE_strand\toverlap_bp\tTE_fraction_of_transcript\tstrand_relation\tcategory\texon_positions\tjunction_supported\n'
    awk -F '\t' -v OFS='\t' -v minov="$MIN_OV" -v afrac="$AUTO_FRAC" "$AWK_ATTR"'
        fileno == 1 { je[$1 SUBSEP $2 SUBSEP $3 SUBSEP $4] = 1; next }
        fileno == 2 { txlen[$1] = $2; next }
        fileno == 3 { if (FNR > 1) gene[$1] = $2; next }
        fileno == 4 {
            # exon: $1-$9 (rank $7, n_exons $8, genomic index $9); TE: $10-$18; overlap bp: $19
            k = $4 SUBSEP $13
            if (!(k in tot)) {
                order[++np] = k
                te[k] = $13 OFS $16 OFS $17 OFS $18 OFS $10 OFS $11 OFS $12 OFS $15
                rel[k] = ($6 == $15) ? "sense" : "antisense"
                sup[k] = "NA"
            }
            tot[k] += $19
            if ($7 == 1) ovf[k] += $19
            if ($7 == $8) ovl[k] += $19
            if ($7 > 1 && $7 < $8 && $19 > ovi[k]) ovi[k] = $19
            # Exon positions (in transcription order) touched by this TE
            if ($7 == 1) pf[k] = 1
            if ($7 == $8) pl[k] = 1
            if ($7 > 1 && $7 < $8) pi[k] = 1
            # Splice sites of this exon: left boundary unless genomically first, right unless last
            if ($8 > 1) {
                ok = 1
                if ($9 > 1) { f = 0; for (d = -5; d <= 5; d++) if (has_end(je, $1, "E", $2 + d, $6)) f = 1; if (!f) ok = 0 }
                if ($9 < $8) { f = 0; for (d = -5; d <= 5; d++) if (has_end(je, $1, "S", $3 + d, $6)) f = 1; if (!f) ok = 0 }
                if (sup[k] == "NA") sup[k] = ok ? "yes" : "no"
                else if (!ok) sup[k] = "no"
            }
        }
        END {
            for (i = 1; i <= np; i++) {
                k = order[i]
                if (tot[k] < minov) continue
                split(k, a, SUBSEP)
                frac = (txlen[a[1]] > 0) ? tot[k] / txlen[a[1]] : 0
                if (frac >= afrac)            cat = "autonomous"
                else if (ovf[k] + 0 >= minov) cat = "TSS"
                else if (ovl[k] + 0 >= minov) cat = "PAS"
                else if (ovi[k] + 0 >= minov) cat = "internal_exonization"
                else                          cat = "other_exonic"
                g = (a[1] in gene) ? gene[a[1]] : "NA"
                pos = ""
                if (k in pf) pos = "first"
                if (k in pi) pos = pos (pos == "" ? "" : ";") "internal"
                if (k in pl) pos = pos (pos == "" ? "" : ";") "last"
                printf "%s\t%s\t%s\t%d\t%.4f\t%s\t%s\t%s\t%s\n", a[1], g, te[k], tot[k], frac, rel[k], cat, pos, sup[k]
            }
        }
    ' fileno=1 junction_ends.txt fileno=2 tx_exonic_length.txt fileno=3 counts.tsv fileno=4 exon_te_overlap.txt | sort
} > isoform_TE_annotation.tsv
echo "TE annotation rows (transcript, TE instance): $(( $(wc -l < isoform_TE_annotation.tsv) - 1 ))"
awk -F '\t' 'NR > 1 { n[$14]++ } END { for (c in n) printf "  %s: %d\n", c, n[c] }' isoform_TE_annotation.tsv | sort

############################################################
# 4. TE family / instance counts
############################################################
printf 'TE_Family\t%s\n' "$SAMPLES" > TE_family_counts.tsv
printf 'TE_id\tTE_chr\tTE_start\tTE_end\t%s\n' "$SAMPLES" > TE_instance_counts.tsv

if [ "$MODE" = span ]; then
    # v0.1.0 counting rule: same-strand overlap of the whole transcript span with a TE; when a
    # transcript overlaps several TEs, only the last overlap reported by bedtools is kept (data[$4]
    # overwrite). Kept as a comparison option; rows are sorted.
    bedtools intersect -s -wa -wb -a isoforms.bed -b te.bed > exon_overlap_general.txt
    awk -F '\t' 'NR == FNR { data[$4] = $1 "\t" $2 "\t" $3 "\t" $10 "\t" $7 "\t" $8 "\t" $9; next }
        ($1 in data) { printf "%s\t%s\t%s", $1, $2, data[$1]; for (i = 3; i <= NF; i++) printf "\t%s", $i; printf "\n" }' \
        exon_overlap_general.txt counts.tsv > all_joined.tsv
    awk -F '\t' -v n="$NS" '{
        fam = $6; sub(/_dup[0-9]+/, "", fam); for (i = 1; i <= n; i++) sum[fam, i] += $(i + 9); ids[fam] = 1
    } END { for (f in ids) { printf "%s", f; for (j = 1; j <= n; j++) printf "\t%s", (sum[f, j] ? sum[f, j] : 0); printf "\n" } }' \
        all_joined.tsv | sort >> TE_family_counts.tsv
    awk -F '\t' -v n="$NS" '{
        inst = $6 "\t" $7 "\t" $8 "\t" $9; for (i = 1; i <= n; i++) inst_sum[inst, i] += $(i + 9); locs[inst] = 1
    } END { for (l in locs) { printf "%s", l; for (j = 1; j <= n; j++) printf "\t%s", (inst_sum[l, j] ? inst_sum[l, j] : 0); printf "\n" } }' \
        all_joined.tsv | sort >> TE_instance_counts.tsv
else
    # Exon mode: (transcript, TE) pairs from the annotation (exonic overlap >= MIN_OVERLAP). Each
    # transcript contributes its reads once per subfamily and once per TE instance.
    awk -F '\t' -v OFS='\t' -v n="$NS" '
        FNR == NR { if (FNR > 1) { row[$1] = $0 } next }
        FNR == 1 { next }
        ($1 in row) {
            split(row[$1], c, "\t")
            if (!(($1 SUBSEP $4) in seen_fam)) {
                seen_fam[$1 SUBSEP $4] = 1; fams[$4] = 1
                for (i = 1; i <= n; i++) fsum[$4, i] += c[i + 2]
            }
            inst = $3 OFS $7 OFS $8 OFS $9
            if (!(($1 SUBSEP inst) in seen_inst)) {
                seen_inst[$1 SUBSEP inst] = 1; insts[inst] = 1
                for (i = 1; i <= n; i++) isum[inst, i] += c[i + 2]
            }
        }
        END {
            for (f in fams)  { line = f; for (i = 1; i <= n; i++) line = line OFS sprintf("%.10g", fsum[f, i] + 0); print line > "fam.tmp" }
            for (l in insts) { line = l; for (i = 1; i <= n; i++) line = line OFS sprintf("%.10g", isum[l, i] + 0); print line > "inst.tmp" }
        }
    ' counts.tsv isoform_TE_annotation.tsv
    touch fam.tmp inst.tmp
    sort fam.tmp >> TE_family_counts.tsv
    sort inst.tmp >> TE_instance_counts.tsv
    rm -f fam.tmp inst.tmp
fi
echo "TE families with counts: $(( $(wc -l < TE_family_counts.tsv) - 1 )); TE instances: $(( $(wc -l < TE_instance_counts.tsv) - 1 ))"
