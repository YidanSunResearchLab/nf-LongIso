#!/bin/sh
# host_pairs.sh <extended_annotations.gtf> <te.gtf> <out.tsv>: (gene_id, TE subfamily) for genes whose
# span (first to last exon of all its isoforms) overlaps a TE instance on the same strand
set -eu
awk -F '\t' -v OFS='\t' '$3 == "exon" { match($9, /gene_id "[^"]+"/); g = substr($9, RSTART + 9, RLENGTH - 10)
    k = g; if (!(k in s) || $4 - 1 < s[k]) s[k] = $4 - 1; if ($5 > e[k]) e[k] = $5; c[k] = $1; st[k] = $7 }
    END { for (k in s) print c[k], s[k], e[k], k, ".", st[k] }' "$1" | sort -k1,1 -k2,2n > genes_span.bed
awk -F '\t' -v OFS='\t' '$3 == "exon" { match($9, /gene_id "[^"]+"/); sf = substr($9, RSTART + 9, RLENGTH - 10); print $1, $4 - 1, $5, sf, ".", $7 }' "$2" > te_sf.bed
bedtools intersect -s -wa -wb -a genes_span.bed -b te_sf.bed | awk -F '\t' -v OFS='\t' '{ print $4, $10 }' | sort -u > "$3"
