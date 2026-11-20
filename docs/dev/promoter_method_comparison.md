# Promoter usage: proActiv vs transcript-5'-end (tss) method

Status: method description complete; demo numbers are filled in after the demo run of review
phase 12 (`demo_v0.2.0_r2` with the default `--promoter_method tss`; proActiv results from
`demo_v0.2.0`).

## Methods

| | `--promoter_method proactiv` | `--promoter_method tss` (default) |
|---|---|---|
| Promoter definition | proActiv annotation from the GTF (first exons of transcripts) | clusters of Bambu transcript 5' ends per gene and strand (`--tss_cluster_window` 50 bp, maximal width `--tss_max_cluster_width` 100 bp) |
| Reads used | splice junctions converted to STAR `SJ.out.tab` files with placeholder `motif`/`annotated` columns | Bambu read counts of the transcripts starting in the cluster |
| Single-exon transcripts | annotated but never quantified (no first-intron junction): their promoters get NA | included |
| Absolute activity | log2(normalised junction count + 1), DESeq2 size factors from promoter counts | log2(normalised count + 1), DESeq2 size factors from gene-level counts |
| Relative activity | absolute activity / sum of the gene's absolute activities | promoter count / gene count |
| Classes | Major / Minor / Inactive (threshold 0.25); internal promoters not classified | Major / Minor / Inactive (`--promoter_min_activity` 0.25); internal promoters classified and flagged |
| Differential usage | linear models on absolute and relative activity (`getAlternativePromoters`), no covariates | DEXSeq (`bin/usage_dexseq.R`, same as APA): `~ sample + exon + <covariates>:exon + condition:exon`, promoter padj, gene q |
| Annotation columns | proActiv rowData | annotated, internal_promoter, dist_to_annotated_tss, n_single_exon |

Caveat (README): ONT cDNA reads are often truncated at the 5' end, which moves or splits TSS clusters;
CAGE data, where available, should be used to confirm TSS-level results.

## Pre-run figures (demo Bambu annotation, `tss_sites.R` only)

199,302 promoters in 78,936 genes; 28,993 genes with ≥ 2 promoters; 26,741 promoters made only of
single-exon transcripts (proActiv: 25,638 single-exon-only promoters, all NA, see
`docs/dev/review_v0.2.0_verification.md`); 120,364 promoters lie inside the body of another transcript
of the same gene (internal promoters).

## Demo comparison

To be completed after the demo run: promoters detected and active per method, contribution of
single-exon transcripts, overlap of differential promoters (proActiv `alternative_promoters.tsv` vs tss
`alternative_promoters.tsv`), and agreement with the IsoformSwitchAnalyzeR ATSS events
(`isoswitch/atss_results.csv`).
