# Promoter usage: proActiv vs transcript-5'-end (tss) method

Demo comparison: `demo_v0.2.0_r2` (default `--promoter_method tss`) vs the proActiv results of
`demo_v0.2.0` (proActiv code unchanged since then apart from the move to `bin/`; the run with the
current code was not repeated).

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

| | proActiv (`demo_v0.2.0/proactiv/`) | tss (`demo_v0.2.0_r2/promoter/`) |
|---|---|---|
| Promoters | 160,497 | 199,302 in 78,936 genes |
| NA in all samples | 25,638 (single-exon-only promoters) | 0 |
| Class NA (internal or unquantified) | 67,931 | 0 (120,364 internal promoters are classified and flagged) |
| Active (Major + Minor), control / treatment | 890 / 1,208 | 2,591 / 3,241 |
| Active promoters made only of single-exon transcripts | 0 (always NA) | 112 |
| Differential usage | no alternative promoters (`getAlternativePromoters`) | 72 promoters in 35 genes tested; **3 genes** with gene q < 0.05 |

The three genes with differential promoter usage (ENSG00000144381 HSPD1, ENSG00000188846 RPL14,
ENSG00000169764 UGP2) all have an IsoformSwitchAnalyzeR ATSS event (842 genes with ATSS events).
RPL14 is also one of the three significant APA genes.

### Notes on the proActiv input

`bin/run_proactiv.R` writes STAR `SJ.out.tab` files with constant `motif` and `annotated` columns.
proActiv reads them with `GenomicAlignments::readSTARJunctions()` and counts, per intron, the unique
reads of the junction with exactly the same coordinates and a compatible strand; the `motif` and
`annotated` columns are not used, so the placeholders do not change its counts. The junctions come
from `summarizeJunctions()` without a genome and therefore have strand `*`, which matches introns on
either strand; this only matters for identical intron coordinates on both strands. Feeding proActiv
the stranded junction tables of `EXTRACT_JUNCTIONS` would remove this; it is not done in v0.2.0.
