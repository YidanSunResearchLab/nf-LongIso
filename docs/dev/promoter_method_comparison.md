# Promoter usage: proActiv (default) and transcript-5'-end (tss) method

proActiv is the default promoter method (`--promoter_method proactiv`, `results/proactiv/`, as in
v0.1.0). The transcript-5'-end method is optional (`--promoter_method tss`, `results/promoter_tss/`) and
complements proActiv where its junction-based design cannot see promoters.

## Methods

| | `--promoter_method proactiv` (default) | `--promoter_method tss` (optional) |
|---|---|---|
| Promoter definition | proActiv annotation from the GTF (first exons of transcripts) | clusters of Bambu transcript 5' ends per gene and strand (`--tss_cluster_window` 50 bp, maximal width `--tss_max_cluster_width` 100 bp) |
| Reads used | junction reads of the first intron: stranded junction tables of `EXTRACT_JUNCTIONS` converted to STAR `SJ.out.tab` (strand `+` = 1, `-` = 2, real motif and annotation codes; junctions of unknown strand not used) | Bambu read counts of the transcripts starting in the cluster |
| Single-exon transcripts | annotated but not quantified (no first-intron junction): their promoters get NA | included |
| Absolute activity | log2(normalised junction count + 1), DESeq2 size factors from promoter counts | log2(normalised count + 1), DESeq2 size factors from gene-level counts |
| Relative activity | absolute activity / sum of the gene's absolute activities | promoter count / gene count |
| Classes | Major / Minor / Inactive (threshold 0.25); internal promoters not classified | Major / Minor / Inactive (`--promoter_min_activity` 0.25); internal promoters classified and flagged |
| Internal promoter | the intron after the promoter's first exon is a later (not the first) intron of another transcript, i.e. the first exon is an internal exon there (`MaxMergedIntronRank > 1`) | TSS in a non-first exon (transcription order) of another transcript of the same gene (same idea, defined on the TSS position because single-exon transcripts have no intron) |
| Differential usage | linear models on absolute and relative activity (`getAlternativePromoters`), no covariates | DEXSeq (`bin/usage_dexseq.R`, same as APA): `~ sample + exon + <covariates>:exon + condition:exon`, promoter padj, gene q |
| Annotation columns | proActiv rowData | annotated, internal_promoter, dist_to_annotated_tss, n_single_exon |

Limitations of proActiv for long reads (README): promoter activity is inferred from junctions, so
promoters of single-exon transcripts are not detected, and the transcription start sites contained in
full-length reads are not used. Limitation of the tss method: ONT cDNA reads are often truncated at the
5' end, which moves or splits TSS clusters; CAGE data, where available, should be used to confirm
TSS-level results.

## proActiv strand fix

### Cause

Up to and including the v0.2.0 development builds (and in v0.1.0), the junctions were taken from the
BAM files with `summarizeJunctions()` without a genome. Their strand is `*`, which was written as STAR
strand code 0. proActiv reads the files with `GenomicAlignments::readSTARJunctions()`, which sets the
strand with `strand(column4 == 2)`: code 0 becomes `+`. proActiv then counts, per intron, the reads of
the junction with the same coordinates on the same strand (`findOverlaps(type = "equal")`), so
minus-strand introns never received reads. (An earlier version of this document stated that `*`
matches introns on either strand; that was wrong.) The `motif` and `annotated` columns are not used by
proActiv for counting; they now carry the real codes.

### Effect

| Data set | Run | Minus-strand promoters | with counts | Active (control / treatment) | Plus-strand promoters with counts | Alternative promoters |
|---|---|---|---|---|---|---|
| test data | before (`tiny_p7_proactiv`) | 505 | 0 | 0 / 0 | 188 | 0 |
| test data | after (`tiny_r3_default`) | 505 | 172 | 33 / 55 | 188 | 0 |
| demo | v0.1.0 baseline and `demo_v0.2.0` | 79,200 | 0 | 0 / 0 | 4,539 (active 890 / 1,208) | 0 |
| demo | after the fix | to be measured (demo re-run pending) | | | | |

On the test data, 34 of the 517 plus-strand promoters lose some reads (52 of 1,800 reads in those
promoters, no promoter gains reads): reads of minus-strand junctions with the same coordinates as a
plus-strand intron, of supplementary alignments, and of junctions of unknown strand are no longer
counted for them. The BAM-input run gives byte-identical proActiv counts.

## tss method on the demo (`demo_v0.2.0_r2`, before the changes of this round)

Figures from `demo_v0.2.0_r2/promoter/` (then the default output folder). Since then, internal promoters
are defined as above (TSS in a non-first exon of another transcript; on the test data 254 instead of
822 of 1,260 promoters), and annotated transcripts are recognised by the transcript IDs of `--gtf`.

| | proActiv (`demo_v0.2.0/proactiv/`, before the strand fix) | tss (`demo_v0.2.0_r2/promoter/`) |
|---|---|---|
| Promoters | 160,497 | 199,302 in 78,936 genes |
| NA in all samples | 25,638 (single-exon-only promoters) | 0 |
| Active (Major + Minor), control / treatment | 890 / 1,208 (plus strand only) | 2,591 / 3,241 |
| Active promoters made only of single-exon transcripts | 0 (always NA) | 112 |
| Differential usage | no alternative promoters (`getAlternativePromoters`) | 72 promoters in 35 genes tested; **3 genes** with gene q < 0.05 |

The three genes with differential promoter usage (ENSG00000144381 HSPD1, ENSG00000188846 RPL14,
ENSG00000169764 UGP2) all have an IsoformSwitchAnalyzeR ATSS event (842 genes with ATSS events).
RPL14 is also one of the three significant APA genes.
