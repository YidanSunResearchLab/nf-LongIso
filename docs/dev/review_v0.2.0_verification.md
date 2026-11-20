# Verification of the methodological review of v0.2.0

Checked on 2026-11-20 against the v0.2.0 code (`834ded3`) and the demo run
`nf-LongIso_v2_runs/demo_v0.2.0` (GM12878 ×2 = control, H1 ×2 = treatment, 10,000 reads per sample,
GRCh38 / GENCODE / TEtranscripts TE GTF). No code was changed for these checks.

| Point | Review statement | Verdict |
|---|---|---|
| A1 | Junction strand is wrong: minimap2 writes no `XS` tag, so the read orientation is used | **Confirmed** |
| A2 | APA uses pseudo-counts built from condition means (IF × gene_value) and a Fisher test, not sample-level counts | **Confirmed** |
| A3 | TE DESeq2 estimates size factors from transcript and TE counts combined | **Confirmed** |
| A4 | proActiv is fed fake STAR SJ files; promoters of single-exon transcripts cannot be quantified | **Confirmed** |

## A1: junction strand

`EXTRACT_JUNCTIONS` takes the strand from `XS:A` and otherwise from the read orientation (FLAG 16).
minimap2 does not write `XS`; it writes `ts:A`, the transcript strand relative to the read
(`ts:A:+`: transcript on the read's alignment strand; `ts:A:-`: on the opposite strand).

All primary, mapped, spliced reads (CIGAR with `N`) of each demo BAM were examined (each BAM has fewer
than 10,000 spliced reads, so "the first 10,000" is all of them).

| Sample | Spliced reads | `XS` | `ts` | `ts:+` | `ts:-` | Reads whose current strand ≠ ts-derived strand |
|---|---|---|---|---|---|---|
| GM12878_demo_rep1 | 8,883 | 0 | 8,879 | 8,468 | 411 | 411 (4.6 %) |
| GM12878_demo_rep2 | 9,248 | 0 | 9,244 | 9,006 | 238 | 238 (2.6 %) |
| H1_demo_rep1 | 8,922 | 0 | 8,918 | 8,894 | 24 | 24 (0.3 %) |
| H1_demo_rep2 | 9,041 | 0 | 9,035 | 9,003 | 32 | 32 (0.4 %) |
| **Total** | 36,094 | **0** | 36,076 (99.95 %) | 35,371 | 705 | **705 (2.0 %)** |

Junction level (junctions seen in these reads; ts-derived strand = majority over their reads):

| Sample | Junctions | BED strand differs from the ts-derived strand | Junction written on both strands |
|---|---|---|---|
| GM12878_demo_rep1 | 12,887 | 3,231 (25.1 %) | 3,071 (23.8 %) |
| GM12878_demo_rep2 | 15,345 | 2,165 (14.1 %) | 2,036 (13.3 %) |
| H1_demo_rep1 | 15,972 | 61 (0.4 %) | 39 (0.2 %) |
| H1_demo_rep2 | 24,506 | 76 (0.3 %) | 17 (0.1 %) |

The demo reads are mostly already in transcript orientation (98 % `ts:+`), which keeps the read-level
error at 2 %. Because every read on the wrong strand creates a second copy of its junction, up to a
quarter of the GM12878 junctions carry a wrong or duplicated strand. For unoriented ONT cDNA reads
(about half of them reverse-complemented) the error would approach 50 %.

## A2: APA statistics

`bin/isoswitch.R` (APA part):

* `pm$c1 <- pm$IF1 * pm$gene_value_1` and `pm$c2 <- pm$IF2 * pm$gene_value_2` (lines 167-168): IF and
  gene_value come from `isoformFeatures`, i.e. **condition-level means**;
* these values are summed per PAS cluster, rounded and tested with `fisher.test(..., simulate.p.value = TRUE)`;
* the sample-level count matrix (`counts_data`) is only passed to `importRdata`; it is never used for APA.

The test therefore treats condition means as counts, ignores replicate variability and the
library-size scale of `gene_value`, and its p-values are not calibrated.

## A3: TE normalisation

`bin/te_deseq2.R` line 62: `combined <- rbind(as.matrix(tx_mat), as.matrix(te_mat))`; the DESeq2 object
(`DESeqDataSetFromMatrix(round(combined))`, line 72) estimates size factors from transcripts and TEs
together. TE counts are themselves sums of transcript counts (each read is counted again for every TE
its transcript overlaps), so TE rows are not independent of the transcript rows that set the size factors.

## A4: proActiv input

`bin/proactiv.R` converts long-read junctions to STAR `SJ.out.tab` rows with constant placeholder columns:
`motif = 1L` (every junction GT/AG), `annotated = 1L` (every junction annotated), `multi = 0L`,
`overhang = 10L` (lines 57-61); only the unique-read count is real.

Single-exon transcripts in the Bambu extended annotation:

| Item | Value |
|---|---|
| Transcripts | 388,075 |
| Single-exon transcripts | 28,988 (7.5 %; none novel) |
| Expressed single-exon transcripts (count > 0) | 143 of 6,205 expressed transcripts; 693 of 37,473 reads (1.8 %) |
| proActiv promoters made only of single-exon transcripts | 25,638 (101 with an expressed transcript, 610 reads) |
| proActiv counts of those 101 promoters | **NA in all 4 samples** |

proActiv does list promoters of single-exon transcripts, but it measures promoter activity from
first-intron junction reads, which these transcripts do not have, so their activity is always NA. In
addition, 3,183 promoters are shared by single- and multi-exon transcripts; only the multi-exon reads
contribute to their counts.
