# TE quantification: exon mode vs span mode

NextLongIso v0.2.0 counts TE expression in two ways (`--te_count_mode`):

* **exon** (default): an isoform is counted for a TE when one of its **exons** overlaps the TE by at
  least `--te_min_overlap` bp (20) on the same strand. Each isoform is counted once per subfamily and
  once per TE instance it contains.
* **span** (v0.1.0): an isoform is counted for a TE that overlaps its **span** (first to last base,
  including introns) on the same strand; when several TEs overlap, only the last one reported by
  bedtools receives the counts.

## Data and commands

Demo data (GM12878 ×2 vs H1 ×2, 10,000 long reads per sample), GRCh38, GENCODE annotation and
`TEtranscript_hg38_rmsk_TE.gtf`; runs `demo_v0.2.0` (exon) and `demo_v0.2.0_span` (span, `-resume`), see
`docs/dev/v0.2.0_vs_v0.1.0.md`. The span run's TE counts and TE DESeq2 results are identical to v0.1.0.

```bash
# host genes of each subfamily (bedtools container)
sh tests/te_host_pairs.sh demo_v0.2.0/bambu/extended_annotations.gtf TEtranscript_hg38_rmsk_TE.gtf host_pairs.tsv
# report (the span work directory is that of the span run's TE_ANALYSIS task)
python3 tests/compare_te_modes.py --exon demo_v0.2.0 --span demo_v0.2.0_span \
    --host-pairs host_pairs.tsv --span-workdir <work dir of TE_ANALYSIS in the span run>
```

## Main observations

1. **Span counts mostly measure host transcripts, not TE-containing RNA.** In span mode, 79.1 % of the
   counted isoform-TE assignments, carrying **89.3 % of all span-mode counts**, involve a TE that does not
   overlap any exon of the isoform (intronic TE). Span-mode subfamily counts also correlate more strongly
   with the expression of their host genes (Spearman 0.67 vs 0.58 over counted subfamilies).
2. **The two modes measure different things.** Subfamily counts correlate only moderately between modes
   (Spearman 0.57 over the 321 subfamilies counted in both); 146 subfamilies are counted only in exon mode
   and 88 only in span mode. Exon mode counts fewer reads in total (19,562 vs 33,730) because intronic
   overlaps are no longer counted, but more subfamilies (467 vs 409) because an isoform containing
   several TEs is now counted for each of them.
3. **Differential expression.** Exon mode finds HERVH-int (log2FC +10.0, padj 0.004) and LTR7 (log2FC +8.6,
   padj 0.019) up in H1. HERVH elements and their LTR7 promoters are known to be highly transcribed in
   human embryonic stem cells, so this is the expected signal. Span mode reports no significant subfamily
   (its top-20 list consists of ties at padj 0.925), and the two top-20 lists share no subfamily.
4. **Categories** (`isoform_TE_annotation.tsv`): most TE-containing isoforms carry the TE in their last
   exon (PAS, 51.8 % of isoform-TE pairs), followed by internal exons (23.7 %), first exons (TSS, 20.7 %) and
   TE-dominated (autonomous, 3.8 %) isoforms.

Caveat: the demo has only 10,000 reads per sample, so absolute numbers are small; the comparison should
be repeated on a full dataset before drawing biological conclusions.

## Full report
### 1. Counted features and totals

| Metric | exon | span |
|---|---|---|
| subfamilies with counts > 0 | 467 | 409 |
| total counts over subfamilies (all samples) | 19,562 | 33,730 |
| TE instances with counts > 0 | 3,353 | 3,561 |
| total counts over TE instances (all samples) | 22,472 | 33,730 |
| subfamilies counted in both modes | 321 | |
| subfamilies counted in one mode only | 146 (exon only) | 88 (span only) |

### 2. Spearman correlation between modes

Counts summed over samples. 'Union' sets missing features to 0.

| Level | Features (union) | rho (union) | Features (both) | rho (both) |
|---|---|---|---|---|
| subfamily | 1,127 | 0.652 | 321 | 0.574 |
| instance | 177,726 | 0.155 | 595 | 0.897 |

### 3. TE counts vs host-gene expression

Host genes of a subfamily: genes whose span overlaps one of its instances on the same strand (defined independently of both counting modes). Columns 2-3: all subfamilies with host-gene expression > 0 (TE count 0 included). Columns 4-6: only subfamilies counted (> 0) in that mode; ratio = TE count / host-gene expression (a ratio near 1 means the TE count is essentially the expression of its host genes).

| Mode | Subfamilies (host > 0) | rho | Counted subfamilies | rho | median ratio |
|---|---|---|---|---|---|
| exon | 1,076 | 0.669 | 467 | 0.579 | 0.022 |
| span | 1,076 | 0.727 | 409 | 0.672 | 0.025 |

Span mode assigns each transcript to one TE overlapping its span. For these assignments:

| | Assignments | Share of span-mode counts |
|---|---|---|
| TE overlaps an exon of the transcript | 67,393 (20.9%) | 10.7% |
| TE only in introns of the transcript | 255,203 (79.1%) | 89.3% |

### 4. TE differential expression (padj < 0.05)

| Level | exon | span | shared |
|---|---|---|---|
| Family | 2 | 0 | 0 |
| Instance | 0 | 1 | 0 |

Top 20 subfamilies by padj:

| Rank | exon | padj | span | padj |
|---|---|---|---|---|
| 1 | HERVH-int | 4.30e-03 | 7SK | 9.25e-01 |
| 2 | LTR7 | 1.91e-02 | Alu | 9.25e-01 |
| 3 | MER67C | 9.45e-02 | AluJb | 9.25e-01 |
| 4 | Ricksha_a | 9.45e-02 | AluSc | 9.25e-01 |
| 5 | THE1A | 1.55e-01 | AluSc5 | 9.25e-01 |
| 6 | Charlie10a | 1.69e-01 | AluSg4 | 9.25e-01 |
| 7 | Ricksha_c | 2.15e-01 | AluSg7 | 9.25e-01 |
| 8 | LTR2B | 3.51e-01 | AluSq10 | 9.25e-01 |
| 9 | L1MD3 | 4.71e-01 | AluSq2 | 9.25e-01 |
| 10 | L1PA4 | 4.71e-01 | AluSq4 | 9.25e-01 |
| 11 | L1MD2 | 4.77e-01 | AluSx4 | 9.25e-01 |
| 12 | MER110A | 5.28e-01 | AluSz6 | 9.25e-01 |
| 13 | MamRep1894 | 5.28e-01 | AluYa5 | 9.25e-01 |
| 14 | HAL1ME | 5.67e-01 | AluYb8 | 9.25e-01 |
| 15 | MER41A | 6.67e-01 | AluYb9 | 9.25e-01 |
| 16 | L1MB5 | 7.07e-01 | AluYj4 | 9.25e-01 |
| 17 | MER61A | 7.07e-01 | AluYk11 | 9.25e-01 |
| 18 | MamTip2 | 7.07e-01 | AluYk3 | 9.25e-01 |
| 19 | LTR48B | 7.44e-01 | AluYk4 | 9.25e-01 |
| 20 | FLAM_C | 7.91e-01 | AmnSINE1 | 9.25e-01 |

Subfamilies in both top-20 lists: 0

### 5. isoform_TE_annotation.tsv categories (exon-overlap definition; identical in both runs)

| Category | (isoform, TE) pairs | % |
|---|---|---|
| autonomous | 9,761 | 3.8 |
| TSS | 52,862 | 20.7 |
| PAS | 132,392 | 51.8 |
| internal_exonization | 60,456 | 23.7 |
| other_exonic | 9 | 0.0 |
| **total** | 255,480 (149,263 isoforms) | 100 |
