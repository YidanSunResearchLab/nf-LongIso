# APA: old method (v0.2.0 draft, Fisher) vs new method (DEXSeq)

Demo comparison: `demo_v0.2.0` (old method) vs `demo_v0.2.0_r2` (new method), both on the bundled demo
(GM12878 ×2 = control, H1 ×2 = treatment, 10,000 reads per sample, GRCh38 / GENCODE).

## Methods

| | Old (until review phase 6; inside ISOSWITCH) | New (APA process, `results/apa/`) |
|---|---|---|
| PAS of a transcript | 3' end of its last exon | same (strand-aware) |
| Clustering | ±12.5 bp resize + `reduce()` over all transcripts (no width limit; single linkage can chain sites) | per gene and strand, window `--apa_cluster_window` (25 bp), maximal width `--apa_max_cluster_width` (50 bp) |
| Internal priming | not considered | novel 3' ends with ≥ `--apa_ip_a_frac` (0.7) A or a run of ≥ 6 A in the 20 nt downstream, unless an annotated 3' end is ≤ 50 bp away or a poly(A) signal lies 10-40 nt upstream; `--apa_internal_priming filter` (default) or `flag` |
| Counts | condition-level pseudo-counts IF × gene_value (isoformFeatures) | sample-level Bambu transcript counts summed per PAS cluster |
| Test | Fisher's exact test with simulated p-values on the rounded pseudo-counts, BH per comparison | DEXSeq per comparison (`~ sample + exon + <covariates>:exon + condition:exon`), PAS-level padj and gene-level q (perGeneQValue) |
| Replicates | not used | required (≥ 2 per condition), dispersion estimated from them |
| Covariates | no | `--design` |
| Filters | none | PAS with < `--apa_min_count` (10) reads removed; genes need ≥ 2 PAS |
| Gene summary | site with the smallest \|dPAU\| | significant site with the largest \|dPAU\|, gene q, PAU-weighted shift of the PAS position (distal = longer 3' ends) |

## Synthetic checks (tests/apa_unit)

Internal-priming rules on a synthetic contig: an A-rich novel end without exception is flagged and
removed (`filter`) or marked (`flag`); A-rich novel ends with a poly(A) signal 25 nt upstream or an
annotated end 30 bp away are kept; annotated ends are never flagged.

DEXSeq usage test (`bin/usage_dexseq.R`) on 60 simulated genes with 2 sites, 3 vs 3 replicates and a
usage shift (0.5 -> 0.85) in 6 genes: exactly the 6 shifted genes have q < 0.05, with and without a
batch covariate.

## Demo comparison

| | Old (`demo_v0.2.0/isoswitch/`) | New (`demo_v0.2.0_r2/apa/`) |
|---|---|---|
| Transcript 3' ends | from switch-analysis isoforms only | 388,075 (131 novel) |
| Ends flagged as internal priming | - | 0 (1 A-rich novel end, rescued by a PAS signal / nearby annotated end) |
| PAS clusters | 1,564 in 874 genes (expressed isoforms of tested genes) | 234,252 in 78,936 genes |
| Tested | 1,564 clusters | 75 PAS in 37 genes (≥ 10 reads, ≥ 2 PAS per gene) |
| Significant | 1,041 clusters (q < 0.05) in 450 genes | 6 PAS (padj < 0.05); **3 genes** (gene q < 0.05) |
| Gene summary | 450 genes | 37 genes, 3 significant (2 distal, 1 proximal shift) |

The three significant genes (ENSG00000188846 RPL14, PAS shift +5.3 kb distal, dPAU −0.85;
ENSG00000132780 NASP, +448 bp distal; ENSG00000164587 RPS14, −305 bp proximal) are all among the 450
old genes. 28 of the 450 old genes reach the new test (the others have < 10 reads per PAS or a
single PAS with reads); none of the other 25 is significant.

The old method called 58 % of its tested clusters significant with 2 × 10,000 reads per condition: the
Fisher test treated condition-mean pseudo-counts as independent reads, so any usage difference in a
moderately expressed gene became "significant" regardless of replicate agreement. With DEXSeq the
replicate variability enters the test, and only consistent shifts remain. The demo is too shallow
for the new test to have much power; the numbers show the change in calibration, not sensitivity.

Pre-run figures on the demo Bambu annotation (`demo_v0.2.0/bambu`, `apa_sites.R` only): 388,075
transcript 3' ends (131 novel); 1 novel end with an A-rich downstream sequence, rescued by an
exception, so 0 ends were flagged; 117 of the 131 novel ends carry a poly(A) signal 10-40 nt upstream;
234,252 PAS clusters in 78,936 genes, 29,941 genes with ≥ 2 clusters before count filtering.
