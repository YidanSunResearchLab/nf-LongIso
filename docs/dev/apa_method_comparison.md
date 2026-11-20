# APA: old method (v0.2.0 draft, Fisher) vs new method (DEXSeq)

Status: method description complete; demo numbers are filled in after the demo run of review
phase 12 (`demo_v0.2.0_r2`), compared with `demo_v0.2.0` (old method).

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

To be completed after the demo run (review phase 12): PAS clusters, genes tested, significant genes
(old: PAS-cluster q < 0.05 in `clustered_base_level_PAS_usage.csv` / genes in `APA_gene_summary.csv`;
new: gene q < 0.05), overlap, and the number of 3' ends removed by the internal-priming filter.

Pre-run figures on the demo Bambu annotation (`demo_v0.2.0/bambu`, `apa_sites.R` only): 388,075
transcript 3' ends (131 novel); 1 novel end with an A-rich downstream sequence, rescued by an
exception, so 0 ends were flagged; 117 of the 131 novel ends carry a poly(A) signal 10-40 nt upstream;
234,252 PAS clusters in 78,936 genes, 29,941 genes with ≥ 2 clusters before count filtering.
