# Changelog

All notable changes to NextLongIso (nf-LongIso) are documented here.
Issue numbers refer to the v0.2.0 review list (issues 1-22). A1-A4, B1-B5 and C refer to the
methodological review of v0.2.0 (`docs/dev/review_v0.2.0_verification.md`); changes without a label
were found during the work on v0.2.0.

## [0.2.0] - 2026-11-20

### Breaking changes

- **Sample and condition names** must match `^[A-Za-z][A-Za-z0-9_]*$` (letters, digits, underscores;
  start with a letter). Names with `.` or `-` were silently rewritten by R/DESeq2 in v0.1.0 and are now
  rejected at start-up. A missing condition is an error (v0.1.0 silently used `control`) (issue 18).
- **The samplesheet is validated against a schema** (`assets/schema_input.json`): required columns
  `sample`, `condition` and `fastq` or `bam`; unique sample names; files must exist (issue 18). All
  parameters are validated against `nextflow_schema.json`; invalid values (e.g. a misspelt
  `--input_type`) stop the run (issue 20).
- **SAM files are no longer written** (`results/MINIMAP2/*.sam` is gone). Alignment streams directly into
  a sorted BAM; `results/minimap2/` holds BAM, BAI, `.stats` and `.flagstat` (issue 14).
- **`DE_genes/` is now gene level** (from `counts_gene.txt`); v0.1.0 wrote transcript-level results there.
  Transcript-level results moved to the new **`DE_transcripts/`** (issue 6).
- **TE counting now defaults to the exon mode** (`--te_count_mode exon`). `--te_count_mode span` keeps the
  v0.1.0 counting rule (overlap with the whole transcript span) as a comparison option. To reproduce the
  published v0.1.0 results exactly, run the pipeline with `-r v0.1.0` (issue 12).
- **`isoform_TE_exonization.tsv` uses a stricter junction test** (in both counting modes): each exon
  boundary must match an observed junction endpoint. Novel *terminal* exons (one splice site, e.g. a
  TE-derived first or last exon) are therefore no longer listed there; they remain visible in
  `isoform_TE_annotation.tsv` (categories `TSS`, `PAS` or `autonomous`). Columns 3-5 are named `exon_chr`,
  `exon_start`, `exon_end` (v0.1.0 labelled the same exon coordinates `tx_*`) (issue 12).
- **Default resources changed**: requests are capped by `--max_cpus 16 --max_memory 128.GB --max_time 72.h`
  (v0.1.0 requested up to 64 CPUs for minimap2 and 200 GB for Bambu) (issue 1); memory and time grow on
  retry (issue 15).
- With more than two conditions, DESeq2 results are written per comparison
  (`DE_results_<condition>_vs_<reference>.csv`, `DE_TE_<type>_results_<condition>_vs_<reference>.csv`)
  against `--reference_condition` (issue 7).
- **APA moved to `results/apa/`** and is now a sample-level DEXSeq test (A2): `isoswitch/APA_gene_summary.csv`
  and `isoswitch/clustered_base_level_PAS_usage.csv` are replaced by `apa/pas_transcripts.tsv`,
  `pas_annotation.tsv`, `pas_counts.tsv`, `apa_dexseq_results.tsv` and a new `apa/APA_gene_summary.csv`
  (gene q-value, top PAS, PAS shift and direction). PAS clusters are bounded (25 bp window, 50 bp maximal
  width), and novel 3' ends that look like internal priming are removed by default
  (`--apa_internal_priming filter`). `--seed` is removed (no simulated p-values any more).
- **Junction files changed** (A1): `junctions/<sample>.junctions.bed` now uses the transcript strand
  (minimap2 `ts` tag, or `XS`) instead of the read orientation, counts primary alignments only, and has the
  read count as score (capped at 1,000); a new `junctions/<sample>.junctions.tsv` adds motif and
  annotation status.
- Nextflow **≥ 25.10.0** is required (nf-schema plugin 2.7.3) (issue 20).
- The test data are no longer part of the repository; `-profile test` downloads them from Zenodo (issue 21).

### Added

- Parameter validation with the nf-schema plugin; `--help`, `--helpFull`, `--help <parameter>`;
  parameter summary in the log (issue 20).
- Samplesheet schema with explanatory error messages; warning for conditions with fewer than 2
  replicates (issue 18).
- `MINIMAP2_INDEX` (built once; `--minimap2_index` to reuse an `.mmi`, `--minimap2_index_args`) and
  `JUNC_BED` (built once) (issue 13).
- `SAMTOOLS_STATS` (stats and flagstat) for both input types.
- BAM input: automatic indexing, coordinate-sort check, staging as `<sample>.sorted.bam`, and the same
  downstream steps as FASTQ input (statistics, junctions, bigWigs, TE junction validation) (issue 17).
- `--reference_condition` (default `control` if present, else the alphabetically first condition), used
  by DESeq2, TE DESeq2, IsoformSwitchAnalyzeR and proActiv; every other condition is compared with the
  reference; the condition used is printed in the log (issue 7).
- Gene-level and transcript-level DESeq2 (`DE_genes/`, `DE_transcripts/`) (issue 6).
- Tool arguments: `--minimap2_preset`, `--minimap2_args`, `--minimap2_junc_bed`, `--bamcoverage_args`,
  `--bambu_args`, `--suppa2_events`, `--proactiv_altprom_args`, `--isoswitch_alpha`,
  `--isoswitch_dif_cutoff`; all can also be set with `ext.args` in a custom config.
- `--proactiv_species` (issue 9).
- TE: `--te_count_mode exon|span`, `--te_min_overlap`, `--te_strand sense|both`, `--te_autonomous_frac`
  and the new `TE/isoform_TE_annotation.tsv`: one row per isoform-TE pair with subfamily/family/class,
  overlap, strand relation, category (autonomous > TSS > PAS > internal_exonization > other_exonic),
  `exon_positions` (first/internal/last exons touched by the TE) and junction support (issue 12).
- `versions.yml` for every process, merged into `pipeline_info/software_versions.yml` and shown in
  MultiQC (issue 19).
- Profiles `slurm` and `lsf` (`--hpc_queue`, `--hpc_cluster_options`; LSF uses `perJobMemLimit`)
  (issue 16) and `apptainer`.
- Profile `test` (tiny dataset on Zenodo, reproducible with `tests/make_test_data.sh`), `stub:` blocks for
  every process, and GitHub Actions CI (schema check, `--help`, lint, TE unit test, stub run; real
  test-profile run on manual dispatch) (issue 21).
- Execution report, timeline, trace and DAG in `pipeline_info/` (timestamped).
- `tests/te_unit/` (synthetic TE test), `tests/summarize_results.py`, `tests/compare_te_modes.py`.
- `docs/dev/`: baseline, verification of suspected issues, v0.2.0 vs v0.1.0 comparison;
  `docs/te_count_mode_comparison.md`.
- **Experimental design with covariates** (B5): `--design` (e.g. `'~ batch + condition'`) over samplesheet
  columns, used by DESeq2 (genes, transcripts, TEs), IsoformSwitchAnalyzeR and the DEXSeq usage tests;
  the design matrix is checked before alignment (confounded designs stop with the columns involved).
- **Isoform switch consequences** (B3): transcript sequences from the Bambu annotation (`GFFREAD`,
  gffread 0.12.7); ORFs from the CDS of `--gtf` for annotated transcripts (`addORFfromGTF`) and predicted
  for novel transcripts (`analyzeNovelIsoformORF`; all predicted when `--gtf` has no CDS), with the new
  column `orf_source` (annotated/predicted) in the switch tables; `analyzeSwitchConsequences` (intron
  retention, ORF similarity, NMD, TSS, TTS, last exon, length, exon number) with `switch_consequences.csv`,
  consequence summary and enrichment tables/plots and `switch_plots.pdf` (`--isoswitch_consequences`,
  `--isoswitch_n_plots`).
- **SUPPA2 differential splicing** (B4): `SUPPA2_DIFFSPLICE` runs `diffSplice` (`--suppa2_diff_args`,
  default `-m empirical -gc`) for each condition with ≥ 2 replicates vs the reference, on the events of
  all types merged (as in the SUPPA2 workflow), with `suppa2/diffsplice/<condition>_vs_<reference>.dpsi/.psivec`,
  `suppa2_diffsplice_results.tsv` (all events with a genome-wide Benjamini-Hochberg `padj`; the `-gc`
  p-values are only corrected within genes) and `suppa2_diffsplice_significant.tsv` (padj <
  `--suppa2_alpha`, |dPSI| > `--suppa2_dpsi_cutoff`; `--suppa2_diff false` turns it off).
- **Read QC** (B1): NanoPlot 1.46.2 and `seqkit stats` 2.10.1 per sample for FASTQ input
  (`results/read_qc/`, parsed by MultiQC; `--skip_read_qc`).
- **SQANTI3 QC** (B2, optional): `--run_sqanti3` classifies the Bambu transcripts with `sqanti3_qc.py`
  (official image `anaconesalab/sqanti3:v6.0.1`); classification, junctions, HTML report and
  `sqanti3_category_counts.tsv` (FSM, ISM, NIC, NNC, ... for annotated and novel transcripts);
  `--sqanti3_args`.
- Alternative polyadenylation (`APA` process, `results/apa/`, A2): novel 3' ends (transcripts not in
  `--gtf`) with an A-rich downstream sequence are treated as internal priming unless an annotated end is
  within 50 bp or AATAAA/ATTAAA lies 10-35 nt upstream; parameters `--apa_*`, `--usage_alpha`.
- **Transcript-based promoter usage** (`--promoter_method tss`, optional; `PROMOTER_TSS`,
  `results/promoter_tss/`, A4): promoters are clusters of Bambu transcript 5' ends (single-exon transcripts
  included), with proActiv-style activity and classes (internal promoter = TSS in a non-first exon of
  another transcript of the gene) and a DEXSeq usage test with the design covariates; parameters
  `--tss_*`, `--promoter_min_activity`, `--promoter_min_count`. proActiv remains the default.
- Shared DEXSeq usage script `bin/usage_dexseq.R` (APA, `tss` promoters): each comparison is tested
  separately and skipped with a warning when it cannot be fitted.
- Junction TSV with motif and annotation, `--junction_min_mapq`, `SAMTOOLS_FAIDX` (genome index when no
  `.fai` exists) and `GTF_INTRONS` (A1).
- `tests/apa_unit/` (internal priming on both strands and clustering); `tests/summarize_results.py` reads
  the new `apa/` and `promoter_tss/` folders; CI stub run with `--run_sqanti3` and a design covariate, and
  real Docker runs of the test profile with BAM input, a single condition and three conditions
  (`tests/samplesheet_test_single.csv`, `tests/samplesheet_test_3cond.csv`; manual dispatch);
  `docs/dev/review_v0.2.0_verification.md`, `docs/dev/apa_method_comparison.md`,
  `docs/dev/promoter_method_comparison.md`, `docs/dev/v0.3_backlog.md`.
- README section "Known issues in v0.1.0": problems of v0.1.0 that affect its outputs, with the files
  concerned.

### Changed

- Resource requests: base values per process sized for a 16-core / 64 GB machine, capped by
  `process.resourceLimits` built from `--max_cpus/--max_memory/--max_time` (v0.1.0 defined
  `max_memory` but never used it) (issue 1).
- Memory and time requests scale with `task.attempt`; retry on exit codes 104, 130 (LSF memory-limit
  kill), 134, 137, 139, 140, 143, 247 (up to 2 retries) (issue 15).
- Containers referenced without the `docker://` prefix (issue 2). The pipeline no longer sets a
  Singularity/Apptainer cache directory: images go to `NXF_SINGULARITY_CACHEDIR` / `NXF_APPTAINER_CACHEDIR`
  when set (recommended on HPC; see README), otherwise to `<work directory>/singularity`.
- IsoformSwitchAnalyzeR imports the data with `detectUnwantedEffects = FALSE`, so the switch test uses the
  `--design` formula as given (by default it adds surrogate variables, which made its model differ from
  DESeq2 and the usage tests).
- Under `-profile docker`, `SQANTI3_QC` disables the image's `conda run` entrypoint
  (`--entrypoint ""`); the process activates the environment itself.
- Alignment is a single pipe `zcat | awk | minimap2 | samtools sort` (no SAM or truncated FASTQ on disk)
  in the mulled minimap2 2.28 + samtools 1.20 image (same minimap2 version as v0.1.0); per-process
  `publishDir` replaces the global default that published the SAM files (issue 14).
- MultiQC reads its inputs from channels (alignment statistics, DESeq2 summaries, software versions)
  instead of scanning the output directory (issue 5).
- R code of DESeq2, TE DESeq2, IsoformSwitchAnalyzeR and proActiv moved to `bin/` scripts.
- Inline R/awk code of all processes moved to `bin/` scripts with `--key=value` arguments parsed in base R
  (`bin/nli_design.R` holds the shared design helpers); outputs on the test data are unchanged (C).
- `TE_ANALYSIS` rewritten in POSIX awk (no gawk extensions such as `ARGIND`; tested with BusyBox awk). In
  the shipped container `ARGIND` happened to work, but other awks (mawk, POSIX mode) silently returned an
  empty exonization table. Junction extraction no longer uses gawk's `and()` (output unchanged)
  (issue 10, `docs/dev/verification.md` V1).
- Sample names are handled identically in every step: only the exact `.sorted(.bam)` suffix is removed
  (no unescaped `gsub`, `make.names` or prefix matching) (issue 18).
- Naming unified as NextLongIso; log header `NextLongIso v<version>` (v0.1.0 printed `IsoLong Pipeline
  (v1.0)`); completion message reports success or failure with the work directory; spelling
  `IsoformSwitchAnalyzeR`, `proActiv`, `deepTools` (issue 22).
- Removed the ineffective `nxf_tmp` environment block (its temporary directory was never visible inside
  the containers) and its `mkdirs()` statement, which is invalid in Nextflow's strict syntax.

### Fixed

- **SUPPA2 PSI values were almost all `nan`**: transcript versions were removed from the expression file
  but not from the events, so only 2 of 691,201 demo events could be computed. IDs are now kept as in the
  Bambu GTF; on the demo data every event with expression has a PSI value (issue 11,
  `docs/dev/verification.md` V2).
- IsoformSwitchAnalyzeR: a single condition made the process fail (required outputs missing), as did
  conditions without replicates (DEXSeq) and `aggregate()` on empty APA tables; placeholder outputs with a
  `Message` column are written instead, unreplicated conditions are excluded, and all summaries check for
  empty tables first (issue 3).
- IsoformSwitchAnalyzeR read sample and condition by column position (issue 4).
- `DE_genes/` contained transcript-level results (issue 6).
- No way to choose the reference condition; with more than two conditions only one comparison was
  reported (issue 7).
- APA: `APA_gene_summary.csv` aggregated each column separately and reported the site with the *smallest*
  |dPAU| per gene, without its coordinates; APA was computed for one comparison only (issue 8). The APA
  analysis was then replaced as a whole (A2, below).
- proActiv: species and alternative-promoter thresholds were hard-coded; samples were matched to
  conditions by prefix/substring (issue 9).
- TE: a transcript overlapping several TEs kept only the last one (`data[$4]` overwrite) in
  the exon-mode counts and in `isoform_TE_exonization.tsv`; junction validation counted junction hits per
  transcript (transcripts with two TE exons were dropped) and accepted exons only spanned by an
  exon-skipping junction (`docs/dev/verification.md` V3); TE GTF attributes are parsed by name instead of
  a greedy `sed` (issue 12).
- With BAM input, Bambu column names were the BAM file names instead of the sample names, and no
  junctions were extracted, so TE exonization could not be validated (issue 17).
- `deseq2_de_mqc.tsv` was declared but never produced; the declaration is removed and DESeq2 writes a real
  MultiQC table instead (issue 22).
- `significant_switches.csv` could contain all-NA rows (filter without `which()`); it now also applies
  the |dIF| cutoff.
- SUPPA2 "PSI in ≥ 2 samples" filter ignored the last sample.
- **Junction strand** (A1): minimap2 writes no `XS` tag, so junctions took the read orientation; on the
  demo data 2.0 % of spliced reads and up to 25 % of a sample's junctions had a wrong or duplicated
  strand (for unoriented cDNA reads this approaches 50 %). The strand now comes from the `ts` tag (or
  `XS`), only primary alignments are used, and the TE junction support only accepts junctions on the
  transcript's strand (or of unknown strand).
- **APA statistics** (A2): the APA test used condition-level pseudo-counts (IF × gene mean) and a Fisher
  test that ignored replicates, library size and the sample-level counts. APA now uses sample-level
  Bambu counts per PAS cluster and DEXSeq with the design covariates (see Breaking changes). On the demo
  data the old test called 1,041 of 1,564 clusters (450 genes) significant; the new one finds 3 genes,
  all among the old 450 (`docs/dev/apa_method_comparison.md`).
- **TE normalisation** (A3): TE DESeq2 estimated size factors from transcript and TE counts combined
  (TE counts are sums of transcript counts); size factors now come from the gene-level counts alone.
  On the demo data the H1 rep1 size factor changes from 1.31 to 1.49 and TE-family log2 fold changes
  shift by a median of −0.12.
- **proActiv counted junctions on the plus strand only** (also in v0.1.0): the junctions passed to
  proActiv had STAR strand code 0, which proActiv (`readSTARJunctions`) reads as `+`, so promoters of
  minus-strand genes had no counts (demo: 79,200 promoters). proActiv now reads the stranded junction
  tables (`+` = 1, `-` = 2; junctions of unknown strand are not used) with their real motif and annotation
  codes. On the test data, minus-strand promoters with counts go from 0 to 172 (active: 0 to 33/55 in the
  two conditions); plus-strand counts lose a few reads that came from antisense junctions and
  supplementary alignments.
- **Promoter quantification** (A4): proActiv infers promoter activity from first-intron junctions, so
  promoters of single-exon transcripts are not quantified (25,638 promoters in the demo annotation); the
  optional `tss` method covers them (see Added).
- DESeq2: the fallback used when no dispersion trend can be fitted failed ("first call
  estimateSizeFactors"); it now estimates size factors first.
- Designs without residual degrees of freedom (e.g. one sample per condition) made DESeq2 and TE DESeq2
  fail after alignment and Bambu; they now write a `Message` instead of results. A TE table without any
  read (e.g. chromosome names `1` vs `chr1`) no longer makes TE DESeq2 fail, and `TE_ANALYSIS` warns
  when the TE GTF and the isoform GTF share no chromosome name.
- DEXSeq usage tests (APA, `tss` promoters): one condition with a single replicate cancelled all
  comparisons, and a comparison whose design was not of full rank stopped the process; comparisons are
  now checked one by one.
- Junction motifs were matched to junctions by line position: a junction on a contig missing from the
  genome (or past its end) shifted the motifs of all following junctions. Sequences are now looked up by
  region, such junctions get motif `NA` with a warning, and junction extraction fails on a failing
  `samtools`/`awk` (pipefail).
- SUPPA2 diffSplice failed when an event type had no events; such types are left out.
- Novel transcripts were recognised by the `BambuTx` prefix, which changes when the Bambu prefix is set
  (e.g. `--bambu_args='opt.discovery = list(prefix = "Tx")'`); APA and `tss` promoters now use the
  transcript IDs of `--gtf`.
- Annotated introns (`annotated` column of the junction tables) are built per transcript and chromosome
  (a transcript ID reused on another chromosome produced a false intron), without zero-length introns.
- Design tables are read as text, so that values such as `T`, `F` or `NA` stay valid names.
- DESeq2 `results()` uses `alpha = 0.05`, the cutoff used for significance.
- `bambu/versions.yml` was published by two processes; only the declared result files are published.

## [0.1.0] - 2026-08-03

Version described in Tan J. et al. (2026), *Bioinformatics*.
