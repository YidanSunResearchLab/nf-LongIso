# Changelog

All notable changes to NextLongIso (nf-LongIso) are documented here.
Issue numbers refer to the v0.2.0 review list (issues 1-22); changes without a number were found
during the work on v0.2.0.

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
- **TE counting now defaults to the exon mode** (`--te_count_mode exon`). Use `--te_count_mode span` to
  reproduce all v0.1.0 TE outputs (issue 12).
- **In exon mode, `isoform_TE_exonization.tsv` uses a stricter junction test**: each exon boundary must
  match an observed junction endpoint. Novel *terminal* exons (one splice site, e.g. a TE-derived first or
  last exon) are therefore no longer listed there; they remain visible in `isoform_TE_annotation.tsv`
  (categories `TSS`, `PAS` or `autonomous`). Columns 3-5 are named `exon_chr`, `exon_start`, `exon_end`
  (v0.1.0 labelled the same exon coordinates `tx_*`; span mode keeps the v0.1.0 names) (issue 12).
- **Default resources changed**: requests are capped by `--max_cpus 16 --max_memory 128.GB --max_time 72.h`
  (v0.1.0 requested up to 64 CPUs for minimap2 and 200 GB for Bambu) (issue 1); memory and time grow on
  retry (issue 15).
- With more than two conditions, DESeq2 results are written per comparison
  (`DE_results_<condition>_vs_<reference>.csv`, `DE_TE_<type>_results_<condition>_vs_<reference>.csv`)
  against `--reference_condition` (issue 7).
- `APA_gene_summary.csv` and `clustered_base_level_PAS_usage.csv` have `condition_1`/`condition_2` columns
  (one block per comparison); the summary also reports the cluster ID/coordinate and p-value (issue 8).
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
- `--proactiv_species` (issue 9) and `--seed` for the APA Fisher tests (issue 8).
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

### Changed

- Resource requests: base values per process sized for a 16-core / 64 GB machine, capped by
  `process.resourceLimits` built from `--max_cpus/--max_memory/--max_time` (v0.1.0 defined
  `max_memory` but never used it) (issue 1).
- Memory and time requests scale with `task.attempt`; retry on exit codes 104, 134, 137, 139, 140, 143,
  247 (up to 2 retries) (issue 15).
- Containers referenced without the `docker://` prefix; existing Singularity caches are reused (issue 2).
- Alignment is a single pipe `zcat | awk | minimap2 | samtools sort` (no SAM or truncated FASTQ on disk)
  in the mulled minimap2 2.28 + samtools 1.20 image (same minimap2 version as v0.1.0); per-process
  `publishDir` replaces the global default that published the SAM files (issue 14).
- MultiQC reads its inputs from channels (alignment statistics, DESeq2 summaries, software versions)
  instead of scanning the output directory (issue 5).
- R code of DESeq2, TE DESeq2, IsoformSwitchAnalyzeR and proActiv moved to `bin/` scripts.
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
  |dPAU| per gene, without its coordinates; it now reports the significant site with the largest |dPAU|
  with its cluster ID, coordinate and q-value, all from that site (on the demo data |dPAU| is now larger
  or equal for 450/450 genes). APA is computed per comparison. Fisher tests with simulated p-values are
  seeded (`--seed`) (issue 8).
- proActiv: species and alternative-promoter thresholds were hard-coded; samples were matched to
  conditions by prefix/substring (issue 9).
- TE (exon mode): a transcript overlapping several TEs kept only the last one (`data[$4]` overwrite) in
  the counts and in `isoform_TE_exonization.tsv`; junction validation counted junction hits per
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

## [0.1.0] - 2026-08-03

Version described in Tan J. et al. (2026), *Bioinformatics*.
