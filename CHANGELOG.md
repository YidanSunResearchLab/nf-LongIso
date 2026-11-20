# Changelog

All notable changes to NextLongIso (nf-LongIso) are documented here.
Issue numbers refer to the v0.2.0 review list (issues 1-22).

## [0.2.0] - 2026-11-20

### Breaking changes

- **Sample and condition names** must match `^[A-Za-z][A-Za-z0-9_]*$` (letters, digits, underscores;
  start with a letter). Names with `.` or `-` were silently rewritten by R/DESeq2 in v0.1.0 and are now
  rejected at start-up. Missing conditions are an error (v0.1.0 silently used `control`).
- **The samplesheet is validated against a schema** (`assets/schema_input.json`): required columns
  `sample`, `condition` and `fastq` or `bam`; unique sample names; files must exist. All parameters are
  validated against `nextflow_schema.json`; unknown values (e.g. a misspelt `--input_type`) stop the run.
- **SAM files are no longer written** (`results/MINIMAP2/*.sam` is gone). Alignment streams directly into
  a sorted BAM; `results/minimap2/` holds BAM, BAI, `.stats` and `.flagstat`.
- **`DE_genes/` is now gene-level** (from `counts_gene.txt`); v0.1.0 wrote transcript-level results there.
  Transcript-level results moved to the new **`DE_transcripts/`**.
- **TE counting now defaults to the exon mode** (`--te_count_mode exon`). Use `--te_count_mode span` to
  reproduce v0.1.0 TE counts and TE differential expression.
- **Default resources changed**: requests are capped by `--max_cpus 16 --max_memory 128.GB --max_time 72.h`
  (v0.1.0 requested up to 64 CPUs for minimap2 and 200 GB for Bambu). Memory and time grow on retry.
- With more than two conditions, DESeq2 results are written per comparison
  (`DE_results_<condition>_vs_<reference>.csv`, `DE_TE_<type>_results_<condition>_vs_<reference>.csv`).
- `isoform_TE_exonization.tsv`: columns 3-5 are now named `exon_chr`, `exon_start`, `exon_end`
  (they always contained the exon coordinates; v0.1.0 labelled them `tx_*`).
- `APA_gene_summary.csv` and `clustered_base_level_PAS_usage.csv` have `condition_1`/`condition_2` columns
  (one block per comparison); the summary also reports the cluster ID/coordinate and p-value.
- Nextflow **≥ 25.10.0** is required (nf-schema plugin 2.7.3).

### Added

- Parameter and samplesheet validation with the nf-schema plugin; `--help`, `--helpFull`,
  `--help <parameter>`; parameter summary and a run-settings block in the log (issue 20).
- `MINIMAP2_INDEX` (built once; `--minimap2_index` to reuse an `.mmi`, `--minimap2_index_args`) and
  `JUNC_BED` (built once) (issue 13).
- `SAMTOOLS_STATS` (stats and flagstat) for both input types.
- BAM input: automatic indexing, coordinate-sort check, staging as `<sample>.sorted.bam`, and the same
  downstream steps as FASTQ input (statistics, junctions, bigWigs, TE junction validation) (issue 17).
- `--reference_condition` (default `control` if present, else the alphabetically first condition),
  used by DESeq2, TE DESeq2, IsoformSwitchAnalyzeR and proActiv; comparisons of every condition against
  the reference.
- Gene-level and transcript-level DESeq2 (`DE_genes/`, `DE_transcripts/`) (issue 6).
- Tool arguments: `--minimap2_preset`, `--minimap2_args`, `--minimap2_junc_bed`, `--bamcoverage_args`,
  `--bambu_args`, `--suppa2_events`, `--proactiv_species`, `--proactiv_altprom_args`,
  `--isoswitch_alpha`, `--isoswitch_dif_cutoff`, `--seed`; all can also be set with `ext.args` in a
  custom config.
- TE: `--te_count_mode exon|span`, `--te_min_overlap`, `--te_strand sense|both`, `--te_autonomous_frac`
  and the new `TE/isoform_TE_annotation.tsv` (one row per isoform-TE pair with subfamily/family/class,
  overlap, strand relation, category TSS/internal_exonization/PAS/autonomous/other_exonic and junction
  support) (issue 12).
- `versions.yml` for every process, merged into `pipeline_info/software_versions.yml` and shown in
  MultiQC (issue 19).
- Execution report, timeline, trace and DAG in `pipeline_info/` (timestamped) (issues 15, 16).
- Profiles `apptainer`, `slurm`, `lsf` (`--hpc_queue`, `--hpc_cluster_options`; LSF uses
  `perJobMemLimit`) and `test` (tiny dataset in `tests/data`, built by `tests/make_test_data.sh`)
  (issues 15, 16, 21).
- `stub:` blocks for every process; GitHub Actions CI (schema check, `--help`, lint, TE unit test,
  stub run; real test-profile run on manual dispatch) (issue 21).
- `tests/te_unit/`: synthetic TE test with hand-derived expected outputs;
  `tests/summarize_results.py` and `tests/compare_te_modes.py` for result comparisons.
- `docs/dev/`: baseline, verification of suspected issues, v0.2.0 vs v0.1.0 comparison.

### Changed

- Per-process `publishDir` instead of a global default (issue 1).
- Containers referenced without the `docker://` prefix; existing Singularity caches are reused (issue 2).
- Resource requests: base values per process, capped by `process.resourceLimits`; memory/time scale with
  `task.attempt`; retry on exit codes 104, 134, 137, 139, 140, 143, 247 (up to 2 retries) (issues 15, 16).
- Alignment uses the mulled minimap2 2.28 + samtools 1.20 image (same minimap2 version as v0.1.0) and a
  single pipe `zcat | awk | minimap2 | samtools sort` (issue 14).
- R code of DESeq2, TE DESeq2, IsoformSwitchAnalyzeR and proActiv moved to `bin/` scripts.
- MultiQC reads its inputs from channels (alignment statistics, DESeq2 summaries, software versions)
  instead of scanning the output directory (issue 5).
- `TE_ANALYSIS` rewritten in POSIX awk (no gawk extensions such as `ARGIND`; tested with BusyBox awk).
  In the shipped container `ARGIND` happened to work, but other awks (mawk, POSIX mode) silently
  returned an empty exonization table (issue 10, `docs/dev/verification.md` V1).
- Sample names are handled identically everywhere: only the exact `.sorted(.bam)` suffix is removed
  (issue 18).
- Naming unified as NextLongIso (log header, MultiQC title); `IsoformSwitchAnalyzeR` spelling (issue 22).
- Removed the ineffective `nxf_tmp` environment block (its temporary directory was never visible inside
  the containers) and its `mkdirs()` statement, which is invalid in Nextflow's strict syntax.

### Fixed

- **SUPPA2 PSI values were almost all `nan`**: transcript versions were removed from the expression file
  but not from the events, so only 2 of 691,201 demo events could be computed. IDs are now kept as in
  the Bambu GTF (issue 11, `docs/dev/verification.md` V2).
- SUPPA2 "PSI in ≥ 2 samples" filter ignored the last sample.
- DESeq2: names were cleaned with unescaped regular expressions (`gsub(".bam")`, `gsub("^X")`), which
  could alter valid sample names (issue 18).
- `deseq2_de_mqc.tsv` was declared but never produced; DESeq2 now writes a real MultiQC table (issue 7).
- IsoformSwitchAnalyzeR: samplesheet columns read by position; fewer than two conditions made the
  process fail; conditions without replicates made the DEXSeq test fail (issues 3, 4).
- APA: `isoformFeatures` has one row per comparison but the first matching row was used for every
  comparison; APA is now computed per comparison (issue 8).
- APA: `aggregate()` on empty tables failed; `APA_gene_summary.csv` aggregated each column separately
  and reported the site with the *smallest* |dPAU| per gene, without its coordinates; it now reports the
  significant site with the largest |dPAU| together with its cluster ID, coordinate and q-value, all from
  that site (on the demo data |dPAU| is now larger or equal for 450/450 genes) (issue 8).
- APA Fisher tests used unseeded simulated p-values (now `--seed`) (issue 8).
- `significant_switches.csv` could contain all-NA rows (filter without `which()`); it now also applies
  the |dIF| cutoff (issue 3).
- proActiv: species and alternative-promoter thresholds were hard-coded; samples were matched to
  conditions by prefix/substring (issue 9).
- TE: a transcript overlapping several TEs kept only the last one (`data[$4]` overwrite) in
  `isoform_TE_exonization.tsv`; junction validation counted junction hits per transcript (transcripts
  with two TE exons were dropped) and accepted exons that were only spanned by an exon-skipping junction
  (`docs/dev/verification.md` V3) (issue 12).
- TE GTF attributes are parsed by name instead of a greedy `sed` (issue 12).
- Junction extraction no longer uses gawk's `and()` (output unchanged) (issue 10).
- With BAM input, Bambu column names were the BAM file names instead of the sample names (issue 17).

## [0.1.0] - 2026-08-03

Version described in Tan J. et al. (2026), *Bioinformatics*.
