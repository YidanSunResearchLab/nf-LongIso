# Verification of suspected issues (before any code change)

Both checks were run on 2026-11-20 against v0.1.0 (`6fcc6f1`) and the v0.1.0 demo
baseline (`docs/dev/baseline_v0.1.0.md`).

## V1 (issue 10): `TE_ANALYSIS` uses `ARGIND`

`modules/te_analysis.nf` builds `isoform_TE_exonization.tsv` with

```awk
ARGIND==1 { ok[$1]=1; next }
ARGIND==2 { data[$4] = ...; next }
ARGIND==3 { if ($1 in data && $1 in ok) { ... } }
```

`ARGIND` is not POSIX. If an awk does not implement it, `ARGIND` is an
uninitialised variable (0), none of the three rules matches, and the step writes a
header-only file **without any error** (`set -euo pipefail` does not help).

### Which awk runs in the pipeline

Container `quay.io/biocontainers/bedtools:2.31.1--h13024bc_3`, tested with Apptainer 1.5.4:

| Check | Result |
|---|---|
| `command -v awk` | `/usr/bin/awk` → BusyBox multi-call binary |
| `awk --version` | `BusyBox v1.36.1 (2024-06-02 11:42:27 UTC)` |
| `awk '{print ARGIND, $1}' f1 f2 f3` | `1 a` / `2 x` / `3 z` (ARGIND supported) |
| v0.1.0 three-file pattern | routes file1 / file2 / file3 correctly |

### Other awk implementations (same minimal test, host system)

| awk | ARGIND pattern output |
|---|---|
| GNU awk 5.1.0 (`gawk`, also host `awk`/`nawk`) | 3 lines (works) |
| `gawk --posix` | **0 lines, no error** |
| `mawk` (default `awk` on Debian/Ubuntu) | **0 lines, no error** |

### Baseline

`TE/isoform_TE_exonization.tsv` has 2 data rows, which is consistent with `ARGIND`
working in the BusyBox container.

### Conclusion

**Not a bug in the shipped configuration.** With the pinned bedtools container, the
`ARGIND` code produces correct output. It is a **portability risk**: the same code
silently returns an empty table under mawk or a POSIX awk (for example a container
rebuilt on a different base image, or a run without containers). Rewriting with the
POSIX idioms `FNR==NR` / `FILENAME` removes the risk without changing results.

Independently of `ARGIND`, the same awk block has a real defect: `data[$4]` is keyed
by transcript ID, so when a transcript overlaps several TEs only the last TE is kept
(this is the "one TE per transcript" issue listed in phase 6).

## V2 (issue 11): SUPPA2 transcript-ID versions

`modules/suppa2.nf` strips the version from every transcript ID in the expression
file (`parts[0].split('.')[0]`), but `suppa.py generateEvents` reads the Bambu GTF,
which keeps versions. The `.ioe` events therefore reference IDs such as
`ENST00000704336.1`, while the expression file contains `ENST00000704336`.

### Measurements (baseline run)

| Item | ENST with version | ENST without version | BambuTx (no version) |
|---|---|---|---|
| `extended_annotations.gtf` transcript IDs (388,075) | 387,944 | 0 | 131 |
| `cleaned_counts.txt` (SUPPA2 expression input, 388,075 rows) | 0 | 387,944 | 131 |
| Transcript IDs referenced by `events_*_strict.ioe` (304,452) | 304,332 | 0 | 120 |

| Item | Value |
|---|---|
| `.ioe` transcript IDs present in the expression file | 120 of 304,452 (only the BambuTx IDs) |
| Events whose transcripts are all present in the expression file | 2 of 691,201 |
| PSI cells that are `nan` (all `psi/*_samples.psi`) | 2,764,799 of 2,764,804 = **100.000 %** |
| Events with at least one non-nan PSI | 2 (both AF events made only of matching IDs) |
| Events kept by the "≥ 2 non-nan" filter | 1 (AF) |

PSI event IDs keep the gene version (`ENSG00000000003.16;A5:chrX:...`), so the
mismatch is only between transcript IDs in the expression file and in the events.

### Conclusion

**Confirmed bug.** Stripping versions only from the expression file makes virtually
every SUPPA2 PSI value `nan`. The fix is to keep transcript IDs exactly as written in
the Bambu GTF (no version stripping), so that the expression file and the `.ioe`
files use identical IDs. After the fix, the PSI nan fraction should drop sharply; the
remaining nan values should come only from events whose transcripts have zero total
expression in a sample.

## V3 (found while fixing issue 10/12): junction validation of `isoform_TE_exonization.tsv`

The v0.1.0 exonization step validates a TE-derived novel exon by intersecting 10-bp windows
around both exon ends with the junction BEDs and then requiring that the transcript ID occurs
exactly twice in the combined "left ok" + "right ok" list. On the synthetic locus in
`tests/te_unit/` (bedtools container, v0.1.0 code rendered from the `v0.1.0` tag):

| Transcript | Situation | Expected | v0.1.0 | v0.2.0 |
|---|---|---|---|---|
| N1 | two novel intronic exons, each fully covered by a TE (AluY, L1M), both ends of both exons matched by observed junctions | 2 rows | **0 rows** (the ID occurs 4 times, not 2) | 2 rows |
| N3 | novel intronic exon fully covered by MER5A; the only junction skips the exon (spans it) | 0 rows | **1 row** (a spanning junction overlaps both windows) | 0 rows |

v0.2.0 validates each exon separately and requires an observed junction **endpoint** within
±5 bp of each exon boundary (a junction ending at the exon start and one starting at the exon
end), and keeps every (transcript, exon, TE) combination instead of one TE per transcript.

The same test confirms that `--te_count_mode span` reproduces the v0.1.0 family and instance
counts (including the v0.1.0 behaviour that a transcript overlapping several TEs contributes
only to the last one reported by bedtools). On the tiny test dataset, span-mode
`TE_family_counts.tsv` (116 rows) and `TE_instance_counts.tsv` (610 rows) are identical to the
v0.1.0 module's output.
