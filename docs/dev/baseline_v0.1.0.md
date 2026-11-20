# v0.1.0 baseline (demo data)

This file records the outputs of the v0.1.0 pipeline on the bundled demo data
(`demo/samplesheet.csv`: GM12878 ×2 = `control`, H1 ×2 = `treatment`, 10,000 reads each).
It is the reference point for the v0.2.0 comparison (`docs/dev/v0.2.0_vs_v0.1.0.md`).

## Provenance

| Item | Value |
|---|---|
| Results directory | `/mnt/citadel2/research/syidan/Tools/nf-LongIso_v2_runs/baseline_v0.1.0_demo` |
| Code | `v0.1.0` (`6fcc6f1`) plus an argument-passing layer (`13a62db`, branch `feature/custom-tool-args`) whose default commands were verified to be identical to v0.1.0 (see below) |
| Run date | 2026-11-20, Nextflow 25.10.0, Apptainer 1.5.4, `-profile singularity`, local executor |
| References | `GRCh38.genome.fa`, `gencode.annotation.gtf`, `TEtranscript_hg38_rmsk_TE.gtf` |

The run was cross-checked against an independent v0.1.0 run of the same reads
(`scripts_Jiang/IsoLong/demo_results`, June 2026, sample names `GM12878_SRR10838647` etc.).
After mapping sample names, alignments (SAM records, samtools stats), junction BEDs and
bigWigs were byte-identical; Bambu output differed only in the numbering of the 131 novel
`BambuTx` IDs (same exon structures, genes and counts); SUPPA2 events differed only in the
order of transcript IDs inside each list; DESeq2 statistics agreed to ~5e-13 relative
error. Every metric below is identical between the two runs.

## Key metrics

| Metric | Value |
|---|---|
| minimap2: GM12878_demo_rep1 reads total / mapped | 10000 / 10000 |
| minimap2: GM12878_demo_rep2 reads total / mapped | 10000 / 10000 |
| minimap2: H1_demo_rep1 reads total / mapped | 10000 / 10000 |
| minimap2: H1_demo_rep2 reads total / mapped | 10000 / 10000 |
| minimap2: SAM files published (MINIMAP2/*.sam) | 4 |
| bigwig files | 4 |
| junction BED files | 4 |
| bambu: transcripts | 388075 |
| bambu: novel transcripts (BambuTx*) | 131 |
| bambu: genes | 78936 |
| DE_genes/DE_results.csv: rows | 388075 |
| DE_genes/DE_results.csv: padj < 0.05 | 0 |
| isoswitch: isoformFeatures rows | 2270 |
| isoswitch: q < 0.05 (any dIF) | 0 |
| isoswitch: q < 0.05 & |dIF| > 0.1 | 0 |
| isoswitch: significant_switches.csv rows | 0 |
| isoswitch: significant_switches.csv all-NA rows | 0 |
| APA: PAS clusters | 1564 |
| APA: PAS clusters q < 0.05 | 1041 |
| APA: APA_gene_summary.csv genes | 450 |
| isoswitch: atss_results.csv rows | 1317 |
| isoswitch: atts_results.csv rows | 1130 |
| isoswitch: all_AS_events.csv rows | 2270 |
| proactiv: promoters | 160497 |
| proactiv: alternative promoters | - (not written) |
| suppa2 A3: events / PSI nan % / events kept after filter | 42591 / 100.0% / 0 |
| suppa2 A5: events / PSI nan % / events kept after filter | 39553 / 100.0% / 0 |
| suppa2 AF: events / PSI nan % / events kept after filter | 309402 / 100.0% / 1 |
| suppa2 AL: events / PSI nan % / events kept after filter | 166493 / 100.0% / 0 |
| suppa2 MX: events / PSI nan % / events kept after filter | 24423 / 100.0% / 0 |
| suppa2 RI: events / PSI nan % / events kept after filter | 15134 / 100.0% / 0 |
| suppa2 SE: events / PSI nan % / events kept after filter | 93605 / 100.0% / 0 |
| suppa2: overall PSI nan % | 100.0% |
| TE: TE_family_counts.tsv data rows | 1080 |
| TE: TE_instance_counts.tsv data rows | 105548 |
| TE: isoform_TE_exonization.tsv data rows | 2 |
| TE: isoform_TE_annotation.tsv data rows | - |
| TE DE Family: padj < 0.05 | 0 |
| TE DE Instance: padj < 0.05 | 1 |

### Files

| File | Lines / size |
|---|---|
| `DE_genes/DE_results.csv` | 388076 lines |
| `DE_genes/normalized_counts.csv` | 388076 lines |
| `DE_genes/pca_plot.pdf` | 4957 bytes |
| `MINIMAP2/GM12878_demo_rep1.sam` | 35192052 bytes |
| `MINIMAP2/GM12878_demo_rep2.sam` | 37187089 bytes |
| `MINIMAP2/H1_demo_rep1.sam` | 53864415 bytes |
| `MINIMAP2/H1_demo_rep2.sam` | 59730911 bytes |
| `TE/DE/DE_TE_Family_results.csv` | 1081 lines |
| `TE/DE/DE_TE_Instance_results.csv` | 105549 lines |
| `TE/DE/norm_counts_Family.csv` | 1081 lines |
| `TE/DE/norm_counts_Instance.csv` | 105549 lines |
| `TE/DE/pca_Family.pdf` | 5079 bytes |
| `TE/DE/pca_Instance.pdf` | 5098 bytes |
| `TE/DE/volcano_Family.pdf` | 24081 bytes |
| `TE/DE/volcano_Instance.pdf` | 108858 bytes |
| `TE/TE_family_counts.tsv` | 1081 lines |
| `TE/TE_instance_counts.tsv` | 105549 lines |
| `TE/isoform_TE_exonization.tsv` | 3 lines |
| `bambu/CPM_transcript.txt` | 388076 lines |
| `bambu/bambu_annotations.rds` | 17016736 bytes |
| `bambu/bambu_results.rds` | 17576614 bytes |
| `bambu/counts_gene.txt` | 78937 lines |
| `bambu/counts_transcript.txt` | 388076 lines |
| `bambu/extended_annotations.gtf` | 2553913 lines |
| `bigwig/GM12878_demo_rep1.bw` | 435270 bytes |
| `bigwig/GM12878_demo_rep2.bw` | 504362 bytes |
| `bigwig/H1_demo_rep1.bw` | 430122 bytes |
| `bigwig/H1_demo_rep2.bw` | 654955 bytes |
| `isoswitch/APA_gene_summary.csv` | 451 lines |
| `isoswitch/all_AS_events.csv` | 2271 lines |
| `isoswitch/atss_results.csv` | 1318 lines |
| `isoswitch/atts_results.csv` | 1131 lines |
| `isoswitch/clustered_base_level_PAS_usage.csv` | 1565 lines |
| `isoswitch/isoform_switch_results.csv` | 2271 lines |
| `isoswitch/significant_switches.csv` | 1 lines |
| `isoswitch/switchlist.rds` | 551682 bytes |
| `junctions/GM12878_demo_rep1.junctions.bed` | 16050 lines |
| `junctions/GM12878_demo_rep2.junctions.bed` | 17540 lines |
| `junctions/H1_demo_rep1.junctions.bed` | 16059 lines |
| `junctions/H1_demo_rep2.junctions.bed` | 24559 lines |
| `minimap2/GM12878_demo_rep1.sorted.bam.bai` | 982664 bytes |
| `minimap2/GM12878_demo_rep1.sorted.bam.stats` | 55325 lines |
| `minimap2/GM12878_demo_rep1.sorted.bam` | 4982054 bytes |
| `minimap2/GM12878_demo_rep2.sorted.bam.bai` | 1121424 bytes |
| `minimap2/GM12878_demo_rep2.sorted.bam.stats` | 79414 lines |
| `minimap2/GM12878_demo_rep2.sorted.bam` | 5341450 bytes |
| `minimap2/H1_demo_rep1.sorted.bam.bai` | 546376 bytes |
| `minimap2/H1_demo_rep1.sorted.bam.stats` | 63543 lines |
| `minimap2/H1_demo_rep1.sorted.bam` | 5379852 bytes |
| `minimap2/H1_demo_rep2.sorted.bam.bai` | 975144 bytes |
| `minimap2/H1_demo_rep2.sorted.bam.stats` | 66247 lines |
| `minimap2/H1_demo_rep2.sorted.bam` | 6620038 bytes |
| `proactiv/promoter_activity.tsv` | 160498 lines |
| `proactiv/promoter_counts.tsv` | 160498 lines |
| `resource_report_demo.html` | 3054364 bytes |
| `suppa2/events_A3_strict.ioe` | 42592 lines |
| `suppa2/events_A5_strict.ioe` | 39554 lines |
| `suppa2/events_AF_strict.ioe` | 309403 lines |
| `suppa2/events_AL_strict.ioe` | 166494 lines |
| `suppa2/events_MX_strict.ioe` | 24424 lines |
| `suppa2/events_RI_strict.ioe` | 15135 lines |
| `suppa2/events_SE_strict.ioe` | 93606 lines |
| `suppa2/psi/events_A3_samples.psi` | 42592 lines |
| `suppa2/psi/events_A3_samples_filtered.psi` | 1 lines |
| `suppa2/psi/events_A5_samples.psi` | 39554 lines |
| `suppa2/psi/events_A5_samples_filtered.psi` | 1 lines |
| `suppa2/psi/events_AF_samples.psi` | 309403 lines |
| `suppa2/psi/events_AF_samples_filtered.psi` | 2 lines |
| `suppa2/psi/events_AL_samples.psi` | 166494 lines |
| `suppa2/psi/events_AL_samples_filtered.psi` | 1 lines |
| `suppa2/psi/events_MX_samples.psi` | 24424 lines |
| `suppa2/psi/events_MX_samples_filtered.psi` | 1 lines |
| `suppa2/psi/events_RI_samples.psi` | 15135 lines |
| `suppa2/psi/events_RI_samples_filtered.psi` | 1 lines |
| `suppa2/psi/events_SE_samples.psi` | 93606 lines |
| `suppa2/psi/events_SE_samples_filtered.psi` | 1 lines |
