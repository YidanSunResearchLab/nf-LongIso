# 🧬 NextLongIso: a comprehensive Nextflow pipeline for multi-dimensional long-read RNA-seq analysis

NextLongIso is a modular Nextflow pipeline for streamlined, reproducible, and scalable analysis of long-read RNA-seq data from PacBio Iso-Seq and Oxford Nanopore platforms. It integrates all major analysis steps,from read alignment and bigWig generation to isoform quantification, alternative splicing, alternative polyadenylation (APA), alternative promoter usage, and repeated sequences (simple repeats and TE) analysis, by using state-of-the-art tools and custom scripts.

🧭 **Workflow Overview**

Below is the schematic overview of the NextLongIso pipeline:

![Workflow Overview](data/figure.png)

Figure 1. Overview of the NextLongIso pipeline, showing read alignment, isoform quantification, alternative promoter usage, APA, alternative splicing, and repeats (simple repeats and TE) exonization/differential expression analysis.

✳️ **Key Features**

* **End-to-end automation** of long-read RNA-seq analysis from FASTQ to results.
* **Supports multiple platforms:** Fully compatible with both **PacBio** and **Nanopore** long-read data.
* **Integrated analysis tools:** Utilizes minimap2, samtools, deeptools, Bambu, DESeq2, IsoformSwitchAnalyzer, SUPPA2, ProActiv, bedtools and custom scripts.
* **Containerized environment:** Uses Singularity/Apptainer with pre-built images from Docker Hub.
* **Fully modular, reproducible, and scalable**—features overridable reference paths, making it ideal for large datasets and HPC environments.

🧩 **Installation**

To run NextLongIso, you need **Nextflow (≥ 25.10.0)** and **Apptainer** or **Singularity**. Internet access is required on the first run to pull containers (~several GB total).

**1. Environment Setup**: Use Conda to install the dependencies required to run the pipeline:

```bash
# Create and activate environment
conda create -n nf-LongIso -c conda-forge -c bioconda \
    nextflow=25.10.0 \
    apptainer \
    openjdk=17 \
    graphviz \
    -y
```

**2. Pipeline Installation**

```bash
conda activate nf-LongIso

git clone https://github.com/YidanSunResearchLab/nf-LongIso.git

cd nf-LongIso
```

🚀 **Quick Start for demo test**

**1. Prepare Reference Files (GRCh38/hg38 as example)**
```bash
mkdir -p data/genome
cd data/genome

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz

gunzip GRCh38.primary_assembly.genome.fa.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz

gunzip gencode.v47.annotation.gtf.gz

wget https://zenodo.org/records/21074209/files/TEtranscript_hg38_rmsk_TE.gtf?download=1
```
You will need the Genome FASTA, Gene annotation GTF, and TE annotation GTF. Ensure your fasta and gtf are from the same assembly.
* **Genome FASTA:** Download from Ensembl or UCSC (primary assembly).
* **Gene GTF:** Download from Gencode (v47 or newer).
* **TE GTF:** Download from the RepeatMasker (https://www.repeatmasker.org/) or generate from the UCSC Table Browser (https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=4105313869_RcAnJvPCP00P9TQImIQTYZR0qr7U&db=hg38&hgta_group=rep&hgta_track=knownGene&hgta_table=0&hgta_regionType=genome&position=chr7%3A155%2C799%2C529-155%2C812%2C871&hgta_outputType=gff&hgta_outFileName=). Ensure the assembly version matches your Gene GTF and that chromosome nomenclature (e.g., chr1 vs 1) is consistent across all files.

*Important Note: Inconsistent chromosome naming (e.g., "chr1" in one file and "1" in another) is a common cause of pipeline failure. Always verify the first few lines of your files to ensure they match before starting.*

**2. Start**

Navigate to the project directory and execute the pipeline using the command below. Note that the file paths are relative to the project root; ensure your reference files are correctly placed in the **data/genome/** directory.

```bash
cd nf-LongIso

nextflow run main.nf \
  -profile singularity \
  --input_type fastq \
  -with-report resource_report.html \
  --genome data/genome/GRCh38.primary_assembly.genome.fa \
  --gtf data/genome/gencode.v47.annotation.gtf \
  --te_gtf data/genome/TEtranscript_hg38_rmsk_TE.gtf \
  --samplesheet demo/samplesheet.csv \
  -resume
```

Understanding the Flags:
**-profile singularity**: Uses Singularity/Apptainer containers to ensure software reproducibility \
**--input_type**: Defines the data format (fastq) \ 
**-with-report**: Generates an HTML report showing execution metrics (CPU/RAM) \
**-resume**: This allows the pipeline to pick up exactly where it left off if a run is interrupted \
**--genome**, **--gtf**, **--te_gtf**: Paths to the required reference files \ 
**--samplesheet**: Path to the input CSV file that specifies the location of your input files (e.g., fastq) \

**Note**: During the first run, the pipeline will download required containers, which may take some time depending on your internet connection. We recommend using the -resume flag for all future executions, it intelligently tracks progress and reuses successful results, drastically reducing runtime.

**3. Output**

After a successful run, NextLongIso generates the following output directories and representative files.

| Output folder  | Output files                                                                                                                                                                                          | Description                                                                                                                                                         |
| -------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **bambu/**     | `extended_annotations.gtf`, `counts_gene.txt`, `counts_transcript.txt`, `CPM_transcript.txt`, `bambu_annotations.rds`, `bambu_results.rds`                                                            | Transcript annotation, gene/transcript expression quantification, and Bambu intermediate objects.                                                                   |
| **MINIMAP2/**  | `*.sam`                                                                                                                                                                                               | Raw read alignments generated by Minimap2.                                                                                                                          |
| **minimap2/**  | `*.sorted.bam`, `*.sorted.bam.bai`, `*.sorted.bam.stats`                                                                                                                                              | Sorted and indexed BAM files together with alignment statistics.                                                                                                    |
| **bigwig/**    | `*.bw`                                                                                                                                                                                                | BigWig coverage tracks for genome browser visualization (e.g., IGV, UCSC Genome Browser).                                                                           |
| **junctions/** | `*.junctions.bed`                                                                                                                                                                                     | Splice junction BED files for visualization and downstream analyses.                                                                                                |
| **DE_genes/**  | `DE_results.csv`, `normalized_counts.csv`, `pca_plot.pdf`                                                                                                                                             | Differential gene expression results, normalized expression matrix, and PCA plot.                                                                                   |
| **isoswitch/** | `significant_switches.csv`, `isoform_switch_results.csv`, `atts_results.csv`, `atss_results.csv`, `all_AS_events.csv`, `APA_gene_summary.csv`, `clustered_base_level_PAS_usage.csv`, `switchlist.rds` | Isoform switch analysis, alternative transcription start/termination sites (ATSS/ATTS), alternative splicing events, and alternative polyadenylation (APA) results. |
| **suppa2/**    | `events_A3/A5/AF/AL/MX/RI/SE_*.ioe`, `psi/`                                                                                                                                                           | Alternative splicing event definitions and PSI (Percent Spliced In) matrices generated by SUPPA2.                                                                   |
| **proactiv/**  | `alternative_promoters.tsv`, `promoter_activity.tsv`, `promoter_counts.tsv`                                                                                                                           | Alternative promoter identification, promoter activity estimation, and promoter-level quantification.                                                               |
| **TE/**        | `isoform_TE_exonization.tsv`, `TE_family_counts.tsv`, `TE_instance_counts.tsv`                                                                                                                        | Transposable element exonization events and TE expression quantification.                                                                                           |
| **TE/DE/**     | `DE_TE_Family_results.csv`, `DE_TE_Instance_results.csv`                                                                                                                               | Differential expression analysis of TE families and individual TE instances.                                                                                        |
| **multiqc/**   | `multiqc_report.html`                                                                                                                                                                                 | Integrated quality control report summarizing the entire workflow.                                                                                                  |

**Note**:Some output files (e.g., isoform_TE_exonization.tsv) may not be generated if the corresponding biological events are not detected in your dataset.

💡 **Troubleshooting**

* **Pull fails?** Check your internet connection and ensure your Singularity/Apptainer version is up to date (≥ 3.8).
* **Reference mismatch?** Ensure your FASTA and GTF are from the exact same assembly build.

🛠 **Maintenance and Support**

NextLongIso is actively maintained by the Sun Lab. We are committed to maintaining and improving the software.
Maintenance activities include:
* Updating software dependencies
* Fixing reported bugs
* Improving compatibility with new systems and sequencing technologies
* Implementing usability improvements when appropriate

We also welcome community contributions. Users are encouraged to:
* Report issues via the GitHub **Issues** page
* Suggest new features
* Contribute improvements through **Pull Requests**

While we aim to maintain long-term support for NextLongIso, development priorities may evolve as the software and research needs grow.

🧾 **Citation**

If you use NextLongIso in your work, please cite:
* **This repository:** [https://github.com/YidanSunResearchLab/nf-LongIso.git](https://github.com/YidanSunResearchLab/nf-LongIso.git)
* **The individual tools used within the pipeline:** minimap2, samtools, deeptools, Bambu, SUPPA2, IsoformSwitchAnalyzer, ProActiv, DESeq2, bedtools, etc.

⚖️ **License**

MIT LICENSE.
