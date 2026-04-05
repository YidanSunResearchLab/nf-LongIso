#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ------------------------------------------------------------
// PARAMETERS (DEFAULTS)
// ------------------------------------------------------------
params.samplesheet = null
params.input_type  = 'fastq'   
params.genome      = null
params.gtf         = null
params.te_gtf      = null      
params.outdir      = "results"

// ------------------------------------------------------------
// BASIC VALIDATION
// ------------------------------------------------------------
if (!params.samplesheet)
    error "Please provide --samplesheet (CSV with sample,condition,fastq/bam)"

if (!params.gtf)
    error "Please provide --gtf (Reference annotation)"

if (params.input_type == 'fastq' && !params.genome)
    error "FASTQ input requires --genome (Fasta reference)"

log.info """
============================================================
             IsoLong Pipeline (v1.0)
============================================================
Input type  : ${params.input_type}
Sample Sheet: ${params.samplesheet}
Output Dir  : ${params.outdir}
TE Analysis : ${params.te_gtf ? 'Enabled' : 'Disabled'}
============================================================
"""

// ------------------------------------------------------------
// MODULE INCLUDES
// ------------------------------------------------------------
include { MINIMAP2; EXTRACT_JUNCTIONS; SAMTOOLS_SORT; MAKE_BIGWIG } from './modules/align.nf'
include { BAMBU_PREP; BAMBU }                                       from './modules/bambu.nf'
include { SUPPA2 }                                                  from './modules/suppa2.nf'
include { TE_ANALYSIS }                                             from './modules/te_analysis.nf'
include { DESEQ2; TE_DESEQ2 as TE_DESEQ2_FAMILY; TE_DESEQ2 as TE_DESEQ2_INSTANCE } from './modules/deseq2.nf'
include { ISOSWITCH }                                               from './modules/isoswitch.nf'
include { PROACTIV }                                                from './modules/proactiv.nf'
include { MULTIQC }                                                 from './modules/multiqc.nf'

// ------------------------------------------------------------
// WORKFLOW DEFINITION
// ------------------------------------------------------------
workflow {

    // 1. INPUT HANDLING
    ch_samplesheet = Channel.fromPath(params.samplesheet, checkIfExists: true)

    ch_input = ch_samplesheet
        .splitCsv(header:true)
        .map { row ->
            def meta = [ id: row.sample, condition: row.condition ?: 'control' ]
            def f = params.input_type == 'bam' ? row.bam : row.fastq
            if (!f) error "Column '${params.input_type}' not found in samplesheet for sample: ${row.sample}"
            tuple(meta, file(f))
        }

    // 2. ALIGNMENT STRATEGY
    if (params.input_type == 'fastq') {
        
        MINIMAP2(ch_input, file(params.genome), file(params.gtf))
        SAMTOOLS_SORT(MINIMAP2.out.sam)
        
        ch_bams_indexed = SAMTOOLS_SORT.out.bam_bai 

        EXTRACT_JUNCTIONS(ch_bams_indexed, file(params.gtf))
        MAKE_BIGWIG(ch_bams_indexed)

        ch_junctions = EXTRACT_JUNCTIONS.out.junctions

    } else {
        // Start from BAM
        ch_bams_indexed = ch_input.map { meta, bam ->
            def bai = file("${bam}.bai")
            tuple(meta, bam, bai)
        }
        ch_junctions = Channel.empty()
    }

    // 3. COHORT COLLECTION
    ch_bams_for_cohort = ch_bams_indexed.map { meta, bam, bai -> bam }.collect()
    ch_bais_for_cohort = ch_bams_indexed.map { meta, bam, bai -> bai }.collect()

    // 4. BAMBU ISOFORM DISCOVERY
    BAMBU_PREP(file(params.gtf))

    BAMBU(
        ch_bams_for_cohort,
        ch_bais_for_cohort,
        BAMBU_PREP.out.annot_rds,
        file(params.genome)
    )

    ch_extended_gtf         = BAMBU.out.extended_gtf
    ch_transcript_counts    = BAMBU.out.transcript_counts
    ch_transcript_abundance = BAMBU.out.transcript_abundance

    // 5. TE ANALYSIS (Optional)
    if (params.te_gtf) {
        
        ch_all_junctions = ch_junctions
            .map { it[1] }
            .collect()
            .ifEmpty { Channel.value([]) }

        TE_ANALYSIS(
            ch_extended_gtf,
            file(params.te_gtf),
            ch_transcript_counts,
            ch_all_junctions,
            file(params.gtf)
        )

        TE_DESEQ2_FAMILY(
            TE_ANALYSIS.out.te_counts, 
            ch_transcript_counts, 
            ch_samplesheet, 
            "Family"
        )
        
        TE_DESEQ2_INSTANCE(
            TE_ANALYSIS.out.te_instances, 
            ch_transcript_counts, 
            ch_samplesheet, 
            "Instance"
        )
    }

    // 6. CORE ISOFORM AND APA ANALYSES
    ISOSWITCH(
        ch_transcript_counts,
        ch_transcript_abundance,
        ch_extended_gtf,
        ch_samplesheet
    )

    DESEQ2(ch_transcript_counts, ch_samplesheet, ch_extended_gtf)
    
    SUPPA2(ch_extended_gtf, ch_transcript_counts)

    // 7. PROMOTER ACTIVITY (PROACTIV)
    PROACTIV(
        ch_extended_gtf,
        ch_bams_indexed.map { meta, bam, bai -> meta }.collect(),
        ch_bams_for_cohort,
        ch_bais_for_cohort
    )

    // 8. MULTIQC
    MULTIQC(
        Channel
            .fromPath("${params.outdir}/**/*")
            .filter { it.name =~ /(log|stats|metrics|summary)/ }
            .collect()
            .ifEmpty([])
    )
}

// ------------------------------------------------------------
// COMPLETION HANDLER
// ------------------------------------------------------------
workflow.onComplete {
    println """
    ============================================================
    IsoLong pipeline finished successfully
    Output directory: ${params.outdir}
    Duration: ${workflow.duration}
    ============================================================
    """
}