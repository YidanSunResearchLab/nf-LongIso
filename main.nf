#!/usr/bin/env nextflow
/*
 * NextLongIso (nf-LongIso): multi-dimensional long-read RNA-seq analysis
 * https://github.com/YidanSunResearchLab/nf-LongIso
 *
 * All parameter defaults live in nextflow.config; nextflow_schema.json documents and validates them.
 */

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

include { MINIMAP2_INDEX; JUNC_BED; MINIMAP2_ALIGN; SAMTOOLS_INDEX;
          SAMTOOLS_STATS; EXTRACT_JUNCTIONS; MAKE_BIGWIG }                    from './modules/align.nf'
include { BAMBU_PREP; BAMBU }                                                  from './modules/bambu.nf'
include { SUPPA2 }                                                             from './modules/suppa2.nf'
include { TE_ANALYSIS }                                                        from './modules/te_analysis.nf'
include { DESEQ2 as DESEQ2_GENE; DESEQ2 as DESEQ2_TRANSCRIPT;
          TE_DESEQ2 as TE_DESEQ2_FAMILY; TE_DESEQ2 as TE_DESEQ2_INSTANCE }   from './modules/deseq2.nf'
include { ISOSWITCH }                                                          from './modules/isoswitch.nf'
include { PROACTIV }                                                           from './modules/proactiv.nf'
include { MULTIQC }                                                            from './modules/multiqc.nf'

/*
 * Read and check the samplesheet. Returns a list of [meta, file] where file is the
 * FASTQ or BAM column, depending on --input_type.
 */
def readSamplesheet(samplesheet, input_type) {
    def rows = samplesheetToList(samplesheet, "${projectDir}/assets/schema_input.json")
    def samples = rows.collect { row ->
        def meta = [ id: row[0].id.toString(), condition: row[0].condition.toString() ]
        def f = input_type == 'bam' ? row[2] : row[1]
        if (!f) {
            error "Sample '${meta.id}' has no '${input_type}' value, but --input_type ${input_type} was given. " +
                  "Expected samplesheet columns: sample,${input_type},condition"
        }
        [ meta, f ]
    }
    samples.groupBy { s -> s[0].condition }.each { cond, s ->
        if (s.size() < 2) {
            log.warn "Condition '${cond}' has only ${s.size()} replicate(s). Differential analyses need at least 2 replicates per condition for reliable results."
        }
    }
    return samples
}

/*
 * Reference condition for all comparisons: --reference_condition if given, otherwise
 * 'control' if present, otherwise the alphabetically first condition.
 */
def resolveReference(samples, requested) {
    def conditions = samples.collect { s -> s[0].condition }.unique().sort()
    if (requested) {
        if (!conditions.contains(requested)) {
            error "--reference_condition '${requested}' is not a condition in the samplesheet (conditions: ${conditions.join(', ')})"
        }
        return requested
    }
    return conditions.contains('control') ? 'control' : conditions[0]
}

workflow {

    validateParameters()
    log.info "NextLongIso v${workflow.manifest.version}"
    log.info paramsSummaryLog(workflow)

    // 1. INPUT
    def samples   = readSamplesheet(params.samplesheet, params.input_type)
    def reference = resolveReference(samples, params.reference_condition)
    log.info "Reference condition: ${reference}" +
             (params.reference_condition ? '' : ' (not set with --reference_condition; chosen automatically)')
    ch_input = channel.fromList(samples)

    // Normalised design table (sample,condition) used by all R steps
    ch_design = channel.fromList(samples)
        .map { meta, _f -> "${meta.id},${meta.condition}" }
        .collectFile(name: 'design.csv', seed: 'sample,condition', newLine: true, sort: false)
        .first()

    // 2. ALIGNMENT (FASTQ) OR BAM PREPARATION
    if (params.input_type == 'fastq') {
        ch_index = params.minimap2_index
            ? channel.value(file(params.minimap2_index, checkIfExists: true))
            : MINIMAP2_INDEX(file(params.genome)).index
        ch_junc  = params.minimap2_junc_bed
            ? JUNC_BED(file(params.gtf)).bed
            : channel.value([])
        MINIMAP2_ALIGN(ch_input, ch_index, ch_junc)
        ch_bams_indexed = MINIMAP2_ALIGN.out.bam_bai
    } else {
        ch_bam_in = ch_input.map { meta, bam ->
            def bai = [ file("${bam}.bai"), file("${bam.parent}/${bam.baseName}.bai") ].find { f -> f.exists() }
            [ meta, bam, bai ?: [] ]
        }
        SAMTOOLS_INDEX(ch_bam_in)
        ch_bams_indexed = SAMTOOLS_INDEX.out.bam_bai
    }

    SAMTOOLS_STATS(ch_bams_indexed)
    EXTRACT_JUNCTIONS(ch_bams_indexed)
    MAKE_BIGWIG(ch_bams_indexed)
    ch_junctions = EXTRACT_JUNCTIONS.out.junctions

    // 3. COHORT COLLECTION
    ch_bams_for_cohort = ch_bams_indexed.map { _meta, bam, _bai -> bam }.collect()
    ch_bais_for_cohort = ch_bams_indexed.map { _meta, _bam, bai -> bai }.collect()

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

    // 5. TE ANALYSIS (optional)
    if (params.te_gtf) {
        ch_all_junctions = ch_junctions.map { _meta, bed -> bed }.collect()
        TE_ANALYSIS(
            ch_extended_gtf,
            file(params.te_gtf),
            ch_transcript_counts,
            ch_all_junctions,
            file(params.gtf)
        )
        TE_DESEQ2_FAMILY(TE_ANALYSIS.out.te_counts, ch_transcript_counts, ch_design, "Family", reference)
        TE_DESEQ2_INSTANCE(TE_ANALYSIS.out.te_instances, ch_transcript_counts, ch_design, "Instance", reference)
    }

    // 6. ISOFORM SWITCHING, DE AND SPLICING
    ISOSWITCH(ch_transcript_counts, ch_transcript_abundance, ch_extended_gtf, ch_design, reference)
    DESEQ2_GENE(BAMBU.out.gene_counts, ch_design, 'gene', reference)
    DESEQ2_TRANSCRIPT(ch_transcript_counts, ch_design, 'transcript', reference)
    SUPPA2(ch_extended_gtf, ch_transcript_counts)

    // 7. PROMOTER ACTIVITY
    PROACTIV(ch_extended_gtf, ch_design, ch_bams_for_cohort, ch_bais_for_cohort, reference)

    // 8. MULTIQC
    MULTIQC(
        channel
            .fromPath("${params.outdir}/**/*")
            .filter { f -> f.name =~ /(log|stats|metrics|summary)/ }
            .collect()
            .ifEmpty([])
    )

    // Completion message. `workflow` and `params` are captured in local variables because
    // they do not resolve inside a handler closure registered in the entry workflow.
    def wf     = workflow
    def outdir = params.outdir
    wf.onComplete {
        def line = '=' * 60
        if (wf.success) {
            log.info "\n${line}\nNextLongIso v${wf.manifest.version} finished successfully\n" +
                     "Output directory: ${outdir}\nDuration: ${wf.duration}\n${line}"
        } else {
            log.error "\n${line}\nNextLongIso v${wf.manifest.version} failed (exit status ${wf.exitStatus})\n" +
                      "Error: ${wf.errorMessage ?: 'see the failed task above'}\n" +
                      "${wf.errorReport ? wf.errorReport + '\n' : ''}" +
                      "Work directory: ${wf.workDir}\n" +
                      "Inspect the failed task's .command.err/.command.log, fix the cause and rerun with -resume.\n${line}"
        }
    }
}
