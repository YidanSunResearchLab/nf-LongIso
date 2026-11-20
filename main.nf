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
include { SOFTWARE_VERSIONS; MULTIQC }                                         from './modules/multiqc.nf'
include { MAKE_DESIGN }                                                        from './modules/design.nf'

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
    def conditions = samples.collect { s -> s[0].condition }.unique().sort()
    log.info """
        Run settings
          Input type          : ${params.input_type} (${samples.size()} samples)
          Conditions          : ${conditions.join(', ')}
          Reference condition : ${reference}${params.reference_condition ? '' : ' (chosen automatically)'}
          Genome / GTF        : ${params.genome} / ${params.gtf}
          TE analysis         : ${params.te_gtf ? "on (${params.te_gtf}); te_count_mode=${params.te_count_mode}, te_min_overlap=${params.te_min_overlap}, te_strand=${params.te_strand}" : 'off (no --te_gtf)'}
          minimap2            : ${params.input_type == 'fastq' ? "preset=${params.minimap2_preset}, junc_bed=${params.minimap2_junc_bed}, index=${params.minimap2_index ?: 'built by the pipeline'}, args='${params.minimap2_args}'" : 'not used (BAM input)'}
          Tool arguments      : bamCoverage='${params.bamcoverage_args}', bambu='${params.bambu_args}', SUPPA2 events='${params.suppa2_events}'
          proActiv            : species=${params.proactiv_species}; ${params.proactiv_altprom_args}
          Output directory    : ${params.outdir}
        """.stripIndent()
    ch_input = channel.fromList(samples)
    ch_versions = channel.empty()

    // Normalised design table (sample,condition) used by all R steps
    ch_design = MAKE_DESIGN(samples.collect { s -> "${s[0].id},${s[0].condition}" }).csv

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
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
        if (!params.minimap2_index) ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
        if (params.minimap2_junc_bed) ch_versions = ch_versions.mix(JUNC_BED.out.versions)
    } else {
        ch_bam_in = ch_input.map { meta, bam ->
            def bai = [ file("${bam}.bai"), file("${bam.parent}/${bam.baseName}.bai") ].find { f -> f.exists() }
            [ meta, bam, bai ?: [] ]
        }
        SAMTOOLS_INDEX(ch_bam_in)
        ch_bams_indexed = SAMTOOLS_INDEX.out.bam_bai
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    }

    SAMTOOLS_STATS(ch_bams_indexed)
    EXTRACT_JUNCTIONS(ch_bams_indexed)
    MAKE_BIGWIG(ch_bams_indexed)
    ch_junctions = EXTRACT_JUNCTIONS.out.junctions
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions, EXTRACT_JUNCTIONS.out.versions, MAKE_BIGWIG.out.versions)

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
    ch_versions = ch_versions.mix(BAMBU_PREP.out.versions, BAMBU.out.versions)

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
        ch_versions = ch_versions.mix(TE_ANALYSIS.out.versions, TE_DESEQ2_FAMILY.out.versions, TE_DESEQ2_INSTANCE.out.versions)
    }

    // 6. ISOFORM SWITCHING, DE AND SPLICING
    ISOSWITCH(ch_transcript_counts, ch_transcript_abundance, ch_extended_gtf, ch_design, reference)
    DESEQ2_GENE(BAMBU.out.gene_counts, ch_design, 'gene', reference)
    DESEQ2_TRANSCRIPT(ch_transcript_counts, ch_design, 'transcript', reference)
    SUPPA2(ch_extended_gtf, ch_transcript_counts)
    ch_versions = ch_versions.mix(ISOSWITCH.out.versions, DESEQ2_GENE.out.versions,
                                  DESEQ2_TRANSCRIPT.out.versions, SUPPA2.out.versions)

    // 7. PROMOTER ACTIVITY
    PROACTIV(ch_extended_gtf, ch_design, ch_bams_for_cohort, ch_bais_for_cohort, reference)

    ch_versions = ch_versions.mix(PROACTIV.out.versions)

    // 8. SOFTWARE VERSIONS AND MULTIQC (inputs come from channels, not from the output directory)
    SOFTWARE_VERSIONS(ch_versions.collect())
    ch_multiqc = channel.empty().mix(
        SAMTOOLS_STATS.out.stats.map { _meta, f -> f },
        SAMTOOLS_STATS.out.flagstat.map { _meta, f -> f },
        DESEQ2_GENE.out.mqc,
        DESEQ2_TRANSCRIPT.out.mqc,
        SOFTWARE_VERSIONS.out.mqc
    )
    MULTIQC(ch_multiqc.collect())

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
            log.error "\n${line}\nNextLongIso v${wf.manifest.version} failed\n" +
                      "Error: ${wf.errorMessage ?: 'see the failed task above'}\n" +
                      "${wf.errorReport ? wf.errorReport + '\n' : ''}" +
                      "Work directory: ${wf.workDir}\n" +
                      "Inspect the failed task's .command.err/.command.log, fix the cause and rerun with -resume.\n${line}"
        }
    }
}
