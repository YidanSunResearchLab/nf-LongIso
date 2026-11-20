#!/usr/bin/env nextflow
/*
 * NextLongIso (nf-LongIso): multi-dimensional long-read RNA-seq analysis
 * https://github.com/YidanSunResearchLab/nf-LongIso
 *
 * All parameter defaults live in nextflow.config; nextflow_schema.json documents and validates them.
 */

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

include { MINIMAP2_INDEX; JUNC_BED; MINIMAP2_ALIGN; SAMTOOLS_INDEX; SAMTOOLS_FAIDX; GTF_INTRONS;
          SAMTOOLS_STATS; EXTRACT_JUNCTIONS; MAKE_BIGWIG }                    from './modules/align.nf'
include { BAMBU_PREP; BAMBU }                                                  from './modules/bambu.nf'
include { SUPPA2; SUPPA2_DIFFSPLICE }                                          from './modules/suppa2.nf'
include { TE_ANALYSIS }                                                        from './modules/te_analysis.nf'
include { DESEQ2 as DESEQ2_GENE; DESEQ2 as DESEQ2_TRANSCRIPT;
          TE_DESEQ2 as TE_DESEQ2_FAMILY; TE_DESEQ2 as TE_DESEQ2_INSTANCE }   from './modules/deseq2.nf'
include { ISOSWITCH; GFFREAD }                                                 from './modules/isoswitch.nf'
include { PROACTIV }                                                           from './modules/proactiv.nf'
include { SOFTWARE_VERSIONS; MULTIQC }                                         from './modules/multiqc.nf'
include { MAKE_DESIGN }                                                        from './modules/design.nf'
include { APA }                                                                from './modules/apa.nf'
include { PROMOTER_TSS }                                                       from './modules/promoter.nf'
include { NANOPLOT; SEQKIT_STATS }                                             from './modules/read_qc.nf'
include { SQANTI3_QC }                                                         from './modules/sqanti3.nf'

/*
 * Covariates of --design ('~ condition' or '~ cov1 + cov2 + condition'), in formula order.
 */
def designCovariates(design) {
    def d = design.toString().trim()
    if (!(d ==~ /^~\s*([A-Za-z][A-Za-z0-9_]*\s*\+\s*)*condition$/)) {
        error "--design '${design}' must look like '~ condition' or '~ batch + condition': samplesheet columns joined by '+', with condition as the last term"
    }
    def covs = d.replaceFirst(/^~/, '').split(/\+/).collect { t -> t.trim() }.findAll { t -> t != 'condition' }
    def reserved = covs.findAll { c -> c in ['sample', 'fastq', 'bam'] }
    if (reserved) error "--design: column(s) ${reserved.join(', ')} cannot be used as covariates"
    if (covs.size() != covs.unique(false).size()) error "--design '${design}' lists a covariate twice"
    return covs
}

/*
 * Read and check the samplesheet. Returns a list of [meta, file] where file is the
 * FASTQ or BAM column, depending on --input_type; meta.covariates holds the design covariates.
 */
def readSamplesheet(samplesheet, input_type, covariates) {
    def rows = samplesheetToList(samplesheet, "${projectDir}/assets/schema_input.json")
    // Covariate columns are not part of the schema: read them from the CSV itself
    def table = covariates ? file(samplesheet).splitCsv(header: true) : []
    def name_rule = /^[A-Za-z][A-Za-z0-9_]*$/
    covariates.each { c ->
        if (!(c ==~ name_rule)) error "Covariate column name '${c}' must match ${name_rule}"
        if (!table || !table[0].containsKey(c)) error "--design uses column '${c}', which is not in the samplesheet"
        def values = table.collect { r -> r[c]?.toString()?.trim() }
        def bad = values.findAll { v -> !(v ==~ name_rule) }
        if (bad) {
            error "Samplesheet column '${c}' (used in --design) has invalid value(s) ${bad.unique().collect { v -> "'${v}'" }.join(', ')}: values must start with a letter and contain only letters, digits and underscores"
        }
        if (values.unique(false).size() < 2) error "Samplesheet column '${c}' (used in --design) has a single value; remove it from --design"
    }
    def cov_by_sample = table.collectEntries { r -> [ (r.sample.toString().trim()): covariates.collectEntries { c -> [ (c): r[c].toString().trim() ] } ] }
    def samples = rows.collect { row ->
        def meta = [ id: row[0].id.toString(), condition: row[0].condition.toString() ]
        meta.covariates = covariates ? cov_by_sample[meta.id] : [:]
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

    // --help / --helpFull: nf-schema has already printed the help text
    if (params.help || params.helpFull) {
        return
    }

    validateParameters()
    log.info "NextLongIso v${workflow.manifest.version}"
    log.info paramsSummaryLog(workflow)

    // 1. INPUT
    def covariates = designCovariates(params.design)
    def samples   = readSamplesheet(params.samplesheet, params.input_type, covariates)
    def reference = resolveReference(samples, params.reference_condition)
    def conditions = samples.collect { s -> s[0].condition }.unique().sort()
    def n_reps = conditions.collectEntries { c -> [ (c): samples.count { s -> s[0].condition == c } ] }
    // SUPPA2 diffSplice: each non-reference condition with >= 2 replicates vs a reference with >= 2 replicates
    def diff_conditions = n_reps[reference] >= 2 ? conditions.findAll { c -> c != reference && n_reps[c] >= 2 } : []
    log.info """
        Run settings
          Input type          : ${params.input_type} (${samples.size()} samples)
          Conditions          : ${conditions.join(', ')}
          Reference condition : ${reference}${params.reference_condition ? '' : ' (chosen automatically)'}
          Design              : ${params.design}${covariates ? " (covariates: ${covariates.join(', ')})" : ''}
          Genome / GTF        : ${params.genome} / ${params.gtf}
          TE analysis         : ${params.te_gtf ? "on (${params.te_gtf}); te_count_mode=${params.te_count_mode}, te_min_overlap=${params.te_min_overlap}, te_strand=${params.te_strand}" : 'off (no --te_gtf)'}
          minimap2            : ${params.input_type == 'fastq' ? "preset=${params.minimap2_preset}, junc_bed=${params.minimap2_junc_bed}, index=${params.minimap2_index ?: 'built by the pipeline'}, args='${params.minimap2_args}'" : 'not used (BAM input)'}
          Tool arguments      : bamCoverage='${params.bamcoverage_args}', bambu='${params.bambu_args}', SUPPA2 events='${params.suppa2_events}'
          SUPPA2 diffSplice   : ${!params.suppa2_diff ? 'off' : diff_conditions ? "${diff_conditions.join(', ')} vs ${reference}; ${params.suppa2_diff_args}" : 'skipped (needs >= 2 replicates in the reference and another condition)'}
          Read QC             : ${params.input_type != 'fastq' ? 'not run (BAM input)' : params.skip_read_qc ? 'skipped (--skip_read_qc)' : 'NanoPlot, seqkit stats'}
          SQANTI3 QC          : ${params.run_sqanti3 ? 'on' + (params.sqanti3_args ? " (${params.sqanti3_args})" : '') : 'off (--run_sqanti3 to enable)'}
          proActiv            : species=${params.proactiv_species}; ${params.proactiv_altprom_args}
          Output directory    : ${params.outdir}
        """.stripIndent()
    ch_input = channel.fromList(samples)
    ch_versions = channel.empty()
    ch_read_qc  = channel.empty()

    // Normalised design table (sample,condition) used by all R steps
    def design_lines = [ (['sample', 'condition'] + covariates).join(',') ] +
        samples.collect { s -> ([ s[0].id, s[0].condition ] + covariates.collect { c -> s[0].covariates[c] }).join(',') }
    MAKE_DESIGN(design_lines, reference)
    ch_design   = MAKE_DESIGN.out.csv
    ch_versions = ch_versions.mix(MAKE_DESIGN.out.versions)

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

        // Read QC (FASTQ input only)
        if (!params.skip_read_qc) {
            NANOPLOT(ch_input)
            SEQKIT_STATS(ch_input)
            ch_read_qc  = ch_read_qc.mix(NANOPLOT.out.stats.map { _meta, f -> f }, SEQKIT_STATS.out.stats.map { _meta, f -> f })
            ch_versions = ch_versions.mix(NANOPLOT.out.versions, SEQKIT_STATS.out.versions)
        }
    } else {
        ch_bam_in = ch_input.map { meta, bam ->
            def bai = [ file("${bam}.bai"), file("${bam.parent}/${bam.baseName}.bai") ].find { f -> f.exists() }
            [ meta, bam, bai ?: [] ]
        }
        SAMTOOLS_INDEX(ch_bam_in)
        ch_bams_indexed = SAMTOOLS_INDEX.out.bam_bai
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    }

    // Genome index (.fai): use <genome>.fai when present, otherwise build it
    def fai_file = file("${params.genome}.fai")
    if (fai_file.exists()) {
        ch_fai = channel.value(fai_file)
    } else {
        ch_fai = SAMTOOLS_FAIDX(file(params.genome)).fai
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }
    GTF_INTRONS(file(params.gtf))

    SAMTOOLS_STATS(ch_bams_indexed)
    EXTRACT_JUNCTIONS(ch_bams_indexed, file(params.genome), ch_fai, GTF_INTRONS.out.introns)
    MAKE_BIGWIG(ch_bams_indexed)
    ch_junctions = EXTRACT_JUNCTIONS.out.junctions
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions, EXTRACT_JUNCTIONS.out.versions, MAKE_BIGWIG.out.versions,
                                  GTF_INTRONS.out.versions)

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
        // TE size factors come from the gene-level counts
        TE_DESEQ2_FAMILY(TE_ANALYSIS.out.te_counts, BAMBU.out.gene_counts, ch_design, "Family", reference)
        TE_DESEQ2_INSTANCE(TE_ANALYSIS.out.te_instances, BAMBU.out.gene_counts, ch_design, "Instance", reference)
        ch_versions = ch_versions.mix(TE_ANALYSIS.out.versions, TE_DESEQ2_FAMILY.out.versions, TE_DESEQ2_INSTANCE.out.versions)
    }

    // 6. ISOFORM SWITCHING, DE AND SPLICING
    // Transcript sequences for ORF prediction and switch consequences
    if (params.isoswitch_consequences) {
        GFFREAD(ch_extended_gtf, file(params.genome), ch_fai)
        ch_nt_fasta = GFFREAD.out.fasta
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    } else {
        ch_nt_fasta = channel.value([])
    }
    ISOSWITCH(ch_transcript_counts, ch_transcript_abundance, ch_extended_gtf, ch_design, reference, ch_nt_fasta)
    DESEQ2_GENE(BAMBU.out.gene_counts, ch_design, 'gene', reference)
    DESEQ2_TRANSCRIPT(ch_transcript_counts, ch_design, 'transcript', reference)
    SUPPA2(ch_extended_gtf, ch_transcript_counts)
    if (params.suppa2_diff && diff_conditions) {
        def no_reps = conditions.findAll { c -> c != reference && !(c in diff_conditions) }
        if (no_reps) {
            log.warn "SUPPA2 diffSplice: condition(s) ${no_reps.join(', ')} have fewer than 2 replicates and are not compared"
        }
        if (covariates) {
            log.warn "SUPPA2 diffSplice does not model the design covariates (${covariates.join(', ')}); it compares conditions only"
        }
        SUPPA2_DIFFSPLICE(SUPPA2.out.ioe, SUPPA2.out.psi_per_sample, ch_transcript_abundance, ch_design, reference, diff_conditions)
        ch_versions = ch_versions.mix(SUPPA2_DIFFSPLICE.out.versions)
    } else if (params.suppa2_diff) {
        log.warn "SUPPA2 diffSplice skipped: it needs >= 2 replicates in the reference condition (${reference}) and in at least one other condition"
    }
    // Alternative polyadenylation (DEXSeq on PAS clusters)
    APA(ch_extended_gtf, ch_transcript_counts, ch_design, reference, file(params.genome), ch_fai)
    ch_versions = ch_versions.mix(APA.out.versions)

    ch_versions = ch_versions.mix(ISOSWITCH.out.versions, DESEQ2_GENE.out.versions,
                                  DESEQ2_TRANSCRIPT.out.versions, SUPPA2.out.versions)

    // 7. PROMOTER USAGE: transcript 5' ends (default) or proActiv (junction-based)
    if (params.promoter_method == 'tss') {
        PROMOTER_TSS(ch_extended_gtf, ch_transcript_counts, BAMBU.out.gene_counts, ch_design, reference)
        ch_promoter_versions = PROMOTER_TSS.out.versions
    } else {
        PROACTIV(ch_extended_gtf, ch_design, ch_bams_for_cohort, ch_bais_for_cohort, reference)
        ch_promoter_versions = PROACTIV.out.versions
    }

    ch_versions = ch_versions.mix(ch_promoter_versions)

    // 8. TRANSCRIPT CLASSIFICATION (optional): SQANTI3 QC of the Bambu extended annotation
    if (params.run_sqanti3) {
        SQANTI3_QC(ch_extended_gtf, file(params.gtf), file(params.genome), ch_fai)
        ch_versions = ch_versions.mix(SQANTI3_QC.out.versions)
    }

    // 9. SOFTWARE VERSIONS AND MULTIQC (inputs come from channels, not from the output directory)
    SOFTWARE_VERSIONS(ch_versions.collect())
    ch_multiqc = channel.empty().mix(
        SAMTOOLS_STATS.out.stats.map { _meta, f -> f },
        SAMTOOLS_STATS.out.flagstat.map { _meta, f -> f },
        DESEQ2_GENE.out.mqc,
        DESEQ2_TRANSCRIPT.out.mqc,
        ch_read_qc,
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
