/*
 * Isoform switching and alternative TSS/TTS (IsoformSwitchAnalyzeR). R code: bin/run_isoswitch.R
 */

process ISOSWITCH {
    tag "isoform_switch"
    publishDir "${params.outdir}/isoswitch", mode: params.publish_dir_mode, pattern: "*.{csv,rds,pdf}"

    input:
    path counts
    path abundance
    path iso_gtf
    path design
    val  reference
    path nt_fasta       // transcript sequences (GFFREAD); [] when --isoswitch_consequences false
    path ref_gtf        // reference annotation (--gtf): annotated CDS for reference transcripts

    output:
    path "isoform_switch_results.csv"        , emit: switch_results
    path "significant_switches.csv"          , emit: sig_results
    path "atss_results.csv"                  , emit: atss_results, optional: true
    path "atts_results.csv"                  , emit: atts_results, optional: true
    path "all_AS_events.csv"                 , emit: all_as_events, optional: true
    path "switchlist.rds"                    , emit: switchlist, optional: true
    path "switch_consequences.csv"           , emit: consequences, optional: true
    path "consequence_summary.{csv,pdf}"     , emit: consequence_summary, optional: true
    path "consequence_enrichment.{csv,pdf}"  , emit: consequence_enrichment, optional: true
    path "switch_plots.pdf"                  , emit: switch_plots, optional: true
    path "versions.yml"                      , emit: versions

    script:
    def consequences = nt_fasta ? "--nt-fasta=${nt_fasta} --ref-gtf=${ref_gtf} --n-plots=${params.isoswitch_n_plots}" : ''
    """
    run_isoswitch.R \\
        --counts=${counts} \\
        --abundance=${abundance} \\
        --gtf=${iso_gtf} \\
        --design=${design} \\
        --design-formula='${params.design}' \\
        --reference=${reference} \\
        --alpha=${params.isoswitch_alpha} \\
        --dif-cutoff=${params.isoswitch_dif_cutoff} ${consequences}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-isoformswitchanalyzer: \$(Rscript -e 'cat(as.character(packageVersion("IsoformSwitchAnalyzeR")))')
    END_VERSIONS
    """

    stub:
    """
    touch isoform_switch_results.csv significant_switches.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}

process GFFREAD {
    tag "${gtf.name}"

    input:
    path gtf            // Bambu extended_annotations.gtf
    path fasta
    path fai

    output:
    path "transcripts.fa", emit: fasta
    path "versions.yml"  , emit: versions

    script:
    """
    gffread -w transcripts.fa -g ${fasta} ${gtf}
    echo "Transcript sequences: \$(grep -c '>' transcripts.fa)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch transcripts.fa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}
