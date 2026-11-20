/*
 * Isoform switching, alternative TSS/TTS and APA (IsoformSwitchAnalyzeR). R code: bin/run_isoswitch.R
 */

process ISOSWITCH {
    tag "isoform_switch"
    publishDir "${params.outdir}/isoswitch", mode: params.publish_dir_mode, pattern: "*.{csv,rds}"

    input:
    path counts
    path abundance
    path iso_gtf
    path design
    val  reference

    output:
    path "isoform_switch_results.csv"        , emit: switch_results
    path "significant_switches.csv"          , emit: sig_results
    path "atss_results.csv"                  , emit: atss_results, optional: true
    path "atts_results.csv"                  , emit: atts_results, optional: true
    path "all_AS_events.csv"                 , emit: all_as_events, optional: true
    path "clustered_base_level_PAS_usage.csv", emit: pas_usage, optional: true
    path "APA_gene_summary.csv"              , emit: apa_summary, optional: true
    path "switchlist.rds"                    , emit: switchlist, optional: true
    path "versions.yml"                      , emit: versions

    script:
    """
    run_isoswitch.R \\
        --counts=${counts} \\
        --abundance=${abundance} \\
        --gtf=${iso_gtf} \\
        --design=${design} \\
        --reference=${reference} \\
        --alpha=${params.isoswitch_alpha} \\
        --dif-cutoff=${params.isoswitch_dif_cutoff} \\
        --seed=${params.seed}

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
