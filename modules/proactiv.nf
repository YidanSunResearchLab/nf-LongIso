/*
 * Promoter activity and alternative promoters (proActiv). R code: bin/proactiv.R
 */

process PROACTIV {
    tag "proactiv"
    publishDir "${params.outdir}/proactiv", mode: params.publish_dir_mode, pattern: "*.tsv"

    input:
    path gtf
    path design
    path bams
    path bais
    val  reference

    output:
    path "promoter_activity.tsv"    , emit: activity
    path "promoter_counts.tsv"      , emit: counts
    path "alternative_promoters.tsv", emit: alt_promoters, optional: true
    path "versions.yml"             , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # Keep only stranded records (makeTxDbFromGFF rejects strand '.')
    awk -F '\\t' '\$7 == "+" || \$7 == "-"' ${gtf} > cleaned_annotations.gtf

    # getAlternativePromoters() arguments, written verbatim (parsed as R by proactiv.R)
    cat <<'ALTPROM_ARGS' > altprom_args.txt
    ${args}
    ALTPROM_ARGS

    proactiv.R \\
        --gtf cleaned_annotations.gtf \\
        --design ${design} \\
        --species ${params.proactiv_species} \\
        --reference ${reference} \\
        --altprom-args-file altprom_args.txt \\
        --cpus ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-proactiv: \$(Rscript -e 'cat(as.character(packageVersion("proActiv")))')
    END_VERSIONS
    """

    stub:
    """
    touch promoter_activity.tsv promoter_counts.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}
