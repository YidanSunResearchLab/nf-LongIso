/*
 * Software versions and MultiQC report
 */

process SOFTWARE_VERSIONS {
    tag "versions"
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode, pattern: "software_versions.yml"

    input:
    path versions, stageAs: "versions/versions????.yml"   // versions.yml of every task

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc

    script:
    """
    # Version of the tools in this container (also used by MULTIQC)
    cat <<-END_VERSIONS > versions/this_process.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/^multiqc, version //')
        python: \$(python3 --version | sed 's/^Python //')
    END_VERSIONS

    dump_software_versions.py \\
        --pipeline-version ${workflow.manifest.version} \\
        --nextflow-version ${workflow.nextflow.version} \\
        versions/*.yml
    """

    stub:
    """
    touch software_versions.yml software_versions_mqc.yml
    """
}

process MULTIQC {
    tag "multiqc"
    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode, pattern: "multiqc_*"

    input:
    path multiqc_files, stageAs: "?/*"     // samtools stats/flagstat, DESeq2 tables, software versions

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data
    path "versions.yml"       , emit: versions

    script:
    """
    multiqc . \\
        --filename multiqc_report.html \\
        --title "NextLongIso report" \\
        --force \\
        --comment "NextLongIso v${workflow.manifest.version}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/^multiqc, version //')
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_report_data
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}
