process MULTIQC {
    label 'process_medium'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path 'inputs/*' 

    output:
    path "multiqc_report.html", emit: report, optional: true
    path "multiqc_data",         emit: data,   optional: true

    script:
    """
    multiqc . \
        --filename multiqc_report.html \
        --title "IsoLong Pipeline Report" \
        --force \
        --comment "Report generated for: ${workflow.runName}"
    """
}