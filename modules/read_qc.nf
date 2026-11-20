/*
 * Read-level QC for FASTQ input (skipped for BAM input or with --skip_read_qc):
 *   NANOPLOT      read length / quality distributions per sample (NanoPlot)
 *   SEQKIT_STATS  read count, length and quality summary per sample (seqkit stats)
 * Both summaries are parsed by MultiQC.
 */

process NANOPLOT {
    tag "${meta.id}"
    publishDir "${params.outdir}/read_qc/nanoplot", mode: params.publish_dir_mode, pattern: "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}")                    , emit: report
    tuple val(meta), path("${meta.id}/${meta.id}_NanoStats.txt"), emit: stats
    path "versions.yml"                                    , emit: versions

    script:
    // Static images need a headless browser (kaleido) that is not in the container: interactive HTML only
    """
    NanoPlot \\
        --threads ${task.cpus} \\
        --fastq ${reads} \\
        --outdir ${meta.id} \\
        --prefix ${meta.id}_ \\
        --tsv_stats \\
        --no_static

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(NanoPlot --version | sed 's/^NanoPlot //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}
    printf 'Metrics\\tdataset\\nnumber_of_reads\\t0\\n' > ${meta.id}/${meta.id}_NanoStats.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}

process SEQKIT_STATS {
    tag "${meta.id}"
    publishDir "${params.outdir}/read_qc/seqkit", mode: params.publish_dir_mode, pattern: "*.tsv"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.seqkit_stats.tsv"), emit: stats
    path "versions.yml"                                 , emit: versions

    script:
    // The file column is set to the sample name so that MultiQC labels the row with it
    """
    seqkit stats --all --tabular --threads ${task.cpus} ${reads} \\
        | awk -v id='${meta.id}' 'BEGIN { FS = OFS = "\\t" } NR > 1 { \$1 = id } { print }' \\
        > ${meta.id}.seqkit_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/^seqkit v//')
    END_VERSIONS
    """

    stub:
    """
    printf 'file\\tformat\\ttype\\tnum_seqs\\tsum_len\\tmin_len\\tavg_len\\tmax_len\\n${meta.id}\\tFASTQ\\tDNA\\t0\\t0\\t0\\t0\\t0\\n' > ${meta.id}.seqkit_stats.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}
