/*
 * Differential expression (DESeq2). R code: bin/deseq2.R and bin/te_deseq2.R
 */

process DESEQ2 {
    tag "${level}"
    publishDir "${params.outdir}/${level == 'gene' ? 'DE_genes' : 'DE_transcripts'}",
        mode: params.publish_dir_mode, pattern: "*.{csv,pdf}"

    input:
    path counts            // counts_gene.txt (gene) or counts_transcript.txt (transcript)
    path design            // sample,condition
    val  level             // 'gene' or 'transcript'
    val  reference         // reference condition

    output:
    path "DE_results*.csv"        , emit: de_results
    path "normalized_counts.csv"  , emit: norm_counts
    path "pca_plot.pdf"           , emit: pca_pdf
    path "deseq2_${level}_mqc.tsv", emit: mqc
    path "versions.yml"           , emit: versions

    script:
    """
    deseq2.R \\
        --counts ${counts} \\
        --design ${design} \\
        --level ${level} \\
        --reference ${reference}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-deseq2: \$(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
    END_VERSIONS
    """

    stub:
    """
    touch DE_results.csv normalized_counts.csv pca_plot.pdf deseq2_${level}_mqc.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}

process TE_DESEQ2 {
    tag "DE_TE_${type}"
    publishDir "${params.outdir}/TE/DE", mode: params.publish_dir_mode, pattern: "*.{csv,pdf}"

    input:
    path te_counts
    path tx_counts
    path design
    val  type              // 'Family' or 'Instance'
    val  reference

    output:
    path "DE_TE_${type}_results*.csv", optional: true, emit: de_results
    path "norm_counts_${type}.csv"   , emit: norm_counts
    path "volcano_${type}.pdf"       , optional: true, emit: volcano
    path "pca_${type}.pdf"           , emit: pca
    path "versions.yml"              , emit: versions

    script:
    """
    te_deseq2.R \\
        --te-counts ${te_counts} \\
        --tx-counts ${tx_counts} \\
        --design ${design} \\
        --type ${type} \\
        --reference ${reference}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-deseq2: \$(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
    END_VERSIONS
    """

    stub:
    """
    touch DE_TE_${type}_results.csv norm_counts_${type}.csv volcano_${type}.pdf pca_${type}.pdf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}
