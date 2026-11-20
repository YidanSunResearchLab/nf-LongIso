/*
 * Promoter usage from transcript 5' ends (--promoter_method tss, default).
 * R code: bin/tss_sites.R (promoters, counts, activity), bin/usage_dexseq.R (differential usage)
 * The alternative proActiv implementation is in modules/proactiv.nf (--promoter_method proactiv).
 */

process PROMOTER_TSS {
    tag "promoter_tss"
    publishDir "${params.outdir}/promoter", mode: params.publish_dir_mode, pattern: "*.tsv"

    input:
    path iso_gtf        // Bambu extended_annotations.gtf
    path counts         // Bambu counts_transcript.txt
    path gene_counts    // Bambu counts_gene.txt (size factors)
    path design
    val  reference

    output:
    path "promoter_annotation.tsv"     , emit: annotation
    path "promoter_counts.tsv"         , emit: counts
    path "promoter_activity.tsv"       , emit: activity
    path "promoter_dexseq_results.tsv" , emit: results
    path "alternative_promoters.tsv"   , emit: alt_promoters
    path "versions.yml"                , emit: versions

    script:
    """
    # Exons of the extended annotation: transcript_id, gene_id, chr, start, end, strand
    awk -F '\\t' -v OFS='\\t' '\$3 == "exon" {
        match(\$9, /transcript_id "[^"]+"/); t = substr(\$9, RSTART + 15, RLENGTH - 16)
        match(\$9, /gene_id "[^"]+"/);       g = substr(\$9, RSTART + 9, RLENGTH - 10)
        print t, g, \$1, \$4, \$5, \$7 }' ${iso_gtf} > transcript_exons.tsv

    tss_sites.R \\
        --exons=transcript_exons.tsv \\
        --counts=${counts} \\
        --gene-counts=${gene_counts} \\
        --design=${design} \\
        --design-formula='${params.design}' \\
        --reference=${reference} \\
        --window=${params.tss_cluster_window} \\
        --max-width=${params.tss_max_cluster_width} \\
        --min-activity=${params.promoter_min_activity}

    usage_dexseq.R \\
        --counts=promoter_counts.tsv \\
        --design=${design} \\
        --design-formula='${params.design}' \\
        --reference=${reference} \\
        --min-count=${params.promoter_min_count} \\
        --prefix=promoter \\
        --usage-label=usage \\
        --cpus=${task.cpus}
    mv promoter_usage.tsv promoter_dexseq_results.tsv

    # Alternative promoters: promoters of genes with differential promoter usage (gene q < usage_alpha)
    awk -F '\\t' -v a=${params.usage_alpha} 'NR == 1 { print; for (i = 1; i <= NF; i++) if (\$i == "gene_qvalue") q = i; next }
        q && \$q != "NA" && \$q + 0 < a' promoter_dexseq_results.tsv > alternative_promoters.tsv
    echo "Alternative promoters (gene q < ${params.usage_alpha}): \$(( \$(wc -l < alternative_promoters.tsv) - 1 ))"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-dexseq: \$(Rscript -e 'cat(as.character(packageVersion("DEXSeq")))')
        bioconductor-deseq2: \$(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
    END_VERSIONS
    """

    stub:
    """
    touch promoter_annotation.tsv promoter_counts.tsv promoter_activity.tsv promoter_dexseq_results.tsv alternative_promoters.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}
