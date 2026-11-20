/*
 * Alternative polyadenylation: PAS clusters from transcript 3' ends, internal-priming check,
 * PAS usage tested with DEXSeq. R code: bin/apa_sites.R, bin/usage_dexseq.R, bin/apa_summary.R
 */

process APA {
    tag "apa"
    publishDir "${params.outdir}/apa", mode: params.publish_dir_mode, pattern: "*.{tsv,csv}"

    input:
    path iso_gtf        // Bambu extended_annotations.gtf
    path counts         // Bambu counts_transcript.txt
    path design
    val  reference
    path fasta
    path fai

    output:
    path "pas_annotation.tsv"     , emit: annotation
    path "pas_transcripts.tsv"    , emit: transcripts
    path "pas_counts.tsv"         , emit: counts
    path "apa_dexseq_results.tsv" , emit: results
    path "APA_gene_summary.csv"   , emit: summary
    path "versions.yml"           , emit: versions

    script:
    """
    # Exons of the extended annotation: transcript_id, gene_id, chr, start, end, strand
    awk -F '\\t' -v OFS='\\t' '\$3 == "exon" {
        match(\$9, /transcript_id "[^"]+"/); t = substr(\$9, RSTART + 15, RLENGTH - 16)
        match(\$9, /gene_id "[^"]+"/);       g = substr(\$9, RSTART + 9, RLENGTH - 10)
        print t, g, \$1, \$4, \$5, \$7 }' ${iso_gtf} > transcript_exons.tsv

    apa_sites.R \\
        --exons=transcript_exons.tsv \\
        --counts=${counts} \\
        --genome=${fasta} \\
        --design=${design} \\
        --window=${params.apa_cluster_window} \\
        --max-width=${params.apa_max_cluster_width} \\
        --ip-a-frac=${params.apa_ip_a_frac} \\
        --ip-mode=${params.apa_internal_priming}

    usage_dexseq.R \\
        --counts=pas_counts.tsv \\
        --design=${design} \\
        --design-formula='${params.design}' \\
        --reference=${reference} \\
        --min-count=${params.apa_min_count} \\
        --prefix=apa \\
        --usage-label=PAU \\
        --cpus=${task.cpus}
    mv apa_usage.tsv apa_dexseq_results.tsv

    apa_summary.R --usage=apa_dexseq_results.tsv --annotation=pas_annotation.tsv --alpha=${params.usage_alpha}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-dexseq: \$(Rscript -e 'cat(as.character(packageVersion("DEXSeq")))')
        bioconductor-rsamtools: \$(Rscript -e 'cat(as.character(packageVersion("Rsamtools")))')
    END_VERSIONS
    """

    stub:
    """
    touch pas_annotation.tsv pas_transcripts.tsv pas_counts.tsv apa_dexseq_results.tsv APA_gene_summary.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}
