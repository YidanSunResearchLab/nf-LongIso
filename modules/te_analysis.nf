/*
 * Transposable-element analysis of Bambu isoforms. Shell/awk code: bin/te_analysis.sh
 *
 * --te_count_mode exon (default): counts from exon-TE overlaps >= --te_min_overlap bp,
 *     one count per (transcript, subfamily) and per (transcript, TE instance)
 * --te_count_mode span: v0.1.0 counting (whole transcript span), to reproduce published results
 */

process TE_ANALYSIS {
    tag "TE_analysis"
    publishDir "${params.outdir}/TE", mode: params.publish_dir_mode, pattern: "*.tsv"

    input:
    path iso_gtf            // Bambu extended_annotations.gtf
    path te_gtf             // TE annotation (TEtranscripts format)
    path counts             // Bambu counts_transcript.txt
    path junctions_bed      // per-sample junction BEDs
    path ref_gtf            // reference annotation

    output:
    path "TE_family_counts.tsv"      , emit: te_counts
    path "TE_instance_counts.tsv"    , emit: te_instances
    path "isoform_TE_exonization.tsv", emit: te_exonization
    path "isoform_TE_annotation.tsv" , emit: te_annotation
    path "versions.yml"              , emit: versions

    script:
    """
    te_analysis.sh \\
        ${params.te_count_mode} \\
        ${params.te_min_overlap} \\
        ${params.te_strand} \\
        ${params.te_autonomous_frac} \\
        ${iso_gtf} \\
        ${te_gtf} \\
        ${counts} \\
        ${ref_gtf} \\
        ${junctions_bed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/^bedtools v//')
    END_VERSIONS
    """

    stub:
    """
    touch TE_family_counts.tsv TE_instance_counts.tsv isoform_TE_exonization.tsv isoform_TE_annotation.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}
