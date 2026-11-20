process BAMBU_PREP {
    tag "bambu_prep"
    publishDir "${params.outdir}/bambu", mode: params.publish_dir_mode

    input:
    path gtf

    output:
    path "bambu_annotations.rds", emit: annot_rds
    path "versions.yml"         , emit: versions

    script:
    """
    run_bambu_prep.R --gtf=${gtf} --out=bambu_annotations.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-bambu: \$(Rscript -e 'cat(as.character(packageVersion("bambu")))')
    END_VERSIONS
    """

    stub:
    """
    touch bambu_annotations.rds
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}

process BAMBU {
    tag "bambu"
    publishDir "${params.outdir}/bambu", mode: params.publish_dir_mode

    input:
    path bams
    path bais
    path annot_rds
    path fasta

    output:
    path "bambu_results.rds",        emit: bambu_rds
    path "extended_annotations.gtf", emit: extended_gtf
    path "counts_transcript.txt",    emit: transcript_counts
    path "counts_gene.txt",          emit: gene_counts
    path "CPM_transcript.txt",       emit: transcript_abundance
    path "versions.yml",             emit: versions

    script:
    // Extra arguments are appended to the bambu() call by bin/run_bambu.R. Single quotes are
    // rejected (as in the parameter schema): use double quotes for R strings.
    def args = task.ext.args ?: ''
    if (args.contains("'")) {
        error "--bambu_args must not contain single quotes; use double quotes for R strings (got: ${args})"
    }
    """
    # Extra bambu() arguments, written verbatim (appended to the call by run_bambu.R)
    cat <<'BAMBU_ARGS' > bambu_args.txt
    ${args}
    BAMBU_ARGS

    run_bambu.R \\
        --annotations=${annot_rds} \\
        --genome=${fasta} \\
        --ncore=${task.cpus} \\
        --extra-args-file=bambu_args.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-bambu: \$(Rscript -e 'cat(as.character(packageVersion("bambu")))')
    END_VERSIONS
    """

    stub:
    """
    touch bambu_results.rds extended_annotations.gtf counts_transcript.txt counts_gene.txt CPM_transcript.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}
