/*
 * Structural classification of the Bambu transcripts against the reference annotation
 * (SQANTI3 QC module only; optional, --run_sqanti3).
 */

process SQANTI3_QC {
    tag "sqanti3"
    publishDir "${params.outdir}/sqanti3", mode: params.publish_dir_mode,
        pattern: "{nextlongiso_classification.txt,nextlongiso_junctions.txt,nextlongiso_SQANTI3_report.html,sqanti3_category_counts.tsv}"

    input:
    path gtf            // Bambu extended_annotations.gtf (the transcripts to classify)
    path ref_gtf        // reference annotation (--gtf)
    path fasta          // genome FASTA
    path fai            // genome index

    output:
    path "nextlongiso_classification.txt"   , emit: classification
    path "nextlongiso_junctions.txt"        , emit: junctions
    path "nextlongiso_SQANTI3_report.html"  , emit: report, optional: true
    path "sqanti3_category_counts.tsv"      , emit: categories
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # The official image activates its conda environment in the Docker ENTRYPOINT, which
    # Singularity/Apptainer exec does not run: put the environment on PATH explicitly.
    if [ -d /conda/miniconda3/envs/sqanti3/bin ]; then export PATH="/conda/miniconda3/envs/sqanti3/bin:\$PATH"; fi
    export HOME="\$PWD"

    sqanti3_qc.py \\
        --isoforms ${gtf} \\
        --refGTF ${ref_gtf} \\
        --refFasta ${fasta} \\
        --output nextlongiso \\
        --dir . \\
        --cpus ${task.cpus} \\
        ${args}

    awk -f "\$(command -v sqanti3_categories.awk)" fileno=1 ${ref_gtf} fileno=2 nextlongiso_classification.txt \\
        > sqanti3_category_counts.tsv
    cat sqanti3_category_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \$(sqanti3_qc.py --version 2>&1 | grep -o 'SQANTI3 [0-9][0-9.]*' | sed 's/^SQANTI3 //' | head -n 1)
    END_VERSIONS
    """

    stub:
    """
    printf 'isoform\\tstructural_category\\nBambuTx1\\tnovel_not_in_catalog\\n' > nextlongiso_classification.txt
    touch nextlongiso_junctions.txt nextlongiso_SQANTI3_report.html
    awk -f "\$(command -v sqanti3_categories.awk)" fileno=1 /dev/null fileno=2 nextlongiso_classification.txt \\
        > sqanti3_category_counts.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}
