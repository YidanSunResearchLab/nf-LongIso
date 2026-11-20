/*
 * Alignment, BAM preparation and per-sample BAM-derived outputs.
 *
 * FASTQ input: MINIMAP2_INDEX (once) + JUNC_BED (once) -> MINIMAP2_ALIGN (per sample)
 * BAM input:   SAMTOOLS_INDEX (per sample; stages the BAM as <sample>.sorted.bam)
 * Both:        SAMTOOLS_STATS, EXTRACT_JUNCTIONS, MAKE_BIGWIG
 *
 * All awk code is POSIX (no gawk extensions) and runs with BusyBox awk.
 */

process MINIMAP2_INDEX {
    tag "${fasta.name}"

    input:
    path fasta

    output:
    path "genome.mmi"   , emit: index
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    minimap2 -x ${params.minimap2_preset} ${args} -t ${task.cpus} -d genome.mmi ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    """
    touch genome.mmi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}

process JUNC_BED {
    tag "${gtf.name}"

    input:
    path gtf

    output:
    path "junctions_guide.bed", emit: bed
    path "versions.yml"       , emit: versions

    script:
    """
    paftools.js gff2bed ${gtf} > junctions_guide.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    """
    touch junctions_guide.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}

process MINIMAP2_ALIGN {
    tag "${meta.id}"
    publishDir "${params.outdir}/minimap2", mode: params.publish_dir_mode, pattern: "*.sorted.bam*"

    input:
    tuple val(meta), path(reads)
    path index
    path junc_bed      // [] when --minimap2_junc_bed false

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam_bai
    path "versions.yml", emit: versions

    script:
    def args         = task.ext.args ?: ''
    def junc_opt     = junc_bed ? "--junc-bed ${junc_bed}" : ''
    def reader       = reads.name.endsWith('.gz') ? 'zcat' : 'cat'
    def sort_threads = Math.max(1, (task.cpus as int).intdiv(4))
    def minimap2_cmd = "minimap2 -ax ${params.minimap2_preset} -t ${task.cpus} ${junc_opt} --secondary=no ${args} ${index} -"
    """
    set -o pipefail

    echo "Running alignment command:"
    cat <<'ALIGN_CMD'
    ${reader} ${reads} | awk '<truncate FASTQ headers to 250 characters>' | ${minimap2_cmd} | samtools sort -@ ${sort_threads} -o ${meta.id}.sorted.bam -
    ALIGN_CMD

    # Header lines (including '@') are cut to 250 characters: samtools rejects read names > 254
    ${reader} ${reads} \\
        | awk 'NR % 4 == 1 && length(\$0) > 250 { \$0 = substr(\$0, 1, 250) } { print }' \\
        | ${minimap2_cmd} \\
        | samtools sort -@ ${sort_threads} -T ${meta.id}.sorttmp -o ${meta.id}.sorted.bam -
    samtools index ${meta.id}.sorted.bam

    MAX_LEN=\$(samtools view ${meta.id}.sorted.bam | awk '{ if (length(\$1) > m) m = length(\$1) } END { print m + 0 }')
    echo "Max read name length in BAM: \$MAX_LEN"
    if [ "\$MAX_LEN" -gt 251 ]; then
        echo "ERROR: read names longer than 251 characters in ${meta.id}.sorted.bam" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -n 1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}

process SAMTOOLS_INDEX {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam, stageAs: 'input/*'), path(bai, stageAs: 'input_index/*')   // bai is [] if absent

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam_bai
    path "versions.yml", emit: versions

    script:
    """
    samtools view -H ${bam} > header.sam
    if ! grep -q '^@HD.*SO:coordinate' header.sam; then
        echo "ERROR: BAM for sample '${meta.id}' (${bam.name}) is not coordinate-sorted (no SO:coordinate in the @HD line). Sort it with 'samtools sort' first." >&2
        exit 1
    fi

    # Stage under the sample name so that Bambu and all downstream tables use <sample> as column name
    ln -s ${bam} ${meta.id}.sorted.bam
    if [ -n "${bai}" ]; then
        ln -s ${bai} ${meta.id}.sorted.bam.bai
    else
        echo "No index found for ${bam.name}; running samtools index"
        samtools index -@ ${task.cpus} ${meta.id}.sorted.bam
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}

process SAMTOOLS_STATS {
    tag "${meta.id}"
    publishDir "${params.outdir}/minimap2", mode: params.publish_dir_mode, pattern: "*.{stats,flagstat}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam.stats")   , emit: stats
    tuple val(meta), path("${meta.id}.sorted.bam.flagstat"), emit: flagstat
    path "versions.yml", emit: versions

    script:
    """
    samtools stats    -@ ${task.cpus} ${bam} > ${meta.id}.sorted.bam.stats
    samtools flagstat -@ ${task.cpus} ${bam} > ${meta.id}.sorted.bam.flagstat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam.stats ${meta.id}.sorted.bam.flagstat
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}

process EXTRACT_JUNCTIONS {
    tag "${meta.id}"
    publishDir "${params.outdir}/junctions", mode: params.publish_dir_mode, pattern: "*.junctions.bed"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.junctions.bed"), emit: junctions
    path "versions.yml", emit: versions

    script:
    """
    set -o pipefail
    # Junctions from the CIGAR strings (bin/extract_junctions.awk)
    samtools view -F 4 ${bam} \\
        | awk -f "\$(command -v extract_junctions.awk)" \\
        | sort -k1,1 -k2,2n -k3,3n | uniq > ${meta.id}.junctions.bed

    if [ ! -s ${meta.id}.junctions.bed ]; then
        echo "Warning: no junctions extracted from ${bam}" >&2
    else
        echo "Extracted \$(wc -l < ${meta.id}.junctions.bed) junctions"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.junctions.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}

process MAKE_BIGWIG {
    tag "${meta.id}"
    publishDir "${params.outdir}/bigwig", mode: params.publish_dir_mode, pattern: "*.bw"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.bw"), emit: bigwig
    path "versions.yml", emit: versions

    script:
    // Extra options are appended after the defaults; argparse keeps the last value of a
    // repeated option, so e.g. --binSize 10 overrides --binSize 5.
    def args = task.ext.args ?: ''
    """
    bamCoverage \\
        -b ${bam} \\
        -o ${meta.id}.bw \\
        --binSize 5 \\
        --minMappingQuality 5 \\
        --normalizeUsing CPM \\
        --numberOfProcessors ${task.cpus} \\
        --skipNonCoveredRegions \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bamCoverage --version | sed 's/^bamCoverage //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.bw
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}
