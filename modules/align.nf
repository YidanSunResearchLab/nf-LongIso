process MINIMAP2 {
    tag "${meta.id}"
    input:
    tuple val(meta), path(reads)
    path fasta
    path gtf 
    output:
    tuple val(meta), path("${meta.id}.sam"), emit: sam
    script:
    """
    # Create the junction guide from the GTF
    paftools.js gff2bed ${gtf} > junctions_guide.bed
    
    # Truncate read names in FASTQ (use \$0 for whole line, not \$1)
    echo "Truncating FASTQ read names..."
    
    if [[ ${reads} == *.gz ]]; then
        zcat ${reads} | awk '
        {
            if (NR % 4 == 1) {
                # Header line - truncate ENTIRE line to 250 chars
                if (length(\$0) > 250) {
                    print substr(\$0, 1, 250)
                } else {
                    print \$0
                }
            } else {
                # Other lines - keep as is
                print \$0
            }
        }' > reads_truncated.fastq
    else
        awk '
        {
            if (NR % 4 == 1) {
                if (length(\$0) > 250) {
                    print substr(\$0, 1, 250)
                } else {
                    print \$0
                }
            } else {
                print \$0
            }
        }' ${reads} > reads_truncated.fastq
    fi
    
    # Verify truncation worked
    echo "First 3 read names from truncated FASTQ:"
    head -12 reads_truncated.fastq | awk 'NR % 4 == 1 {print "Length:", length(\$0), "Name:", \$0}'
    
    # Check max length
    MAX_NAME_LEN=\$(head -1000 reads_truncated.fastq | awk 'NR % 4 == 1 {print length(\$0)}' | sort -rn | head -1)
    echo "Max FASTQ header length: \$MAX_NAME_LEN"
    
    if [ "\$MAX_NAME_LEN" -gt 250 ]; then
        echo "ERROR: FASTQ headers still too long!"
        exit 1
    fi
    
    # Align
    minimap2 -ax splice \
        -t ${task.cpus} \
        --junc-bed junctions_guide.bed \
        --secondary=no \
        ${fasta} reads_truncated.fastq > ${meta.id}.sam
    
    # Verify SAM read names
    echo "First 3 SAM read names:"
    grep -v "^@" ${meta.id}.sam | head -3 | cut -f1 | while read name; do
        echo "Length: \${#name} Name: \$name"
    done
    
    # Check SAM max length
    SAM_MAX_LEN=\$(grep -v "^@" ${meta.id}.sam | head -1000 | cut -f1 | awk '{print length(\$0)}' | sort -rn | head -1)
    echo "Max SAM read name length: \$SAM_MAX_LEN"
    
    if [ "\$SAM_MAX_LEN" -gt 251 ]; then
        echo "ERROR: SAM read names still too long!"
        exit 1
    fi
    
    rm reads_truncated.fastq
    """
}

process SAMTOOLS_SORT {
    tag "${meta.id}"
    publishDir "${params.outdir}/minimap2", mode: 'copy'
    input:
    tuple val(meta), path(sam)
    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam_bai
    path "*.stats", emit: log
    
    script:
    """
    echo "Final check before sorting..."
    MAX_LEN=\$(grep -v "^@" ${sam} | head -1000 | cut -f1 | awk '{print length(\$0)}' | sort -rn | head -1)
    echo "Max read name length in SAM: \$MAX_LEN"
    
    if [ "\$MAX_LEN" -gt 251 ]; then
        echo "FATAL ERROR: Read names too long for samtools!"
        echo "Showing problematic read names:"
        grep -v "^@" ${sam} | cut -f1 | awk 'length(\$0) > 251 {print "Length:", length(\$0), "Name:", \$0}' | head -10
        exit 1
    fi
    
    samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam ${sam}
    samtools index ${meta.id}.sorted.bam
    samtools quickcheck -v ${meta.id}.sorted.bam || exit 1
    samtools stats ${meta.id}.sorted.bam > ${meta.id}.sorted.bam.stats
    """
}

process EXTRACT_JUNCTIONS {
    tag "${meta.id}"
    publishDir "${params.outdir}/junctions", mode: 'copy'
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    output:
    tuple val(meta), path("${meta.id}.junctions.bed"), emit: junctions
    script:
    """
    samtools view -h ${bam} | \
    awk -v OFS="\\t" '
    /^@/ { next }
    {
        if (and(\$2, 4)) next
        
        chr = \$3
        pos = \$4
        cigar = \$6
        flag = \$2
        
        strand = "."
        for (i = 12; i <= NF; i++) {
            if (\$i ~ /^XS:A:/) {
                strand = substr(\$i, 6, 1)
                break
            }
        }
        
        if (strand == ".") {
            strand = (and(flag, 16)) ? "-" : "+"
        }
        
        ref_pos = pos
        cigar_copy = cigar
        while (match(cigar_copy, /([0-9]+)([MIDNSHPX=])/)) {
            len = substr(cigar_copy, RSTART, RLENGTH-1)
            op = substr(cigar_copy, RSTART+RLENGTH-1, 1)
            cigar_copy = substr(cigar_copy, RSTART+RLENGTH)
            
            if (op == "M" || op == "=" || op == "X" || op == "D") {
                ref_pos += len
            }
            else if (op == "N") {
                junc_start = ref_pos - 1
                junc_end = ref_pos + len - 1
                junc_id = "junc_" chr "_" junc_start "_" junc_end "_" strand
                print chr, junc_start, junc_end, junc_id, "1000", strand
                ref_pos += len
            }
        }
    }' | sort -k1,1 -k2,2n -k3,3n | uniq > ${meta.id}.junctions.bed
    
    if [ ! -s ${meta.id}.junctions.bed ]; then
        echo "Warning: No junctions extracted" >&2
        touch ${meta.id}.junctions.bed
    else
        echo "Extracted \$(wc -l < ${meta.id}.junctions.bed) junctions"
    fi
    """
}

process MAKE_BIGWIG {
    tag "${meta.id}"
    publishDir "${params.outdir}/bigwig", mode: 'copy'
    input:
    tuple val(meta), path(bam), path(bai)
    output:
    path "${meta.id}.bw", emit: bigwig
    script:
    """
    bamCoverage \
        -b ${bam} \
        -o ${meta.id}.bw \
        --binSize 5 \
        --minMappingQuality 5 \
        --normalizeUsing CPM \
        --numberOfProcessors ${task.cpus} \
        --skipNonCoveredRegions \
    """
}