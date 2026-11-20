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
    Rscript -e '
    library(bambu)
    annot <- prepareAnnotations("${gtf}")
    saveRDS(annot, "bambu_annotations.rds")
    '

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
    // Extra arguments are spliced into the bambu() call inside a single-quoted Rscript -e
    // program, so they must not contain single quotes (use double quotes for R strings).
    def args = task.ext.args ?: ''
    if (args.contains("'")) {
        error "--bambu_args must not contain single quotes; use double quotes for R strings (got: ${args})"
    }
    def extra     = args ? ", ${args}" : ''
    def extra_msg = args ? args.replace('"', '\\"') : 'none'   // escaped for the R print() string
    """
    Rscript -e '
    library(bambu)
    library(Biostrings)
    library(rtracklayer)

    print("Loading genome...")
    genome_seq <- readDNAStringSet("${fasta}")
    names(genome_seq) <- sub(" .*", "", names(genome_seq))

    print("Setting up BAM files...")
    all_files <- list.files(path = ".", full.names = TRUE)
    bam_paths <- all_files[endsWith(all_files, ".bam")]
    if (length(bam_paths) == 0) stop("ERROR: No BAM files found!")
    print(paste("Found", length(bam_paths), "BAM files:"))
    print(bam_paths)

    print("Running Bambu...")
    print("Extra bambu() arguments: ${extra_msg}")
    se <- bambu(
        reads       = bam_paths,
        annotations = readRDS("${annot_rds}"),
        genome      = genome_seq,
        discovery   = TRUE,
        quant       = TRUE,
        ncore       = ${task.cpus},
        lowMemory   = TRUE,
        verbose     = FALSE${extra}
    )

    print("Writing output files...")
    # writeBambuOutput writes:
    #   counts_transcript.txt, counts_gene.txt,
    #   CPM_transcript.txt, CPM_gene.txt  etc.
    writeBambuOutput(se, path = ".")

    print("Saving extended GTF...")
    writeToGTF(rowRanges(se), file = "extended_annotations.gtf")

    print("Verifying CPM_transcript.txt was written...")
    # writeBambuOutput should produce CPM_transcript.txt automatically.
    # If it does not (older bambu versions), write it manually:
    if (!file.exists("CPM_transcript.txt")) {
        print("CPM_transcript.txt not found - writing manually...")
        se_tx      <- se[["transcript"]]
        cpm_mat    <- as.data.frame(assay(se_tx, "CPM"))
        cpm_out    <- cbind(
            data.frame(
                TXNAME = rowData(se_tx)\$TXNAME,
                GENEID = rowData(se_tx)\$GENEID
            ),
            cpm_mat
        )
        write.table(cpm_out, "CPM_transcript.txt",
                    sep = "\\t", quote = FALSE, row.names = FALSE)
        print(paste("Written CPM_transcript.txt:", nrow(cpm_out), "rows"))
    }

    # Verify all expected files exist before finishing
    expected <- c("counts_transcript.txt", "counts_gene.txt",
                  "CPM_transcript.txt", "extended_annotations.gtf")
    for (f in expected) {
        if (file.exists(f)) {
            print(paste("OK:", f))
        } else {
            stop(paste("ERROR: expected output missing:", f))
        }
    }

    saveRDS(se, "bambu_results.rds")
    print("Bambu completed successfully!")
    '

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
