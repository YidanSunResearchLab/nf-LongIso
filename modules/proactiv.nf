process PROACTIV {
    tag "proactiv"
    publishDir "${params.outdir}/proactiv", mode: 'copy'

    input:
    path gtf
    val  metas
    path bams
    path bais

    output:
    path "promoter_activity.tsv"         , emit: activity
    path "promoter_counts.tsv"           , emit: counts
    path "alternative_promoters.tsv"     , emit: alt_promoters, optional: true

    script:
    def sample_ids   = metas.collect { it.id }.join('","')
    def sample_conds = metas.collect { it.condition }.join('","')
    def dollar = '$'

    """
    awk '${dollar}7 == "+" || ${dollar}7 == "-"' ${gtf} > cleaned_annotations.gtf

    cat <<'R_SCRIPT' > run_proactiv.R
    library(proActiv)
    library(GenomicFeatures)
    library(GenomicAlignments)
    library(SummarizedExperiment)

    # -- 1. Promoter annotation ------------------------------------------------
    message("Preparing promoter annotation...")
    txdb           <- makeTxDbFromGFF("cleaned_annotations.gtf")
    promoter_annot <- preparePromoterAnnotation(txdb = txdb, species = "Homo_sapiens")

    # -- 2. Extract junctions from BAM files -----------------------------------
    message("Extracting junctions from BAM files...")
    bam_files <- sort(list.files(pattern = "[.]bam${dollar}", full.names = TRUE))
    bam_files <- bam_files[!grepl("[.]bai${dollar}", bam_files)]
    message(sprintf("Found %d BAM files", length(bam_files)))

    junc_files   <- character(0)
    sample_names <- character(0)

    for (bam_file in bam_files) {
        message(sprintf("Processing: %s", basename(bam_file)))
        tryCatch({
            gal <- readGAlignments(
                bam_file,
                param = ScanBamParam(
                    flag = scanBamFlag(
                        isSecondaryAlignment = FALSE,
                        isUnmappedQuery      = FALSE
                    )
                )
            )
            juncs <- summarizeJunctions(gal)
            if (length(juncs) == 0) {
                warning(sprintf("No junctions found in %s - skipping", basename(bam_file)))
                next
            }
            df         <- as.data.frame(juncs)
            strand_num <- ifelse(as.character(df${dollar}strand) == "+", 1L,
                          ifelse(as.character(df${dollar}strand) == "-", 2L, 0L))
            star_df <- data.frame(
                chrom     = as.character(df${dollar}seqnames),
                start     = df${dollar}start,
                end       = df${dollar}end,
                strand    = strand_num,
                motif     = 1L,
                annotated = 1L,
                unique    = pmax(as.integer(df${dollar}score), 1L),
                multi     = 0L,
                overhang  = 10L
            )
            base_name <- sub("[.]sorted[.]bam${dollar}", "", basename(bam_file))
            base_name <- sub("[.]bam${dollar}",          "", base_name)
            out_name  <- paste0(base_name, ".SJ.out.tab")
            write.table(star_df, file = out_name, sep = "\\t",
                        quote = FALSE, row.names = FALSE, col.names = FALSE)
            junc_files   <- c(junc_files,   out_name)
            sample_names <- c(sample_names, base_name)
            message(sprintf("  Created: %s (%d junctions)", out_name, nrow(star_df)))
        }, error = function(e) {
            warning(sprintf("ERROR processing %s: %s", basename(bam_file), e${dollar}message))
        })
    }

    if (length(junc_files) == 0) stop("No junction files created - aborting.")

    # -- 3. Sample → condition mapping -----------------------------------------
    sample_ids        <- c("${sample_ids}")
    sample_conditions <- c("${sample_conds}")

    conditions <- sapply(sample_names, function(sname) {
        idx <- which(sample_ids == sname)
        if (length(idx) == 0)
            idx <- which(sapply(sample_ids, function(id) startsWith(sname, id)))
        if (length(idx) == 0)
            idx <- which(sapply(sample_ids, function(id) grepl(id, sname, fixed = TRUE)))
        if (length(idx) == 0) { warning(sprintf("No condition for '%s'", sname)); return("unknown") }
        sample_conditions[idx[1]]
    })

    message("Sample to condition mapping:")
    print(data.frame(Sample = sample_names, Condition = conditions, row.names = NULL))

    # -- 4. Run proActiv -------------------------------------------------------
    message("Running proActiv...")
    result <- proActiv(
        files              = junc_files,
        promoterAnnotation = promoter_annot,
        condition          = conditions,
        ncores             = ${task.cpus}
    )

    # -- 5. Export activity + counts -------------------------------------------
    message("Exporting results...")
    pos_info <- as.data.frame(rowData(result))

    list_cols <- which(sapply(pos_info, is.list))
    for (col in names(list_cols)) {
        pos_info[[col]] <- sapply(pos_info[[col]], function(x) {
            if (is.null(x) || length(x) == 0) return(NA_character_)
            paste(x, collapse = ",")
        })
    }

    assay_list      <- assays(result)
    activity_matrix <- as.data.frame(assay_list[["absolutePromoterActivity"]])
    counts_matrix   <- as.data.frame(assay_list[["promoterCounts"]])

    write.table(cbind(pos_info, activity_matrix), "promoter_activity.tsv",
                sep = "\\t", quote = FALSE, row.names = FALSE)
    write.table(cbind(pos_info, counts_matrix),   "promoter_counts.tsv",
                sep = "\\t", quote = FALSE, row.names = FALSE)

    # -- 6. getAlternativePromoters — the correct built-in method -------------
    # Requires >1 condition and uses both absolute + relative linear models
    # with a gene expression stability filter to exclude DE-driven changes.
    unique_conds <- unique(conditions[conditions != "unknown"])

    if (length(unique_conds) >= 2) {
        # Use control as reference so fold changes are treatment-vs-control
        ref_cond <- if ("control" %in% unique_conds) "control" else unique_conds[1]
        message(sprintf("Running getAlternativePromoters (reference: %s)...", ref_cond))

        tryCatch({
            alt <- getAlternativePromoters(
                result             = result,
                referenceCondition = ref_cond,
                minAbs             = 0.25,   # promoter must have abs activity >= 0.25
                minRel             = 0.05,   # promoter must have rel activity >= 0.05
                maxPval            = 0.05,   # BH-adjusted p-value threshold
                promoterFC         = 2.0,    # >= 2x fold change in promoter activity
                geneFC             = 1.5     # gene expression stable within 1.5x
                                             # <- this is what excludes DE-driven hits
            )

            up   <- alt${dollar}upReg
            down <- alt${dollar}downReg

            message(sprintf("Alternative promoters UP in reference (%s):   %d", ref_cond, nrow(up)))
            message(sprintf("Alternative promoters DOWN in reference (%s): %d", ref_cond, nrow(down)))

            if ((nrow(up) + nrow(down)) > 0) {
                up${dollar}direction   <- "up_in_reference"
                down${dollar}direction <- "down_in_reference"
                combined <- rbind(up, down)

                # Merge in coordinates from rowData for easier downstream use
                rdata    <- as.data.frame(rowData(result))
                rdata${dollar}promoterId <- rdata${dollar}promoterId
                combined <- merge(combined, 
                                  rdata[, c("promoterId", "seqnames", "start", 
                                            "strand", "internalPromoter", 
                                            "promoterPosition")],
                                  by = "promoterId", all.x = TRUE)

                combined <- combined[order(combined${dollar}padjAbs), ]
                write.table(combined, "alternative_promoters.tsv",
                            sep = "\\t", quote = FALSE, row.names = FALSE)
                message(sprintf("Written alternative_promoters.tsv (%d rows)", nrow(combined)))
            } else {
                message("No alternative promoters found - consider relaxing thresholds")
                message("  Try: minAbs=0.1, promoterFC=1.5, geneFC=2.0")
            }
        }, error = function(e) {
            message(sprintf("getAlternativePromoters error: %s", e${dollar}message))
        })
    } else {
        message("Only one condition found - skipping alternative promoter detection")
    }

    message("proActiv completed successfully")
R_SCRIPT

    Rscript run_proactiv.R
    """
}