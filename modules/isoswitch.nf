process ISOSWITCH {
    tag "isoform_switch"
    publishDir "${params.outdir}/isoswitch", mode: 'copy'

    input:
    path counts
    path abundance
    path iso_gtf
    path samplesheet

    output:
    path "isoform_switch_results.csv", emit: switch_results
    path "significant_switches.csv", emit: sig_results
    path "atss_results.csv", emit: atss_results, optional: true
    path "atts_results.csv", emit: atts_results, optional: true
    path "all_AS_events.csv", emit: all_as_events, optional: true
    path "clustered_base_level_PAS_usage.csv", emit: pas_usage, optional: true
    path "APA_gene_summary.csv", emit: apa_summary, optional: true
    path "switchlist.rds", emit: switchlist, optional: true

    script:
    """
    cat <<'R_SCRIPT' > run_isoswitch.R
    .libPaths(unique(c(.libPaths(), "/usr/local/lib/R/site-library", "/usr/lib/R/library")))

    suppressPackageStartupMessages({
        library(IsoformSwitchAnalyzeR)
        library(ggplot2)
        library(GenomicRanges)
    })

    # [Keep all your existing data import code exactly the same up to analyzeAlternativeSplicing]
    
    counts_data <- read.table("${counts}", header=TRUE, check.names=FALSE, sep="\\t")
    abund_data  <- read.table("${abundance}", header=TRUE, check.names=FALSE, sep="\\t")
    colnames(counts_data)[1] <- "isoform_id"
    colnames(abund_data)[1]  <- "isoform_id"

    clean_names <- function(x) {
        x <- sub("\\\\..*", "", x) 
        x <- make.names(x)
        return(x)
    }
    colnames(counts_data) <- clean_names(colnames(counts_data))
    colnames(abund_data)  <- clean_names(colnames(abund_data))

    design_full <- read.csv("${samplesheet}", stringsAsFactors=FALSE)
    design_full[,1] <- clean_names(design_full[,1])
    common_samples <- intersect(design_full[,1], colnames(counts_data))

    minimal_design <- data.frame(
        sampleID  = design_full[design_full[,1] %in% common_samples, 1],
        condition = make.names(design_full[design_full[,1] %in% common_samples, 3])
    )
    
    counts_data <- counts_data[, c("isoform_id", common_samples)]
    abund_data  <- abund_data[, c("isoform_id", common_samples)]
    counts_data[, common_samples] <- round(counts_data[, common_samples])

    if (length(unique(minimal_design\$condition)) >= 2) {
        switchList <- importRdata(
            isoformCountMatrix   = counts_data,
            isoformRepExpression = abund_data,
            designMatrix         = minimal_design,
            isoformExonAnnoation = "${iso_gtf}",
            showProgress         = FALSE,
            fixStringTieAnnotationProblem = TRUE
        )

        switchList <- preFilter(switchList)
        switchList <- isoformSwitchTestDEXSeq(switchList, reduceToSwitchingGenes = FALSE)
        switchList <- analyzeAlternativeSplicing(switchList, onlySwitchingGenes = FALSE)

        # === SAVE ALL AS EVENTS ===
        as_analysis <- switchList\$AlternativeSplicingAnalysis
        if (!is.null(as_analysis) && nrow(as_analysis) > 0) {
            message("Total AS events found: ", nrow(as_analysis))
            write.csv(as_analysis, "all_AS_events.csv", row.names=FALSE)
            
            # ATSS - any gene with ATSS > 0
            atss_all <- as_analysis[!is.na(as_analysis\$ATSS) & as_analysis\$ATSS > 0, ]
            message("Total ATSS events: ", nrow(atss_all))
            if (nrow(atss_all) > 0) {
                write.csv(atss_all, "atss_results.csv", row.names=FALSE)
            }
            
            # ATTS - any gene with ATTS > 0
            atts_all <- as_analysis[!is.na(as_analysis\$ATTS) & as_analysis\$ATTS > 0, ]
            message("Total ATTS events: ", nrow(atts_all))
            if (nrow(atts_all) > 0) {
                write.csv(atts_all, "atts_results.csv", row.names=FALSE)
            }
        } else {
            message("No AS analysis results found")
        }

        # [Keep all your PAS clustering code - starting from here]
        
        message("--- Clustering 3' Ends ---")
        exons_df <- as.data.frame(switchList\$exons)
        
        pas_raw <- do.call(rbind, lapply(split(exons_df, exons_df\$isoform_id), function(x) {
            data.frame(
                isoform_id = as.character(x\$isoform_id[1]),
                seqnames   = as.character(x\$seqnames[1]),
                strand     = as.character(x\$strand[1]),
                pas_coord  = if(x\$strand[1] == "+") max(x\$end) else min(x\$start),
                stringsAsFactors = FALSE
            )
        }))

        pas_gr <- makeGRangesFromDataFrame(pas_raw, start.field="pas_coord", end.field="pas_coord")
        pas_clusters <- reduce(resize(pas_gr, width=25, fix="center"))
        pas_clusters <- resize(pas_clusters, width=1, fix="center") 
        
        hits_df <- as.data.frame(findOverlaps(pas_gr, pas_clusters, maxgap=25))
        hits_df <- hits_df[!duplicated(hits_df\$queryHits), ]
        
        mapping_df <- data.frame(
            row_idx        = hits_df\$queryHits,
            pas_cluster_id = paste0("PAS_", hits_df\$subjectHits),
            cluster_coord  = start(pas_clusters[hits_df\$subjectHits]),
            stringsAsFactors = FALSE
        )

        pas_raw\$row_idx <- 1:nrow(pas_raw)
        pas_mapped <- merge(pas_raw, mapping_df, by="row_idx", all.x=FALSE)

        feat <- switchList\$isoformFeatures
        clean_id <- function(x) trimws(gsub('"', '', sub(".*:", "", as.character(x))))
        
        pas_mapped\$iso_clean <- clean_id(pas_mapped\$isoform_id)
        feat\$iso_clean       <- clean_id(feat\$isoform_id)

        idx <- match(pas_mapped\$iso_clean, feat\$iso_clean)
        pas_merged <- cbind(pas_mapped, feat[idx, ])
        pas_merged <- pas_merged[!is.na(idx), ]

        if_cols <- grep("^IF[0-9]+", colnames(feat), value=TRUE)
        val_cols <- grep("gene_value_[0-9]", colnames(feat), value=TRUE)
        if(length(val_cols) < 2) val_cols <- grep("condition_[0-9]", colnames(feat), value=TRUE)
        
        message("Using IF columns: ", paste(if_cols, collapse=", "))
        message("Using Expression columns: ", paste(val_cols, collapse=", "))

        pas_merged\$c1 <- pas_merged[[if_cols[1]]] * pas_merged[[val_cols[1]]]
        pas_merged\$c2 <- pas_merged[[if_cols[2]]] * pas_merged[[val_cols[2]]]
        
        group_vars <- c("gene_id", "seqnames", "cluster_coord", "strand", "pas_cluster_id")
        if(any(is.na(pas_merged\$gene_id))) pas_merged\$gene_id <- pas_merged\$iso_clean
        
        formula_str <- paste0("cbind(", paste(if_cols, collapse=","), ", c1, c2) ~ ", paste(group_vars, collapse=" + "))
        pas_final <- aggregate(as.formula(formula_str), data = pas_merged, FUN = sum, na.rm = TRUE)
        
        pas_final\$dPAU <- pas_final[[if_cols[2]]] - pas_final[[if_cols[1]]]

        message("Running statistical tests...")
        gene_ids <- unique(pas_final\$gene_id)
        pas_final\$pas_p_value <- NA_real_
        pas_final\$n_sites <- NA_integer_
        
        for (g in gene_ids) {
            sub <- pas_final[pas_final\$gene_id == g, ]
            pas_final\$n_sites[pas_final\$gene_id == g] <- nrow(sub)
            
            if (nrow(sub) > 1) {
                mat <- as.matrix(sub[, c("c1", "c2")])
                
                if (all(mat >= 0) && sum(mat) > 10) {
                    p <- tryCatch({
                        mat_int <- round(mat)
                        mat_int <- mat_int[rowSums(mat_int) > 0, , drop=FALSE]
                        
                        if (nrow(mat_int) > 1) {
                            fisher.test(mat_int, simulate.p.value = TRUE, B = 10000)\$p.value
                        } else {
                            NA_real_
                        }
                    }, error = function(e) NA_real_)
                    pas_final\$pas_p_value[pas_final\$gene_id == g] <- p
                }
            }
        }
        
        valid_p <- !is.na(pas_final\$pas_p_value)
        pas_final\$pas_q_value <- NA_real_
        if (sum(valid_p) > 0) {
            pas_final\$pas_q_value[valid_p] <- p.adjust(pas_final\$pas_p_value[valid_p], method = "BH")
        }
        
        message("Genes tested: ", sum(pas_final\$n_sites > 1, na.rm=TRUE))
        message("Significant APA events (q < 0.05): ", sum(pas_final\$pas_q_value < 0.05, na.rm=TRUE))
        
        write.csv(pas_final[order(pas_final\$pas_q_value, na.last=TRUE), ], 
                  "clustered_base_level_PAS_usage.csv", row.names=FALSE)
        
        # Create APA gene summary
        apa_genes <- aggregate(cbind(dPAU, pas_q_value) ~ gene_id + seqnames + strand, 
                               data=pas_final[pas_final\$pas_q_value < 0.05 & !is.na(pas_final\$pas_q_value), ],
                               FUN=function(x) x[which.min(abs(x))[1]])
        if (nrow(apa_genes) > 0) {
            write.csv(apa_genes, "APA_gene_summary.csv", row.names=FALSE)
        }

        saveRDS(switchList, "switchlist.rds")
        write.csv(feat, "isoform_switch_results.csv", row.names=FALSE)
        write.csv(feat[feat\$isoform_switch_q_value < 0.05, ], "significant_switches.csv", row.names=FALSE)
    }
    R_SCRIPT

    Rscript run_isoswitch.R
    """
}