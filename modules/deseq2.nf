process DESEQ2 {
    tag "deseq2"
    publishDir "${params.outdir}/DE_genes", mode: 'copy'

    input:
    path counts
    path samplesheet
    path gtf

    output:
    path "DE_results.csv"           , optional: true, emit: de_results
    path "deseq2_de_mqc.tsv"        , optional: true, emit: mqc_pca
    path "normalized_counts.csv"    , emit: norm_counts
    path "pca_plot.pdf"             , emit: pca_pdf

    script:
    """
    cat <<'R_SCRIPT' > run_deseq2.R
    library(DESeq2)
    library(ggplot2)

    # 1. Load Data
    counts_data <- read.table("${counts}", header=TRUE, row.names=1, check.names=FALSE)
    meta <- read.table("${samplesheet}", header=TRUE, sep=",", row.names=1)

    colnames(counts_data) <- gsub("^X", "", colnames(counts_data))
    colnames(counts_data) <- gsub(".sorted", "", colnames(counts_data))
    colnames(counts_data) <- gsub(".bam", "", colnames(counts_data))
    
    common_samples <- intersect(colnames(counts_data), rownames(meta))
    counts_data <- counts_data[, common_samples]
    meta <- meta[common_samples, , drop=FALSE]

    meta\$condition <- as.factor(meta\$condition)
    num_conditions <- length(unique(meta\$condition))

    # 2. ID Mapping
    gtf_df <- read.table("${gtf}", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
    gtf_tx <- gtf_df[gtf_df[[3]] == "transcript", ]
    get_attr <- function(attr_string, key) {
        val <- sub(paste0(".*", key, " \\"([^\\"]+)\\".*"), "\\\\1", attr_string)
        return(ifelse(grepl(key, attr_string), val, NA))
    }
    id_map <- unique(data.frame(
        id = get_attr(gtf_tx[[9]], "transcript_id"),
        gene_id = get_attr(gtf_tx[[9]], "gene_id"),
        stringsAsFactors = FALSE
    ))

    # 3. Initialize Object (Fix for single condition)
    design_formula <- if(num_conditions > 1) ~ condition else ~ 1
    dds <- DESeqDataSetFromMatrix(countData = round(counts_data), colData = meta, design = design_formula)

    # 4. PCA
    vsd <- vst(dds, blind=TRUE) 
    pdf("pca_plot.pdf", width=7, height=5)
    if(ncol(dds) > 1){
        print(plotPCA(vsd, intgroup="condition") + theme_minimal())
    } else { plot.new(); text(0.5,0.5,"PCA requires >1 sample") }
    dev.off()

    # 5. Analysis
    if (num_conditions > 1) {
        dds <- DESeq(dds)
        res <- as.data.frame(results(dds))
        res\$id <- rownames(res)
        final_res <- merge(id_map, res, by="id", all.y=TRUE)
        final_res <- final_res[order(final_res\$padj), ]
        write.csv(final_res, "DE_results.csv", row.names=FALSE)

    } else {
        dds <- estimateSizeFactors(dds)
        write.csv(data.frame(Message="Single condition - no DE"), "DE_results.csv", row.names=FALSE)
    }

    # 6. Export Normalized counts (includes all replicates)
    norm_counts_df <- as.data.frame(counts(dds, normalized=TRUE))
    norm_counts_df\$id <- rownames(norm_counts_df)
    write.csv(merge(id_map, norm_counts_df, by="id", all.y=TRUE), "normalized_counts.csv", row.names=FALSE)
R_SCRIPT
    Rscript run_deseq2.R
    """
}

process TE_DESEQ2 {
    tag "DE_TE_${type}"
    publishDir "${params.outdir}/TE/DE", mode: 'copy'

    input:
    path te_counts
    path gene_counts
    path samplesheet
    val type

    output:
    path "DE_TE_${type}_results.csv"    , optional: true, emit: de_results
    path "norm_counts_${type}.csv"      , emit: norm_counts
    path "volcano_${type}.pdf"          , optional: true, emit: volcano
    path "pca_${type}.pdf"              , emit: pca

    script:
    """
    cat <<'R_SCRIPT' > run_te_de.R
    library(DESeq2)
    library(ggplot2)

    # 1. LOAD & CLEAN
    te_data   <- read.table("${te_counts}", header=TRUE, sep="\t", check.names=FALSE)
    gene_data <- read.table("${gene_counts}", header=TRUE, row.names=1, check.names=FALSE)
    meta      <- read.table("${samplesheet}", header=TRUE, sep=",", row.names=1)

    if ("${type}" == "Family") {
        rownames(te_data) <- make.unique(as.character(te_data[,1]))
        te_meta_cols <- te_data[, 1, drop=FALSE]
        te_counts_mat <- te_data[, 2:ncol(te_data)]
    } else {
        unique_ids <- paste(te_data[,1], te_data[,2], te_data[,3], te_data[,4], sep="_")
        rownames(te_data) <- make.unique(unique_ids)
        te_meta_cols <- te_data[, 1:4]
        te_counts_mat <- te_data[, 5:ncol(te_data)]
    }
    te_ids <- rownames(te_counts_mat)
    te_meta_cols\$ID_JOIN <- te_ids

    clean_nm <- function(x) gsub(".bam|.sorted|^X", "", x)
    colnames(te_counts_mat) <- clean_nm(colnames(te_counts_mat))
    colnames(gene_data)     <- clean_nm(colnames(gene_data))

    common <- intersect(intersect(colnames(te_counts_mat), colnames(gene_data)), rownames(meta))
    te_counts_mat <- te_counts_mat[, common]; gene_data <- gene_data[, common]; meta <- meta[common,,drop=FALSE]

    # 2. NORMALIZATION
    combined <- rbind(as.matrix(gene_data), as.matrix(te_counts_mat))
    meta\$condition <- as.factor(meta\$condition)
    num_conds <- length(unique(meta\$condition))
    
    design_f <- if(num_conds > 1) ~ condition else ~ 1
    dds <- DESeqDataSetFromMatrix(countData = round(combined), colData = meta, design = design_f)

    if (num_conds > 1) {
        dds <- DESeq(dds)
        res <- as.data.frame(results(dds))[te_ids, ]
        res\$ID_JOIN <- rownames(res)
        write.csv(merge(te_meta_cols, res, by="ID_JOIN"), "DE_TE_${type}_results.csv", row.names=FALSE)

        pdf("volcano_${type}.pdf")
        print(ggplot(res, aes(x=log2FoldChange, y=-log10(pvalue), color=padj<0.05)) + geom_point(alpha=0.4) + theme_minimal())
        dev.off()
    } else {
        dds <- estimateSizeFactors(dds)
    }

    # 3. EXPORT ALL REPLICATES
    norm_mat <- counts(dds, normalized=TRUE)[te_ids, ]
    te_nc <- as.data.frame(norm_mat)
    te_nc\$ID_JOIN <- rownames(te_nc)
    write.csv(merge(te_meta_cols, te_nc, by="ID_JOIN"), "norm_counts_${type}.csv", row.names=FALSE)

    # 4. PCA
    vsd <- vst(dds, blind=TRUE)[te_ids, ]
    pdf("pca_${type}.pdf")
    if(ncol(dds) > 1) print(plotPCA(vsd, intgroup="condition")) else plot.new()
    dev.off()
R_SCRIPT
    Rscript run_te_de.R
    """
}