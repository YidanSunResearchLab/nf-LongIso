#!/usr/bin/env Rscript
# NextLongIso: isoform switching (IsoformSwitchAnalyzeR) and alternative TSS/TTS usage from Bambu
# isoforms. (Poly(A)-site usage is analysed by the APA process: apa_sites.R, usage_dexseq.R.)
#
# Usage:
#   run_isoswitch.R --counts=counts_transcript.txt --abundance=CPM_transcript.txt \
#               --gtf=extended_annotations.gtf --design=design.csv --design-formula="~ condition" \
#               --reference=control \
#               --alpha=0.05 --dif-cutoff=0.1 [--nt-fasta=transcripts.fa --n-plots=10]
#
# With --nt-fasta (transcript sequences, e.g. from gffread), ORFs are predicted (analyzeORF, longest
# ORF, no BSgenome needed) and switch consequences are analysed (types that need no external tools).
#
# Every non-reference condition is compared with the reference (condition_1 = reference,
# condition_2 = other condition; dIF and dPAU are condition_2 - condition_1).
#
# Outputs
#   isoform_switch_results.csv        isoformFeatures for all comparisons (required)
#   significant_switches.csv          q < alpha and |dIF| > dif-cutoff (required)
#   all_AS_events.csv, atss_results.csv, atts_results.csv
#   switchlist.rds
#   With --nt-fasta: switch_consequences.csv, consequence_summary.csv, consequence_enrichment.csv,
#   consequence_summary.pdf, consequence_enrichment.pdf, switch_plots.pdf (top --n-plots switches);
#   files contain a 'Message' column (PDFs a text page) when there is nothing to report.
# With fewer than 2 conditions the two required files contain a single 'Message' column.

.libPaths(unique(c(.libPaths(), "/usr/local/lib/R/site-library", "/usr/lib/R/library")))
suppressPackageStartupMessages({
    library(IsoformSwitchAnalyzeR)
    library(GenomicRanges)
})

# Shared helpers (parse_args, design handling) from bin/nli_design.R
.script <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "nli_design.R"))
opt <- parse_args(commandArgs(trailingOnly = TRUE))
for (k in c("counts", "abundance", "gtf", "design", "design-formula", "reference", "alpha", "dif-cutoff")) {
    if (is.null(opt[[k]])) stop(sprintf("Missing required argument --%s", k))
}
alpha      <- as.numeric(opt$alpha)
dif_cutoff <- as.numeric(opt[["dif-cutoff"]])

strip_suffix  <- function(x) sub("[.]sorted([.]bam)?$", "", x)
write_message <- function(file, msg) write.csv(data.frame(Message = msg), file, row.names = FALSE)

# ---- 1. Load ---------------------------------------------------------------------------------
counts_data <- read.delim(opt$counts,    check.names = FALSE, stringsAsFactors = FALSE)
abund_data  <- read.delim(opt$abundance, check.names = FALSE, stringsAsFactors = FALSE)
colnames(counts_data)[1] <- "isoform_id"
colnames(abund_data)[1]  <- "isoform_id"
colnames(counts_data) <- strip_suffix(colnames(counts_data))
colnames(abund_data)  <- strip_suffix(colnames(abund_data))

design  <- read.csv(opt$design, stringsAsFactors = FALSE)
missing <- setdiff(design$sample, intersect(colnames(counts_data), colnames(abund_data)))
if (length(missing) > 0) stop(sprintf("Samples missing from the count tables: %s", paste(missing, collapse = ", ")))

conditions <- unique(design$condition)
if (length(conditions) < 2) {
    msg <- sprintf("Isoform switch analysis skipped: the samplesheet has %d condition(s); at least 2 are required.",
                   length(conditions))
    message(msg)
    write_message("isoform_switch_results.csv", msg)
    write_message("significant_switches.csv", msg)
    quit(save = "no", status = 0)
}
reference <- opt$reference
if (!reference %in% conditions) stop(sprintf("Reference condition '%s' not found in the design", reference))

# The DEXSeq-based switch test needs biological replicates: only conditions with >= 2 samples
# are compared, and the reference itself must have >= 2 samples.
n_rep  <- table(design$condition)
single <- names(n_rep)[n_rep < 2]
others <- sort(setdiff(names(n_rep)[n_rep >= 2], reference))
if (n_rep[[reference]] < 2 || length(others) == 0) {
    msg <- sprintf(paste0("Isoform switch analysis skipped: it needs >= 2 replicates in the reference ",
                          "condition ('%s', %d sample(s)) and in at least one other condition. ",
                          "Conditions with < 2 replicates: %s."),
                   reference, n_rep[[reference]], if (length(single)) paste(single, collapse = ", ") else "none")
    message(msg)
    write_message("isoform_switch_results.csv", msg)
    write_message("significant_switches.csv", msg)
    quit(save = "no", status = 0)
}
if (length(single) > 0) {
    message("Conditions with < 2 replicates are excluded from isoform switching: ", paste(single, collapse = ", "))
}
design <- design[design$condition %in% c(reference, others), ]
covariates <- design_covariates(opt[["design-formula"]])
if (length(covariates) > 0) {
    # Covariates are passed as extra designMatrix columns (used by the DEXSeq-based switch test)
    cd <- design_coldata(design, opt[["design-formula"]], reference)
    check_full_rank(design_formula(opt[["design-formula"]], cd), cd)
    message("Covariates: ", paste(covariates, collapse = ", "))
}
comparisons <- data.frame(condition_1 = reference, condition_2 = others, stringsAsFactors = FALSE)
message("Comparisons: ", paste(comparisons$condition_2, "vs", comparisons$condition_1, collapse = "; "))

counts_data <- counts_data[, c("isoform_id", design$sample)]
abund_data  <- abund_data[,  c("isoform_id", design$sample)]
counts_data[, design$sample] <- round(counts_data[, design$sample])

# ---- 2. Isoform switching and alternative splicing ---------------------------------------------
switchList <- importRdata(
    isoformCountMatrix   = counts_data,
    isoformRepExpression = abund_data,
    designMatrix         = data.frame(sampleID = design$sample, condition = design$condition,
                                      design[, covariates, drop = FALSE], stringsAsFactors = FALSE),
    isoformExonAnnoation = opt$gtf,
    isoformNtFasta       = opt[["nt-fasta"]],
    comparisonsToMake    = comparisons,
    showProgress         = FALSE,
    fixStringTieAnnotationProblem = TRUE
)
switchList <- preFilter(switchList)
switchList <- isoformSwitchTestDEXSeq(switchList, reduceToSwitchingGenes = FALSE)
switchList <- analyzeAlternativeSplicing(switchList, onlySwitchingGenes = FALSE)

as_analysis <- switchList$AlternativeSplicingAnalysis
if (!is.null(as_analysis) && nrow(as_analysis) > 0) {
    message("Total AS events: ", nrow(as_analysis))
    write.csv(as_analysis, "all_AS_events.csv", row.names = FALSE)
    atss <- as_analysis[!is.na(as_analysis$ATSS) & as_analysis$ATSS > 0, ]
    atts <- as_analysis[!is.na(as_analysis$ATTS) & as_analysis$ATTS > 0, ]
    message("ATSS events: ", nrow(atss), "; ATTS events: ", nrow(atts))
    if (nrow(atss) > 0) write.csv(atss, "atss_results.csv", row.names = FALSE)
    if (nrow(atts) > 0) write.csv(atts, "atts_results.csv", row.names = FALSE)
} else {
    message("No alternative splicing results")
}

feat <- switchList$isoformFeatures
clean_id <- function(x) trimws(gsub('"', '', sub(".*:", "", as.character(x))))
feat$iso_clean <- clean_id(feat$isoform_id)

# ---- 3. Switch consequences (with transcript sequences) ------------------------------------------
CONSEQUENCES <- c("intron_retention", "ORF_seq_similarity", "NMD_status", "tss", "tts",
                  "last_exon", "isoform_length", "exon_number")
message_pdf <- function(file, msg) { pdf(file, width = 8, height = 4); plot.new(); text(0.5, 0.5, msg, cex = 0.9); invisible(dev.off()) }
if (!is.null(opt[["nt-fasta"]])) {
    n_plots <- as.integer(if (is.null(opt[["n-plots"]])) 10 else opt[["n-plots"]])
    switchList <- analyzeORF(switchList, orfMethod = "longest", showProgress = FALSE)
    message("ORFs predicted: ", sum(!is.na(switchList$orfAnalysis$orfTransciptStart)))
    # ORF amino-acid sequences (needed for ORF_seq_similarity); uses the transcript sequences already
    # in the switch list, no genome object
    switchList <- tryCatch(
        extractSequence(switchList, onlySwitchingGenes = FALSE, alpha = alpha, dIFcutoff = dif_cutoff,
                        extractNTseq = TRUE, extractAAseq = TRUE, addToSwitchAnalyzeRlist = TRUE,
                        writeToFile = FALSE, quiet = TRUE),
        error = function(e) { message("extractSequence: ", conditionMessage(e)); switchList })
    switchList <- tryCatch(
        analyzeSwitchConsequences(switchList, consequencesToAnalyze = CONSEQUENCES, alpha = alpha,
                                  dIFcutoff = dif_cutoff, showProgress = FALSE),
        error = function(e) { message("analyzeSwitchConsequences: ", conditionMessage(e)); switchList })
    cons <- switchList$switchConsequence
    if (!is.null(cons) && nrow(cons) > 0) {
        write.csv(cons, "switch_consequences.csv", row.names = FALSE)
        summ <- extractConsequenceSummary(switchList, consequencesToAnalyze = "all", alpha = alpha, dIFcutoff = dif_cutoff,
                                          plot = FALSE, returnResult = TRUE)
        write.csv(summ, "consequence_summary.csv", row.names = FALSE)
        pdf("consequence_summary.pdf", width = 10, height = 6)
        print(extractConsequenceSummary(switchList, consequencesToAnalyze = "all", alpha = alpha, dIFcutoff = dif_cutoff,
                                        plot = TRUE, returnResult = FALSE))
        invisible(dev.off())
        enr <- tryCatch(extractConsequenceEnrichment(switchList, consequencesToAnalyze = "all", alpha = alpha,
                                                     dIFcutoff = dif_cutoff, plot = FALSE, returnResult = TRUE,
                                                     returnSummary = TRUE),
                        error = function(e) data.frame(Message = conditionMessage(e)))
        write.csv(enr, "consequence_enrichment.csv", row.names = FALSE)
        pdf("consequence_enrichment.pdf", width = 10, height = 6)
        tryCatch(print(extractConsequenceEnrichment(switchList, consequencesToAnalyze = "all", alpha = alpha,
                                                    dIFcutoff = dif_cutoff, plot = TRUE, returnResult = FALSE)),
                 error = function(e) { plot.new(); text(0.5, 0.5, conditionMessage(e)) })
        invisible(dev.off())
        message("Switch consequences: ", nrow(cons), " rows")
    } else {
        msg <- sprintf("No isoform switches with q < %g and |dIF| > %g: no switch consequences to report.", alpha, dif_cutoff)
        message(msg)
        for (f in c("switch_consequences.csv", "consequence_summary.csv", "consequence_enrichment.csv")) {
            write.csv(data.frame(Message = msg), f, row.names = FALSE)
        }
        message_pdf("consequence_summary.pdf", msg); message_pdf("consequence_enrichment.pdf", msg)
    }

    # switchPlot for the top switching genes (smallest isoform switch q-value)
    f2 <- switchList$isoformFeatures
    top <- f2[which(f2$isoform_switch_q_value < alpha & abs(f2$dIF) > dif_cutoff), ]
    top <- top[order(top$isoform_switch_q_value, -abs(top$dIF)), ]
    top <- top[!duplicated(paste(top$gene_id, top$condition_1, top$condition_2)), ]
    top <- head(top, n_plots)
    pdf("switch_plots.pdf", width = 12, height = 7)
    if (nrow(top) == 0) {
        plot.new(); text(0.5, 0.5, sprintf("No isoform switches with q < %g and |dIF| > %g", alpha, dif_cutoff))
    } else {
        # switchPlot measures text before its own grid.newpage(), which leaves a blank first page on a
        # fresh device; capture each plot off-screen and draw it on its own page instead.
        for (i in seq_len(nrow(top))) {
            tryCatch({
                g <- grid::grid.grabExpr(switchPlot(switchList, gene = top$gene_id[i], condition1 = top$condition_1[i],
                                                    condition2 = top$condition_2[i], alphas = c(alpha, alpha / 50),
                                                    dIFcutoff = dif_cutoff, plotTopology = FALSE),
                                         wrap.grobs = TRUE, width = 12, height = 7)
                grid::grid.newpage(); grid::grid.draw(g)
            }, error = function(e) { plot.new(); text(0.5, 0.5, paste(top$gene_id[i], conditionMessage(e))) })
        }
    }
    invisible(dev.off())
    message("switch_plots.pdf: ", nrow(top), " gene(s)")
}

# ---- 4. Switch tables --------------------------------------------------------------------------
saveRDS(switchList, "switchlist.rds")
write.csv(feat, "isoform_switch_results.csv", row.names = FALSE)
sig_switch <- feat[which(feat$isoform_switch_q_value < alpha & abs(feat$dIF) > dif_cutoff), ]
write.csv(sig_switch, "significant_switches.csv", row.names = FALSE)
message(sprintf("Significant switches (q < %g, |dIF| > %g): %d", alpha, dif_cutoff, nrow(sig_switch)))
