#!/usr/bin/env Rscript
# NextLongIso: transcript discovery and quantification with Bambu.
#
# Usage:
#   run_bambu.R --annotations=bambu_annotations.rds --genome=genome.fa --ncore=16 \
#               --extra-args-file=bambu_args.txt
#
# BAM files (*.bam) are read from the working directory. --extra-args-file contains extra R
# arguments for bambu() (e.g. NDR = 0.1); they are appended to the call exactly as written.
#
# Outputs: counts_transcript.txt, counts_gene.txt, CPM_transcript.txt (and the other
# writeBambuOutput() tables), extended_annotations.gtf, bambu_results.rds

parse_args <- function(args) {
    # Base-R parser for --key=value arguments (no optparse dependency)
    out <- list()
    for (a in args) {
        if (!grepl("^--[A-Za-z0-9_-]+=", a)) stop(sprintf("Arguments must be given as --key=value (got '%s')", a))
        out[[sub("^--([^=]+)=.*$", "\\1", a)]] <- sub("^--[^=]+=", "", a)
    }
    out
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
for (k in c("annotations", "genome", "ncore", "extra-args-file")) {
    if (is.null(opt[[k]])) stop(sprintf("Missing required argument --%s", k))
}

library(bambu)
library(Biostrings)
library(rtracklayer)

print("Loading genome...")
genome_seq <- readDNAStringSet(opt$genome)
names(genome_seq) <- sub(" .*", "", names(genome_seq))

print("Setting up BAM files...")
all_files <- list.files(path = ".", full.names = TRUE)
bam_paths <- all_files[endsWith(all_files, ".bam")]
if (length(bam_paths) == 0) stop("ERROR: No BAM files found!")
print(paste("Found", length(bam_paths), "BAM files:"))
print(bam_paths)

extra <- trimws(paste(readLines(opt[["extra-args-file"]], warn = FALSE), collapse = " "))
print(paste("Extra bambu() arguments:", if (nzchar(extra)) extra else "none"))

print("Running Bambu...")
annotations <- readRDS(opt$annotations)
ncore <- as.integer(opt$ncore)
# The call is built as text so that the extra arguments are evaluated exactly as written by the user
bambu_call <- sprintf(paste0(
    "bambu(reads = bam_paths, annotations = annotations, genome = genome_seq, discovery = TRUE, ",
    "quant = TRUE, ncore = ncore, lowMemory = TRUE, verbose = FALSE%s)"),
    if (nzchar(extra)) paste0(", ", extra) else "")
se <- eval(parse(text = bambu_call))

print("Writing output files...")
# writeBambuOutput writes counts_transcript.txt, counts_gene.txt, CPM_transcript.txt, CPM_gene.txt, ...
writeBambuOutput(se, path = ".")

print("Saving extended GTF...")
writeToGTF(rowRanges(se), file = "extended_annotations.gtf")

# writeBambuOutput should produce CPM_transcript.txt; older Bambu versions may not
if (!file.exists("CPM_transcript.txt")) {
    print("CPM_transcript.txt not found - writing manually...")
    se_tx   <- se[["transcript"]]
    cpm_out <- cbind(data.frame(TXNAME = rowData(se_tx)$TXNAME, GENEID = rowData(se_tx)$GENEID),
                     as.data.frame(assay(se_tx, "CPM")))
    write.table(cpm_out, "CPM_transcript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
    print(paste("Written CPM_transcript.txt:", nrow(cpm_out), "rows"))
}

for (f in c("counts_transcript.txt", "counts_gene.txt", "CPM_transcript.txt", "extended_annotations.gtf")) {
    if (file.exists(f)) print(paste("OK:", f)) else stop(paste("ERROR: expected output missing:", f))
}

saveRDS(se, "bambu_results.rds")
print("Bambu completed successfully!")
