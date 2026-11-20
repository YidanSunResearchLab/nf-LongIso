# NextLongIso: convert Bambu counts_transcript.txt (TXNAME, GENEID, samples...) into a SUPPA2
# expression file: header = sample names only; rows = transcript ID + values.
# Transcript IDs are kept exactly as in the Bambu GTF (including versions such as
# ENST00000456328.2), so that they match the transcript IDs in the .ioe event files.
# Usage: awk -f suppa2_expression.awk counts_transcript.txt > expression.tsv   (POSIX awk)
BEGIN { FS = OFS = "\t" }
NR == 1 {
    out = ""
    for (i = 3; i <= NF; i++) { n = $i; sub(/[.]sorted([.]bam)?$/, "", n); out = out (i > 3 ? OFS : "") n }
    print out
    next
}
{ out = $1; for (i = 3; i <= NF; i++) out = out OFS $i; print out }
