# NextLongIso: count SQANTI3 structural categories (FSM, ISM, NIC, NNC, ...) of the classified
# transcripts, for all transcripts and separately for annotated and novel ones. A transcript is
# annotated when its ID is a transcript_id of the reference GTF (Bambu keeps reference IDs).
# Usage: awk -f sqanti3_categories.awk fileno=1 reference.gtf fileno=2 <prefix>_classification.txt   (POSIX awk)
BEGIN {
    FS = "\t"; OFS = "\t"
    n = split("full-splice_match incomplete-splice_match novel_in_catalog novel_not_in_catalog genic antisense fusion intergenic genic_intron", order, " ")
    split("FSM ISM NIC NNC Genic_genomic Antisense Fusion Intergenic Genic_intron", abbr, " ")
    for (i = 1; i <= n; i++) short[order[i]] = abbr[i]
}
fileno == 1 {
    if ($0 ~ /^#/ || $3 != "transcript" && $3 != "exon") next
    if (match($9, /transcript_id "[^"]+"/)) ref[substr($9, RSTART + 15, RLENGTH - 16)] = 1
    next
}
fileno == 2 && FNR == 1 {
    for (i = 1; i <= NF; i++) { if ($i == "isoform") ic = i; if ($i == "structural_category") cc = i }
    if (!ic || !cc) { print "sqanti3_categories.awk: isoform/structural_category columns not found" > "/dev/stderr"; bad = 1; exit 1 }
    next
}
fileno == 2 {
    c = $cc
    if (!(c in all)) { extra[++ne] = c; all[c] = 0; annotated[c] = 0; novel[c] = 0 }
    all[c]++; tot++
    if ($ic in ref) { annotated[c]++; tot_a++ } else { novel[c]++; tot_n++ }
}
END {
    if (bad) exit 1
    print "structural_category", "abbreviation", "all_transcripts", "annotated_transcripts", "novel_transcripts"
    for (i = 1; i <= n; i++) {
        c = order[i]
        print c, short[c], all[c] + 0, annotated[c] + 0, novel[c] + 0
        done[c] = 1
    }
    for (i = 1; i <= ne; i++) if (!(extra[i] in done)) print extra[i], extra[i], all[extra[i]], annotated[extra[i]], novel[extra[i]]
    print "total", "total", tot + 0, tot_a + 0, tot_n + 0
}
