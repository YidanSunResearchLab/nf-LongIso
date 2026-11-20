# NextLongIso: keep the columns of one condition from a SUPPA2-style matrix (PSI or expression file:
# header = sample names only; data rows = ID + one value per sample). Sample-to-condition assignments
# come from the design table (sample,condition[,covariates]; comma-separated, with a header).
# Columns keep the order of the matrix, so PSI and expression files split from matrices with the same
# sample order stay aligned.
# Usage: awk -v cond=treatment -f split_condition.awk fileno=1 design.csv fileno=2 matrix.tsv   (POSIX awk)
BEGIN { FS = OFS = "\t" }
fileno == 1 {
    n = split($0, f, ",")
    if (FNR == 1) { for (i = 1; i <= n; i++) { if (f[i] == "sample") sc = i; if (f[i] == "condition") cc = i }; next }
    if (f[cc] == cond) { want[f[sc]] = 1; nwant++ }
    next
}
fileno == 2 && FNR == 1 {
    if (nwant == 0) { print "split_condition.awk: no samples for condition " cond " in the design" > "/dev/stderr"; bad = 1; exit 1 }
    out = ""; k = 0
    for (i = 1; i <= NF; i++) if ($i in want) { col[++k] = i + 1; out = out (k > 1 ? OFS : "") $i; seen[$i] = 1 }
    for (s in want) if (!(s in seen)) { print "split_condition.awk: sample " s " (condition " cond ") not found in the matrix header" > "/dev/stderr"; bad = 1 }
    if (bad) exit 1
    print out
    next
}
fileno == 2 { out = $1; for (j = 1; j <= k; j++) out = out OFS $(col[j]); print out }
END { if (bad) exit 1 }
