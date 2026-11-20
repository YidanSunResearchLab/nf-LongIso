# NextLongIso: combine SUPPA2 diffSplice .dpsi files into one table of significant events
# (p-value < alpha and |dPSI| > dpsi). A .dpsi file has a two-field header (<cond1>-<cond2>_dPSI,
# <cond1>-<cond2>_p-val) and data rows event_id, dPSI (cond2 - cond1), p-value; nan values are skipped.
# The comparison and event type are passed per file with command-line assignments.
# Usage: awk -v alpha=0.05 -v dpsi=0.1 -f suppa2_significant.awk \
#            cmp=treatment_vs_control type=SE a.dpsi cmp=... type=... b.dpsi ...   (POSIX awk)
BEGIN {
    FS = OFS = "\t"
    if (alpha == "") alpha = 0.05
    if (dpsi == "") dpsi = 0.1
    print "comparison", "event_type", "event_id", "gene_id", "dPSI", "p_value"
}
FNR == 1 { next }
$2 == "nan" || $3 == "nan" || $2 == "" || $3 == "" { next }
{
    d = $2 + 0; p = $3 + 0
    if (p < alpha && (d > dpsi || -d > dpsi)) {
        gene = $1; sub(/;.*$/, "", gene)
        print cmp, type, $1, gene, $2, $3
    }
}
