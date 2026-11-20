# NextLongIso: keep SUPPA2 events with a PSI value (not nan/NA) in at least `min` samples.
# PSI files have a header of sample names only; data rows are event_id + one column per sample.
# Usage: awk -v min=2 -f psi_min_samples.awk events_SE_samples.psi   (POSIX awk)
BEGIN { FS = "\t"; if (min == "") min = 2 }
NR == 1 { print; next }
{ n = 0; for (i = 2; i <= NF; i++) if ($i != "nan" && $i != "NA" && $i != "") n++; if (n >= min) print }
