#!/usr/bin/env python3
"""Compare TE quantification modes (exon vs span) of two NextLongIso runs on the same data.

Usage:
    compare_te_modes.py --exon <outdir_exon> --span <outdir_span> --host-pairs host_pairs.tsv

--span-workdir (optional): work directory of the span-mode TE_ANALYSIS task; used to report how
much of the span-mode count comes from TEs that do not overlap any exon of the counted transcript.

--host-pairs: two columns (gene_id, subfamily), one row per gene whose genomic span overlaps an
instance of the subfamily on the same strand. It defines "host genes" independently of both
counting modes; docs/te_count_mode_comparison.md shows the bedtools command used to create it.

Prints a Markdown report.
"""
import argparse
import csv
import math
from collections import Counter, defaultdict


def read_counts(path, n_id):
    """TE count table -> {id: total count over samples}, plus number of samples."""
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out["\t".join(f[:n_id])] = sum(float(x) for x in f[n_id:])
    return out, len(header) - n_id


def ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    r = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def spearman(x, y):
    if len(x) < 3:
        return float("nan")
    rx, ry = ranks(x), ranks(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    vx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    vy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return cov / (vx * vy) if vx and vy else float("nan")


def read_de(path, id_cols):
    rows = list(csv.DictReader(open(path)))
    for r in rows:
        try:
            r["_padj"] = float(r["padj"])
        except (ValueError, KeyError):
            r["_padj"] = float("nan")
        r["_id"] = "_".join(r[c] for c in id_cols)
    sig = [r for r in rows if not math.isnan(r["_padj"]) and r["_padj"] < 0.05]
    top = sorted([r for r in rows if not math.isnan(r["_padj"])], key=lambda r: (r["_padj"], r["_id"]))[:20]
    return rows, sig, top


def fmt(v, d=3):
    return "NA" if v is None or (isinstance(v, float) and math.isnan(v)) else f"{v:.{d}f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exon", required=True)
    ap.add_argument("--span", required=True)
    ap.add_argument("--host-pairs", required=True)
    ap.add_argument("--span-workdir", default=None)
    a = ap.parse_args()

    fam = {m: read_counts(f"{d}/TE/TE_family_counts.tsv", 1)[0] for m, d in (("exon", a.exon), ("span", a.span))}
    inst = {m: read_counts(f"{d}/TE/TE_instance_counts.tsv", 4)[0] for m, d in (("exon", a.exon), ("span", a.span))}

    print("## 1. Counted features and totals\n")
    print("| Metric | exon | span |\n|---|---|---|")
    for label, tab in (("subfamilies", fam), ("TE instances", inst)):
        print(f"| {label} with counts > 0 | {sum(v > 0 for v in tab['exon'].values()):,} | {sum(v > 0 for v in tab['span'].values()):,} |")
        print(f"| total counts over {label} (all samples) | {sum(tab['exon'].values()):,.0f} | {sum(tab['span'].values()):,.0f} |")
    both_f = set(k for k, v in fam["exon"].items() if v > 0) & set(k for k, v in fam["span"].items() if v > 0)
    print(f"| subfamilies counted in both modes | {len(both_f):,} | |")
    only_e = sum(1 for k, v in fam["exon"].items() if v > 0 and fam["span"].get(k, 0) == 0)
    only_s = sum(1 for k, v in fam["span"].items() if v > 0 and fam["exon"].get(k, 0) == 0)
    print(f"| subfamilies counted in one mode only | {only_e:,} (exon only) | {only_s:,} (span only) |")

    print("\n## 2. Spearman correlation between modes\n")
    print("Counts summed over samples. 'Union' sets missing features to 0.\n")
    print("| Level | Features (union) | rho (union) | Features (both) | rho (both) |\n|---|---|---|---|---|")
    for label, tab in (("subfamily", fam), ("instance", inst)):
        union = sorted(set(tab["exon"]) | set(tab["span"]))
        both = sorted(k for k in union if tab["exon"].get(k, 0) > 0 and tab["span"].get(k, 0) > 0)
        ru = spearman([tab["exon"].get(k, 0) for k in union], [tab["span"].get(k, 0) for k in union])
        rb = spearman([tab["exon"][k] for k in both], [tab["span"][k] for k in both])
        print(f"| {label} | {len(union):,} | {fmt(ru)} | {len(both):,} | {fmt(rb)} |")

    # Host-gene expression per subfamily (gene-level Bambu counts; same in both runs)
    gene_expr = {}
    with open(f"{a.exon}/bambu/counts_gene.txt") as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            gene_expr[f[0]] = sum(float(x) for x in f[1:])
    host = defaultdict(set)
    for line in open(a.host_pairs):
        g, s = line.rstrip("\n").split("\t")[:2]
        host[s].add(g)
    host_expr = {s: sum(gene_expr.get(g, 0.0) for g in gs) for s, gs in host.items()}

    print("\n## 3. TE counts vs host-gene expression\n")
    print("Host genes of a subfamily: genes whose span overlaps one of its instances on the same strand "
          "(defined independently of both counting modes). Columns 2-3: all subfamilies with host-gene "
          "expression > 0 (TE count 0 included). Columns 4-6: only subfamilies counted (> 0) in that mode; "
          "ratio = TE count / host-gene expression (a ratio near 1 means the TE count is essentially the "
          "expression of its host genes).\n")
    print("| Mode | Subfamilies (host > 0) | rho | Counted subfamilies | rho | median ratio |\n|---|---|---|---|---|---|")
    for m in ("exon", "span"):
        keys = [k for k in host_expr if host_expr[k] > 0]
        rho_all = spearman([fam[m].get(k, 0.0) for k in keys], [host_expr[k] for k in keys])
        cnt = [k for k in keys if fam[m].get(k, 0.0) > 0]
        rho_cnt = spearman([fam[m][k] for k in cnt], [host_expr[k] for k in cnt])
        ratios = sorted(fam[m][k] / host_expr[k] for k in cnt)
        med = ratios[len(ratios) // 2] if ratios else float("nan")
        print(f"| {m} | {len(keys):,} | {fmt(rho_all)} | {len(cnt):,} | {fmt(rho_cnt)} | {fmt(med)} |")

    if a.span_workdir:
        exonic = set()
        for line in open(f"{a.span_workdir}/exon_te_overlap.txt"):
            f = line.rstrip("\n").split("\t")
            exonic.add((f[3], f[12]))
        n = n_ex = 0
        w = w_ex = 0.0
        for line in open(f"{a.span_workdir}/all_joined.tsv"):
            f = line.rstrip("\n").split("\t")
            c = sum(float(x) for x in f[9:])
            n += 1
            w += c
            if (f[0], f[5]) in exonic:
                n_ex += 1
                w_ex += c
        print("\nSpan mode assigns each transcript to one TE overlapping its span. For these assignments:\n")
        print("| | Assignments | Share of span-mode counts |\n|---|---|---|")
        print(f"| TE overlaps an exon of the transcript | {n_ex:,} ({100 * n_ex / n:.1f}%) | {100 * w_ex / w:.1f}% |")
        print(f"| TE only in introns of the transcript | {n - n_ex:,} ({100 * (n - n_ex) / n:.1f}%) | {100 * (w - w_ex) / w:.1f}% |")

    print("\n## 4. TE differential expression (padj < 0.05)\n")
    print("| Level | exon | span | shared |\n|---|---|---|---|")
    tops = {}
    for level, cols in (("Family", ["TE_Family"]), ("Instance", ["TE_id", "TE_chr", "TE_start", "TE_end"])):
        res = {m: read_de(f"{d}/TE/DE/DE_TE_{level}_results.csv", cols) for m, d in (("exon", a.exon), ("span", a.span))}
        se = {r["_id"] for r in res["exon"][1]}
        ss = {r["_id"] for r in res["span"][1]}
        print(f"| {level} | {len(se)} | {len(ss)} | {len(se & ss)} |")
        tops[level] = res
    print("\nTop 20 subfamilies by padj:\n")
    print("| Rank | exon | padj | span | padj |\n|---|---|---|---|---|")
    te, ts = tops["Family"]["exon"][2], tops["Family"]["span"][2]
    for i in range(max(len(te), len(ts))):
        e = te[i] if i < len(te) else None
        s = ts[i] if i < len(ts) else None
        print(f"| {i + 1} | {e['_id'] if e else ''} | {('%.2e' % e['_padj']) if e else ''} | "
              f"{s['_id'] if s else ''} | {('%.2e' % s['_padj']) if s else ''} |")
    shared = {r["_id"] for r in te} & {r["_id"] for r in ts}
    print(f"\nSubfamilies in both top-20 lists: {len(shared)}")

    print("\n## 5. isoform_TE_annotation.tsv categories (exon-overlap definition; identical in both runs)\n")
    cats = Counter()
    tx = set()
    with open(f"{a.exon}/TE/isoform_TE_annotation.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cats[r["category"]] += 1
            tx.add(r["transcript_id"])
    total = sum(cats.values())
    print("| Category | (isoform, TE) pairs | % |\n|---|---|---|")
    for c in ("autonomous", "TSS", "PAS", "internal_exonization", "other_exonic"):
        print(f"| {c} | {cats.get(c, 0):,} | {100 * cats.get(c, 0) / max(total, 1):.1f} |")
    print(f"| **total** | {total:,} ({len(tx):,} isoforms) | 100 |")


if __name__ == "__main__":
    main()
