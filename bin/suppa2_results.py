#!/usr/bin/env python3
"""NextLongIso: combine SUPPA2 diffSplice results and add a genome-wide FDR.

Usage:
    suppa2_results.py --alpha 0.05 --dpsi 0.1 <comparison>=<file.dpsi> [<comparison>=<file.dpsi> ...]

A .dpsi file has a two-field header (<cond1>-<cond2>_dPSI, <cond1>-<cond2>_p-val) and data rows
event_id, dPSI (cond2 - cond1), p-value. With diffSplice -gc the p-value is already corrected within
each gene (Benjamini-Hochberg per gene), which does not control the false discovery rate over all genes.
padj adds Benjamini-Hochberg over all events with a p-value of the comparison.

Outputs (current directory):
    suppa2_diffsplice_results.tsv      every event with a p-value: comparison, event_type, event_id,
                                       gene_id, dPSI, pvalue (as reported by SUPPA2), padj
    suppa2_diffsplice_significant.tsv  padj < alpha and |dPSI| > dpsi
Standard library only.
"""
import argparse
import math

HEADER = ['comparison', 'event_type', 'event_id', 'gene_id', 'dPSI', 'pvalue', 'padj']


def number(x):
    try:
        v = float(x)
    except ValueError:
        return None
    return None if math.isnan(v) else v


def bh(pvalues):
    """Benjamini-Hochberg adjusted p-values (same order as the input)."""
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    adj = [0.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        running = min(running, pvalues[i] * n / rank)
        adj[i] = running
    return adj


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--alpha', type=float, default=0.05)
    ap.add_argument('--dpsi', type=float, default=0.1)
    ap.add_argument('dpsi_files', nargs='*', help='<comparison>=<file.dpsi>')
    a = ap.parse_args()

    rows = []
    for item in a.dpsi_files:
        comparison, path = item.split('=', 1)
        events = []
        with open(path) as fh:
            next(fh, None)                                  # header
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) < 3:
                    continue
                d, p = number(f[1]), number(f[2])
                if d is None or p is None:
                    continue
                event = f[0]
                gene, _, rest = event.partition(';')
                events.append([comparison, rest.split(':', 1)[0], event, gene, f[1], f[2], p])
        for e, q in zip(events, bh([e[6] for e in events])):
            e[6] = f'{q:.6g}'
        rows.extend(events)
        print(f'{comparison}: {len(events)} events with a p-value')

    sig = [r for r in rows if float(r[6]) < a.alpha and abs(float(r[4])) > a.dpsi]
    for name, data in (('suppa2_diffsplice_results.tsv', rows), ('suppa2_diffsplice_significant.tsv', sig)):
        with open(name, 'w') as out:
            out.write('\t'.join(HEADER) + '\n')
            for r in data:
                out.write('\t'.join(str(x) for x in r) + '\n')
    print(f'Significant events (padj < {a.alpha}, |dPSI| > {a.dpsi}): {len(sig)}')


if __name__ == '__main__':
    main()
