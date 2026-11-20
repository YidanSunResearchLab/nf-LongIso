#!/usr/bin/env python3
"""Summarise a NextLongIso output directory as Markdown.

Used to record the v0.1.0 baseline (docs/dev/baseline_v0.1.0.md) and to compare
v0.2.0 against it. Works on outputs of both versions; missing files are reported
as "-".

Usage:
    python3 tests/summarize_results.py <results_dir> [--title TITLE] [--files]

--files additionally lists every file with its line count (text) or size (binary).
"""
import argparse
import csv
import glob
import math
import os
import re

TEXT_EXT = ('.txt', '.tsv', '.csv', '.gtf', '.ioe', '.psi', '.bed', '.stats', '.flagstat', '.yml', '.yaml')


def n_lines(path):
    with open(path, 'rb') as fh:
        return sum(1 for _ in fh)


def read_csv(path, delim=','):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def to_float(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return None


def count_sig(rows, col, alpha=0.05):
    return sum(1 for r in rows if (v := to_float(r.get(col))) is not None and v < alpha)


def metrics(d):
    m = {}
    # --- alignment ---------------------------------------------------------
    for st in sorted(glob.glob(f'{d}/minimap2/*.stats')):
        sample = re.sub(r'(\.sorted)?(\.bam)?\.stats$', '', os.path.basename(st))
        sn = {}
        for line in open(st):
            if line.startswith('SN\t'):
                f = line.rstrip('\n').split('\t')
                sn[f[1].rstrip(':')] = f[2]
        m[f'minimap2: {sample} reads total / mapped'] = f"{sn.get('raw total sequences', '-')} / {sn.get('reads mapped', '-')}"
    m['minimap2: SAM files published (MINIMAP2/*.sam)'] = len(glob.glob(f'{d}/MINIMAP2/*.sam'))
    m['bigwig files'] = len(glob.glob(f'{d}/bigwig/*.bw'))
    m['junction BED files'] = len(glob.glob(f'{d}/junctions/*.bed'))

    # --- bambu ---------------------------------------------------------------
    ct = f'{d}/bambu/counts_transcript.txt'
    if os.path.exists(ct):
        ids = [l.split('\t', 1)[0] for l in open(ct)][1:]
        m['bambu: transcripts'] = len(ids)
        m['bambu: novel transcripts (BambuTx*)'] = sum(i.startswith('BambuTx') for i in ids)
    cg = f'{d}/bambu/counts_gene.txt'
    if os.path.exists(cg):
        m['bambu: genes'] = n_lines(cg) - 1

    # --- DE -----------------------------------------------------------------
    for sub in ('DE_genes', 'DE_transcripts'):
        for f in sorted(glob.glob(f'{d}/{sub}/DE_results*.csv')):
            rows = read_csv(f)
            key = f'{sub}/{os.path.basename(f)}'
            if rows and 'padj' in rows[0]:
                m[f'{key}: rows'] = len(rows)
                m[f'{key}: padj < 0.05'] = count_sig(rows, 'padj')
            else:
                m[f'{key}'] = 'no DE (message file)'

    # --- isoform switching / APA ---------------------------------------------
    f = f'{d}/isoswitch/isoform_switch_results.csv'
    if os.path.exists(f):
        rows = read_csv(f)
        m['isoswitch: isoformFeatures rows'] = len(rows)
        if rows and 'isoform_switch_q_value' in rows[0]:
            m['isoswitch: q < 0.05 (any dIF)'] = count_sig(rows, 'isoform_switch_q_value')
            m['isoswitch: q < 0.05 & |dIF| > 0.1'] = sum(
                1 for r in rows
                if (q := to_float(r.get('isoform_switch_q_value'))) is not None and q < 0.05
                and (di := to_float(r.get('dIF'))) is not None and abs(di) > 0.1)
    f = f'{d}/isoswitch/significant_switches.csv'
    if os.path.exists(f):
        rows = read_csv(f)
        m['isoswitch: significant_switches.csv rows'] = len(rows)
        m['isoswitch: significant_switches.csv all-NA rows'] = sum(
            1 for r in rows if all(v in ('NA', '') for k, v in r.items() if k != ''))
    f = f'{d}/isoswitch/clustered_base_level_PAS_usage.csv'
    if os.path.exists(f):
        rows = read_csv(f)
        m['APA: PAS clusters'] = len(rows)
        m['APA: PAS clusters q < 0.05'] = count_sig(rows, 'pas_q_value')
    f = f'{d}/isoswitch/APA_gene_summary.csv'
    m['APA: APA_gene_summary.csv genes'] = len(read_csv(f)) if os.path.exists(f) else '- (not written)'
    for name in ('atss_results.csv', 'atts_results.csv', 'all_AS_events.csv'):
        f = f'{d}/isoswitch/{name}'
        m[f'isoswitch: {name} rows'] = n_lines(f) - 1 if os.path.exists(f) else '-'

    # --- proActiv -----------------------------------------------------------
    f = f'{d}/proactiv/promoter_counts.tsv'
    if os.path.exists(f):
        m['proactiv: promoters'] = n_lines(f) - 1
    f = f'{d}/proactiv/alternative_promoters.tsv'
    m['proactiv: alternative promoters'] = n_lines(f) - 1 if os.path.exists(f) else '- (not written)'

    # --- SUPPA2 ---------------------------------------------------------------
    tot_cells = nan_cells = 0
    for psi in sorted(glob.glob(f'{d}/suppa2/psi/*_samples.psi')):
        ev = os.path.basename(psi).split('_')[1]
        lines = open(psi).read().splitlines()[1:]
        cells = [c for l in lines for c in l.split('\t')[1:]]
        nan = sum(c in ('nan', 'NA', '') for c in cells)
        tot_cells += len(cells)
        nan_cells += nan
        flt = psi.replace('_samples.psi', '_samples_filtered.psi')
        nflt = n_lines(flt) - 1 if os.path.exists(flt) else '-'
        m[f'suppa2 {ev}: events / PSI nan % / events kept after filter'] = \
            f'{len(lines)} / {100 * nan / max(len(cells), 1):.1f}% / {nflt}'
    if tot_cells:
        m['suppa2: overall PSI nan %'] = f'{100 * nan_cells / tot_cells:.1f}%'

    # --- TE -------------------------------------------------------------------
    for name in ('TE_family_counts.tsv', 'TE_instance_counts.tsv', 'isoform_TE_exonization.tsv', 'isoform_TE_annotation.tsv'):
        f = f'{d}/TE/{name}'
        m[f'TE: {name} data rows'] = n_lines(f) - 1 if os.path.exists(f) else '-'
    for t in ('Family', 'Instance'):
        f = f'{d}/TE/DE/DE_TE_{t}_results.csv'
        if os.path.exists(f):
            m[f'TE DE {t}: padj < 0.05'] = count_sig(read_csv(f), 'padj')
    return m


def file_table(d):
    rows = []
    for root, _, files in os.walk(d):
        for name in files:
            p = os.path.join(root, name)
            val = f'{n_lines(p)} lines' if name.endswith(TEXT_EXT) else f'{os.path.getsize(p)} bytes'
            rows.append(f'| `{os.path.relpath(p, d)}` | {val} |')
    return '\n'.join(['| File | Lines / size |', '|---|---|'] + sorted(rows))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('results_dir')
    ap.add_argument('--title', default=None)
    ap.add_argument('--files', action='store_true')
    a = ap.parse_args()
    print(f'## {a.title or a.results_dir}\n')
    print('| Metric | Value |\n|---|---|')
    for k, v in metrics(a.results_dir).items():
        print(f'| {k} | {v} |')
    if a.files:
        print('\n### Files\n')
        print(file_table(a.results_dir))


if __name__ == '__main__':
    main()
