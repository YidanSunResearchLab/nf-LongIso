#!/usr/bin/env python3
"""Merge per-task versions.yml files into one software-versions report.

Usage:
    dump_software_versions.py --pipeline-version 0.2.0 --nextflow-version 25.10.0 versions1.yml [versions2.yml ...]

Each input file maps a process name to {tool: version}. Values are read as plain strings
(yaml BaseLoader), so versions such as 1.20 are not turned into the number 1.2.

Outputs
    software_versions.yml      process -> {tool: version}, plus a 'Workflow' entry
    software_versions_mqc.yml  MultiQC custom content (HTML table)
"""
import argparse
import html

import yaml


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipeline-version", required=True)
    ap.add_argument("--nextflow-version", required=True)
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()

    merged = {}
    for path in a.files:
        with open(path) as fh:
            data = yaml.load(fh, Loader=yaml.BaseLoader) or {}
        for process, tools in data.items():
            entry = merged.setdefault(process, {})
            for tool, version in (tools or {}).items():
                entry[tool] = str(version).strip()

    merged = {k: dict(sorted(v.items())) for k, v in sorted(merged.items())}
    merged["Workflow"] = {"NextLongIso": a.pipeline_version, "Nextflow": a.nextflow_version}

    with open("software_versions.yml", "w") as fh:
        yaml.safe_dump(merged, fh, default_flow_style=False, sort_keys=False)

    rows = []
    for process, tools in merged.items():
        for tool, version in tools.items():
            rows.append(
                f"<tr><td>{html.escape(process)}</td><td>{html.escape(tool)}</td>"
                f"<td><samp>{html.escape(version)}</samp></td></tr>"
            )
    table = (
        '<table class="table" style="width:auto">'
        "<thead><tr><th>Process</th><th>Software</th><th>Version</th></tr></thead>"
        "<tbody>" + "".join(rows) + "</tbody></table>"
    )
    mqc = {
        "id": "software_versions",
        "section_name": "NextLongIso software versions",
        "section_href": "https://github.com/YidanSunResearchLab/nf-LongIso",
        "plot_type": "html",
        "description": "Versions of the tools used by each process, collected from the containers at run time.",
        "data": table,
    }
    with open("software_versions_mqc.yml", "w") as fh:
        yaml.safe_dump(mqc, fh, default_flow_style=False, sort_keys=False)


if __name__ == "__main__":
    main()
