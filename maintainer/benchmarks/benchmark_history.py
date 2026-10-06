#!/usr/bin/env python3
#
# Copyright (C) 2018-2026 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""
Append the results of the ESPResSo benchmark runs to a single HTML page.

The page holds one collapsible table per run directory (one execution of
suite.sh), newest on top, each collapsed until its summary line is clicked.
The run directories are selected with ``--prefix`` exactly as in
plot_benchmarks.py. Every run directory that is not yet on the page is added
as a new table; tables already on the page are never rewritten, so calling the
script twice with the same runs leaves the page unchanged. The first call
against an existing history adds every earlier run.

Usage:
    python3 benchmark_history.py --prefix PREFIX_REGEX -o OUTPUT.html
"""

import argparse
import csv
import html
import os
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Iterable, TypedDict, cast

from plot_benchmarks import (
    parse_timestamp,
    resolve_perflogs,
    shorten_label,
    split_metric,
    to_float,
)

# New run tables are inserted right after this marker.
RUNS_MARKER = "<!-- RUNS -->"
RUN_ID_RE = re.compile(r'<details class="run" id="run-([^"]+)"')
UPDATED_RE = re.compile(r"(<p class=\"updated\">Last updated: )[^<]*(</p>)")

# (build_config, descr, metric label) of one table row.
MetricKey = tuple[str, str, str]


class Metric(TypedDict, total=False):
    """One table row: a metric's mean/ci pair and its unit and result."""
    mean: float
    ci: float
    unit: str
    result: str


class RunSummary(TypedDict):
    """The rows and summary line data of one run directory."""
    metrics: dict[MetricKey, Metric]
    n_pass: int
    n_fail: int
    commits: list[str]
    hosts: list[str]
    time: datetime | None


PAGE_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>ESPResSo benchmark history</title>
<style>
body {{ font-family: system-ui, sans-serif; margin: 1.5em; color: #222;
       background: #fff; }}
h1 {{ font-size: 1.4em; }}
.updated {{ color: #666; }}
details.run {{ border: 1px solid #ccc; border-radius: 4px;
              margin: 0.4em 0; padding: 0.3em 0.6em; }}
details.run > summary {{ cursor: pointer; font-family: monospace; }}
details.run[open] > summary {{ margin-bottom: 0.5em; }}
.fail {{ color: #b00020; font-weight: bold; }}
table {{ border-collapse: collapse; font-size: 0.9em; }}
th, td {{ border: 1px solid #ddd; padding: 0.2em 0.6em; text-align: left; }}
th {{ background: #f3f3f3; }}
td.num {{ text-align: right; font-family: monospace; }}
th[scope="row"] {{ background: none; font-weight: normal; }}
summary:focus-visible, .table-wrap:focus-visible {{
    outline: 3px solid #1a5fb4; outline-offset: 2px; }}
.table-wrap {{ overflow-x: auto; }}
/* Read by screen readers, not shown on screen. */
.visually-hidden {{ position: absolute; width: 1px; height: 1px;
                    margin: -1px; padding: 0; overflow: hidden;
                    clip: rect(0 0 0 0); white-space: nowrap; border: 0; }}
</style>
</head>
<body>
<main>
<h1>ESPResSo benchmark history</h1>
<p class="updated">Last updated: {updated}</p>
<p>Every entry corresponds to one benchmark run. Each entry is collapsed;
activate it (click, or Enter/Space on the focused entry) to show the
results table of that run.</p>
{marker}
</main>
</body>
</html>
"""


def read_run(log_paths: Iterable[Path]) -> RunSummary | None:
    """Collect the table rows and summary data of one run directory.

    Returns ``None`` if no perflog of the run holds any timing data.
    """
    # (build, descr, label) -> {"mean", "ci", "unit", "result"}
    metrics: defaultdict[MetricKey, Metric] = defaultdict(Metric)
    # (build, descr) -> set of results, to count passed/failed tests
    tests: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    commits: set[str] = set()
    hosts: set[str] = set()
    times: list[datetime] = []

    for log_path in log_paths:
        try:
            with open(log_path, newline="") as f:
                reader = csv.DictReader(f)
                fields = {name.strip(): name
                          for name in reader.fieldnames or ()}
                if not {"descr", "perf_value", "perf_var"} <= fields.keys():
                    raise ValueError("not a long-format ESPResSo perflog")
                for row in reader:
                    if None in row.values():
                        continue  # truncated final line of a log being written

                    def get(col: str) -> str:
                        if col not in fields:
                            return ""
                        return (row[fields[col]] or "").strip()

                    value = to_float(get("perf_value"))
                    descr = get("descr")
                    if value is None or not descr:
                        continue  # build/dry-run row
                    build = get("build_config") or "unknown"
                    label, kind = split_metric(get("perf_var"))
                    entry = metrics[(build, descr, label)]
                    if kind == "ci":
                        entry["ci"] = value
                    else:
                        entry["mean"] = value
                    if kind == "mean" or "unit" not in entry:
                        entry["unit"] = get("perf_unit")
                    result = get("result") or "unknown"
                    entry["result"] = result
                    tests[(build, descr)].add(result)
                    if get("espresso_commit") not in ("", "null", "unknown"):
                        commits.add(get("espresso_commit"))
                    if get("hostname"):
                        hosts.add(get("hostname"))
                    try:
                        times.append(parse_timestamp(
                            get("job_completion_time")))
                    except ValueError:
                        pass
        except (OSError, ValueError, csv.Error) as exc:
            print(f"Skipping {log_path}: {exc}", file=sys.stderr)

    if not metrics:
        return None
    n_fail = sum(1 for results in tests.values() if results != {"pass"})
    return {
        "metrics": dict(metrics),
        "n_pass": len(tests) - n_fail,
        "n_fail": n_fail,
        "commits": sorted(commits),
        "hosts": sorted(hosts),
        "time": max(times) if times else None,
    }


# Separator between the parts of a summary line: a dot on screen, a plain
# pause for screen readers (which would otherwise read "middle dot").
SEP = ('<span aria-hidden="true"> · </span>'
       '<span class="visually-hidden">, </span>')


def spoken(shown: str, read: str) -> str:
    """Show ``shown`` on screen but have screen readers read ``read``."""
    return (f'<span aria-hidden="true">{shown}</span>'
            f'<span class="visually-hidden">{read}</span>')


NONE_HTML = spoken("–", "none")


def fmt_number(value: float | None) -> str:
    """Format a value; exponents are spelled out for screen readers."""
    if value is None:
        return NONE_HTML
    text = f"{value:.4g}"
    mantissa, has_exp, exponent = text.partition("e")
    if not has_exp:
        return text
    return spoken(text,
                  f"{mantissa} times 10 to the power of {int(exponent)}")


# Spoken names of the units the benchmarks report.
UNIT_NAMES = {"s": "seconds", "ms": "milliseconds", "us": "microseconds"}


def fmt_unit(unit: str) -> str:
    if not unit:
        return NONE_HTML
    if unit in UNIT_NAMES:
        return spoken(html.escape(unit), UNIT_NAMES[unit])
    return html.escape(unit)


def render_run(name: str, run: RunSummary) -> str:
    """Render one run directory as a collapsed ``<details>`` block."""
    esc = html.escape
    when = run["time"].strftime("%Y-%m-%d %H:%M") if run["time"] else "unknown"
    commits = run["commits"] or ["unknown"]
    commit_html = ", ".join(
        f'<span title="{esc(c)}">{esc(c[:7])}</span>' for c in commits)
    fail_class = ' class="fail"' if run["n_fail"] else ""
    summary = SEP.join((
        f"Run {esc(name)}",
        f"completed {when}",
        f"commit {commit_html}",
        f"host {esc(', '.join(run['hosts']) or 'unknown')}",
        f"{run['n_pass']} passed, "
        f"<span{fail_class}>{run['n_fail']} failed</span>",
    ))

    rows: list[str] = []
    for (build, descr, label), m in sorted(run["metrics"].items()):
        result = m.get("result", "unknown")
        short = shorten_label(descr)
        result_class = "" if result == "pass" else ' class="fail"'
        rows.append(
            "<tr>"
            # Underscores are read out one by one; speak them as spaces.
            f'<th scope="row">{
                spoken(esc(short), esc(short.replace("_", " ")))}'
            "</th>"
            f"<td>{esc(build)}</td>"
            f"<td>{esc(label) or NONE_HTML}</td>"
            f'<td class="num">{fmt_number(m.get("mean"))}</td>'
            f'<td class="num">{fmt_number(m.get("ci"))}</td>'
            f"<td>{fmt_unit(m.get('unit', ''))}</td>"
            f"<td{result_class}>{esc(result)}</td>"
            "</tr>"
        )
    caption = f"Results of run {esc(name)}"
    return (
        f'<details class="run" id="run-{esc(name)}">\n'
        f"<summary>{summary}</summary>\n"
        # Focusable so keyboard users can scroll a table wider than the screen.
        f'<div class="table-wrap" role="region" tabindex="0" '
        f'aria-label="{caption}">\n'
        f'<table>\n<caption class="visually-hidden">{caption}</caption>\n'
        '<thead><tr><th scope="col">Benchmark</th><th scope="col">Build</th>'
        '<th scope="col">Metric</th><th scope="col">Mean</th>'
        '<th scope="col">CI<span class="visually-hidden">, confidence '
        'interval</span></th><th scope="col">Unit</th>'
        '<th scope="col">Result</th></tr></thead>\n'
        "<tbody>\n" + "\n".join(rows) + "\n</tbody>\n</table>\n</div>\n"
        "</details>\n"
    )


def load_page(path: Path) -> str:
    """Return the existing page, or a fresh skeleton if there is none."""
    if path.is_file():
        page = path.read_text(encoding="utf-8")
        if RUNS_MARKER not in page:
            raise ValueError(f"{path} lacks the {RUNS_MARKER!r} marker")
        return page
    return PAGE_TEMPLATE.format(updated="", marker=RUNS_MARKER)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-o", "--output", required=True,
                        help="HTML page to create or append to.")
    parser.add_argument(
        "--prefix",
        required=True,
        help="Run directories to include, as a regular expression with the "
        "same meaning as in plot_benchmarks.py.",
    )
    args = parser.parse_args(argv)
    output = Path(args.output)

    log_paths, error = resolve_perflogs(args.prefix)
    if error is not None:
        parser.error(error)
    log_paths = cast(list[Path], log_paths)

    try:
        page = load_page(output)
    except (OSError, ValueError) as exc:
        print(f"Error: cannot read {output}: {exc}", file=sys.stderr)
        return 1
    known = set(html.unescape(m) for m in RUN_ID_RE.findall(page))

    by_run: defaultdict[Path, list[Path]] = defaultdict(list)
    for log in log_paths:
        by_run[log.parents[3]].append(log)

    new_runs: list[tuple[str, RunSummary]] = []
    for run_dir, logs in by_run.items():
        if run_dir.name in known:
            continue
        run = read_run(logs)
        if run is None:
            print(f"No timing data in {run_dir}; not added.", file=sys.stderr)
            continue
        new_runs.append((run_dir.name, run))

    if not new_runs:
        print(f"No new runs to add to {output}.")
        return 0 if output.is_file() else 1

    new_runs.sort(key=lambda item: (item[1]["time"] or datetime.min, item[0]),
                  reverse=True)
    blocks = "".join(render_run(name, run) for name, run in new_runs)
    page = page.replace(RUNS_MARKER, RUNS_MARKER +
                        "\n" + blocks.rstrip("\n"), 1)
    updated = datetime.now().strftime("%Y-%m-%d %H:%M")
    page = UPDATED_RE.sub(lambda m: m.group(1) + updated + m.group(2), page, 1)

    tmp = output.with_name(output.name + ".tmp")
    tmp.write_text(page, encoding="utf-8")
    os.replace(tmp, output)
    print(f"Added {len(new_runs)} run(s) to {output}: "
          + ", ".join(name for name, _ in new_runs))
    return 0


if __name__ == "__main__":
    sys.exit(main())
