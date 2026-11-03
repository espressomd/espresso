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
Read a ReFrame perflog produced by the ESPResSo benchmark suite and render an
interactive Plotly timeline of the recorded benchmark performances across
ESPResSo commits.

The x-axis is ordinal: one tick per ESPResSo commit, ordered chronologically by
the earliest completion time at which the commit was benchmarked. Each tick is
labelled with the short commit hash; the completion time is shown on hover. A
commit that was benchmarked more than once gets one adjacent column per run.
Each benchmark metric (``MC``/``MD``/... or the unnamed default) is drawn as its
own line connecting its mean value across commits, with the confidence interval
as an error bar.

The suite can be run with several ESPResSo builds (the ``build_config`` column,
e.g. ``maxset``/``default``/``empty``); each build gets its own figure covering
its whole history, written as Plotly figure JSON:

    <stem>[_<build>].json     e.g. EspressoBenchmark_maxset.json

A figure opens on the newest ``--window`` commit columns (default 20); the range
slider below the plot reaches the older ones. benchmark_history.py embeds the
figures in the history page.

The run directories to read are selected with ``--prefix``. ``suite.sh`` writes
every invocation to its own ``<prefix>_dd_mm_yyyy_<n>`` directory, so ``--prefix``
takes a *regular expression* (Python ``re``) and the timeline is assembled from
all the run directories it matches:

    --prefix '/path/to/benchmarks_[0-9]{2}_[0-9]{2}_[0-9]{4}_[0-9]+'

The value is split at its last ``/``: everything before it is a literal parent
directory, and the trailing component is matched against the *names* of that
directory's entries with ``re.fullmatch``. An existing directory is taken
verbatim, so passing a plain path still plots exactly that one directory.

In each matched directory the perflog lives at
``perflogs/<system>/<partition>/EspressoBenchmark.log`` -- ReFrame's ``filelog``
handler names those two directories after the system and partition it ran on, so
they are ``local/default`` on a workstation but e.g. ``ant_cluster/debug`` on a
cluster. The path is therefore discovered rather than assumed. Timings from
different systems or partitions are not comparable, so matching a mixture of them
is an error rather than a silently meaningless plot -- narrow the regex (debug
runs of the suite live in ``_DEBUG`` directories) or pass ``--log`` to select a
single perflog explicitly.

Usage:
    python3 plot_benchmarks.py --prefix PREFIX_REGEX -o OUTPUT.json
                               [--window N] [--log PERFLOG]
"""

import argparse
import csv
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import cast

import plotly.graph_objects as go
from plotly.colors import qualitative

# Basename of the perflog ReFrame's ``filelog`` handler writes (named after the
# benchmark's check class) and the glob that locates it below a ``--prefix``.
# The two wildcards are ReFrame's system and partition names.
PERFLOG_NAME = "EspressoBenchmark.log"
PERFLOG_GLOB = f"perflogs/*/*/{PERFLOG_NAME}"

# Sentinel marking the fields a truncated (short) perflog row does not have.
_TRUNCATED = "\x00truncated"


def find_run_dirs(pattern):
    """Return ``(dirs, error)``: the run directories selected by ``--prefix``.

    ``pattern`` is split at its last ``/``: everything before it is a literal
    parent directory, the trailing component is a regular expression matched
    against the *names* of that directory's entries with ``re.fullmatch``. An
    existing directory is used verbatim, so a plain path still selects exactly
    itself even when it contains regex metacharacters.

    ``fullmatch`` rather than ``search`` is what keeps the ``_DEBUG``,
    ``_DRYRUN`` and ``_LIST`` run directories written by suite.sh out of a
    production pattern.
    """
    literal = Path(pattern)
    if literal.is_dir():
        return [literal], None

    parent_str, sep, name_pattern = pattern.rpartition("/")
    if sep and not name_pattern:
        return None, f"--prefix must not end with '/': {pattern!r}"
    parent = Path(parent_str) if parent_str else Path(".")
    try:
        regex = re.compile(name_pattern)
    except re.error as exc:
        return None, (f"Invalid --prefix regular expression {name_pattern!r}: {exc}")
    try:
        entries = list(parent.iterdir())
    except OSError as exc:
        return None, f"Cannot list {parent}: {exc}"

    # resolve() collapses a symlinked or copied duplicate of a run directory,
    # which would otherwise contribute the same data points twice.
    matches = {
        entry.resolve()
        for entry in entries
        if entry.is_dir() and regex.fullmatch(entry.name)
    }
    if not matches:
        return None, (
            f"No directory under {parent} matches {name_pattern!r}. "
            "Has the suite been run with this prefix?"
        )
    return sorted(matches), None


def perflog_group(log_path):
    """The ``<system>/<partition>`` pair a perflog was written for."""
    return f"{log_path.parts[-3]}/{log_path.parts[-2]}"


def resolve_perflogs(pattern, explicit=None):
    """Resolve every perflog to merge into one timeline.

    Returns ``(paths, error)``; exactly one of the two is ``None``. Matched
    directories holding no perflog are skipped silently: listing and dry runs
    leave such directories behind, and so does a run without timing data.
    Perflogs of different systems/partitions are not comparable, so matching a
    mixture of them is reported as an error.
    """
    if explicit is not None:
        path = Path(explicit)
        if not path.is_file():
            return None, f"Log file not found: {path}"
        return [path], None

    run_dirs, error = find_run_dirs(pattern)
    if error is not None:
        return None, error
    run_dirs = cast(list[Path], run_dirs)

    candidates = sorted(
        log for run_dir in run_dirs for log in run_dir.glob(PERFLOG_GLOB)
    )
    if not candidates:
        return None, (
            f"None of the {len(run_dirs)} directories matching {pattern!r} "
            f"contains {PERFLOG_GLOB}."
        )

    groups = defaultdict(list)
    for log in candidates:
        groups[perflog_group(log)].append(log)
    if len(groups) > 1:
        listing = "\n  ".join(
            f"{group}: " + ", ".join(str(log.parents[3]) for log in logs)
            for group, logs in sorted(groups.items())
        )
        return None, (
            "Perflogs of several systems/partitions matched; their timings are "
            "not comparable. Narrow the --prefix regex (debug runs of the suite "
            "live in _DEBUG directories) or select a single log with --log:\n  "
            + listing
        )
    return candidates, None


def parse_timestamp(value):
    """Parse a ReFrame ``job_completion_time`` into a datetime.

    ReFrame has emitted timestamps both as ``2026-05-19T13:00:51`` (ISO, with
    or without a timezone offset) and as ``2026-04-10 20:03:25`` (space
    separated). ``datetime.fromisoformat`` handles both, so try it first and
    fall back to a couple of explicit formats.

    An offset, when present, is converted to local time and dropped: merging
    perflogs written before and after a ReFrame upgrade can mix aware and naive
    timestamps, which cannot be compared with each other.
    """
    value = value.strip()
    try:
        return _drop_timezone(datetime.fromisoformat(value))
    except ValueError:
        pass
    for fmt in ("%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S"):
        try:
            return _drop_timezone(datetime.strptime(value, fmt))
        except ValueError:
            continue
    raise ValueError(f"Unrecognised timestamp format: {value!r}")


def _drop_timezone(stamp):
    """Return ``stamp`` as local time without a timezone."""
    if stamp.tzinfo is None:
        return stamp
    return stamp.astimezone().replace(tzinfo=None)


def to_float(value):
    """Convert a perflog field to float, treating empty/``null`` as missing."""
    if value is None:
        return None
    value = value.strip()
    if value == "" or value.lower() == "null":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def split_metric(perf_var):
    """Split a perf-variable name into a (label, kind) pair.

    ``set_perf_variables`` in espresso_benchmarks.py names variables
    ``<label>_mean`` / ``<label>_ci`` (or plain ``mean`` / ``ci`` when the
    benchmark uses a single, unnamed label). Returns ``(label, "mean")``,
    ``(label, "ci")`` or ``(perf_var, "mean")`` for anything unrecognised (so
    stray metrics still show up rather than being silently dropped).
    """
    name = perf_var.strip()
    for suffix, kind in (("_mean", "mean"), ("_ci", "ci")):
        if name.endswith(suffix):
            return name[: -len(suffix)], kind
    if name in ("mean", "ci"):
        return "", name
    return name, "mean"


class SeriesStore:
    """Accumulates data points keyed by (build, descr, label) then by timestamp."""

    def __init__(self):
        # (build, descr, label) -> {timestamp: {"mean", "ci", "unit", "commit"}}
        self._data = defaultdict(lambda: defaultdict(dict))
        self.units = set()
        self.commits = set()

    def add(self, build, descr, label, time, kind, value, unit=None, commit=None):
        point = self._data[(build, descr, label)][time]
        point[kind] = value
        if unit and unit.lower() != "null":
            point["unit"] = unit
            self.units.add(unit)
        commit = (commit or "").strip()
        if commit and commit.lower() not in ("null", "unknown"):
            point["commit"] = commit
            self.commits.add(commit)
        else:
            point.setdefault("commit", "unknown")

    def series(self):
        """Yield (build, descr, label, sorted_points) with a defined mean value."""
        for (build, descr, label), by_time in self._data.items():
            points = [
                {
                    "time": t,
                    "mean": p.get("mean"),
                    "ci": p.get("ci"),
                    "commit": p.get("commit", "unknown"),
                }
                for t, p in sorted(by_time.items())
                if p.get("mean") is not None
            ]
            if points:
                yield build, descr, label, points

    def builds(self):
        """Return the sorted set of build names present."""
        return sorted({key[0] for key in self._data})

    def is_empty(self):
        return not any(True for _ in self.series())


def read_records(log_path, store=None):
    """Read a perflog into a SeriesStore, handling long and wide layouts.

    Points are added to ``store`` when one is given, so the perflogs of several
    run directories accumulate into a single timeline.
    """
    store = SeriesStore() if store is None else store
    with open(log_path, newline="") as f:
        reader = csv.DictReader(f, restval=_TRUNCATED)
        if reader.fieldnames is None:
            raise ValueError(f"{log_path} is empty")

        fields = {name.strip(): name for name in reader.fieldnames}
        if "descr" not in fields or "job_completion_time" not in fields:
            raise ValueError(
                f"{log_path} is missing the 'descr'/'job_completion_time' "
                "columns; is this a valid ESPResSo perflog?"
            )

        descr_col = fields["descr"]
        time_col = fields["job_completion_time"]
        commit_col = fields.get("espresso_commit")
        build_col = fields.get("build_config")

        # Detect the layout.
        long_format = "perf_var" in fields and "perf_value" in fields
        # Wide layout: every column that ends in "mean_value" defines a metric.
        wide_labels = _wide_metric_labels(fields) if not long_format else {}

        if not long_format and not wide_labels:
            raise ValueError(
                f"{log_path} contains no performance columns "
                "(no 'perf_value' and no '*mean_value'). This looks like a "
                "build-only or dry-run log with no timing data."
            )

        for row in reader:
            if _TRUNCATED in row.values():
                continue  # truncated final line of a log being written
            descr = row[descr_col].strip()
            if not descr:
                continue
            try:
                time = parse_timestamp(row[time_col])
            except ValueError:
                continue
            commit = row[commit_col] if commit_col else None
            build = (row[build_col].strip() if build_col else "") or "unknown"
            if build.lower() in ("null", ""):
                build = "unknown"

            if long_format:
                value = to_float(row[fields["perf_value"]])
                if value is None:
                    continue  # build/dry-run row
                label, kind = split_metric(row[fields["perf_var"]])
                unit = row[fields["perf_unit"]
                           ] if "perf_unit" in fields else None
                store.add(build, descr, label, time, kind, value, unit, commit)
            else:
                for label, cols in wide_labels.items():
                    mean = to_float(row.get(cols["mean"], ""))
                    if mean is not None:
                        unit = row.get(cols.get("unit", ""), None)
                        store.add(build, descr, label, time,
                                  "mean", mean, unit, commit)
                    ci = to_float(row.get(cols.get("ci", ""), ""))
                    if ci is not None:
                        store.add(build, descr, label, time,
                                  "ci", ci, None, commit)

    return store


def read_all_records(log_paths):
    """Merge every perflog in ``log_paths`` into one SeriesStore.

    A single unreadable file (an empty log, or one holding no timing data, as
    left behind by a build-only or dry run) is reported and skipped rather than
    aborting the whole timeline. Returns ``(store, n_read)``.
    """
    store = SeriesStore()
    n_read = 0
    for path in log_paths:
        try:
            read_records(path, store)
        except (ValueError, OSError) as exc:
            print(f"Skipping {path}: {exc}", file=sys.stderr)
        else:
            n_read += 1
    return store, n_read


def _wide_metric_labels(fields):
    """Map metric label -> {mean, ci, unit} column names for a wide log."""
    labels = {}
    for name in fields:
        if name.endswith("mean_value"):
            label = name[: -len("mean_value")].rstrip("_")
            prefix = f"{label}_" if label else ""
            labels[label] = {
                "mean": fields[name],
                "ci": fields.get(f"{prefix}ci_value"),
                "unit": fields.get(f"{prefix}mean_unit"),
            }
    return labels


def shorten_label(descr):
    """Strip the redundant ``ESPRESSO_`` prefix for a compact legend entry."""
    return descr[len("ESPRESSO_"):] if descr.startswith("ESPRESSO_") else descr


def series_legend(descr, label):
    base = shorten_label(descr)
    return f"{base} [{label}]" if label else base


def build_commit_columns(all_series):
    """Lay out the x-axis as one column per (commit, run#).

    Commits are ordered chronologically by the earliest completion time at which
    they were seen. Because each benchmark in a suite invocation finishes at a
    slightly different time, completion time cannot align metrics from the same
    build -- the *commit* is what groups them into a shared column.

    A commit that was benchmarked more than once gets one adjacent column per
    run (the k-th run of a commit, ordered by completion time). This preserves
    every data point instead of averaging repeated runs.

    Returns ``(columns, plotted)`` where:
      * ``columns`` is the ordered list of ``{"commit", "run", "time", "runs"}``
        dicts, one per x position (``time`` is the representative completion
        time, ``runs`` the number of columns this commit spans);
      * ``plotted`` is a list of ``(descr, label, points)`` where each point is
        ``{"x", "mean", "ci", "time", "commit"}`` (``time`` is that point's own
        completion time; ``x`` is the integer column index).
    """
    # Earliest time each commit was seen -> chronological commit order.
    commit_first_time = {}
    for _descr, _label, points in all_series:
        for p in points:
            c, t = p["commit"], p["time"]
            if c not in commit_first_time or t < commit_first_time[c]:
                commit_first_time[c] = t
    commit_order = sorted(
        commit_first_time, key=lambda c: commit_first_time[c])

    # Per series, group points by commit (sorted by time -> run index), and
    # track how many runs each commit has across all series.
    grouped = []  # (descr, label, {commit: [points sorted by time]})
    max_runs = {c: 0 for c in commit_order}
    for descr, label, points in all_series:
        by_commit = defaultdict(list)
        for p in sorted(points, key=lambda p: p["time"]):
            by_commit[p["commit"]].append(p)
        grouped.append((descr, label, by_commit))
        for c, pts in by_commit.items():
            max_runs[c] = max(max_runs[c], len(pts))

    # Ordered columns: (commit, run#).
    columns = [
        {"commit": c, "run": k, "time": None, "runs": max_runs[c]}
        for c in commit_order
        for k in range(max_runs[c])
    ]
    col_index = {(col["commit"], col["run"]): i for i,
                 col in enumerate(columns)}

    # Assign each point to its column and record the representative (earliest)
    # completion time per column.
    plotted = []
    for descr, label, by_commit in grouped:
        pts_out = []
        for c, pts in by_commit.items():
            for k, p in enumerate(pts):
                i = col_index[(c, k)]
                pts_out.append(
                    {
                        "x": i,
                        "mean": p["mean"],
                        "ci": p["ci"],
                        "time": p["time"],
                        "commit": c,
                    }
                )
                col = columns[i]
                if col["time"] is None or p["time"] < col["time"]:
                    col["time"] = p["time"]
        pts_out.sort(key=lambda d: d["x"])
        plotted.append((descr, label, pts_out))

    return columns, plotted


def format_tick_label(column):
    """X tick: short commit hash (completion time now lives in the hover tooltip).

    When a commit spans several runs, a ``#k`` suffix keeps the otherwise
    identical adjacent ticks distinguishable.
    """
    commit = column["commit"]
    commit = commit[:9] if commit != "unknown" else commit
    if column["runs"] > 1:
        return f"{commit} (#{column['run'] + 1})"
    else:
        return f"{commit} (#1)"


# Default number of commit columns (x-axis positions) shown when a figure opens;
# the range slider below the plot reaches the rest of the history.
DEFAULT_WINDOW = 20


def hover_text(descr, label, point, unit):
    """Tooltip of one data point: value, CI, metric, commit and completion time."""
    commit = point["commit"]
    lines = [f"<b>{point['mean']:.4g} {unit}</b>"]
    if point["ci"] is not None:
        lines.append(f"± {point['ci']:.3g} {unit}")
    lines.append(series_legend(descr, label))
    lines.append("commit " + (commit if commit == "unknown" else commit[:12]))
    if point["time"]:
        lines.append(point["time"].strftime("%Y-%m-%d %H:%M:%S"))
    return "<br>".join(lines)


def build_figure(build, columns, plotted, unit, window=DEFAULT_WINDOW):
    """Build the interactive timeline of one build as a Plotly figure.

    The x-axis spans every commit column; it opens on the newest ``window``
    columns, and the range slider below the plot reaches the older ones.
    """
    fig = go.Figure()
    for descr, label, points in plotted:
        fig.add_trace(
            go.Scatter(
                x=[p["x"] for p in points],
                y=[p["mean"] for p in points],
                error_y=dict(
                    type="data",
                    array=[p["ci"] if p["ci"] is not None else 0.0
                           for p in points],
                    thickness=1,
                    width=3,
                ),
                mode="lines+markers",
                marker=dict(size=6),
                line=dict(width=1.2),
                name=series_legend(descr, label),
                text=[hover_text(descr, label, p, unit) for p in points],
                hovertemplate="%{text}<extra></extra>",
            )
        )

    n = len(columns)
    title = "ESPResSo benchmark performance timeline"
    if build and build != "unknown":
        title += f"  —  build: {build}"
    fig.update_layout(
        title=title,
        template="plotly_white",
        legend=dict(title="Benchmark [metric]"),
        hovermode="closest",
        # More colours than the default ten, so series rarely share one.
        colorway=qualitative.Dark24,
    )
    # The x-axis is numeric rather than categorical: the tick labels of a
    # commit's repeated runs may coincide.
    fig.update_xaxes(
        title="ESPResSo commit (chronological)",
        tickvals=list(range(n)),
        ticktext=[format_tick_label(c) for c in columns],
        tickangle=-30,
        range=[max(-0.5, n - window - 0.5), n - 0.5],
        rangeslider=dict(visible=True, thickness=0.08),
    )
    fig.update_yaxes(
        title=f"Mean execution time ({unit})",
        type="log",
        # Label the decades only, as powers of ten.
        dtick=1,
        exponentformat="power",
        showgrid=True,
        minor=dict(showgrid=True),
    )
    return fig


def _sanitize(text):
    """Turn a string into a filesystem-safe filename token."""
    return "".join(c if (c.isalnum() or c in "._-") else "_" for c in text)


def figure_path(base, build):
    """Figure filename ``<stem>[_<build>]<suffix>`` of one build."""
    base = Path(base)
    stem = base.stem
    if build and build != "unknown":
        stem = f"{stem}_{_sanitize(build)}"
    return base.with_name(f"{stem}{base.suffix}")


def parse_figure_name(path, stem):
    """Inverse of ``figure_path``: the (sanitized) build of a figure file.

    The build is ``"unknown"`` if the name has no build token. Returns ``None``
    for other files.
    """
    name = Path(path).stem
    if name == stem:
        return "unknown"
    if name.startswith(stem + "_") and len(name) > len(stem) + 1:
        return name[len(stem) + 1:]
    return None


def plot_all(store, output_path, window=DEFAULT_WINDOW):
    """Write one Plotly figure (as JSON) per build, spanning its whole history.

    Returns a list of ``(path, n_series, n_points)`` for every figure written.
    """
    if len(store.units) > 1:
        print(
            f"Warning: mixed units {sorted(store.units)} in the merged "
            'perflogs; labelling the axis as "s"',
            file=sys.stderr,
        )
    unit = next(iter(store.units)) if len(store.units) == 1 else "s"

    by_build = defaultdict(list)
    for build, descr, label, points in store.series():
        by_build[build].append((descr, label, points))

    outputs = []
    for build in sorted(by_build):
        series_list = sorted(by_build[build], key=lambda s: (s[0], s[1]))
        columns, plotted = build_commit_columns(series_list)
        fig = build_figure(build, columns, plotted, unit, window)
        out = figure_path(output_path, build)
        out.write_text(fig.to_json(), encoding="utf-8")
        n_points = sum(len(pts) for _, _, pts in plotted)
        outputs.append((out, len(plotted), n_points))

    return outputs


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output figure path (JSON); the build is appended to its stem.",
        required=True,
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help="Run directories to plot, as a regular expression: the value is "
        "split at its last '/' into a literal parent directory and a pattern "
        "full-matched against that directory's entries, e.g. "
        "'/path/benchmarks_[0-9]{2}_[0-9]{2}_[0-9]{4}_[0-9]+'.",
    )
    parser.add_argument(
        "--log",
        default=None,
        help="Explicit perflog path, bypassing discovery under --prefix. "
        "Plots that single log instead of the merged history.",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=DEFAULT_WINDOW,
        help="Commit columns (x-axis positions) shown when a figure opens "
        f"(default: {DEFAULT_WINDOW}).",
    )
    args = parser.parse_args(argv)
    if args.window < 1:
        parser.error("--window must be >= 1")

    log_paths, error = resolve_perflogs(args.prefix, args.log)
    if error is not None:
        parser.error(error)
    log_paths = cast(list[Path], log_paths)

    store, n_read = read_all_records(log_paths)
    if n_read == 0:
        print(
            f"None of the {len(log_paths)} matched perflog(s) could be read. "
            "Nothing to plot.",
            file=sys.stderr,
        )
        return 1

    if store.is_empty():
        listing = ", ".join(str(path) for path in log_paths)
        print(
            f"No performance data found in {listing}. Nothing to plot.",
            file=sys.stderr,
        )
        return 1

    print(f"Merged {n_read} perflog(s) into the timeline.")

    outputs = plot_all(store, Path(args.output), window=args.window)
    for path, n_series, n_points in outputs:
        print(
            f"Wrote {path} ({n_series} metric series,"
            f" {n_points} data point(s))."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
