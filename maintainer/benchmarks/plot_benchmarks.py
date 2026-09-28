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
SVG timeline of the recorded benchmark performances across ESPResSo commits.

The x-axis is ordinal: one tick per ESPResSo commit, ordered chronologically by
the earliest completion time at which the commit was benchmarked. Each tick is
labelled with the short commit hash and that build's completion time. A commit
that was benchmarked more than once gets one adjacent column per run. Each
benchmark metric (``MC``/``MD``/... or the unnamed default) is drawn as its own
line connecting its mean value across commits, with the confidence interval as
an error bar.

The suite can be run with several ESPResSo builds (the ``build_config`` column,
e.g. ``maxset``/``default``/``empty``); each build is plotted separately. Every
build's timeline is paginated into pages of at most ``--max-points`` commit
columns (default 20), and each page is written to its own SVG tagged with the
calendar span it covers:

    <stem>[_<build>]_<start>_<end>.svg     e.g. EspressoBenchmark_maxset_2026-03-10_2026-06-29.svg

Full pages carry an immutable ``<start>_<end>`` ISO-date span; the trailing,
not-yet-full page is written as ``<start>_current.svg`` and overwritten on later
runs until it reaches ``--max-points`` columns, at which point it is finalised
under its span name and a new ``_current`` page begins.

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
    python3 plot_benchmarks.py --prefix PREFIX_REGEX -o OUTPUT.svg
                               [--max-points N] [--log PERFLOG]
"""

import argparse
import csv
import io
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import cast

import matplotlib

matplotlib.use("Agg")  # headless: no display needed
import matplotlib.pyplot as plt  # noqa: E402

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
        return None, (
            f"Invalid --prefix regular expression {name_pattern!r}: {exc}"
        )
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


# Marker size of plotted points.
HIT_MARKERSIZE = 22

# Default number of commit columns (x-axis positions) per SVG page.
DEFAULT_MAX_POINTS = 20

# Self-contained hover layer injected into the SVG: styling, an (initially
# hidden) tooltip container, and vanilla JS that fills and positions it from the
# ``data-*`` attributes on each hit target. No external libraries; works offline.
# Interactivity only runs when the SVG is opened in a browser (or inlined into an
# HTML page) -- as a static image it degrades to the plain chart.
_INTERACTIVITY_TEMPLATE = """
<style type="text/css"><![CDATA[
  .benchmark-hit { cursor: pointer; }
  #benchmark-tooltip { pointer-events: none; }
  #benchmark-tooltip rect { fill: #1f2933; fill-opacity: 0.96;
    stroke: #0b0f14; stroke-width: 1; }
  #benchmark-tooltip text { font-family: "DejaVu Sans", sans-serif; }
  .bt-value { fill: #ffffff; font-weight: bold; }
  .bt-line { fill: #cbd2d9; }
]]></style>
<g id="benchmark-tooltip" visibility="hidden"></g>
<script type="text/javascript"><![CDATA[
(function(){
  var script = document.currentScript;
  var svg = (script && script.closest) ? script.closest("svg") : null;
  if (!svg) { var s = document.getElementsByTagName("svg"); svg = s[s.length-1]; }
  if (!svg) return;
  var NS = "http://www.w3.org/2000/svg";
  var tip = svg.querySelector("#benchmark-tooltip");
  if (!tip) return;

  function clear(el){ while (el.firstChild) el.removeChild(el.firstChild); }
  function attr(t, n){ return t.getAttribute("data-" + n) || ""; }

  function toSvg(cx, cy){
    var ctm = svg.getScreenCTM();
    if (!ctm) return {x: cx, y: cy};
    var p = svg.createSVGPoint(); p.x = cx; p.y = cy;
    return p.matrixTransform(ctm.inverse());
  }

  function show(target){
    var lines = [
      {t: attr(target, "value"),            cls: "bt-value"},
      {t: attr(target, "ci"),               cls: "bt-line"},
      {t: attr(target, "metric"),           cls: "bt-line"},
      {t: "commit " + attr(target, "commit"), cls: "bt-line"},
      {t: attr(target, "time"),             cls: "bt-line"}
    ].filter(function(l){ return l.t && l.t.trim() !== ""; });

    clear(tip);
    var PAD = 8, LH = 16, SW = 14, FS = 13, TX = PAD + SW + 6;

    var rect = document.createElementNS(NS, "rect");
    rect.setAttribute("rx", "4"); rect.setAttribute("ry", "4");
    rect.setAttribute("x", "0"); rect.setAttribute("y", "0");
    tip.appendChild(rect);

    var sw = document.createElementNS(NS, "line");
    sw.setAttribute("stroke", attr(target, "color") || "#ffffff");
    sw.setAttribute("stroke-width", "3");
    sw.setAttribute("x1", String(PAD)); sw.setAttribute("x2", String(PAD));
    sw.setAttribute("y1", String(PAD)); sw.setAttribute("y2", String(PAD + FS));
    tip.appendChild(sw);

    var txt = document.createElementNS(NS, "text");
    txt.setAttribute("x", String(TX)); txt.setAttribute("y", String(PAD + FS));
    txt.setAttribute("font-size", String(FS));
    for (var i = 0; i < lines.length; i++){
      var ts = document.createElementNS(NS, "tspan");
      ts.setAttribute("x", String(TX));
      if (i > 0) ts.setAttribute("dy", String(LH));
      ts.setAttribute("class", lines[i].cls);
      ts.appendChild(document.createTextNode(lines[i].t));
      txt.appendChild(ts);
    }
    tip.appendChild(txt);

    var bb = txt.getBBox();
    var w = bb.width + TX + PAD;
    var h = Math.max(LH * lines.length + 2 * PAD, bb.height + 2 * PAD);
    rect.setAttribute("width", String(w));
    rect.setAttribute("height", String(h));

    var r = target.getBoundingClientRect();
    var c = toSvg(r.left + r.width / 2, r.top + r.height / 2);
    var tx = c.x + 12, ty = c.y - h - 12;
    var vb = svg.viewBox && svg.viewBox.baseVal;
    if (vb && vb.width > 0){
      if (tx + w > vb.x + vb.width) tx = c.x - w - 12;
      if (tx < vb.x) tx = vb.x + 2;
      if (ty < vb.y) ty = c.y + 12;
    }
    tip.setAttribute("transform", "translate(" + tx + "," + ty + ")");
    tip.setAttribute("visibility", "visible");
  }

  function hide(){ tip.setAttribute("visibility", "hidden"); }
  function hit(e){ return e.target.closest ? e.target.closest(".benchmark-hit") : null; }

  svg.addEventListener("mouseover", function(e){ var t = hit(e); if (t) show(t); });
  svg.addEventListener("mouseout",  function(e){ if (hit(e)) hide(); });
  svg.addEventListener("focusin",   function(e){ var t = hit(e); if (t) show(t); });
  svg.addEventListener("focusout",  hide);
})();
]]></script>
"""


def _svg_attr_escape(value):
    """Escape a string for safe inclusion in an XML/SVG attribute value."""
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _inject_interactivity(svg, tips):
    """Turn the static matplotlib SVG into a hover-interactive one.

    ``tips`` maps each hit target's gid to its tooltip fields. Each corresponding
    ``<g id="gid">`` is tagged as a hit target (class, ``pointer-events``,
    ``tabindex``) and given ``data-*`` attributes; the CSS/tooltip/JS block is
    then appended as the last children of the root ``<svg>``.
    """
    for gid, info in tips.items():
        data_attrs = " ".join(
            f'data-{key}="{_svg_attr_escape(info[key])}"'
            for key in ("metric", "commit", "time", "value", "ci", "color")
        )
        old = f'<g id="{gid}">'
        new = (
            f'<g id="{gid}" class="benchmark-hit" pointer-events="all" '
            f'tabindex="0" {data_attrs}>'
        )
        svg = svg.replace(old, new, 1)

    insert_at = svg.rfind("</svg>")
    if insert_at == -1:
        return svg
    return svg[:insert_at] + _INTERACTIVITY_TEMPLATE + svg[insert_at:]


def render_page(build, columns, plotted, unit, output_path, span_label=None):
    """Render one page (a slice of commit columns) to an interactive SVG.

    ``columns``/``plotted`` are page-local: ``x`` indices are already 0-based
    within this page. The SVG is interactive: hovering (or keyboard-focusing) a
    data point reveals a tooltip with its completion time, value, CI, metric and
    commit.
    """
    fig, ax = plt.subplots(figsize=(12, 7))

    tips = {}  # gid -> tooltip fields
    hit_id = 0
    for descr, label, points in plotted:
        xs = [p["x"] for p in points]
        means = [p["mean"] for p in points]
        cis = [p["ci"] if p["ci"] is not None else 0.0 for p in points]
        container = ax.errorbar(
            xs,
            means,
            yerr=cis,
            marker="o",
            markersize=5,
            capsize=3,
            linewidth=1.2,
            label=series_legend(descr, label),
        )
        color = container[0].get_color(
        ) if container[0] is not None else "#000000"

        # One transparent, generously sized hit target per point, carrying the
        # tooltip data. Drawn on top (high zorder) so it always receives hovers.
        for p in points:
            gid = f"benchhit-{hit_id}"
            hit_id += 1
            # Pass an explicit colour so this invisible artist does not consume
            # the property cycle and shift the visible series' colours.
            marker, = ax.plot(
                [p["x"]],
                [p["mean"]],
                marker="o",
                markersize=HIT_MARKERSIZE,
                linestyle="none",
                color=color,
                markerfacecolor=(0, 0, 0, 0),
                markeredgecolor="none",
                alpha=0,
                zorder=6,
            )
            marker.set_gid(gid)
            when = p["time"]
            commit = p["commit"]
            tips[gid] = {
                "metric": series_legend(descr, label),
                "commit": commit if commit == "unknown" else commit[:12],
                "time": when.strftime("%Y-%m-%d %H:%M:%S") if when else "",
                "value": f"{p['mean']:.4g} {unit}",
                "ci": f"± {p['ci']:.3g} {unit}" if p["ci"] is not None else "",
                "color": color,
            }

    ax.set_yscale("log")
    ax.set_xlabel("ESPResSo commit (chronological)")
    ax.set_ylabel(f"Mean execution time ({unit})")
    title = "ESPResSo benchmark performance timeline"
    if build and build != "unknown":
        title += f"  —  build: {build}"
    if span_label:
        title += f"\n{span_label}"
    ax.set_title(title)
    ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.6)

    # One tick per commit column, labelled with the commit hash (completion time
    # is shown in the hover tooltip).
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels(
        [format_tick_label(c) for c in columns], rotation=30, ha="right",
        fontsize="small",
    )
    # Pad the x-range so single-column plots and edge markers are not clipped.
    ax.set_xlim(-0.5, len(columns) - 0.5)

    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize="small",
        title="Benchmark [metric]",
        borderaxespad=0.0,
    )

    fig.tight_layout()

    # Render to an in-memory SVG, inject the hover layer, then write it out.
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", bbox_inches="tight")
    plt.close(fig)
    svg = buf.getvalue().decode("utf-8")
    svg = _inject_interactivity(svg, tips)
    Path(output_path).write_text(svg, encoding="utf-8")


def slice_page(columns, plotted, start, end):
    """Restrict ``columns``/``plotted`` to column indices ``[start, end)``.

    ``x`` indices are rebased to 0 within the page; series with no points in the
    range are dropped.
    """
    cols_sub = columns[start:end]
    plotted_sub = []
    for descr, label, points in plotted:
        pts = [{**p, "x": p["x"] - start}
               for p in points if start <= p["x"] < end]
        if pts:
            plotted_sub.append((descr, label, pts))
    return cols_sub, plotted_sub


def _sanitize(text):
    """Turn a string into a filesystem-safe filename token."""
    return "".join(c if (c.isalnum() or c in "._-") else "_" for c in text)


def _page_output_path(base, build, start_dt, end_dt, is_current):
    """Page filename ``<stem>[_<build>]_<start>_<end|current><suffix>`` (ISO dates)."""
    date_fmt = "%Y-%m-%d"
    start = start_dt.strftime(date_fmt)
    end = "current" if is_current else end_dt.strftime(date_fmt)
    stem = base.stem
    if build and build != "unknown":
        stem = f"{stem}_{_sanitize(build)}"
    return base.with_name(f"{stem}_{start}_{end}{base.suffix}")


def plot_all(store, output_path, max_points=DEFAULT_MAX_POINTS):
    """Render paginated interactive SVGs: one per build, chunked into pages of
    ``max_points`` commit columns.

    Full pages get an immutable ``<start>_<end>`` span tag; a trailing page with
    fewer than ``max_points`` columns is written as ``<start>_current`` and
    overwritten on later runs until it fills. Superseded ``_current`` pages that
    this run did not (re)write are removed so nothing stale remains.

    Returns a list of ``(path, n_series, n_points)`` for every page written.
    """
    if len(store.units) > 1:
        print(
            f"Warning: mixed units {sorted(store.units)} in the merged "
            'perflogs; labelling the axis as "s"',
            file=sys.stderr,
        )
    unit = next(iter(store.units)) if len(store.units) == 1 else "s"
    date_fmt = "%Y-%m-%d"

    by_build = defaultdict(list)
    for build, descr, label, points in store.series():
        by_build[build].append((descr, label, points))

    outputs = []
    written = set()
    for build in sorted(by_build):
        series_list = sorted(by_build[build], key=lambda s: (s[0], s[1]))
        columns, plotted = build_commit_columns(series_list)
        n = len(columns)
        for start in range(0, n, max_points):
            end = min(start + max_points, n)
            is_current = end == n and (end - start) < max_points
            cols_sub, plotted_sub = slice_page(columns, plotted, start, end)
            span_start = cols_sub[0]["time"]
            span_end = cols_sub[-1]["time"]
            out = _page_output_path(
                output_path, build, span_start, span_end, is_current
            )
            span_label = f"{span_start.strftime(date_fmt)} … " + (
                "current" if is_current else span_end.strftime(date_fmt)
            )
            render_page(build, cols_sub, plotted_sub, unit, out, span_label)
            n_points = sum(len(pts) for _, _, pts in plotted_sub)
            outputs.append((out, len(plotted_sub), n_points))
            written.add(out.resolve())

    # Remove superseded "_current" pages (e.g. ones that have since filled up and
    # been re-emitted under their final span name) that this run did not write.
    base = Path(output_path)
    for stale in base.parent.glob(f"{base.stem}*_current{base.suffix}"):
        if stale.resolve() not in written:
            stale.unlink()

    return outputs


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output SVG path.",
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
        "--max-points",
        type=int,
        default=DEFAULT_MAX_POINTS,
        help="Commit columns (x-axis positions) per SVG page "
        f"(default: {DEFAULT_MAX_POINTS}).",
    )
    args = parser.parse_args(argv)
    if args.max_points < 1:
        parser.error("--max-points must be >= 1")

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

    outputs = plot_all(store, Path(args.output), max_points=args.max_points)
    for path, n_series, n_points in outputs:
        print(f"Wrote {path} ({n_series} metric series, {
              n_points} data point(s)).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
