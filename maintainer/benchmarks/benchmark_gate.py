#!/usr/bin/env python3
#
# Copyright (C) 2026 The ESPResSo project
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
Benchmark gate for the ParticleStore migration.

Runs the agreed LJ/P3M benchmark configurations, records timings to CSV
(same format as maintainer/benchmarks/benchmarks.py:write_report), and
compares against a committed baseline with a regression threshold.

The machine is shared: 'check-load' refuses to measure while other users'
jobs are running. See docs/superpowers/specs/
2026-07-03-array-based-particle-storage-design.md, section 5.

The busy check measures 'foreign' CPU usage: the CPU time consumed over a
short sampling interval by processes owned by users other than the current
one, expressed as a percentage of a single core. This is a two-sample
measurement of utime+stime from /proc/<pid>/stat, summed over all pids whose
owner uid differs from os.getuid(). We abandoned the previous 1-minute load
average because it self-trips: our own multi-rank (4-way MPI) benchmark
configurations drive the load average to ~4 and it decays over minutes, so
the gate would abort its own run after the first multi-rank configuration.
Foreign-CPU usage ignores our own load and only reacts to other users.

Exit codes:
  0 — success (check-load: quiet; run: all benchmarks completed; compare: no regression)
  1 — comparison failed (compare subcommand only): either a regression beyond
      the budget, or a configuration present in the baseline is missing from
      the current CSV
  2 — machine busy (check-load or run subcommand: foreign CPU usage exceeded threshold)
  3 — benchmark subprocess failed (run subcommand: mpirun returned non-zero)
  4 — output file already exists (run subcommand): refusing to append, because
      mixing rows from separate runs lets the min-of-means comparison pick up
      stale timings; delete or rename the file, or point --output elsewhere
"""

import argparse
import csv
import collections
import os
import pathlib
import subprocess
import sys
import time

BENCHMARK_DIRECTORY = pathlib.Path(__file__).resolve().parent

# The six agreed configurations (spec: success bar).
# fields: script, extra script arguments, MPI ranks, OMP threads
CONFIGURATIONS = [
    ("lj.py", ["--particles_per_core", "1000"], 1, 1),
    ("lj.py", ["--particles_per_core", "1000"], 4, 1),
    ("lj.py", ["--particles_per_core", "4000"], 1, 4),
    ("p3m.py", ["--particles_per_core", "1000"], 1, 1),
    ("p3m.py", ["--particles_per_core", "1000"], 4, 1),
    ("p3m.py", ["--particles_per_core", "4000"], 1, 4),
]


def _foreign_ticks():
    """Sum utime+stime ticks over processes not owned by the current user."""
    own_uid = os.getuid()
    total = 0
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            if os.stat(f"/proc/{entry}").st_uid == own_uid:
                continue
            with open(f"/proc/{entry}/stat") as f:
                line = f.read()
            # The comm field (2) may contain spaces and parentheses, so parse
            # the remaining fields from the last ')' in the line.
            fields = line[line.rindex(")") + 1:].split()
            # After comm, fields start at index 0 == field 3 (state); utime is
            # field 14 (index 11), stime is field 15 (index 12).
            total += int(fields[11]) + int(fields[12])
        except (FileNotFoundError, ProcessLookupError, PermissionError,
                ValueError, IndexError):
            # Process may vanish between listing and reading, or be
            # inaccessible; a pid present in only one sample contributes 0.
            continue
    return total


def foreign_cpu_percent(sample_seconds=1.0):
    """Foreign CPU usage as a percentage of one core over sample_seconds."""
    first = _foreign_ticks()
    time.sleep(sample_seconds)
    second = _foreign_ticks()
    delta = max(0, second - first)
    ticks_per_second = os.sysconf("SC_CLK_TCK")
    return (delta / ticks_per_second) / sample_seconds * 100.0


def check_load(max_foreign_cpu):
    foreign = foreign_cpu_percent()
    if foreign >= max_foreign_cpu:
        print(f"REFUSING to benchmark: foreign CPU usage {foreign:.1f}% "
              f">= threshold {max_foreign_cpu:.1f}% (shared machine).")
        print("Top CPU consumers:")
        subprocess.run(["ps", "-eo", "user,pcpu,pmem,comm", "--sort=-pcpu"],
                       check=False, stdout=sys.stdout)
        return False
    print(f"Machine quiet: foreign CPU usage {foreign:.1f}% < "
          f"{max_foreign_cpu:.1f}%.")
    return True


def run_benchmarks(pypresso, output, repetitions, max_foreign_cpu):
    pypresso = pathlib.Path(pypresso).resolve()
    output = pathlib.Path(output).resolve()
    if output.exists():
        print(f"ERROR: output file {output} already exists. Appending would "
              f"mix rows from separate runs and let the min-of-means "
              f"comparison pick up stale timings. Delete or rename it, or "
              f"choose a different --output.")
        return 4
    output.parent.mkdir(parents=True, exist_ok=True)
    for script, arguments, n_ranks, n_threads in CONFIGURATIONS:
        for repetition in range(repetitions):
            if not check_load(max_foreign_cpu):
                return 2
            command = ["mpirun", "-n", str(n_ranks), str(pypresso),
                       str(BENCHMARK_DIRECTORY / script)] + arguments + \
                [f"--output={output}"]
            environment = dict(os.environ)
            environment["OMP_NUM_THREADS"] = str(n_threads)
            environment["OMP_PROC_BIND"] = "false"
            print(f"[{script} {' '.join(arguments)} ranks={n_ranks} "
                  f"threads={n_threads}] repetition {repetition + 1}"
                  f"/{repetitions}")
            try:
                subprocess.run(command, check=True, env=environment)
            except subprocess.CalledProcessError as exc:
                config_name = (f"{script} {' '.join(arguments)} "
                               f"ranks={n_ranks} threads={n_threads}")
                print(f"ERROR: benchmark failed for configuration "
                      f"{config_name!r} (exit code {exc.returncode}).")
                return 3
    print(f"Wrote timings to {output}")
    return 0


def read_rows(filepath):
    """Group rows by configuration; keep the minimum 'mean' per group."""
    groups = collections.OrderedDict()
    with open(filepath, newline="") as f:
        for row in csv.DictReader(f):
            key = (row["script"], row["arguments"],
                   row["ranks"], row["threads"], row["label"])
            mean = float(row["mean"])
            if key not in groups or mean < groups[key]:
                groups[key] = mean
    return groups


def compare(baseline_path, current_path, max_regression):
    baseline = read_rows(baseline_path)
    current = read_rows(current_path)
    missing = [k for k in baseline if k not in current]
    if missing:
        print(f"FAIL: configurations missing from {current_path}:")
        for key in missing:
            print(f"  {key}")
        return 1
    failed = False
    print(f"{'configuration':70s} {'baseline':>12s} {'current':>12s} "
          f"{'ratio':>8s}")
    for key, baseline_mean in baseline.items():
        current_mean = current[key]
        ratio = current_mean / baseline_mean
        marker = ""
        if ratio > 1.0 + max_regression:
            marker = "  <-- REGRESSION"
            failed = True
        label = key[4]
        name = f"{key[0]} {key[1]} ranks={key[2]} threads={key[3]}"
        if label:
            name += f" label={label}"
        print(f"{name:70s} {baseline_mean:12.3e} {current_mean:12.3e} "
              f"{ratio:8.3f}{marker}")
    if failed:
        print(f"FAIL: regression beyond {100 * max_regression:.0f}% budget.")
        return 1
    print(f"PASS: all configurations within {100 * max_regression:.0f}% "
          f"of baseline.")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    parser_load = subparsers.add_parser("check-load")
    parser_load.add_argument("--max-foreign-cpu", type=float, default=50.0)

    parser_run = subparsers.add_parser("run")
    parser_run.add_argument("--pypresso", required=True)
    parser_run.add_argument("--output", required=True)
    parser_run.add_argument("--repetitions", type=int, default=3)
    parser_run.add_argument("--max-foreign-cpu", type=float, default=50.0)

    parser_compare = subparsers.add_parser("compare")
    parser_compare.add_argument("--baseline", required=True)
    parser_compare.add_argument("--current", required=True)
    parser_compare.add_argument("--max-regression", type=float, default=0.05)

    args = parser.parse_args()
    if args.command == "check-load":
        return 0 if check_load(args.max_foreign_cpu) else 2
    if args.command == "run":
        return run_benchmarks(args.pypresso, args.output, args.repetitions,
                              args.max_foreign_cpu)
    if args.command == "compare":
        return compare(args.baseline, args.current, args.max_regression)
    return 1


if __name__ == "__main__":
    sys.exit(main())
