#!/usr/bin/env bash
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

# Usage function
usage() {
    echo "Usage: $0 -p <prefix> [-n <test_names>] [-l] [--debug] [--dry-run]"
    echo "  -p PREFIX       : Installation prefix for ReFrame benchmarks"
    echo "  -n TESTS        : Optional ReFrame test-name filter; repeatable (selects the union)"
    echo "  -l              : List available test cases (overrides -r/--dry-run)"
    echo "  --debug         : On the ant cluster, run on the debug partition"
    echo "                    instead of the production compute nodes"
    echo "  --dry-run       : Optional flag to perform a dry run"
    exit 1
}

# Defaults
DRY_RUN=false
LIST_MODE=false
USE_DEBUG_PARTITION=false
# ReFrame -n filters; repeatable, selects the union.
N_OPTS=()

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -p)
            PREFIX="$2"
            shift 2
            ;;
        -n)
            N_OPTS+=(-n "$2")
            shift 2
            ;;
        -l)
            LIST_MODE=true
            shift
            ;;
        --debug)
            USE_DEBUG_PARTITION=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check required arguments
if [ -z "$PREFIX" ]; then
    usage
fi

#  Create prefix directory if it does not exist yet
mkdir -p "$PREFIX" || exit 1

# Determine final ReFrame action
if [ "$LIST_MODE" = true ]; then
    RUN_OPTION="-l"
elif [ "$DRY_RUN" = true ]; then
    RUN_OPTION="--dry-run"
else
    RUN_OPTION="-r"
fi

# Ensure ReFrame is available (e.g. the virtual environment is activated)
if ! command -v reframe >/dev/null 2>&1; then
    echo "Error: 'reframe' not found. Activate the ReFrame virtual environment" \
         "before running this script." >&2
    exit 1
fi

# Enable the results database and save it to the prefix directory
export RFM_ENABLE_RESULTS_STORAGE=1
export RFM_SQLITE_DB_FILE="${PREFIX}/results.db"
export RFM_PREFIX="$PREFIX"

# Select the ant_cluster partition. The tests turn this into a "+debug" or
# "+compute" constraint on valid_systems; the local system matches neither
# feature, so workstation runs are unaffected either way.
S_OPTS=(-S "use_debug_partition=${USE_DEBUG_PARTITION}")

# Run ReFrame
reframe -C reframe_config.py \
        -c espresso_benchmarks.py \
        --prefix "$PREFIX" \
        "${N_OPTS[@]}" \
        "${S_OPTS[@]}" \
        --performance-report \
        $RUN_OPTION

# Stop if ReFrame failed, so we do not plot from a missing/stale perflog.
reframe_status=$?
if [ "$reframe_status" -ne 0 ]; then
    echo "ReFrame exited with status ${reframe_status}; skipping plot generation." >&2
    exit "$reframe_status"
fi

# After a real run, render an SVG timeline of the benchmark performances.
# Skipped for list (-l) and dry runs, which produce no timing data.
if [ "$RUN_OPTION" = "-r" ]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

    # ReFrame writes the perflog to perflogs/<system>/<partition>/, naming those
    # directories after the system it ran on -- "local/default" on a workstation,
    # but e.g. "ant_cluster/debug" on a cluster -- so the path cannot be
    # hardcoded. Locate it, and put the SVG pages next to the log they came from.
    shopt -s nullglob
    perflogs=("${PREFIX}"/perflogs/*/*/EspressoBenchmark.log)
    shopt -u nullglob

    if [ "${#perflogs[@]}" -eq 0 ]; then
        echo "Warning: no perflog found under ${PREFIX}/perflogs;" \
             "skipping benchmark timeline plot" >&2
    else
        if [ "${#perflogs[@]}" -gt 1 ]; then
            echo "Warning: several perflogs found under ${PREFIX}/perflogs;" \
                 "plotting the most recently modified one" >&2
        fi
        # Most recently modified first, so the run that just finished wins.
        PERFLOG="$(ls -t "${perflogs[@]}" | head -n 1)"
        PERFLOG_DIR="$(dirname "$PERFLOG")"
        echo "Generating benchmark timeline from ${PERFLOG}"
        python3 "${SCRIPT_DIR}/plot_benchmarks.py" \
            --prefix "$PREFIX" \
            --log "$PERFLOG" \
            -o "${PERFLOG_DIR}/EspressoBenchmark.svg" \
            || echo "Warning: could not generate benchmark timeline plot" >&2
    fi
fi