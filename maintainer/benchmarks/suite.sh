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
    echo "  -p PREFIX       : Installation prefix for ReFrame benchmarks. Each run"
    echo "                    writes to its own <PREFIX>_dd_mm_yyyy_<n> directory"
    echo "                    (n counts the runs of that day); the timeline plot"
    echo "                    is assembled from all of them"
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

# The prefix is a naming stem, not a directory: a trailing slash would put the
# run directory *inside* it and desynchronise it from the pattern derived below.
while [ "${PREFIX%/}" != "$PREFIX" ]; do
    PREFIX="${PREFIX%/}"
done
if [ -z "$PREFIX" ]; then
    echo "Error: -p must be a path, not '/'" >&2
    exit 1
fi

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

# Construct unique benchmark directory name
TODAY="$(date +%d_%m_%Y)"
RUN_SUFFIX=""
if [ "$USE_DEBUG_PARTITION" = true ]; then
    RUN_SUFFIX="${RUN_SUFFIX}_DEBUG"
fi
if [ "$DRY_RUN" = true ]; then
    RUN_SUFFIX="${RUN_SUFFIX}_DRYRUN"
fi
if [ "$LIST_MODE" = true ]; then
    RUN_SUFFIX="${RUN_SUFFIX}_LIST"
fi

mkdir -p "$(dirname "$PREFIX")" || exit 1

# Increment run counter if directory already exists
RUN_INDEX=1
while true; do
    RUN_DIR="${PREFIX}_${TODAY}_${RUN_INDEX}${RUN_SUFFIX}"
    if mkdir "$RUN_DIR" 2>/dev/null; then
        break
    fi
    if [ ! -d "$RUN_DIR" ]; then
        echo "Error: cannot create ${RUN_DIR}" >&2
        exit 1
    fi
    RUN_INDEX=$((RUN_INDEX + 1))
done
echo "Writing benchmark artefacts to ${RUN_DIR}"

# Enable the results database and save it to the run directory
export RFM_ENABLE_RESULTS_STORAGE=1
export RFM_SQLITE_DB_FILE="${RUN_DIR}/results.db"
export RFM_PREFIX="$RUN_DIR"

S_OPTS=(-S "use_debug_partition=${USE_DEBUG_PARTITION}")

# Run ReFrame
reframe -C reframe_config.py \
        -c espresso_benchmarks.py \
        --prefix "$RUN_DIR" \
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


if [ "$RUN_OPTION" = "-r" ]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

    shopt -s nullglob
    perflogs=("${RUN_DIR}"/perflogs/*/*/EspressoBenchmark.log)
    shopt -u nullglob

    if [ "${#perflogs[@]}" -eq 0 ]; then
        # Nothing ran (or everything failed), but the sibling run directories
        # still hold the history worth plotting.
        echo "Warning: no perflog found under ${RUN_DIR}/perflogs;" \
             "plotting the history of the earlier runs only" >&2
        PLOT_DIR="$RUN_DIR"
    else
        # A fresh run directory holds a single log, written by this run.
        PLOT_DIR="$(dirname "${perflogs[0]}")"
        echo "Generating benchmark timeline from ${perflogs[0]}"
    fi

    # Plot the whole history of past runs with shared prefix
    PREFIX_BASE_RE="$(python3 -c 'import re, sys; print(re.escape(sys.argv[1]))' \
                      "$(basename "$PREFIX")")"
    PLOT_PREFIX="$(dirname "$PREFIX")/${PREFIX_BASE_RE}_[0-9]{2}_[0-9]{2}_[0-9]{4}_[0-9]+${RUN_SUFFIX}"

    python3 "${SCRIPT_DIR}/plot_benchmarks.py" \
        --prefix "$PLOT_PREFIX" \
        -o "${PLOT_DIR}/EspressoBenchmark.svg" \
        || { echo "Error: benchmark timeline plot failed" >&2; exit 1; }

    # Create link to newest benchmark plots after every run 
    ln -sfn "$(basename "$RUN_DIR")" "${PREFIX}_latest"
fi