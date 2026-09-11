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
    echo "Usage: $0 -p <prefix> [-n <test_names>] [--commit <ref>]" \
         "[--prebuilt-init <script> --prebuilt-module <name>]" \
         "[-l] [--debug] [--dry-run]"
    echo "  -p PREFIX       : Installation prefix for ReFrame benchmarks. Each run"
    echo "                    writes to its own <PREFIX>_dd_mm_yyyy_<n> directory"
    echo "                    (n counts the runs of that day); the timeline plot"
    echo "                    is assembled from all of them"
    echo "  -n TESTS        : Optional ReFrame test-name filter; repeatable (selects the union)"
    echo "  -l              : List available test cases (overrides -r/--dry-run)"
    echo "  --commit REF    : Benchmark a specific ESPResSo commit, tag or branch"
    echo "                    instead of the newest commit of the default branch."
    echo "  --prebuilt-init SCRIPT and --prebuilt-module NAME :"
    echo "                    Benchmark a pre-built ESPResSo instead of building one."
    echo "                    Use -n to exclude cases the build cannot run."
    echo "  --debug         : On the ant cluster, run on the debug partition"
    echo "                    instead of the production compute nodes"
    echo "  --dry-run       : Optional flag to perform a dry run"
    exit 1
}

# Defaults
DRY_RUN=false
LIST_MODE=false
USE_DEBUG_PARTITION=false
# Optional git ref (commit/tag/branch) to benchmark, empty = default branch.
COMMIT_REF=""
# Optional pre-built ESPResSo: script to source and module to load.
PREBUILT_INIT=""
PREBUILT_MODULE=""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# ReFrame -n filters; repeatable, selects the union.
N_OPTS=()

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -p)
            # "shift 2" does nothing when only one token is left, which would
            # spin this loop forever, so check before consuming a value.
            [ $# -ge 2 ] || { echo "Error: -p requires a value" >&2; usage; }
            PREFIX="$2"
            shift 2
            ;;
        -n)
            [ $# -ge 2 ] || { echo "Error: -n requires a value" >&2; usage; }
            N_OPTS+=(-n "$2")
            shift 2
            ;;
        --commit)
            [ $# -ge 2 ] || { echo "Error: --commit requires a value" >&2; usage; }
            COMMIT_REF="$2"
            if [ -z "$COMMIT_REF" ]; then
                echo "Error: --commit value is empty" >&2
                exit 1
            fi
            shift 2
            ;;
        --prebuilt-init)
            [ $# -ge 2 ] || {
                echo "Error: --prebuilt-init requires a value" >&2; usage; }
            PREBUILT_INIT="$2"
            if [ -z "$PREBUILT_INIT" ]; then
                echo "Error: --prebuilt-init value is empty" >&2
                exit 1
            fi
            shift 2
            ;;
        --prebuilt-module)
            [ $# -ge 2 ] || {
                echo "Error: --prebuilt-module requires a value" >&2; usage; }
            PREBUILT_MODULE="$2"
            if [ -z "$PREBUILT_MODULE" ]; then
                echo "Error: --prebuilt-module value is empty" >&2
                exit 1
            fi
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

# The two pre-built flags describe one installation and are useless apart.
if { [ -n "$PREBUILT_INIT" ] && [ -z "$PREBUILT_MODULE" ]; } ||
   { [ -z "$PREBUILT_INIT" ] && [ -n "$PREBUILT_MODULE" ]; }; then
    echo "Error: --prebuilt-init and --prebuilt-module must be given together" >&2
    exit 1
fi
if [ -n "$PREBUILT_MODULE" ] && [ -n "$COMMIT_REF" ]; then
    echo "Error: --commit selects a source revision to build, which is" \
         "meaningless with a pre-built ESPResSo; use one or the other." >&2
    exit 1
fi

# Catch mistyped refs now (pointless without a build)
if [ -n "$COMMIT_REF" ]; then
    if git -C "$SCRIPT_DIR" rev-parse --verify --quiet \
            "${COMMIT_REF}^{commit}" >/dev/null 2>&1; then
        : # known to this clone
    elif git ls-remote --exit-code \
            https://github.com/espressomd/espresso.git \
            "$COMMIT_REF" >/dev/null 2>&1; then
        : # a branch or tag of the upstream repository
    elif [[ "$COMMIT_REF" =~ ^[0-9a-fA-F]{7,40}$ ]]; then
        echo "Note: cannot verify commit ${COMMIT_REF} before the build;" \
             "the checkout in the build job will resolve it." >&2
    else
        echo "Error: '${COMMIT_REF}' is not a commit, tag or branch of" \
             "espressomd/espresso." >&2
        exit 1
    fi
fi

# Test the pre-built installation
if [ -n "$PREBUILT_MODULE" ]; then
    if ! bash -c 'source "$1" >/dev/null 2>&1 && module load "$2" >/dev/null 2>&1' \
            _ "$PREBUILT_INIT" "$PREBUILT_MODULE"; then
        echo "Error: could not activate the pre-built ESPResSo." \
             "Check that '${PREBUILT_INIT}' is sourceable and that module" \
             "'${PREBUILT_MODULE}' exists:" >&2
        bash -c 'source "$1" && module load "$2"' \
            _ "$PREBUILT_INIT" "$PREBUILT_MODULE" >&2
        exit 1
    fi

    ESPRESSO_ID="$(bash -c 'source "$1" >/dev/null 2>&1 && module load "$2" >/dev/null 2>&1 &&
        python3 -c "import espressomd.version as v
print(v.git_commit() or (v.friendly() + \"-\" + \"$2\"))"' \
        _ "$PREBUILT_INIT" "$PREBUILT_MODULE" 2>/dev/null | tail -n 1)"

    if [ -z "$ESPRESSO_ID" ]; then
        # The module loads but this node cannot import espressomd
        echo "Warning: could not query the ESPResSo version from" \
             "'${PREBUILT_MODULE}'; recording the module name instead." >&2
        ESPRESSO_ID="prebuilt-${PREBUILT_MODULE}"
    fi

    if ! bash -c 'source "$1" >/dev/null 2>&1 && module load "$2" >/dev/null 2>&1 &&
            python3 -c "import pint"' \
            _ "$PREBUILT_INIT" "$PREBUILT_MODULE" >/dev/null 2>&1; then
        echo "Warning: 'pint' is not available in the pre-built Python;" \
             "mc_acid_base_reservoir.py will fail. Exclude it with -n." >&2
    fi
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
if [ -n "$COMMIT_REF" ]; then
    RUN_SUFFIX="${RUN_SUFFIX}_COMMIT"
fi
if [ -n "$PREBUILT_MODULE" ]; then
    RUN_SUFFIX="${RUN_SUFFIX}_PREBUILT"
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
if [ -n "$COMMIT_REF" ]; then
    S_OPTS+=(-S "espresso_ref=${COMMIT_REF}")
fi
if [ -n "$PREBUILT_MODULE" ]; then
    export ESPRESSO_PREBUILT_INIT="$PREBUILT_INIT"
    export ESPRESSO_PREBUILT_MODULE="$PREBUILT_MODULE"
    export RFM_RESOLVE_MODULE_CONFLICTS=0
    S_OPTS+=(-S "espresso_commit=${ESPRESSO_ID}")
    S_OPTS+=(-S "espresso_ref=${PREBUILT_MODULE}")
fi

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
    shopt -s nullglob
    perflogs=("${RUN_DIR}"/perflogs/*/*/EspressoBenchmark.log)
    shopt -u nullglob

    NO_PERFLOG=false
    if [ "${#perflogs[@]}" -eq 0 ]; then
        echo "Warning: no perflog found under ${RUN_DIR}/perflogs;" \
             "plotting the history of the earlier runs only" >&2
        NO_PERFLOG=true
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
    ln -sfn "$(basename "$RUN_DIR")" "${PREFIX}_latest${RUN_SUFFIX}"

    if [ "$NO_PERFLOG" = true ]; then
        echo "Error: this run produced no benchmark results" >&2
        exit 1
    fi
fi