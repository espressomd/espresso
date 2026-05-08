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
    echo "Usage: $0 -p <prefix> [-n <test_names>] [-l] [--dry-run]"
    echo "  -p PREFIX       : Installation prefix for ReFrame benchmarks"
    echo "  -n TESTS        : Optional test case filter(s) for ReFrame (-n option)"
    echo "  -l              : List available test cases (overrides -r/--dry-run)"
    echo "  --dry-run       : Optional flag to perform a dry run"
    exit 1
}

# Defaults
DRY_RUN=false
LIST_MODE=false

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -p)
            PREFIX="$2"
            shift 2
            ;;
        -n)
            TESTS="$2"
            shift 2
            ;;
        -l)
            LIST_MODE=true
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

# Determine final ReFrame action
if [ "$LIST_MODE" = true ]; then
    RUN_OPTION="-l"
elif [ "$DRY_RUN" = true ]; then
    RUN_OPTION="--dry-run"
else
    RUN_OPTION="-r"
fi

# Build optional -n argument
N_OPTION=""
if [ -n "$TESTS" ]; then
    N_OPTION="-n $TESTS"
fi

# Save sqlite storage database to prefix directory 
export RFM_SQLITE_DB_FILE="${PREFIX}/results.db"

# Run ReFrame
reframe -C reframe_config.py \
        -c espresso_benchmarks.py \
        --prefix "$PREFIX" \
        --exec-policy serial \
        $N_OPTION \
        --performance-report \
        $RUN_OPTION