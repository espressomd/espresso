#!/bin/bash

# Usage function
usage() {
    echo "Usage: $0 -p <prefix> -n <test_names> [-l] [--dry-run]"
    echo "  -p PREFIX       : Installation prefix for ReFrame benchmarks"
    echo "  -n TESTS        : Test case filter(s) for ReFrame (-n option)"
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
if [ -z "$PREFIX" ] || ([ -z "$TESTS" ] && [ "$LIST_MODE" = false ]); then
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

# Run ReFrame
reframe -C reframe_config.py \
        -c espresso_benchmarks.py \
        --prefix "$PREFIX" \
        --exec-policy serial \
        -n "$TESTS" \
        --performance-report \
        $RUN_OPTION
