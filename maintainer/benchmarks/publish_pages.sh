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

# Run the benchmark suite and publish its history page to a GitLab Pages
# repository. Meant to be called from cron.
#
# Every call pushes one commit: the newest benchmark_history.html (as
# public/index.html) with its plotly.js, and status/status.txt with the result
# of the suite. The CI of the pages repository deploys public/ and fails a
# status job when the suite failed, so GitLab's pipeline mails report both
# outcomes. The push authenticates with whatever credential helper the clone
# is configured with; git never prompts.

usage() {
    echo "Usage: $0 -p <prefix> --pages-repo <dir> [--branch <name>]" \
         "[-- <suite.sh options>]"
    echo "  -p PREFIX        : Prefix passed to suite.sh; the page is taken from"
    echo "                     <PREFIX>_latest<suffix>/"
    echo "  --pages-repo DIR : Clone of the GitLab Pages repository"
    echo "  --branch NAME    : Branch to publish to (default: main)"
    echo "  -- OPTIONS       : Further options for suite.sh, e.g. --debug or"
    echo "                     --rebuild-latest; -l and --dry-run publish nothing"
    echo "                     and are rejected"
    exit 1
}

PREFIX=""
PAGES_REPO=""
BRANCH="main"
SUITE_ARGS=()
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p)
            [ $# -ge 2 ] || { echo "Error: -p requires a value" >&2; usage; }
            PREFIX="$2"
            shift 2
            ;;
        --pages-repo)
            [ $# -ge 2 ] || {
                echo "Error: --pages-repo requires a value" >&2; usage; }
            PAGES_REPO="$2"
            shift 2
            ;;
        --branch)
            [ $# -ge 2 ] || { echo "Error: --branch requires a value" >&2; usage; }
            BRANCH="$2"
            shift 2
            ;;
        --)
            shift
            SUITE_ARGS=("$@")
            break
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

if [ -z "$PREFIX" ] || [ -z "$PAGES_REPO" ]; then
    usage
fi
while [ "${PREFIX%/}" != "$PREFIX" ]; do
    PREFIX="${PREFIX%/}"
done
if ! git -C "$PAGES_REPO" rev-parse --git-dir >/dev/null 2>&1; then
    echo "Error: $PAGES_REPO is not a git repository" >&2
    exit 1
fi

# suite.sh names the _latest directory after these options; keep this in sync
# with its RUN_SUFFIX.
RUN_SUFFIX=""
debug=false
commit=false
prebuilt=false
for arg in "${SUITE_ARGS[@]}"; do
    case "$arg" in
        -p)
            echo "Error: pass -p to $0, not after --" >&2
            exit 1
            ;;
        -l|--dry-run)
            echo "Error: $arg produces no page to publish" >&2
            exit 1
            ;;
        --debug) debug=true ;;
        --commit) commit=true ;;
        --prebuilt-module) prebuilt=true ;;
    esac
done
[ "$debug" = true ] && RUN_SUFFIX="${RUN_SUFFIX}_DEBUG"
[ "$commit" = true ] && RUN_SUFFIX="${RUN_SUFFIX}_COMMIT"
[ "$prebuilt" = true ] && RUN_SUFFIX="${RUN_SUFFIX}_PREBUILT"
LATEST_DIR="${PREFIX}_latest${RUN_SUFFIX}"

# Cron starts the next run even if the previous one still benchmarks.
exec 9>"$(git -C "$PAGES_REPO" rev-parse --absolute-git-dir)/publish.lock"
if ! flock -n 9; then
    echo "Error: another publish into $PAGES_REPO is still running" >&2
    exit 1
fi

SUITE_LOG="$(mktemp)"
trap 'rm -f "$SUITE_LOG"' EXIT

"${SCRIPT_DIR}/suite.sh" -p "$PREFIX" "${SUITE_ARGS[@]}" 2>&1 | tee "$SUITE_LOG"
suite_status=${PIPESTATUS[0]}
if [ "$suite_status" -eq 0 ]; then
    result=passed
else
    result=failed
fi
echo "Benchmark suite ${result} (exit status ${suite_status})"

# Fails instead of waiting for a password nobody will type.
export GIT_TERMINAL_PROMPT=0
pages() { git -C "$PAGES_REPO" "$@"; }

# The clone only mirrors the remote, so local leftovers are discarded.
if ! { pages fetch --quiet origin "$BRANCH" &&
       pages reset --quiet --hard "origin/${BRANCH}" &&
       pages clean --quiet -fd; }; then
    echo "Error: cannot update $PAGES_REPO from origin/${BRANCH};" \
         "nothing was published" >&2
    exit 1
fi

# Without a page (e.g. a failure before the first run) the published one stays.
if [ -f "${LATEST_DIR}/benchmark_history.html" ]; then
    mkdir -p "${PAGES_REPO}/public"
    cp "${LATEST_DIR}/benchmark_history.html" "${PAGES_REPO}/public/index.html"
    shopt -s nullglob
    rm -f "${PAGES_REPO}"/public/plotly-*.min.js
    for js in "${LATEST_DIR}"/plotly-*.min.js; do
        cp "$js" "${PAGES_REPO}/public/"
    done
    shopt -u nullglob
else
    echo "Warning: no ${LATEST_DIR}/benchmark_history.html; the published" \
         "page is left unchanged" >&2
fi

run_name=""
if [ -L "${LATEST_DIR}/run" ]; then
    run_name="$(basename "$(readlink "${LATEST_DIR}/run")")"
fi
mkdir -p "${PAGES_REPO}/status"
{
    echo "result=${result}"
    echo "exit_code=${suite_status}"
    echo "finished=$(date -Iseconds)"
    echo "run=${run_name}"
    echo "host=$(hostname)"
} > "${PAGES_REPO}/status/status.txt"
# The CI status job prints this tail, so it ends up in the failure mail's log.
if [ "$result" = failed ]; then
    tail -n 200 "$SUITE_LOG" > "${PAGES_REPO}/status/suite-tail.log"
else
    rm -f "${PAGES_REPO}/status/suite-tail.log"
fi

# Use the identity configured in the clone, if any.
identity=()
if [ -z "$(pages config user.email)" ]; then
    identity=(-c "user.name=ESPResSo benchmark bot"
              -c "user.email=$(id -un)@$(hostname)")
fi
message="Benchmarks $(date +%Y-%m-%d): ${result} (exit ${suite_status})"
[ -n "$run_name" ] && message="${message}, run ${run_name}"

if ! { pages add -A &&
       pages "${identity[@]}" commit --quiet --allow-empty -m "$message" &&
       pages push --quiet origin "HEAD:${BRANCH}"; }; then
    echo "Error: publishing to origin/${BRANCH} failed; no pipeline mail" \
         "will be sent for this run" >&2
    exit 1
fi
echo "Published: ${message}"

exit "$suite_status"
