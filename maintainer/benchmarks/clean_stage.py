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
Strip the ESPResSo build trees that ReFrame leaves behind in a run's stage dir.

A directory is cleaned only if all of these hold, otherwise it is reported and
left untouched:

  1. RUN_DIR is a real directory, not a symlink, not / or $HOME, and at least
     three path components deep. It contains real stage/ and output/
     directories (ReFrame always creates both).
  2. The candidate sits at exactly stage/<system>/<partition>/<environ>/<test>.
     Nothing is searched recursively and no symlink is followed on the way.
  3. Its name is BuildEspresso_<8 hex digits>, and it holds a real .git
     directory whose origin is espressomd/espresso, or the marker an earlier
     interrupted clean-up left behind.
  4. Unless --no-activity-check is given, nothing in the run directory was
     written during the last 30 minutes.

Exit status: 0 on success (including "nothing to clean"), 1 if some deletion
did not complete, 2 if a safety check refused a run directory.

Usage:
    python3 clean_stage.py RUN_DIR [RUN_DIR ...]            # report only
    python3 clean_stage.py --delete RUN_DIR [RUN_DIR ...]   # actually delete
"""

import argparse
import re
import shutil
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

# Marker written into a build directory before anything is deleted from it
# to let an interrupted clean-up be resumed.
MARKER = ".clean_stage"

BUILD_DIR_RE = re.compile(r"BuildEspresso_[0-9a-f]{8}")

ESPRESSO_ORIGIN_RE = re.compile(
    r"(https://github\.com/|git@github\.com:|ssh://git@github\.com/)"
    r"espressomd/espresso(\.git)?/?"
)

KEEP_NAMES = frozenset({".rfm_mark", ".rfm_testcase.json", "benchmarks.csv",
                        MARKER})
KEEP_RE = re.compile(r"rfm_[A-Za-z_]+\.(sh|out|err)")

# Depth of a test's stage directory below stage/:
# <system>/<partition>/<environ>/<test>.
STAGE_DEPTH = 4

# "/" plus at least three components, so /tikhome/<user> is refused.
MIN_PATH_PARTS = 4

ACTIVITY_WINDOW_S = 30 * 60


class Refused(Exception):
    """The run directory fails a safety check and is not touched at all."""


class NotABuildTree(Exception):
    """A candidate directory that does not qualify for cleaning."""


@dataclass
class RunResult:
    path: Path
    refused: str | None = None
    note: str | None = None
    cleaned: list[tuple[str, int, int]] = field(default_factory=list)
    skipped: list[tuple[str, str]] = field(default_factory=list)
    failures: list[tuple[str, str, OSError]] = field(default_factory=list)


def is_real_dir(path: Path) -> bool:
    return path.is_dir() and not path.is_symlink()


def is_real_file(path: Path) -> bool:
    return path.is_file() and not path.is_symlink()


def is_kept(name: str) -> bool:
    return name in KEEP_NAMES or KEEP_RE.fullmatch(name) is not None


def human(num_bytes: float) -> str:
    for unit in ("B", "KiB", "MiB"):
        if num_bytes < 1024:
            return f"{num_bytes:.1f} {unit}" if unit != "B" else f"{num_bytes:.0f} B"
        num_bytes /= 1024
    return f"{num_bytes:.1f} GiB"


def inode_bytes(path: Path, seen: set[tuple[int, int]]) -> int:
    """Disk space of one inode, counting a hard-linked file only once."""
    try:
        st = path.lstat()
    except OSError:
        return 0  # vanished or unreadable; it is not ours to count
    if st.st_nlink > 1 and not path.is_dir():
        if (st.st_dev, st.st_ino) in seen:
            return 0
        seen.add((st.st_dev, st.st_ino))
    return st.st_blocks * 512


def disk_usage(path: Path, seen: set[tuple[int, int]]) -> int:
    """Bytes on disk taken by ``path``, never following a symlink."""
    total = inode_bytes(path, seen)
    if is_real_dir(path):
        for root, dirs, files in path.walk():
            total += sum(inode_bytes(root / name, seen)
                         for name in (*dirs, *files))
    return total


def check_run_dir(run_dir: Path) -> bool:
    """
    Returns False when there is no stage/ at all -- an -l run, or one that
    failed before staging anything -- which simply means nothing to clean.
    """
    if run_dir.is_symlink():
        raise Refused("is a symlink")
    if not run_dir.is_dir():
        raise Refused("is not a directory")

    resolved = run_dir.resolve()
    if resolved == Path("/"):
        raise Refused("is the filesystem root")
    if resolved == Path.home().resolve():
        raise Refused("is your home directory")
    if len(resolved.parts) < MIN_PATH_PARTS:
        raise Refused(f"is too close to the filesystem root ({resolved})")

    if not (run_dir / "stage").exists():
        return False
    for name in ("stage", "output"):
        if not is_real_dir(run_dir / name):
            raise Refused(f"{name}/ is missing or not a real directory; "
                          "this does not look like a ReFrame run")
    return True


def test_dirs(stage: Path) -> list[Path]:
    """
    Real directories at exactly stage/<system>/<partition>/<environ>/<test>.
    """
    level = [stage]
    for _ in range(STAGE_DEPTH):
        level = [child for parent in level for child in sorted(parent.iterdir())
                 if is_real_dir(child)]
    return level


def recent_activity(run_dir: Path, tests: list[Path], now: float) -> str | None:
    """The first log written within the activity window, or None."""
    logs = [run_dir / "reframe.log",
            *(log for test in tests for log in test.glob("rfm_*.out"))]
    for log in logs:
        if is_real_file(log) and now - log.stat().st_mtime < ACTIVITY_WINDOW_S:
            return str(log.relative_to(run_dir))
    return None


def origin_url(config: Path) -> str | None:
    """The url of [remote "origin"] in a .git/config, read without running git."""
    if not is_real_file(config):
        return None
    section = None
    for raw in config.read_text(errors="replace").splitlines():
        line = raw.strip()
        if line.startswith("["):
            section = line
        elif section == '[remote "origin"]':
            key, sep, value = line.partition("=")
            if sep and key.strip() == "url":
                return value.strip()
    return None


def origin_of(build_dir: Path) -> str:
    """The verified ESPResSo origin of a build tree; raise NotABuildTree if not."""
    marker = build_dir / MARKER
    if marker.exists() or marker.is_symlink():
        if not is_real_file(marker):
            raise NotABuildTree(f"{MARKER} is not a regular file")
        for line in marker.read_text().splitlines():
            if line.startswith("origin="):
                return line.removeprefix("origin=")
        return "(resumed)"

    git = build_dir / ".git"
    if not (git.exists() or git.is_symlink()):
        raise NotABuildTree("no .git directory")
    if not is_real_dir(git):
        # A symlink, or a worktree-style .git file.
        raise NotABuildTree(".git is not a real directory")

    url = origin_url(git / "config")
    if url is None:
        raise NotABuildTree("no origin remote in .git/config")
    if not ESPRESSO_ORIGIN_RE.fullmatch(url):
        raise NotABuildTree(f"origin={url} is not an ESPResSo remote")
    return url


def strip_build_dir(build_dir: Path, origin: str, rel: str, delete: bool,
                    result: RunResult) -> None:
    doomed = sorted((p for p in build_dir.iterdir() if not is_kept(p.name)),
                    key=lambda p: (p.name == ".git", p.name))
    seen: set[tuple[int, int]] = set()
    before = sum(disk_usage(p, seen) for p in doomed)
    if not delete:
        result.cleaned.append((rel, before, len(doomed)))
        return

    (build_dir / MARKER).write_text(
        f"origin={origin}\ncleaned={time.strftime('%Y-%m-%dT%H:%M:%S')}\n")

    def collect(_func: object, path: str, exc: BaseException) -> None:
        if not isinstance(exc, FileNotFoundError):  # a concurrent clean-up
            what = str(Path(path).relative_to(build_dir))
            result.failures.append((rel, what, exc))  # type: ignore

    for path in doomed:
        if is_real_dir(path):
            shutil.rmtree(path, onexc=collect)
        else:
            try:
                path.unlink(missing_ok=True)
            except OSError as exc:
                collect(path.unlink, str(path), exc)

    left = sum(disk_usage(p, set()) for p in doomed
               if p.exists() or p.is_symlink())
    result.cleaned.append((rel, before - left, len(doomed)))


def process_run_dir(run_dir: Path, delete: bool,
                    check_activity: bool) -> RunResult:
    result = RunResult(run_dir)

    try:
        if not check_run_dir(run_dir):
            result.note = "no stage/ directory: nothing to clean"
            return result
        tests = test_dirs(run_dir / "stage")
        if check_activity and (active := recent_activity(run_dir, tests,
                                                         time.time())):
            raise Refused(f"{active} was written in the last "
                          f"{ACTIVITY_WINDOW_S // 60} minutes; the suite may "
                          "still be running (use --no-activity-check to "
                          "override)")
    except (Refused, OSError) as exc:
        result.refused = str(exc)
        return result

    for build_dir in tests:
        if not BUILD_DIR_RE.fullmatch(build_dir.name):
            continue
        rel = str(build_dir.relative_to(run_dir))
        try:
            origin = origin_of(build_dir)
        except (NotABuildTree, OSError) as exc:
            result.skipped.append((rel, str(exc)))
            continue
        strip_build_dir(build_dir, origin, rel, delete, result)

    return result


def report(results: list[RunResult], delete: bool) -> None:
    verb = "freed" if delete else "would free"
    total = 0
    for res in results:
        print(f"{res.path}:")
        if res.refused:
            print(f"  REFUSED: {res.refused}")
            continue
        if res.note:
            print(f"  {res.note}")
        partial = {rel for rel, _what, _exc in res.failures}
        for rel, size, n in res.cleaned:
            total += size
            if not delete:
                label = "would clean"
            elif rel in partial:
                label = "PARTLY cleaned"  # rerun to finish; see FAILED below
            else:
                label = "cleaned"
            print(f"  {label}  {rel}  {human(size)} {verb} ({n} entries)")
        for rel, why in res.skipped:
            print(f"  skipped  {rel}: {why}")
        for rel, what, exc in res.failures:
            print(f"  FAILED   {rel}: {what}: {exc}")
        if not (res.cleaned or res.skipped or res.failures or res.note):
            print("  no build trees to clean")

    n_dirs = sum(len(r.cleaned) for r in results)
    suffix = "" if delete else "  (report only; pass --delete to act)"
    print(f"total: {human(total)} {verb} in {n_dirs} build "
          f"director{'y' if n_dirs == 1 else 'ies'}{suffix}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("run_dirs", nargs="+", type=Path, metavar="RUN_DIR",
                        help="Benchmark run directory, the one holding "
                        "stage/ and output/.")
    parser.add_argument("--delete", action="store_true",
                        help="Actually delete. Without it, only report.")
    parser.add_argument("--no-activity-check", action="store_true",
                        help="Skip the recent-activity check. For callers that "
                        "know the run has finished, i.e. suite.sh.")
    args = parser.parse_args(argv)

    # Ensure that rmtree unlinks symlinks rather than following them
    if args.delete and not shutil.rmtree.avoids_symlink_attacks:
        print("Error: shutil.rmtree on this platform is not resistant to "
              "symlink attacks; refusing to delete.", file=sys.stderr)
        return 2

    results = [process_run_dir(path, args.delete, not args.no_activity_check)
               for path in args.run_dirs]
    report(results, args.delete)

    if any(r.refused for r in results):
        return 2
    if any(r.failures for r in results):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
