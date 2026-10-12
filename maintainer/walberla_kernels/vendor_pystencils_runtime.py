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
Copy the pystencils runtime headers that the committed kernels need into
the waLBerla bridge. Only headers reachable through ``#include "..."``
directives from the kernels are copied; preprocessor conditions are
ignored, so headers behind platform guards are copied too.
"""

import argparse
import pathlib
import re
import shutil

import pystencils.include

REPO = pathlib.Path(__file__).resolve().parents[2]
BRIDGE = REPO / "src" / "walberla_bridge" / "src"
RUNTIME = "pystencils_runtime"
KERNEL_SUFFIXES = (".h", ".hpp", ".cpp", ".cu", ".cuh")
INCLUDE_RE = re.compile(r'^\s*#\s*include\s*"([^"]+)"', re.MULTILINE)


def quoted_includes(path):
    return INCLUDE_RE.findall(path.read_text())


def runtime_includes(kernel_dirs):
    """Runtime headers included by the kernels, e.g. ``pystencils_runtime/cuda.cuh``."""
    return {
        header
        for kernel_dir in kernel_dirs
        for path in kernel_dir.iterdir()
        if path.suffix in KERNEL_SUFFIXES
        for header in quoted_includes(path)
        if header.startswith(f"{RUNTIME}/")
    }


def needed_headers(include_dir, roots):
    """
    Sorted paths, relative to ``include_dir``, of the headers in ``roots``
    and of the runtime headers they include transitively. Includes that
    don't resolve to a runtime header (system or OpenCL headers) are skipped.
    """
    for root in roots:
        if not (include_dir / root).is_file():
            raise FileNotFoundError(f"{root} not found in {include_dir}")
    runtime_dir = (include_dir / RUNTIME).resolve()
    pending = [(include_dir / root).resolve() for root in roots]
    found = set()
    while pending:
        header = pending.pop()
        if header in found:
            continue
        found.add(header)
        for name in quoted_includes(header):
            candidate = (header.parent / name).resolve()
            if candidate.is_file() and candidate.is_relative_to(runtime_dir):
                pending.append(candidate)
    return sorted(str(h.relative_to(include_dir.resolve())) for h in found)


def vendor(include_dir, kernel_dirs, dest):
    """Replace ``dest/pystencils_runtime`` by the headers the kernels need."""
    headers = needed_headers(include_dir, runtime_includes(kernel_dirs))
    shutil.rmtree(dest / RUNTIME, ignore_errors=True)
    for header in headers:
        (dest / header).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(include_dir / header, dest / header)
    return headers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dest", type=pathlib.Path,
                        default=BRIDGE / "third_party",
                        help="Directory that receives pystencils_runtime/")
    args = parser.parse_args()
    kernel_dirs = sorted(BRIDGE.glob("**/generated_kernels"))
    include_dir = pathlib.Path(
        pystencils.include.get_pystencils_include_path())
    for header in vendor(include_dir, kernel_dirs, args.dest):
        print(f"copied {header}")


if __name__ == "__main__":
    main()
