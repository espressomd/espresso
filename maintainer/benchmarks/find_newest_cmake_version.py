#
# Copyright (C) 2013-2026 The ESPResSo project
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
Print the newest CMake version that the module system can load, as used by
``cmake/<version>`` in reframe_config.py.

Exit codes:
    0 (version printed on stdout)
    1 (an error occurred; the reason is printed on stderr)
    2 (no module system found, e.g. on a workstation)
"""

import re
import subprocess
import sys

# "module" is a shell function, not a program, so it only exists in a bash
# that has loaded the module system's init scripts, i.e. a login shell.
NO_MODULE_SYSTEM = 3
MODULE_COMMAND = (
    f"type module >/dev/null 2>&1 || exit {NO_MODULE_SYSTEM}; "
    "module -t avail cmake"
)
# One "cmake/<version>" per line in terse mode; markers such as "(D)" or
# "(default)" follow the version. Builds with suffixes (e.g. "-rc1") are
# skipped, since the version must name the module exactly.
VERSION_RE = re.compile(r"^cmake/(\d+(?:\.\d+)*)(?![\w.-])")


def version_key(version: str) -> tuple[int, ...]:
    return tuple(int(part) for part in version.split("."))


def parse_versions(listing: str) -> list[str]:
    """The cmake versions in the output of ``module -t avail cmake``."""
    return [
        match[1]
        for line in listing.splitlines()
        if (match := VERSION_RE.match(line.strip()))
    ]


def main() -> int:
    # Module systems write their listings to stderr, so merge both streams.
    p = subprocess.run(
        ["bash", "-lc", MODULE_COMMAND],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if p.returncode == NO_MODULE_SYSTEM:
        return 2
    if p.returncode != 0:
        print(f"'module -t avail cmake' failed (exit {p.returncode}):\n"
              f"{p.stdout}", file=sys.stderr)
        return 1

    versions = parse_versions(p.stdout)
    if not versions:
        print(f"No cmake module found in:\n{p.stdout}", file=sys.stderr)
        return 1

    print(max(versions, key=version_key))
    return 0


if __name__ == "__main__":
    sys.exit(main())
