#
# Copyright (C) 2021-2023 The ESPResSo project
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

import os
import re
import hashlib
import lbmpy
import lbmpy_walberla
import pystencils
import pystencils_walberla


def earmark_generated_kernels():
    """
    Add an earmark at the beginning of generated kernels to document the
    pystencils/lbmpy toolchain that was used to create them.
    """
    walberla_root = lbmpy_walberla.__file__.split("/python/lbmpy_walberla/")[0]
    with open(os.path.join(walberla_root, ".git/HEAD")) as f:
        walberla_commit = f.read()
    if walberla_commit.startswith("ref: refs/heads/master"):
        ref = walberla_commit.split()[1]
        with open(os.path.join(walberla_root, f".git/{ref}")) as f:
            walberla_commit = f.read()
    token = "// kernel generated with"
    earmark = (
        f"{token} pystencils v{pystencils.__version__}, "
        f"lbmpy v{lbmpy.__version__}, "
        f"lbmpy_walberla/pystencils_walberla from "
        f"waLBerla commit {walberla_commit}"
    )
    for filename in os.listdir("."):
        if filename.endswith((".h", ".cpp", ".cu", ".cuh")):
            with open(filename, "r+") as f:
                content = f.read()
                if token not in content:
                    pos = 0
                    if content.startswith("/*"):
                        pos = content.find("*/")
                        pos = content.find("\n", pos) + 1
                    elif content.startswith("//====="):
                        pos = content.find("//=====", 5)
                        pos = content.find("\n", pos) + 1
                    f.seek(pos)
                    f.write(f"\n{earmark}\n{content[pos:].rstrip()}\n")


def guard_generated_kernels_clang_format():
    """
    Some namespaces are too long and will break ``clang-format`` versions
    9 and 10. Replace them with a unique string of reasonable size.
    """
    for filename in os.listdir("."):
        if filename.endswith(".cpp"):
            with open(filename, "r") as f:
                content = f.read()
            all_ns = re.findall(r"^namespace (internal_[a-zA-Z0-9_]{54,}) \{$",
                                content, flags=re.MULTILINE)
            if not all_ns:
                continue
            for ns in all_ns:
                ns_hash = hashlib.md5(ns.encode('utf-8')).hexdigest()
                content = re.sub(f"(?<=[^a-zA-Z0-9_]){ns}(?=[^a-zA-Z0-9_])",
                                 f"internal_{ns_hash}", content)
            with open(filename, "w") as f:
                f.write(content)


class CodeGeneration(pystencils_walberla.CodeGeneration):
    """
    This is a patched version of ``CodeGeneration`` that elides parameters
    passed to the command line when running the argument parser, and then
    restores them. It also patches the Jinja templates and earmarks the
    generated kernels.
    """

    def __init__(self):
        import sys
        old_sys_argv = sys.argv
        sys.argv = sys.argv[:1]
        super().__init__()
        sys.argv = old_sys_argv

    def __exit__(self, *args, **kwargs):
        super().__exit__(*args, **kwargs)
        earmark_generated_kernels()
        guard_generated_kernels_clang_format()
