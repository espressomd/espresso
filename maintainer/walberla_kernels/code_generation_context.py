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
import sympy
import lbmpy
import lbmpy_walberla
import pystencils
import pystencils_walberla


def add_metadata_codegen_toolchain():
    """
    Add metadata at the beginning of generated kernels to document the
    pystencils/lbmpy/sympy toolchain that was used to create them.
    """
    walberla_root = lbmpy_walberla.__file__.split("/python/lbmpy_walberla/")[0]
    with open(os.path.join(walberla_root, ".git/HEAD")) as f:
        walberla_commit = f.read()
    if walberla_commit.startswith("ref: refs/heads/master"):
        ref = walberla_commit.split()[1]
        with open(os.path.join(walberla_root, f".git/{ref}")) as f:
            walberla_commit = f.read()
    token = "// kernel generated with"
    metadata = (
        f"{token} pystencils v{pystencils.__version__}, "
        f"lbmpy v{lbmpy.__version__}, "
        f"sympy v{sympy.__version__}, "
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
                    f.write(f"\n{metadata}\n{content[pos:].rstrip()}\n")


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


def patch_file(class_name, extension, patch, *args, **kwargs):
    with open(f"{class_name}.{extension}", "r+") as f:
        old_content = f.read()
        new_content = patch(old_content, *args, **kwargs)
        if new_content != old_content:
            f.seek(0)
            f.truncate()
            f.write(new_content)


class CodeGeneration(pystencils_walberla.CodeGeneration):
    """
    This is a patched version of ``CodeGeneration`` that elides parameters
    passed to the command line when running the argument parser, and then
    restores them. It also patches the Jinja templates and adds metadata
    to the generated kernels.
    """

    def __init__(self):
        import sys
        old_sys_argv = sys.argv
        sys.argv = sys.argv[:1]
        super().__init__()
        sys.argv = old_sys_argv
        self.context.patch_file = patch_file

    def __exit__(self, *args, **kwargs):
        super().__exit__(*args, **kwargs)
        add_metadata_codegen_toolchain()
        guard_generated_kernels_clang_format()
