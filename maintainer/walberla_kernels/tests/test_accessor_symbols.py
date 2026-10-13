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

import os

"""
The CPU and the CUDA accessors are linked into the same library when
ESPResSo is built with CUDA. The kernels of both backends have the same
kinds of parameters (pointers, sizes, strides and scalars), so a function
of the CPU header and a function of the CUDA source file must never have
the same qualified name and parameter types: the linker would bind the
``inline`` CPU kernel to the host stub of the CUDA kernel (or the other way
round) without a diagnostic.
"""

import re
import unittest

from tests import header_api
from tests import test_field_accessors
from tests import test_field_accessors_cuda


def definitions(text):
    """
    Return the set of tuples (qualified name, parameter types) of all
    functions that are defined in the namespaces of a C++ file, whether they
    are ``inline``, ``__global__`` or neither. The names of the namespaces
    are normalized (``walberla::lbm`` and nested namespaces are equal), the
    types are canonical and ignore ``RESTRICT``.
    """
    text = header_api.strip_comments_and_strings(text)
    text = re.sub(r"\bRESTRICT\b|__restrict__|__launch_bounds__\([^()]*\)",
                  "", text)
    stack = []  # names of the namespaces, or None for any other block
    result = set()
    start = 0
    for i, ch in enumerate(text):
        if ch == "{":
            head = text[start:i]
            match = re.search(r"\bnamespace\s+([\w:]+)\s*$", head)
            if match:
                stack.append(match.group(1))
            else:
                if stack and all(s is not None for s in stack):
                    function = parse_function(head)
                    if function is not None:
                        name, types = function
                        result.add(("::".join(
                            "::".join(stack).split("::") + [name]), types))
                stack.append(None)
            start = i + 1
        elif ch == "}":
            stack.pop()
            start = i + 1
        elif ch == ";":
            start = i + 1
    return result


def parse_function(head):
    """The name and the parameter types of a function definition, whose
    text ends before the opening brace, or None."""
    head = re.sub(r"\s*\b(const|noexcept)\s*$", "", head.strip())
    if not head.endswith(")"):
        return None
    depth = 0
    for i in range(len(head) - 1, -1, -1):
        depth += {")": 1, "(": -1}.get(head[i], 0)
        if depth == 0:
            break
    name = re.search(r"(\w+)\s*$", head[:i])
    if name is None or name.group(1) in ("if", "for", "while", "switch"):
        return None
    params = header_api.split_top_level(head[i + 1:-1])
    return name.group(1), tuple(
        header_api.parse_parameter(p)[0] for p in params)


class AccessorSymbolsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cpu_dir = test_field_accessors.generate_headers()
        cls.cuda_dir = test_field_accessors_cuda.generate_files()
        cls.cpu, cls.cuda = {}, {}
        for p in ("Double", "Single"):
            cls.cpu[p] = definitions(
                (cls.cpu_dir / f"FieldAccessors{p}Precision.h").read_text())
            cls.cuda[p] = definitions(
                (cls.cuda_dir / f"FieldAccessors{p}PrecisionCUDA.cu"
                 ).read_text())

    @classmethod
    def tearDownClass(cls):
        test_field_accessors.remove_directory(cls.cpu_dir)
        test_field_accessors_cuda.remove_directory(cls.cuda_dir)

    def test_extraction(self):
        text = """
        namespace walberla { namespace lbm::accessor {
        namespace Density { namespace kernels {
        __global__ void __launch_bounds__(256) get(
            double * RESTRICT const _data, int64_t const n) { if (n) {} }
        } // namespace kernels
        void get(double x, std::vector<double> const& v) const { }
        }}}
        """
        self.assertEqual(definitions(text), {
            ("walberla::lbm::accessor::Density::kernels::get",
             ("double*", "int64_t")),
            ("walberla::lbm::accessor::Density::get",
             ("double", "std::vector<double> const&"))})

    def test_functions_are_found(self):
        for precision in ("Double", "Single"):
            kernels = {n for n, _ in self.cpu[precision]
                       if "cpu_kernels::" in n}
            self.assertEqual(len(kernels), 18)
            self.assertGreater(len(self.cuda[precision]), 30)
            self.assertTrue(any(
                n.endswith("::kernels::get")
                for n, _ in self.cuda[precision]))

    def test_cpu_and_cuda_functions_are_disjoint(self):
        for precision in ("Double", "Single"):
            shared = self.cpu[precision] & self.cuda[precision]
            self.assertEqual(
                shared, set(),
                f"{precision} precision: the functions have the same "
                f"qualified name and parameter types on both targets")


if __name__ == "__main__":
    unittest.main()
