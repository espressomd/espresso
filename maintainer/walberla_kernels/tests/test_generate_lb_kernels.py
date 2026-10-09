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

import pathlib
import subprocess
import sys
import tempfile
import unittest

SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "generate_lb_kernels.py"


class GenerateLbKernelsTest(unittest.TestCase):

    def run_script(self, *args, check=True):
        tmp = pathlib.Path(self.enterContext(tempfile.TemporaryDirectory()))
        result = subprocess.run(
            [sys.executable, str(SCRIPT), *args, "--output-dir", str(tmp)],
            capture_output=True, text=True, check=False)
        if check:
            self.assertEqual(result.returncode, 0, result.stderr)
        return tmp, result

    def test_init_double_cpu(self):
        tmp, _ = self.run_script("--kernels", "init", "--precision", "double",
                                 "--target", "cpu")
        self.assertEqual(sorted(p.name for p in tmp.iterdir()),
                         ["InitialPDFsSetterDoublePrecision.cpp",
                          "InitialPDFsSetterDoublePrecision.h"])
        header = (tmp / "InitialPDFsSetterDoublePrecision.h").read_text()
        self.assertIn(
            "InitialPDFsSetterDoublePrecision (walberla::BlockDataID forceId, "
            "walberla::BlockDataID pdfsId, walberla::BlockDataID velocityId, "
            "double rho_0)", header)
        self.assertIn("double getRho_0() const", header)

    def test_init_all_variants(self):
        tmp, _ = self.run_script("--kernels", "init")
        self.assertEqual(sorted(p.name for p in tmp.iterdir()), sorted(
            f"InitialPDFsSetter{p}Precision{s}"
            for p in ("Single", "Double")
            for s in (".h", ".cpp", "CUDA.h", "CUDA.cu")))

    def test_unsupported_target_is_skipped(self):
        tmp, result = self.run_script("--kernels", "init", "--target", "avx")
        self.assertEqual(list(tmp.iterdir()), [])
        self.assertIn("skipping init for target avx", result.stdout)

    def test_unknown_kernel(self):
        _, result = self.run_script("--kernels", "bogus", check=False)
        self.assertEqual(result.returncode, 2)
        self.assertIn("invalid choice", result.stderr)

    def test_reproducible_up_to_statement_order(self):
        first, _ = self.run_script("--kernels", "init")
        second, _ = self.run_script("--kernels", "init")
        for path in first.iterdir():
            self.assertEqual(
                sorted(path.read_text().splitlines()),
                sorted((second / path.name).read_text().splitlines()),
                path.name)


if __name__ == "__main__":
    unittest.main()
