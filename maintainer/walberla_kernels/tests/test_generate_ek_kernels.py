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

import os
import pathlib
import subprocess
import sys
import tempfile
import unittest

import sfg_driver

SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "generate_ek_kernels.py"


class GenerateEkKernelsTest(unittest.TestCase):

    def run_script(self, *args, seed=None):
        tmp = pathlib.Path(self.enterContext(tempfile.TemporaryDirectory()))
        env = dict(os.environ)
        if seed is not None:
            env["PYTHONHASHSEED"] = seed
        result = subprocess.run(
            [sys.executable, str(SCRIPT), *args, "--output-dir", str(tmp)],
            capture_output=True, text=True, check=False, env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        return tmp, result

    def test_dirichlet_all_variants(self):
        tmp, result = self.run_script("--kernels", "dirichlet")
        self.assertEqual(sorted(p.name for p in tmp.iterdir()), sorted(
            f"Dirichlet{p}Precision{s}"
            for p in ("Single", "Double")
            for s in (".h", ".cpp", "CUDA.h", "CUDA.cu")))
        self.assertIn("skipping dirichlet for target avx", result.stdout)

    def test_all_is_the_default(self):
        # (generating all the kernels takes minutes)
        import generate_ek_kernels
        args = sfg_driver.parse_args(
            "test", generate_ek_kernels.KERNELS, [])
        self.assertEqual((args.kernels, args.precision, args.target),
                         (["all"], ["all"], ["all"]))

    def test_registered_kernels(self):
        import generate_ek_kernels
        reactions = [f"reaction_{kind}_{num}" for kind in ("bulk", "indexed")
                     for num in range(1, 6)]
        self.assertEqual(sorted(generate_ek_kernels.KERNELS), sorted([
            "advective_flux", "continuity", "diffusive_flux",
            "diffusive_flux_electrostatic",
            "diffusive_flux_electrostatic_thermalized",
            "diffusive_flux_thermalized", "dirichlet", "field_accessors",
            "field_accessors_cuda", "fixed_flux", "friction_coupling",
            *reactions]))
        for kernel in generate_ek_kernels.KERNELS.values():
            self.assertIsInstance(kernel, sfg_driver.Kernel)
            self.assertNotIn("avx", kernel.targets)

    def test_continuity_all_variants(self):
        tmp, result = self.run_script("--kernels", "continuity")
        self.assertEqual(sorted(p.name for p in tmp.iterdir()), sorted(
            f"ContinuityKernel{p}Precision{s}"
            for p in ("Single", "Double")
            for s in (".h", ".cpp", "CUDA.h", "CUDA.cu")))
        self.assertIn("skipping continuity for target avx", result.stdout)

    def test_fixed_flux_all_variants(self):
        tmp, _ = self.run_script("--kernels", "fixed_flux")
        self.assertEqual(sorted(p.name for p in tmp.iterdir()), sorted(
            f"FixedFlux{p}Precision{s}"
            for p in ("Single", "Double")
            for s in (".h", ".cpp", "CUDA.h", "CUDA.cu")))
        header = (tmp / "FixedFluxDoublePrecision.h").read_text()
        self.assertIn("class FixedFluxDoublePrecision : public "
                      "walberla::IndexVectorBoundary<", header)
        self.assertIn("walberla::DomainLinksNoCenter >", header)

    def test_field_accessors(self):
        tmp, result = self.run_script(
            "--kernels", "field_accessors", "field_accessors_cuda",
            "--precision", "double")
        self.assertEqual(sorted(p.name for p in tmp.iterdir()), [
            "EKFieldAccessorsDoublePrecision.h",
            "EKFieldAccessorsDoublePrecisionCUDA.cu",
            "EKFieldAccessorsDoublePrecisionCUDA.cuh"])
        self.assertIn("skipping field_accessors for target cuda",
                      result.stdout)
        self.assertIn("skipping field_accessors_cuda for target cpu",
                      result.stdout)

    def test_dirichlet_double_cpu(self):
        tmp, _ = self.run_script("--kernels", "dirichlet",
                                 "--precision", "double", "--target", "cpu")
        self.assertEqual(sorted(p.name for p in tmp.iterdir()),
                         ["DirichletDoublePrecision.cpp",
                          "DirichletDoublePrecision.h"])
        header = (tmp / "DirichletDoublePrecision.h").read_text()
        self.assertIn("Generated by "
                      "maintainer/walberla_kernels/generate_ek_kernels.py",
                      header)
        self.assertIn("class DirichletDoublePrecision : public "
                      "walberla::IndexVectorBoundary<", header)
        self.assertIn("walberla::IndexVectorMemory::Host", header)

    def test_dirichlet_single_cuda(self):
        tmp, _ = self.run_script("--kernels", "dirichlet",
                                 "--precision", "single", "--target", "cuda")
        self.assertEqual(sorted(p.name for p in tmp.iterdir()),
                         ["DirichletSinglePrecisionCUDA.cu",
                          "DirichletSinglePrecisionCUDA.h"])
        header = (tmp / "DirichletSinglePrecisionCUDA.h").read_text()
        self.assertIn("walberla::IndexVectorMemory::Device", header)
        self.assertIn("cudaStream_t stream", header)

    def test_dirichlet_is_registered_for_cpu_and_cuda(self):
        import generate_ek_kernels
        kernel = generate_ek_kernels.KERNELS["dirichlet"]
        self.assertEqual(
            (kernel.stem, kernel.targets, kernel.header_only,
             kernel.header_extension),
            ("Dirichlet", ("cpu", "cuda"), False, "h"))
        self.assertIsInstance(kernel, sfg_driver.Kernel)

    def test_reaction_kernels(self):
        tmp, result = self.run_script(
            "--kernels", "reaction_bulk_2", "reaction_indexed_5",
            "--precision", "single")
        self.assertEqual(sorted(p.name for p in tmp.iterdir()), sorted(
            f"ReactionKernel{kind}SinglePrecision{s}"
            for kind in ("Bulk2", "Indexed5")
            for s in (".h", ".cpp", "CUDA.h", "CUDA.cu")))
        self.assertIn("skipping reaction_bulk_2 for target avx", result.stdout)

    def test_number_of_reactants_is_limited(self):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--kernels", "reaction_bulk_6"],
            capture_output=True, text=True, check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("invalid choice: 'reaction_bulk_6'", result.stderr)

    def test_output_is_deterministic(self):
        # independent of the hash seed of the interpreter
        outputs = [self.run_script("--kernels", "dirichlet", seed=seed)[0]
                   for seed in ("1", "2")]
        names = sorted(p.name for p in outputs[0].iterdir())
        self.assertEqual(len(names), 8)
        self.assertEqual(names, sorted(p.name for p in outputs[1].iterdir()))
        for name in names:
            self.assertEqual((outputs[0] / name).read_bytes(),
                             (outputs[1] / name).read_bytes(), name)

    def test_fixed_flux_and_accessors_output_is_deterministic(self):
        # independent of the hash seed of the interpreter
        outputs = [self.run_script(
            "--kernels", "fixed_flux", "field_accessors",
            "field_accessors_cuda", seed=seed)[0] for seed in ("1", "2")]
        names = sorted(p.name for p in outputs[0].iterdir())
        self.assertEqual(len(names), 8 + 3 * 2)
        self.assertEqual(names, sorted(p.name for p in outputs[1].iterdir()))
        for name in names:
            first, second = ((output / name).read_text()
                             for output in outputs)
            if name.startswith("FixedFlux") and name.endswith((".cpp", ".cu")):
                # the order of the field extraction statements varies
                self.assertEqual(sorted(first.splitlines()),
                                 sorted(second.splitlines()), name)
            else:
                self.assertEqual(first, second, name)

    def test_reaction_output_is_deterministic(self):
        # The headers are independent of the hash seed. The source files only
        # differ in the order of the statements that extract the fields and
        # properties, which pystencils-sfg 0.3 doesn't keep stable (as for the
        # sweeps of the lattice Boltzmann method).
        reactions = [f"reaction_{kind}_{num}" for kind in ("bulk", "indexed")
                     for num in range(1, 6)]
        outputs = [self.run_script("--kernels", *reactions, seed=seed)[0]
                   for seed in ("1", "2")]
        names = sorted(p.name for p in outputs[0].iterdir())
        self.assertEqual(len(names), 80)
        self.assertEqual(names, sorted(p.name for p in outputs[1].iterdir()))
        for name in names:
            first, second = (
                (output / name).read_text() for output in outputs)
            if name.endswith(".h"):
                self.assertEqual(first, second, name)
            else:
                self.assertEqual(sorted(first.splitlines()),
                                 sorted(second.splitlines()), name)


if __name__ == "__main__":
    unittest.main()
