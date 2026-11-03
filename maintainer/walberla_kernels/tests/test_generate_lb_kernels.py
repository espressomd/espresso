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
import pathlib
import re
import subprocess
import sys
import tempfile
import unittest

import sympy as sp

SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "generate_lb_kernels.py"


class GenerateLbKernelsTest(unittest.TestCase):

    def assert_plain_sweep(self, header, name):
        """The header holds the sweepgen class ``name`` and nothing else."""
        namespaces = re.findall(
            r"^ *namespace ([\w:]+) \{", header, re.MULTILINE)
        self.assertEqual(namespaces[0], "walberla::pystencils")
        self.assertIn("detail", namespaces[1:])
        self.assertEqual(re.findall(r"^ *class (\w+)", header, re.MULTILINE),
                         [name])
        self.assertNotIn("::internal::", header)
        self.assertIsNone(re.search(r"namespace internal\b", header))
        self.assertIsNone(re.search(r"\bimpl_\b", header))
        for legacy in ("getTime_step", "setKt", "configure("):
            self.assertNotIn(legacy, header)

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
        self.assert_plain_sweep(header, "InitialPDFsSetterDoublePrecision")
        self.assertIn(
            "InitialPDFsSetterDoublePrecision (walberla::BlockDataID forceId, "
            "walberla::BlockDataID pdfsId, walberla::BlockDataID velocityId, "
            "double rho_0)", header)
        self.assertIn("double& rho_0()", header)
        source = (tmp / "InitialPDFsSetterDoublePrecision.cpp").read_text()
        self.assertIn("_data_force", source)
        self.assertIn("rho_0", source)

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

    def test_field_accessors_are_header_only(self):
        tmp, _ = self.run_script("--kernels", "field_accessors",
                                 "--precision", "double", "--target", "cpu")
        self.assertEqual([p.name for p in tmp.iterdir()],
                         ["FieldAccessorsDoublePrecision.h"])
        header = (tmp / "FieldAccessorsDoublePrecision.h").read_text()
        self.assertIn("namespace Population", header)
        self.assertIn("inline auto", header)
        # the kernels are defined in the header
        self.assertIn("namespace cpu_kernels", header)
        self.assertIn("pystencils_runtime/generic_cpu.hpp", header)

    def test_field_accessors_cuda_is_skipped(self):
        tmp, result = self.run_script("--kernels", "field_accessors",
                                      "--target", "cuda")
        self.assertEqual(list(tmp.iterdir()), [])
        self.assertIn("skipping field_accessors for target cuda",
                      result.stdout)

    def test_field_accessors_cuda(self):
        tmp, _ = self.run_script("--kernels", "field_accessors_cuda",
                                 "--precision", "double", "--target", "cuda")
        self.assertEqual(sorted(p.name for p in tmp.iterdir()),
                         ["FieldAccessorsDoublePrecisionCUDA.cu",
                          "FieldAccessorsDoublePrecisionCUDA.cuh"])
        header = (tmp / "FieldAccessorsDoublePrecisionCUDA.cuh").read_text()
        self.assertIn("namespace Interpolation", header)
        source = (tmp / "FieldAccessorsDoublePrecisionCUDA.cu").read_text()
        self.assertIn('#include "FieldAccessorsDoublePrecisionCUDA.cuh"',
                      source)

    def test_field_accessors_cuda_cpu_targets_are_skipped(self):
        tmp, result = self.run_script("--kernels", "field_accessors_cuda",
                                      "--target", "cpu", "avx")
        self.assertEqual(list(tmp.iterdir()), [])
        for target in ("cpu", "avx"):
            self.assertIn(
                f"skipping field_accessors_cuda for target {target}",
                result.stdout)

    def test_field_accessors_share_one_generator(self):
        import generate_lb_kernels
        cpu = generate_lb_kernels.KERNELS["field_accessors"]
        cuda = generate_lb_kernels.KERNELS["field_accessors_cuda"]
        self.assertIs(cpu.generate, cuda.generate)
        self.assertEqual((cpu.stem, cpu.targets, cpu.header_only),
                         ("FieldAccessors", ("cpu",), True))
        self.assertEqual((cuda.stem, cuda.targets, cuda.header_only),
                         ("FieldAccessors", ("cuda",), False))

    def test_legacy_accessor_generators_are_removed(self):
        import importlib.util
        for module in ("field_accessors", "field_accessors_cuda"):
            self.assertIsNone(importlib.util.find_spec(module), module)

    def test_header_extension_of_kernels(self):
        import generate_lb_kernels
        extensions = {key: kernel.header_extension for key, kernel in
                      generate_lb_kernels.KERNELS.items()}
        self.assertEqual(extensions.pop("field_accessors_cuda"), "cuh")
        self.assertEqual(set(extensions.values()), {"h"})

    def test_pack_info_pdf_all_variants(self):
        tmp, result = self.run_script("--kernels", "pack_info_pdf")
        self.assertEqual(sorted(p.name for p in tmp.iterdir()), sorted(
            f"PackInfoPdf{p}Precision{s}"
            for p in ("Single", "Double")
            for s in (".h", ".cpp", "CUDA.h", "CUDA.cu")))
        self.assertIn("skipping pack_info_pdf for target avx", result.stdout)

    def test_pack_info_vec_double_cpu(self):
        tmp, _ = self.run_script("--kernels", "pack_info_vec",
                                 "--precision", "double", "--target", "cpu")
        self.assertEqual(sorted(p.name for p in tmp.iterdir()),
                         ["PackInfoVecDoublePrecision.cpp",
                          "PackInfoVecDoublePrecision.h"])
        header = (tmp / "PackInfoVecDoublePrecision.h").read_text()
        self.assertIn("class PackInfoVecDoublePrecision : public "
                      "walberla::UniformCpuFieldPackInfoBase<", header)

    def test_pack_info_vec_single_cuda(self):
        tmp, _ = self.run_script("--kernels", "pack_info_vec",
                                 "--precision", "single", "--target", "cuda")
        self.assertEqual(sorted(p.name for p in tmp.iterdir()),
                         ["PackInfoVecSinglePrecisionCUDA.cu",
                          "PackInfoVecSinglePrecisionCUDA.h"])
        header = (tmp / "PackInfoVecSinglePrecisionCUDA.h").read_text()
        self.assertIn("walberla::sweepgen::UniformGpuFieldPackInfoBase<",
                      header)

    def test_pack_infos_are_registered_for_cpu_and_cuda(self):
        import generate_lb_kernels
        for key, stem in (("pack_info_pdf", "PackInfoPdf"),
                          ("pack_info_vec", "PackInfoVec")):
            kernel = generate_lb_kernels.KERNELS[key]
            self.assertEqual(
                (kernel.stem, kernel.targets, kernel.header_only,
                 kernel.header_extension),
                (stem, ("cpu", "cuda"), False, "h"))

    def test_pack_info_output_is_deterministic(self):
        # independent of the hash seed of the interpreter
        outputs = []
        for seed in ("1", "2"):
            tmp = pathlib.Path(self.enterContext(
                tempfile.TemporaryDirectory()))
            subprocess.run(
                [sys.executable, str(SCRIPT), "--kernels", "pack_info_pdf",
                 "pack_info_vec", "--output-dir", str(tmp)],
                env=dict(os.environ, PYTHONHASHSEED=seed),
                capture_output=True, text=True, check=True)
            outputs.append(tmp)
        names = sorted(p.name for p in outputs[0].iterdir())
        self.assertEqual(len(names), 16)
        self.assertEqual(names, sorted(p.name for p in outputs[1].iterdir()))
        for name in names:
            self.assertEqual((outputs[0] / name).read_bytes(),
                             (outputs[1] / name).read_bytes(), name)

    def test_dynamic_ubb_all_variants(self):
        tmp, result = self.run_script("--kernels", "dynamic_ubb")
        self.assertEqual(sorted(p.name for p in tmp.iterdir()), sorted(
            f"DynamicUBB{p}Precision{s}"
            for p in ("Single", "Double")
            for s in (".h", ".cpp", "CUDA.h", "CUDA.cu")))
        self.assertIn("skipping dynamic_ubb for target avx", result.stdout)

    def test_dynamic_ubb_double_cpu(self):
        tmp, _ = self.run_script("--kernels", "dynamic_ubb",
                                 "--precision", "double", "--target", "cpu")
        self.assertEqual(sorted(p.name for p in tmp.iterdir()),
                         ["DynamicUBBDoublePrecision.cpp",
                          "DynamicUBBDoublePrecision.h"])
        header = (tmp / "DynamicUBBDoublePrecision.h").read_text()
        self.assertIn("class DynamicUBBDoublePrecision : public "
                      "walberla::IndexVectorBoundary<", header)
        self.assertIn("walberla::IndexVectorMemory::Host", header)

    def test_dynamic_ubb_single_cuda(self):
        tmp, _ = self.run_script("--kernels", "dynamic_ubb",
                                 "--precision", "single", "--target", "cuda")
        self.assertEqual(sorted(p.name for p in tmp.iterdir()),
                         ["DynamicUBBSinglePrecisionCUDA.cu",
                          "DynamicUBBSinglePrecisionCUDA.h"])
        header = (tmp / "DynamicUBBSinglePrecisionCUDA.h").read_text()
        self.assertIn("walberla::IndexVectorMemory::Device", header)
        self.assertIn("cudaStream_t stream", header)

    def test_dynamic_ubb_is_registered_for_cpu_and_cuda(self):
        import generate_lb_kernels
        kernel = generate_lb_kernels.KERNELS["dynamic_ubb"]
        self.assertEqual(
            (kernel.stem, kernel.targets, kernel.header_only,
             kernel.header_extension),
            ("DynamicUBB", ("cpu", "cuda"), False, "h"))

    def test_dynamic_ubb_output_is_deterministic(self):
        # independent of the hash seed of the interpreter
        outputs = []
        for seed in ("1", "2"):
            tmp = pathlib.Path(self.enterContext(
                tempfile.TemporaryDirectory()))
            subprocess.run(
                [sys.executable, str(SCRIPT), "--kernels", "dynamic_ubb",
                 "--output-dir", str(tmp)],
                env=dict(os.environ, PYTHONHASHSEED=seed),
                capture_output=True, text=True, check=True)
            outputs.append(tmp)
        names = sorted(p.name for p in outputs[0].iterdir())
        self.assertEqual(len(names), 8)
        self.assertEqual(names, sorted(p.name for p in outputs[1].iterdir()))
        for name in names:
            self.assertEqual((outputs[0] / name).read_bytes(),
                             (outputs[1] / name).read_bytes(), name)

    def test_stream_collide_kernels_are_registered(self):
        import generate_lb_kernels
        for key, stem in (
                ("stream_collide_thermalized", "StreamCollideSweepThermalized"),
                ("stream_collide_athermal", "StreamCollideSweepAthermal"),
                ("stream_collide_lees_edwards",
                 "StreamCollideSweepLeesEdwards")):
            kernel = generate_lb_kernels.KERNELS[key]
            self.assertEqual(
                (kernel.stem, kernel.targets, kernel.header_only,
                 kernel.header_extension),
                (stem, ("cpu", "avx", "cuda"), False, "h"))

    def test_sweep_output_is_deterministic(self):
        # independent of the hash seed of the interpreter
        import generate_lb_kernels
        import sfg_driver
        kernels = ("stream_collide_thermalized", "stream_collide_athermal")
        outputs = []
        for seed in ("1", "2"):
            output_dir = pathlib.Path(self.enterContext(
                tempfile.TemporaryDirectory()))
            subprocess.run(
                [sys.executable, str(SCRIPT), "--kernels", *kernels,
                 "--output-dir", str(output_dir)],
                env=dict(os.environ, PYTHONHASHSEED=seed),
                capture_output=True, text=True, check=True)
            outputs.append(output_dir)
        names = sorted(p.name for p in outputs[0].iterdir())
        self.assertEqual(names, sorted(p.name for p in outputs[1].iterdir()))
        # header and source of each precision and target
        self.assertEqual(len(names), sum(
            2 * len(sfg_driver.PRECISIONS)
            * len(generate_lb_kernels.KERNELS[k].targets)
            for k in kernels))
        for name in names:
            self.assertEqual((outputs[0] / name).read_bytes(),
                             (outputs[1] / name).read_bytes(), name)

    def test_duplicate_selections_generate_once(self):
        _, result = self.run_script(
            "--kernels", "init", "init", "--precision", "double", "double",
            "--target", "cpu", "cpu")
        self.assertEqual(
            result.stdout.count("generated InitialPDFsSetterDoublePrecision"),
            1)

    def test_unknown_kernel(self):
        _, result = self.run_script("--kernels", "bogus", check=False)
        self.assertEqual(result.returncode, 2)
        self.assertIn("invalid choice", result.stderr)

    def generate_one(self, kernel, stem, target="cpu", precision="double"):
        tmp, _ = self.run_script("--kernels", kernel, "--precision", precision,
                                 "--target", target)
        suffix = {"cpu": "", "avx": "AVX", "cuda": "CUDA"}[target]
        name = f"{stem}{precision.capitalize()}Precision{suffix}"
        impl = "cu" if target == "cuda" else "cpp"
        return ((tmp / f"{name}.h").read_text(),
                (tmp / f"{name}.{impl}").read_text())

    def test_stream_collide_thermalized(self):
        header, source = self.generate_one(
            "stream_collide_thermalized", "StreamCollideSweepThermalized")
        self.assert_plain_sweep(
            header, "StreamCollideSweepThermalizedDoublePrecision")
        self.assertIn(
            "StreamCollideSweepThermalizedDoublePrecision ("
            "walberla::BlockDataID forceId, walberla::BlockDataID pdfsId, "
            "walberla::Vector3< int64_t > block_offset, double kT, "
            "double omega_bulk, double omega_even, double omega_odd, "
            "double omega_shear, int64_t seed, int64_t time_step)", header)
        self.assertIn("walberla::Vector3< int64_t >& block_offset()", header)
        self.assertIn("double& kT()", header)
        self.assertIn("int64_t& time_step()", header)
        self.assertIn("void swapShadowBuffers(", header)
        self.assertIn("::pystencils::runtime::random::philox_float64", source)
        self.assertNotIn("Athermal", header + source)

    def test_stream_collide_thermalized_cuda(self):
        header, _ = self.generate_one(
            "stream_collide_thermalized", "StreamCollideSweepThermalized",
            "cuda")
        self.assert_plain_sweep(
            header, "StreamCollideSweepThermalizedDoublePrecisionCUDA")
        self.assertIn("walberla::v8::sweep::exectag::GPU exTag", header)
        self.assertNotIn("gpuStream_t stream", header)

    def test_stream_collide_thermalized_avx(self):
        _, source = self.generate_one(
            "stream_collide_thermalized", "StreamCollideSweepThermalized",
            "avx")
        self.assertIn("_mm256_", source)
        self.assertIn("philox_float64<int64_t>(", source)

    def test_stream_collide_thermalized_single_avx(self):
        _, source = self.generate_one(
            "stream_collide_thermalized", "StreamCollideSweepThermalized",
            "avx", "single")
        self.assertIn("int32_t seed", source)
        self.assertIn("_mm256_", source)

    def test_stream_collide_athermal(self):
        header, source = self.generate_one(
            "stream_collide_athermal", "StreamCollideSweepAthermal")
        self.assert_plain_sweep(
            header, "StreamCollideSweepAthermalDoublePrecision")
        # no fluctuations: no temperature, seed, time step or block offset
        self.assertIn(
            "StreamCollideSweepAthermalDoublePrecision ("
            "walberla::BlockDataID forceId, walberla::BlockDataID pdfsId, "
            "double omega_bulk, double omega_even, double omega_odd, "
            "double omega_shear)", header)
        self.assertIn("void swapShadowBuffers(", header)
        for name in ("kT", "seed", "time_step", "block_offset"):
            self.assertNotIn(name, header + source)
        self.assertNotIn("philox", source)

    def test_stream_collide_athermal_avx_and_cuda(self):
        header, source = self.generate_one(
            "stream_collide_athermal", "StreamCollideSweepAthermal", "avx",
            "single")
        self.assert_plain_sweep(
            header, "StreamCollideSweepAthermalSinglePrecisionAVX")
        self.assertIn("_mm256_", source)
        header, source = self.generate_one(
            "stream_collide_athermal", "StreamCollideSweepAthermal", "cuda",
            "single")
        self.assert_plain_sweep(
            header, "StreamCollideSweepAthermalSinglePrecisionCUDA")
        self.assertIn("walberla::v8::sweep::exectag::GPU exTag", header)
        self.assertNotIn("philox", source)

    def test_stream_collide_lees_edwards(self):
        header, source = self.generate_one(
            "stream_collide_lees_edwards", "StreamCollideSweepLeesEdwards")
        self.assertIn(
            "StreamCollideSweepLeesEdwardsDoublePrecision ("
            "walberla::BlockDataID forceId, walberla::BlockDataID pdfsId, "
            "int64_t lebc_bot_index, int64_t lebc_top_index, "
            "double omega_bulk, double omega_even, double omega_odd, "
            "double omega_shear, double v_s)", header)
        self.assert_plain_sweep(
            header, "StreamCollideSweepLeesEdwardsDoublePrecision")
        self.assertIn("double& v_s()", header)
        self.assertIn("le_top", source)

    def test_stream_collide_lees_edwards_avx(self):
        _, source = self.generate_one(
            "stream_collide_lees_edwards", "StreamCollideSweepLeesEdwards",
            "avx")
        self.assertIn("_mm256_", source)
        self.assertIn("std::fmin(", source)

    def test_vel_update(self):
        header, _ = self.generate_one("vel_update", "UpdateVelFromPDF")
        self.assertIn(
            "UpdateVelFromPDFDoublePrecision (walberla::BlockDataID forceId, "
            "walberla::BlockDataID pdfsId, walberla::BlockDataID velocityId)",
            header)
        self.assert_plain_sweep(header, "UpdateVelFromPDFDoublePrecision")
        self.assertIn("void runOnCellInterval(", header)
        self.assertNotIn("swapShadowBuffers", header)

    def test_vel_update_uses_force(self):
        _, source = self.generate_one("vel_update", "UpdateVelFromPDF")
        self.assertIn("_data_force", source)

    def test_reproducible(self):
        first, _ = self.run_script("--kernels", "init")
        second, _ = self.run_script("--kernels", "init")
        for path in first.iterdir():
            self.assertEqual(path.read_bytes(),
                             (second / path.name).read_bytes(), path.name)


class UpdateRuleTest(unittest.TestCase):

    @staticmethod
    def flatten(rule):
        """Main assignments with all subexpressions substituted."""
        subexpressions = {a.lhs: a.rhs for a in rule.subexpressions}
        result = []
        for assignment in rule.main_assignments:
            expr = assignment.rhs
            for _ in range(200):
                new_expr = expr.xreplace(subexpressions)
                if new_expr == expr:
                    break
                expr = new_expr
            result.append(expr)
        return result

    def test_athermal_rule_is_thermalized_rule_at_zero_kT(self):
        import generate_lb_kernels as gen
        import pystencils_espresso
        import sfg_driver
        ctx = sfg_driver.CodegenContext("double", "cpu")
        fields = pystencils_espresso.generate_fields(
            gen.STENCIL, ctx.data_type)
        method = gen.create_lb_method(fields)
        thermal = gen.create_thermalized_update_rule(ctx, fields, method)
        athermal = gen.create_athermal_update_rule(fields, method)
        kT = sp.Symbol("kT")
        thermal_exprs = self.flatten(thermal)
        athermal_exprs = self.flatten(athermal)
        self.assertEqual(len(thermal_exprs), 19)
        self.assertEqual(len(athermal_exprs), 19)
        for a, b in zip(thermal_exprs, athermal_exprs):
            self.assertEqual(sp.simplify(a.subs(kT, 0) - b), 0)


if __name__ == "__main__":
    unittest.main()
