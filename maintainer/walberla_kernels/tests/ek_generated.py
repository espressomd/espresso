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
The generated sweeps of the electrokinetics for
the tests. Generating the flux sweeps takes about a minute each, so the tests
of several modules share the kernels, which are generated once for each test
process, in parallel and for the CPU and the CUDA target.
"""

import atexit
import concurrent.futures
import functools
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "generate_ek_kernels.py"
SWEEPS = ("advective_flux", "continuity", "diffusive_flux",
          "diffusive_flux_thermalized", "diffusive_flux_electrostatic",
          "diffusive_flux_electrostatic_thermalized", "friction_coupling")
# the keys of the registry that these tests generate, and the class names
# (the accessors are cheap: their tests generate them themselves)
BOUNDARIES = ()
CLASS_NAMES = {
    "advective_flux": "AdvectiveFluxKernel",
    "continuity": "ContinuityKernel",
    "diffusive_flux": "DiffusiveFluxKernel",
    "diffusive_flux_thermalized": "DiffusiveFluxKernelThermalized",
    "diffusive_flux_electrostatic": "DiffusiveFluxKernelWithElectrostatic",
    "diffusive_flux_electrostatic_thermalized":
        "DiffusiveFluxKernelWithElectrostaticThermalized",
    "friction_coupling": "FrictionCouplingKernel",
}
# kernels that are generated in single precision too
SINGLE_PRECISION = ("diffusive_flux_electrostatic_thermalized",
                    "advective_flux")


def run_generator(kernel, precision, directory):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--kernels", kernel,
         "--precision", precision, "--output-dir", str(directory)],
        capture_output=True, text=True, check=False,
        env=dict(os.environ, PYTHONHASHSEED="0"))
    if result.returncode:
        raise RuntimeError(
            f"generation of {kernel} failed: {result.stderr[-2000:]}")


@functools.lru_cache(maxsize=None)
def generated(precision):
    """
    The directory that has the files of all the sweeps in ``precision``
    (single precision: only those of `SINGLE_PRECISION`), for all targets. It is removed when the interpreter
    exits.
    """
    directory = pathlib.Path(tempfile.mkdtemp())
    atexit.register(shutil.rmtree, directory, ignore_errors=True)
    kernels = (*SWEEPS, *BOUNDARIES)
    if precision == "single":
        kernels = tuple(k for k in kernels if k in SINGLE_PRECISION)
    workers = min(len(kernels), os.cpu_count() or 1, 8)
    with concurrent.futures.ThreadPoolExecutor(workers) as pool:
        futures = [pool.submit(run_generator, kernel, precision, directory)
                   for kernel in kernels]
        for future in futures:
            future.result()
    return directory
