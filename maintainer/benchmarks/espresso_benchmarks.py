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

from pathlib import Path
from benchmark_utils import generate_test_parameters, CONFIGS
import csv
import subprocess
import reframe as rfm
import reframe.utility.sanity as sn
import reframe.utility.typecheck as typ
from reframe.core.builtins import (
    parameter,
    run_after,
    run_before,
    sanity_function,
    variable,
)


# Feature flags declared by the ant_cluster partitions in reframe_config.py.
# A test selects exactly one of them, so a case is never generated for both.
DEBUG_PARTITION_FEATURE = "debug"
COMPUTE_PARTITION_FEATURE = "compute"


def partition_constraints(use_debug):
    """``valid_systems`` entries selecting one ant_cluster partition.

    ``suite.sh --debug`` sets ``use_debug_partition`` on every test (through
    ReFrame's ``-S``), routing the ant_cluster cases to the debug partition;
    otherwise they run on the production compute nodes. The bare ``local``
    entry keeps workstation runs working: the local system declares neither
    feature, so neither ``+debug`` nor ``+compute`` would match it.

    Both classes must agree, because a benchmark depends on its build within
    the same partition and environment.
    """
    feature = (
        DEBUG_PARTITION_FEATURE if use_debug else COMPUTE_PARTITION_FEATURE
    )
    return [f"+{feature}", "local"]


@rfm.simple_test
class BuildEspresso(rfm.CompileOnlyRegressionTest):
    """Compiles ESPResSo for a specific configuration file."""

    # This forces reframe to schedule a build job instead
    # of building espresso on the login node
    build_locally = False
    build_params = parameter(CONFIGS)

    # Overridden per run in select_partition(); the production compute nodes
    # are the default so an unflagged run never lands on the debug queue.
    valid_systems = ["+compute", "local"]
    valid_prog_environs = ["espresso-env", "local-env"]

    # Route the ant_cluster cases to the debug partition instead of the
    # production compute nodes. Set for every test by ``suite.sh --debug``
    # via ReFrame's ``-S use_debug_partition=true``.
    # typ.Bool (not plain bool) so that "-S use_debug_partition=false" parses
    # as False; bool("false") would be True.
    use_debug_partition = variable(typ.Bool, value=False)

    sourcesdir = "https://github.com/espressomd/espresso.git"
    build_system = "CMake"

    # Commit hash of the ESPResSo checkout
    espresso_commit = variable(str, value="unknown", loggable=True)

    @run_after("init")
    def select_partition(self):
        self.valid_systems = partition_constraints(self.use_debug_partition)

    @run_after("init")
    def set_build_attributes(self):
        self.config_name = self.build_params["config"]  # type: ignore
        self.build_system.builddir = "build"  # type: ignore

        self.descr = f"Build Espresso ({self.config_name})"

    @run_after("setup")
    def set_resources(self):
        if not self.is_local():
            self.build_job.num_cpus_per_task = 64  # type: ignore

    def skip_unsupported_local_configs(self):
        if self.is_local():
            supported_configs = ["maxset"]
            if self.config_name not in supported_configs:
                self.skip(
                    f"Local execution only supports {
                        supported_configs} configs "
                    f"(tried to use {self.config_name})"
                )

    @run_before("compile")
    def set_build_instructions(self):
        config_name = self.build_params["config"]  # type: ignore
        config_dir = Path(__file__).parent.parent / "configs"

        CUDA_VER = "12.8"
        GCC_VER = "13"
        CUDAARCHS = r"75;86"

        self.prebuild_cmds = [
            f'cp {config_dir / "empty.hpp"} .',
            f'cp {config_dir / "default.hpp"} .',
            f'cp {config_dir / "maxset.hpp"} .',
            rf'sed -i "1 i\\#define ELECTROSTATICS\\n#define LENNARD_JONES\\n#define MASS\\n" {
                config_name}.hpp',
            rf'sed -ri "/#define\s+ADDITIONAL_CHECKS/d" {config_name}.hpp',
            rf"cp {config_name}.hpp myconfig.hpp",
        ]
        self.build_system.max_concurrency = 16  # type: ignore

        if not self.is_local():
            self.prebuild_cmds += [
                r"python3 -m venv .reframe_venv",
                r"source .reframe_venv/bin/activate",
                r"pip install --upgrade pip",
                rf"pip install -c {(Path(__file__).parents[2] / 'requirements.txt').resolve(
                )} numpy scipy setuptools cython==3.0.8",
            ]

            self.build_system.max_concurrency = 64  # type: ignore

        self.build_system.config_opts = [  # type: ignore
            "..",
            "-D CMAKE_BUILD_TYPE=Release",
            "-D ESPRESSO_BUILD_BENCHMARKS=ON",
            "-D ESPRESSO_TEST_TIMEOUT=1200",
            "-D ESPRESSO_BUILD_WITH_CUDA=ON",
            f"-D ESPRESSO_CMAKE_CUDA_ARCHITECTURES='{CUDAARCHS}'",
            "-D ESPRESSO_BUILD_WITH_WALBERLA=ON",
            "-D ESPRESSO_BUILD_WITH_CCACHE=OFF",
        ]

        # In local runs reframe does not detect default CUDA architecture
        if self.is_local():
            self.build_system.config_opts += [  # type: ignore
                f"-D CMAKE_C_COMPILER=gcc-{GCC_VER}",
                f"-D CMAKE_CUDA_COMPILER=/usr/local/cuda-{CUDA_VER}/bin/nvcc",
                f"-D CUDAToolkit_ROOT=/usr/local/cuda-{CUDA_VER}",
                f"-D CMAKE_CUDA_FLAGS='--compiler-bindir=/usr/bin/g++-{
                    GCC_VER}'",
            ]
            self.skip_unsupported_local_configs()

    @run_after("compile")
    def record_commit_hash(self):
        """
        Record the commit hash of the ESPResSo checkout that was compiled.
        """
        try:
            self.espresso_commit = subprocess.check_output(
                ["git", "-C", self.stagedir, "rev-parse", "HEAD"],
                stderr=subprocess.DEVNULL,
                text=True,
            ).strip()
        except (subprocess.CalledProcessError, OSError):
            self.espresso_commit = "unknown"

    @sanity_function
    def assert_sanity(self):
        return sn.assert_found(r"Built target pypresso", self.stdout)


@rfm.simple_test
class EspressoBenchmark(rfm.RunOnlyRegressionTest):
    """Executes benchmark tests with different parameters."""

    build_params = parameter(CONFIGS)
    test_case = parameter(generate_test_parameters())

    # Overridden per run in select_partition(); must match BuildEspresso, since
    # the dependency below resolves within one partition and environment.
    valid_systems = ["+compute", "local"]
    valid_prog_environs = ["espresso-env", "local-env"]

    # Route the ant_cluster cases to the debug partition instead of the
    # production compute nodes. Set for every test by ``suite.sh --debug``
    # via ReFrame's ``-S use_debug_partition=true``.
    # typ.Bool (not plain bool) so that "-S use_debug_partition=false" parses
    # as False; bool("false") would be True.
    use_debug_partition = variable(typ.Bool, value=False)

    # This will set the slurm option --exlusive for scheduled jobs
    exclusive_access = True

    # Commit hash of the benchmarked ESPResSo checkout
    espresso_commit = variable(str, value="unknown", loggable=True)

    # Name of the ESPResSo build configuration (maxset/default/empty). Different
    # builds share the same commit, so this is what separates them in the log
    # (and lets the plot script emit one SVG per build). Marked loggable so it is
    # written into the perflog as %(check_build_config)s.
    build_config = variable(str, value="unknown", loggable=True)

    @run_after("init")
    def select_partition(self):
        self.valid_systems = partition_constraints(self.use_debug_partition)

    @run_after("init")
    def setup_test(self):
        mpi_enabled = self.build_params["mpi"]  # type: ignore
        self.build_config = self.build_params["config"]  # type: ignore
        self.script_filename, self.script_args, self.num_mpi_ranks = self.test_case  # type: ignore
        self.use_gpu = "--gpu" in self.script_args

        self.variants = BuildEspresso.get_variant_nums(
            build_params=self.build_params)

        assert (
            len(self.variants) == 1
        ), "Benchmark test should depend on exactly one build test."
        self.depends_on(BuildEspresso.variant_name(self.variants[0]))

        if not mpi_enabled:
            self.num_mpi_ranks = 1

        self.num_tasks = self.num_mpi_ranks
        self.num_tasks_per_node = self.num_mpi_ranks
        self.num_cpus_per_task = 1  # For MPI

        if self.use_gpu:
            self.num_gpus_per_node = self.num_mpi_ranks
        else:
            self.num_gpus_per_node = 0

        args_str = "_".join(
            [a.replace("--", "").replace("=", "_") for a in self.script_args]
        )
        self.descr = f"ESPRESSO_{self.script_filename.replace('.py', '')}_{
            args_str}_cores_{self.num_mpi_ranks}"

    @run_before("run")
    def skip_unsupported_test_configs(self):
        tests_to_valid_mpi_ranks_map = {
            "lb.py": {"local": [1], "cluster": [1, 2]}}

        for test_case, valid_mpi_ranks in tests_to_valid_mpi_ranks_map.items():
            if self.script_filename == test_case:
                if (
                    self.is_local()
                    and self.use_gpu
                    and self.num_mpi_ranks not in valid_mpi_ranks["local"]
                ):
                    self.skip(
                        f"Local execution of test case {
                            self.script_filename} with argument '--gpu' only supports"
                        f" {valid_mpi_ranks['local']} cores (tried to use {
                            self.num_mpi_ranks})"
                    )
                elif (
                    self.use_gpu
                    and self.num_mpi_ranks not in valid_mpi_ranks["cluster"]
                ):
                    self.skip(
                        f"Execution of test case {
                            self.script_filename} with argument '--gpu' only supports"
                        f" {valid_mpi_ranks['cluster']} cores (tried to use {
                            self.num_mpi_ranks})"
                    )

    @run_before("run")
    def skip_unsupported_local_configs(self):
        if self.is_local():
            supported_cores = (1, 4)
            if self.num_mpi_ranks not in supported_cores:
                self.skip(
                    f"Local execution only supports {supported_cores} cores "
                    f"(tried to use {self.num_mpi_ranks})"
                )

    @run_before("run")
    def prepare_execution(self):
        build_target = self.getdep(
            BuildEspresso.variant_name(self.variants[0]))
        self.espresso_commit = build_target.espresso_commit
        build_dir = f"{build_target.stagedir}/build"
        script_path = f"{
            build_dir}/maintainer/benchmarks/{self.script_filename}"
        self.benchmark_file_path = f"{self.stagedir}/benchmarks.csv"

        if self.current_system.name == "local":
            self.executable = f"mpiexec -n {
                self.num_mpi_ranks} {build_dir}/pypresso"
        else:
            self.executable = f"{build_dir}/pypresso"

        self.executable_opts = [script_path, *self.script_args]

    @sanity_function
    def check_test_ran_without_errors(self):
        """
        Check if no errors are found in self.stderr and benchmarks.csv
        file was created.
        """

        return sn.and_(
            sn.assert_not_found(r"(?i)error", self.stderr),
            sn.path_isfile(self.benchmark_file_path),
        )

    def _make_perf_extractor(self, label: str, field: str):
        """
        Factory function that creates extraction functions, which read values
        from generated benchamrk file.
        """

        @sn.deferrable
        def _extract():
            with open(self.benchmark_file_path) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if row["label"] == label:
                        return float(row[field])

            raise ValueError(f"Label '{label}' not found in CSV")

        return _extract

    @run_after("run")
    def set_perf_variables(self):

        if self.is_dry_run():
            return

        perf_vars = {}
        csv_path = f"{self.stagedir}/benchmarks.csv"

        with open(csv_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        for row in rows:
            label = row["label"]
            logging_label = label if label == "" else label + "_"

            perf_vars[f"{logging_label}mean"] = sn.make_performance_function(
                self._make_perf_extractor(label, "mean"), "s"
            )

            perf_vars[f"{logging_label}ci"] = sn.make_performance_function(
                self._make_perf_extractor(label, "ci"), "s"
            )

        self.perf_variables = perf_vars
