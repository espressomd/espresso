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

site_configuration = {
    "systems": [
        {
            "name": "ant_cluster",
            "descr": "ICP university compute cluster",
            "hostnames": ["ant"],
            "modules_system": "lmod",
            "partitions": [
                {
                    "name": "debug",
                    "descr": "Debug compute nodes.",
                    "scheduler": "squeue",
                    "launcher": "srun",
                    "time_limit": "0d0h20m0s",
                    "access": ["--partition=debug"],
                    "environs": ["espresso-env"],
                    "features": ["debug"],
                    "resources": [
                        {
                            "name": "_rfm_gpu",
                            "options": ["--gres=gpu:{num_gpus_per_node}"],
                        },
                    ],
                },
                {
                    "name": "compute",
                    "descr": "Production compute nodes with GPU access.",
                    "scheduler": "squeue",
                    "launcher": "srun",
                    "time_limit": "0d2h0m0s",
                    "access": ["--partition=compute"],
                    "environs": ["espresso-env"],
                    "features": ["compute"],
                    "resources": [
                        {
                            "name": "_rfm_gpu",
                            "options": ["--gres=gpu:{num_gpus_per_node}"],
                        },
                    ],
                },
            ],
        },
        {
            "name": "local",
            "descr": "Local system",
            "hostnames": [".*"],
            "partitions": [
                {
                    "name": "default",
                    "descr": "Default partition",
                    "scheduler": "local",
                    "launcher": "local",
                    "environs": ["local-env"],
                    "max_jobs": 1,
                }
            ],
        },
    ],
    "environments": [
        {
            "name": "local-env",
            "cc": "gcc",
            "cxx": "g++",
        },
        {
            "name": "espresso-env",
            "cc": "gcc",
            "cxx": "g++",
            "modules": [
                "spack/default",
                "gcc/12.5.0",
                "cuda/12.3.0",
                "openmpi/4.1.6",
                "fftw/3.3.10",
                "boost/1.84.0",
                "cmake/4.1.2",
                "python/3.12.9",
            ],
        },
    ],
    "logging": [
        {
            "handlers": [
                {
                    "type": "file",
                    "name": "${RFM_PREFIX}/reframe.log",
                    "level": "info",
                    "format": "[%(asctime)s] %(levelname)s: "
                    "%(check_info)s: %(message)s",
                    "append": False,
                }
            ],
            "handlers_perflog": [
                {
                    "type": "filelog",
                    "prefix": "%(check_system)s/%(check_partition)s",
                    "level": "info",
                    "format": "%(check_result)s,"
                    "%(check_descr)s,"
                    "%(check_job_completion_time)s,"
                    "%(hostname)s,"
                    "%(check_system)s:%(check_partition)s,"
                    "%(check_environ)s,"
                    "%(check_num_tasks)s,"
                    "%(check_num_tasks_per_core)s,"
                    "%(check_num_tasks_per_node)s,"
                    "%(check_perf_value)s,"
                    "%(check_perf_unit)s,"
                    "%(check_perf_var)s,"
                    "%(check_espresso_commit)s,"
                    "%(check_build_config)s",
                    "append": True,
                }
            ],
            "perflog_multiline": True,
        }
    ],
}
