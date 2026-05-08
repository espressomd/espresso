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
    'systems': [
        {
            'name': 'local',
            'descr': 'Local system',
            'hostnames': [
                'cip1', 'cip2', 'cip3', 'cip4', 'cip5', 'cip6', 'cip7', 'cip8', 'cip9', 'cip10', 'cip11', 'cip12', 'cip13',
                'cip14', 'cip15', 'cip16', 'cip17', 'cip18', 'cip19', 'cip20', 'cip21', 'cip22', 'cip23', 'cip24', 'kob'
            ],
            'partitions': [
                {
                    'name': 'default',
                    'descr': 'Default partition',
                    'scheduler': 'local',
                    'launcher': 'local',
                    'environs': ['local-env']
                }
            ]
        },
        {
            'name': 'ant_cluster',
            'descr': "ICP university compute cluster",
            'hostnames': ['ant'],
            'modules_system': 'lmod',
            'partitions': [
                {
                    'name': 'debug',
                    'descr': 'Debug compute nodes.',
                    'scheduler': 'squeue',
                    'launcher': 'srun',
                    'time_limit': '0d0h20m0s',
                    'access': ['--partition=debug'],
                    'environs': ['espresso-env'],
                },
            ]
        }
    ],
    'environments': [
        {
            'name': 'local-env',
            'cc': 'gcc',
            'cxx': 'g++',
        },
        {
            'name': 'espresso-env',
            'cc': 'gcc',
            'cxx': 'g++',
            'modules': [
                'spack/default',
                'gcc/12.5.0',
                'cuda/12.3.0',
                'openmpi/4.1.6',
                'fftw/3.3.10',
                'boost/1.84.0',
                'cmake/4.1.2',
                'python/3.12.9',
            ]
        }
    ],
    'logging': [
        {
            'handlers_perflog': [
                {
                    'type': 'filelog',
                    'prefix': '%(check_system)s/%(check_partition)s',
                    'level': 'info',
                    'format': '%(check_result)s,'
                              '%(check_descr)s,'
                              '%(check_job_completion_time)s,'
                              '%(hostname)s,'
                              '%(check_system)s:%(check_partition)s,'
                              '%(check_environ)s,'
                              '%(check_num_tasks)s,'
                              '%(check_num_tasks_per_core)s,'
                              '%(check_num_tasks_per_node)s,'
                              '%(check_perfvalues)s',
                    'format_perfvars': ('%(check_perf_value)s,'
                                        '%(check_perf_unit)s,'),
                    'append': True
                }
            ]
        }
    ]
}
