#
# Copyright (C) 2018-2020 The ESPResSo project
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
import sys
import time
import numpy as np


def minimize(system, energy_target):
    '''
    Run a minimization loop with the steepest descent algorithm.
    Exit with a non-zero error code if the target energy cannot be reached.

    Parameters
    ----------
    system: :class:`espressomd.system.System`
        System to minimize.
    energy_target: :obj:`float`
        Energy threshold.

    '''
    system.integrator.set_steepest_descent(
        f_max=0,
        gamma=0.001,
        max_displacement=0.01)
    for _ in range(20):
        energy = system.analysis.energy()["total"]
        print(f"Minimization: {energy:.1f}")
        if energy < energy_target:
            break
        system.integrator.run(20)
    else:
        print(f"Minimization failed to converge to {energy_target:.1f}")
        exit(1)


def get_timings(system, n_steps, n_iterations, verbose=True):
    '''
    Time the integration loop and write the state of the system to stdout.

    Parameters
    ----------
    system: :class:`espressomd.system.System`
        System to integrate and propagate.
    n_steps: :obj:`int`
        Number of integration steps per timing.
    n_iterations: :obj:`int`
        Number of timings.
    verbose: :obj:`bool`
        Whether to print the state of the system during timing.

    Returns
    -------
    :obj:`ndarray` of :obj:`float`
        Timings.

    '''
    if verbose:
        print(f"Timing every {n_steps} steps")
    timings = []
    for i in range(n_iterations):
        tick = time.time()
        system.integrator.run(n_steps)
        tock = time.time()
        t = (tock - tick) / n_steps
        timings.append(t)
        if verbose:
            energy = system.analysis.energy()["total"]
            verlet = system.cell_system.get_state()["verlet_reuse"]
            print(
                f"step {i}, time: {t:.2e}, verlet: {verlet:.2f}, energy: {energy:.2e}")
    return np.array(timings)


def get_average_time(timings):
    '''
    Calculate the average and 95% confidence interval of the timings.

    Parameters
    ----------
    timings: :obj:`ndarray` of :obj:`float`
        Timings.

    Returns
    -------
    (2,) array_like of :obj:`float`
        Average and confidence interval.

    '''
    avg = np.average(timings)
    ci = 1.96 * np.std(timings) / np.sqrt(len(timings) - 1)
    return (avg, ci)


def write_report(filepath, n_proc, timings, n_steps, label=''):
    '''
    Append timing data to a CSV file. If it doesn't exist, it is created
    with a header.

    Parameters
    ----------
    filepath: :obj:`str`
        Path to the CSV file.
    n_proc: :obj:`int`
        Number of MPI ranks.
    timings: :obj:`ndarray` of :obj:`float`
        Timings.
    n_steps: :obj:`int`
        Number of integration steps per timing.
    label: :obj:`str`, optional
        Label to distinguish e.g. MD from MC or LB steps.

    '''
    script = os.path.basename(sys.argv[0])
    cmd = " ".join(x for x in sys.argv[1:] if not x.startswith("--output"))
    avg, ci = get_average_time(timings)
    header = '"script","arguments","cores","mean","ci","nsteps","duration","label"\n'
    report = f'"{script}","{cmd}",{n_proc},{avg:.3e},{ci:.3e},{n_steps},{np.sum(timings):.1f},"{label}"\n'
    if os.path.isfile(filepath):
        header = ''
    with open(filepath, "a") as f:
        f.write(header + report)
