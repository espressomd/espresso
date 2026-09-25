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

import os
import sys
import time
import pathlib
import numpy as np


def minimize(system, force_target):
    """
    Run a minimization loop with the steepest descent algorithm.
    Exit with a non-zero error code if the target force cannot be reached.

    Parameters
    ----------
    system: :class:`espressomd.system.System`
        System to minimize.
    force_target: :obj:`float`
        Force threshold.

    """
    n_steps = 100
    for i in range(20):
        system.integrator.set_steepest_descent(
            f_max=force_target,
            gamma=0.001 * i,
            max_displacement=0.01)
        steps = system.integrator.run(n_steps)
        print(f"Max force {np.amax(np.abs(system.part.all().f)):.2e}")
        if steps < n_steps:
            break
    if np.amax(np.abs(system.part.all().f)) > force_target:
        print(f"Minimization failed to converge to a force of "
              f"{force_target:.2f}")
        exit(1)


def get_timings(system, n_steps, n_iterations, verbose=True,
                retune_skin_after_steps=None):
    """
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
    retune_skin_after_steps: :obj:`int`, optional
        If provided, retune the skin every this many iterations to within 10%
        of the current skin value.

    Returns
    -------
    :obj:`ndarray` of :obj:`float`
        Timings.

    """
    if verbose:
        print(f"Timing every {n_steps} steps")
    timings = []
    for i in range(n_iterations):
        # Retune skin if requested
        if retune_skin_after_steps is not None and i % retune_skin_after_steps == 0:
            current_skin = system.cell_system.skin
            min_skin_retune = current_skin / 1.1
            max_skin_retune = current_skin * 1.1
            if verbose:
                print(f"Retuning skin at iteration {i} "
                      f"(current: {current_skin:.3f}, "
                      f"range: [{min_skin_retune:.3f}, "
                      f"{max_skin_retune:.3f}])")
            new_skin = system.cell_system.tune_skin(
                min_skin=min_skin_retune, max_skin=max_skin_retune, tol=current_skin * 0.0125, int_steps=n_steps // 4)
            if verbose:
                print(f"New skin: {new_skin:.3f}")

        tick = time.time()
        system.integrator.run(n_steps)
        tock = time.time()
        t = (tock - tick) / n_steps
        timings.append(t)
        if verbose:
            energy = system.analysis.energy()["total"]
            verlet = system.cell_system.get_state()["verlet_reuse"]
            print(
                f"step {i}, time: {1000 * t:.4f} ms, verlet: {verlet:.2f}, energy: {energy:.2e}")
    return np.array(timings)


def get_average_time(timings):
    """
    Calculate the average and 95% confidence interval of the timings.

    Parameters
    ----------
    timings: :obj:`ndarray` of :obj:`float`
        Timings.

    Returns
    -------
    (2,) array_like of :obj:`float`
        Average and confidence interval.

    """
    avg = np.average(timings)
    ci = 1.96 * np.std(timings) / np.sqrt(len(timings) - 1)
    return (avg, ci)


def get_omp_num_threads(system):
    """Number of OpenMP threads per MPI rank (from the script interface)."""
    return int(system.cell_system.get_state()["omp_num_threads"])


def n_cores(system):
    """Total cores = MPI ranks * OpenMP threads."""
    state = system.cell_system.get_state()
    return int(state["n_nodes"]) * int(state["omp_num_threads"])


def add_common_args(parser, default_particles_per_core):
    """
    Register the arguments shared by all benchmark scripts: the optional
    ``tune``/``run`` subcommand, ``--state_file``, ``--skin``, and the
    mutually-exclusive particle-count group. The default particles-per-core
    is stored on the namespace as ``_default_ppc`` for :func:`resolve_n_part`.
    """
    parser.add_argument(
        "mode", nargs="?", choices=["tune", "run"], default=None,
        help="'tune': build and tune, then save --state_file (no timing); "
             "'run': load --state_file and time it (no tuning); "
             "omitted: do both in one go")
    parser.add_argument("--state_file", metavar="PATH", action="store",
                        type=str, default=None, required=False,
                        help="Path to the .npz state file")
    parser.add_argument("--skin", metavar="SKIN", action="store",
                        type=float, default=None, required=False,
                        help="Fix the Verlet skin (disables skin tuning)")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--particles_per_core", metavar="N", action="store",
                       type=int, default=None, required=False,
                       help="Particles per core (core = MPI rank * OMP thread; "
                            f"default: {default_particles_per_core})")
    group.add_argument("--n_particles", metavar="N", action="store",
                       type=int, default=None, required=False,
                       help="Total number of particles, independent of the "
                            "number of cores")
    parser.set_defaults(_default_ppc=default_particles_per_core)


def validate_mode(args):
    """
    Validate and normalize the mode arguments: exit with an error if
    ``tune``/``run`` was requested without a state file, and append the
    ``.npz`` suffix to ``--state_file`` when it is missing.
    """
    if args.mode in ("tune", "run") and not args.state_file:
        raise SystemExit(f"error: '{args.mode}' mode requires --state_file")
    state_file = getattr(args, "state_file", None)
    if state_file and not state_file.endswith(".npz"):
        args.state_file += ".npz"


def resolve_n_part(system, args):
    """
    Total particle count from ``--n_particles`` (fixed) or
    ``--particles_per_core`` (times the number of cores). Falls back to the
    per-script default particles-per-core when neither is given.
    """
    if getattr(args, "n_particles", None) is not None:
        return int(args.n_particles)
    ppc = args.particles_per_core
    if ppc is None:
        ppc = args._default_ppc
    return int(ppc) * n_cores(system)


def topology_meta(system):
    """Parallel topology to embed in a state file."""
    state = system.cell_system.get_state()
    return {
        "n_nodes": int(state["n_nodes"]),
        "node_grid": [int(x) for x in state["node_grid"]],
        "omp_num_threads": int(state["omp_num_threads"]),
    }


def verify_topology(system, meta):
    """
    Refuse to run a state file under a different parallel topology than it was
    tuned with, then restore the exact MPI node grid.
    """
    state = system.cell_system.get_state()
    current_threads = int(state["omp_num_threads"])
    current_nodes = int(state["n_nodes"])
    if current_threads != int(meta["omp_num_threads"]):
        raise SystemExit(
            f"error: state was tuned with {meta['omp_num_threads']} OpenMP "
            f"thread(s), this run has {current_threads}")
    if current_nodes != int(meta["n_nodes"]):
        raise SystemExit(
            f"error: state was tuned with {meta['n_nodes']} MPI rank(s), "
            f"this run has {current_nodes}")
    system.cell_system.node_grid = [int(x) for x in meta["node_grid"]]


def save_state(path, meta, **arrays):
    """
    Save benchmark state to a single ``.npz`` archive: ``meta`` (a dict of
    scalars/short lists) plus named numpy arrays (positions, velocities, ...).
    """
    if not path.endswith(".npz"):
        path = path + ".npz"
    np.savez(path, meta=np.array(meta, dtype=object), **arrays)


def load_state(path):
    """Load a state file. Returns ``(meta_dict, npz_handle)``."""
    if not os.path.exists(path) and not path.endswith(".npz"):
        path = path + ".npz"
    handle = np.load(path, allow_pickle=True)
    meta = handle["meta"].item()
    return meta, handle


def tune_skin_unless_fixed(system, args, min_skin, max_skin, **tune_kwargs):
    """
    Tune the skin, unless the user fixed it via ``--skin`` (in which case set
    that value and skip tuning). Returns the resulting skin.
    """
    if args.skin is not None:
        system.cell_system.skin = args.skin
        return args.skin
    return system.cell_system.tune_skin(
        min_skin=min_skin, max_skin=max_skin, **tune_kwargs)


def write_report(filepath, n_ranks, timings,
                 n_steps, label="", n_threads=None):
    """
    Append timing data to a CSV file. If it doesn't exist, it is created
    with a header.

    Parameters
    ----------
    filepath: :obj:`str`
        Path to the CSV file.
    n_ranks: :obj:`int`
        Number of MPI ranks.
    timings: :obj:`ndarray` of :obj:`float`
        Timings.
    n_steps: :obj:`int`
        Number of integration steps per timing.
    label: :obj:`str`, optional
        Label to distinguish e.g. MD from MC or LB steps.

    """
    if n_threads is None:
        n_threads = int(os.environ.get("OMP_NUM_THREADS", 1))
    script = pathlib.Path(sys.argv[0]).name
    cmd = " ".join(x for x in sys.argv[1:] if not x.startswith("--output"))
    avg, ci = get_average_time(timings)
    header = '"script","arguments","ranks","threads","mean","ci","nsteps","duration","label"\n'
    report = f'"{script}","{cmd}",{n_ranks},{n_threads},{avg:.3e},{ci:.3e},{n_steps},{np.sum(timings):.1f},"{label}"\n'  # nopep8
    if pathlib.Path(filepath).is_file():
        header = ''
    with open(filepath, "a") as f:
        f.write(header + report)
