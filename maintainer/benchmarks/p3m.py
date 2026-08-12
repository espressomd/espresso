#
# Copyright (C) 2013-2026 The ESPResSo project
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

import sys
import espressomd
import espressomd.electrostatics
import benchmarks
import numpy as np
import argparse

DEFAULT_PARTICLES_PER_CORE = 1000
INITIAL_SKIN = 0.5
ACCURACY = 1e-3

# Simulation parameters. The importlib-based benchmark tests (in
# testsuite/benchmarks) override these module-level globals to run a fast
# smoke test. ``measurement_steps = None`` auto-computes the step count;
# ``p3m_params = None`` tunes P3M, otherwise the given pre-tuned parameters
# are used verbatim (skipping P3M tuning).
measurement_steps = None
n_iterations = 30
min_skin = 0.2
max_skin = 1.6
p3m_params = None
KT = 1.0
GAMMA = 1.0
SEED = 42
SPECIES = ["anion", "cation"]
CHARGES = {"anion": -1.0, "cation": 1.0}
TYPES = {"anion": 0, "cation": 0}
LJ_SIGMAS = {"anion": 1.0, "cation": 1.0}
LJ_EPSILONS = {"anion": 1.0, "cation": 1.0}
WCA_CUT = 2.**(1. / 6.)
LJ_CUTS = {"anion": WCA_CUT * LJ_SIGMAS["anion"],
           "cation": WCA_CUT * LJ_SIGMAS["cation"]}

parser = argparse.ArgumentParser(description="Benchmark P3M simulations. "
                                 "Save the results to a CSV file.")
benchmarks.add_common_args(parser, DEFAULT_PARTICLES_PER_CORE)
parser.add_argument("--volume_fraction", metavar="FRAC", action="store",
                    type=float, default=0.25, required=False,
                    help="Fraction of the simulation box volume occupied by "
                    "particles (range: [0.01-0.74], default: 0.25)")
parser.add_argument("--prefactor", metavar="PREFACTOR", action="store",
                    type=float, default=1., required=False,
                    help="P3M prefactor (default: 1)")
parser.add_argument("--gpu", action=argparse.BooleanOptionalAction,
                    default=False, required=False, help="Use GPU implementation")
mesh_group = parser.add_mutually_exclusive_group()
mesh_group.add_argument("--mesh", metavar="M", action="store", type=int,
                        default=None, required=False,
                        help="Fixed P3M mesh size (disables mesh tuning)")
mesh_group.add_argument("--lowest_mesh", metavar="M", action="store", type=int,
                        default=None, required=False,
                        help="Lower limit for the mesh size during tuning")
parser.add_argument("--highest_mesh", metavar="M", action="store", type=int,
                    default=None, required=False,
                    help="Upper limit for the mesh size during tuning")
group = parser.add_mutually_exclusive_group()
group.add_argument("--output", metavar="FILEPATH", action="store",
                   type=str, required=False, default="benchmarks.csv",
                   help="Output file (default: benchmarks.csv)")
group.add_argument("--visualizer", action="store_true",
                   help="Starts the visualizer (for debugging purposes)")

args = parser.parse_args()
benchmarks.validate_mode(args)
if args.mesh is not None and args.highest_mesh is not None:
    parser.error("--mesh cannot be combined with --highest_mesh")

required_features = ["P3M", "LENNARD_JONES"]
if args.gpu:
    required_features.append("CUDA")
espressomd.assert_features(required_features)
np.random.seed(SEED)

system = espressomd.System(box_l=[1, 1, 1])
system.time_step = 0.01
system.cell_system.set_regular_decomposition(use_verlet_lists=True)


def configure_lj(system):
    for i in range(len(SPECIES)):
        ion1 = SPECIES[i]
        for j in range(i, len(SPECIES)):
            ion2 = SPECIES[j]
            lj_sig = (LJ_SIGMAS[ion1] + LJ_SIGMAS[ion2]) / 2
            lj_cut = (LJ_CUTS[ion1] + LJ_CUTS[ion2]) / 2
            lj_eps = (LJ_EPSILONS[ion1] * LJ_EPSILONS[ion2])**(1. / 2.)
            system.non_bonded_inter[TYPES[ion1], TYPES[ion2]].lennard_jones.set_params(
                epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")


def p3m_tune_kwargs(args, mesh_min):
    kwargs = {"prefactor": args.prefactor, "accuracy": ACCURACY,
              "timings": 15, "gpu": args.gpu}
    if args.mesh is not None:
        kwargs["mesh"] = args.mesh
    else:
        low = args.lowest_mesh if args.lowest_mesh is not None else mesh_min
        high = args.highest_mesh
        kwargs["tune_limits"] = [low, high]
    return kwargs


def make_p3m(args, mesh_min):
    """
    Build the P3M solver. Use the pre-tuned ``p3m_params`` when provided
    (skipping tuning, e.g. in the benchmark smoke tests), otherwise tune.
    """
    if p3m_params is not None:
        return espressomd.electrostatics.P3M(**p3m_params)
    return espressomd.electrostatics.P3M(**p3m_tune_kwargs(args, mesh_min))


def build_and_tune(system, args):
    assert args.prefactor > 0, "prefactor must be a positive number"
    assert args.volume_fraction > 0, "volume_fraction must be a positive number"
    assert args.volume_fraction < np.pi / (3 * np.sqrt(2)), \
        "volume_fraction exceeds the physical limit of sphere packing (~0.74)"

    n_part = benchmarks.resolve_n_part(system, args)
    steps = measurement_steps
    if steps is None:
        steps = max(10, int(np.round(
            5e5 / n_part * system.cell_system.get_state()["n_nodes"], -1)))
    if not args.visualizer:
        assert steps >= 10, \
            f"{steps} steps per tick are too short"

    lj_sig = (LJ_SIGMAS["cation"] + LJ_SIGMAS["anion"]) / 2
    box_l = (n_part * 4. / 3. * np.pi * (lj_sig / 2.)**3
             / args.volume_fraction)**(1. / 3.)

    # Lower limit for P3M mesh tuning.
    # A mesh that is too coarse forces a large real-space cutoff, which builds a
    # huge Verlet list: this runs out of memory at large N and makes tuning slow.
    # Anchor the limit at mesh 16 for the reference case (1000 particles at
    # volume_fraction 0.25) and scale it with the box length so the 1D mesh
    # density (mesh points per unit length) stays roughly constant across
    # particle numbers and volume fractions.
    ref_box_l = (1000 * 4. / 3. * np.pi * (lj_sig / 2.)**3 / 0.25)**(1. / 3.)
    mesh_min = max(4, int(np.round(16 * (box_l / ref_box_l)**.5)))
    if mesh_min % 2 == 1:
        mesh_min += 1

    system.box_l = 3 * (box_l,)
    system.cell_system.skin = INITIAL_SKIN if args.skin is None else args.skin
    configure_lj(system)

    pid = 0
    for _ in range(0, n_part, len(SPECIES)):
        for t in SPECIES:
            system.part.add(pos=np.random.random(3) * system.box_l,
                            id=pid, q=CHARGES[t], type=TYPES[t])
            pid += 1

    benchmarks.minimize(system, 1 / system.time_step)
    system.integrator.set_vv()
    system.thermostat.set_langevin(kT=KT, gamma=GAMMA, seed=SEED)

    p3m = make_p3m(args, mesh_min)
    print("Quick equilibration")
    system.time_step /= 10.
    system.integrator.run(100)
    system.time_step *= 10.
    system.integrator.run(min(3 * steps, 1000))
    print("Tune skin: {:.3f}".format(benchmarks.tune_skin_unless_fixed(
        system, args, min_skin, max_skin, tol=0.05, int_steps=100,
        adjust_max_skin=True)))
    print("Equilibration")
    system.integrator.run(min(3 * steps, 3000))
    print("Tune p3m")
    system.electrostatics.solver = p3m
    print("Equilibration")
    system.integrator.run(min(3 * steps, 3000))
    print("Tune skin: {:.3f}".format(benchmarks.tune_skin_unless_fixed(
        system, args, min_skin, max_skin, tol=0.05, int_steps=100,
        adjust_max_skin=True)))
    print("Re-tune p3m")
    p3m = make_p3m(args, mesh_min)
    system.electrostatics.solver = p3m
    print("Tune skin: {:.3f}".format(benchmarks.tune_skin_unless_fixed(
        system, args, min_skin, max_skin, tol=0.05, int_steps=100,
        adjust_max_skin=True)))

    return {"n_part": n_part, "measurement_steps": steps, "p3m": p3m}


def save_p3m_state(system, args, ctx):
    p = system.part.all()
    tuned = ctx["p3m"].get_params()
    meta = {
        "skin": float(system.cell_system.skin),
        "box_l": [float(x) for x in system.box_l],
        "n_part": int(ctx["n_part"]),
        "time_step": float(system.time_step),
        "measurement_steps": int(ctx["measurement_steps"]),
        "n_iterations": n_iterations,
        "volume_fraction": float(args.volume_fraction),
        "prefactor": float(args.prefactor),
        "accuracy": float(ACCURACY),
        "mesh": [int(x) for x in tuned["mesh"]],
        "cao": int(tuned["cao"]),
        "alpha": float(tuned["alpha"]),
        "r_cut": float(tuned["r_cut"]),
        "gpu": bool(args.gpu),
        "kT": KT, "gamma": GAMMA, "seed": SEED,
    }
    meta.update(benchmarks.topology_meta(system))
    benchmarks.save_state(
        args.state_file, meta,
        pos=np.copy(p.pos), vel=np.copy(p.v),
        charge=np.copy(p.q), type=np.copy(p.type))
    print(f"Saved state to {args.state_file}")


def run_from_state(system, args):
    meta, handle = benchmarks.load_state(args.state_file)
    benchmarks.verify_topology(system, meta)
    system.box_l = meta["box_l"]
    system.time_step = meta["time_step"]
    system.cell_system.skin = args.skin if args.skin is not None else meta["skin"]
    configure_lj(system)
    system.part.add(pos=handle["pos"], v=handle["vel"],
                    q=handle["charge"], type=handle["type"].astype(int))
    system.integrator.set_vv()
    system.thermostat.set_langevin(
        kT=meta["kT"], gamma=meta["gamma"], seed=int(meta["seed"]))
    p3m = espressomd.electrostatics.P3M(
        prefactor=meta["prefactor"], accuracy=meta["accuracy"],
        mesh=[int(x) for x in meta["mesh"]], cao=int(meta["cao"]),
        alpha=float(meta["alpha"]), r_cut=float(meta["r_cut"]),
        tune=False, gpu=bool(meta["gpu"]), timings=15)
    system.electrostatics.solver = p3m
    return meta["measurement_steps"], meta["n_iterations"]


def time_and_report(system, args, measurement_steps, n_iterations):
    timings = benchmarks.get_timings(system, measurement_steps, n_iterations)
    avg, ci = benchmarks.get_average_time(timings)
    print(f"average: {avg:.3e} +/- {ci:.3e} (95% C.I.)")
    n_proc = system.cell_system.get_state()["n_nodes"]
    benchmarks.write_report(args.output, n_proc, timings, measurement_steps,
                            n_threads=benchmarks.get_omp_num_threads(system))


if args.mode == "run":
    measurement_steps, n_iterations = run_from_state(system, args)
    time_and_report(system, args, measurement_steps, n_iterations)
    sys.exit(0)

ctx = build_and_tune(system, args)

if args.visualizer:
    import espressomd.visualization
    visualizer = espressomd.visualization.openGLLive(system)
    visualizer.run(1)

if args.state_file:
    save_p3m_state(system, args, ctx)

if args.mode == "tune":
    sys.exit(0)

time_and_report(system, args, ctx["measurement_steps"], n_iterations)
