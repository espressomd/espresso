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

"""
Benchmark Lattice-Boltzmann fluid + Lennard-Jones particles.
"""
import sys
import espressomd
import espressomd.lb
import benchmarks
import numpy as np
import argparse

LJ_EPS = 1.0
LJ_SIG = 1.0
LJ_CUT = LJ_SIG * 2**(1. / 6.)
DEFAULT_PARTICLES_PER_CORE = 125
INITIAL_SKIN = 0.5
KT = 1.0
GAMMA = 1.0
SEED = 42

# Simulation parameters. The importlib-based benchmark tests (in
# testsuite/benchmarks) override these module-level globals to run a fast
# smoke test. ``measurement_steps = None`` means auto-compute.
measurement_steps = None
n_iterations = 30
min_skin = 0.2
max_skin = 1.0

parser = argparse.ArgumentParser(description="Benchmark LB simulations. "
                                 "Save the results to a CSV file.")
benchmarks.add_common_args(parser, DEFAULT_PARTICLES_PER_CORE)
parser.add_argument("--box_l", action="store", nargs="+",
                    type=int, default=argparse.SUPPRESS, required=False,
                    help="Box length (cubic box); requires 0 particles")
parser.add_argument("--lb_sites_per_particle", metavar="N_LB", action="store",
                    type=float, default=28, required=False,
                    help="Number of LB sites per particle")
parser.add_argument("--volume_fraction", metavar="FRAC", action="store",
                    type=float, default=0.03, required=False,
                    help="Fraction of the simulation box volume occupied by "
                    "particles (range: [0.01-0.74], default: 0.03)")
parser.add_argument("--single_precision", action="store_true", required=False,
                    help="Using single-precision floating point accuracy")
parser.add_argument("--gpu", action=argparse.BooleanOptionalAction,
                    default=False, required=False, help="Use GPU implementation")
parser.add_argument("--multi-gpu", action=argparse.BooleanOptionalAction,
                    default=False, required=False, help="Use multi-GPU implementation")
parser.add_argument("--output", metavar="FILEPATH", action="store",
                    type=str, required=False, default="benchmarks.csv",
                    help="Output file (default: benchmarks.csv)")
parser.add_argument("--node_grid", action="store", nargs=3,
                    type=int, default=None, required=False, help="MPI topology")
parser.add_argument("--blocks_per_mpi_rank", action="store", nargs=3,
                    type=int, default=[1, 1, 1], required=False,
                    help="blocks per mpi rank")
parser.add_argument("--weak_scaling", action="store_true", required=False,
                    help="The measurement of weak scaling")

args = parser.parse_args()
benchmarks.validate_mode(args)

required_features = ["WALBERLA"]
if args.gpu:
    required_features.append("CUDA")
espressomd.assert_features(required_features)
np.random.seed(SEED)

system = espressomd.System(box_l=[1, 1, 1])
system.time_step = 0.01


def make_lb(system, meta):
    if meta["multi_gpu"]:
        system.cuda_init_handle.call_method("set_device_id_per_rank")
    return espressomd.lb.LBFluid(
        agrid=meta["agrid"], tau=system.time_step,
        kinematic_viscosity=meta["kinematic_viscosity"],
        density=meta["density"], single_precision=meta["single_precision"],
        gpu=meta["gpu"] or meta["multi_gpu"],
        blocks_per_mpi_rank=list(meta["blocks_per_mpi_rank"]))


def build_and_tune(system, args):
    assert args.volume_fraction > 0, "--volume_fraction must be a positive number"
    assert args.volume_fraction < np.pi / (3 * np.sqrt(2)), \
        "--volume_fraction exceeds the physical limit of sphere packing (~0.74)"
    has_box_l = hasattr(args, "box_l")

    if args.node_grid is not None:
        system.cell_system.node_grid = args.node_grid

    n_part = benchmarks.resolve_n_part(system, args)
    assert not (has_box_l and n_part != 0), \
        "Argument --box_l requires 0 particles (--particles_per_core=0 or --n_particles=0)"

    if n_part == 0:
        box_l = [32] if not has_box_l else args.box_l
        box_l = 3 * box_l if len(box_l) == 1 else box_l
        agrid = 1.
        lb_grid = list(box_l)
        auto_steps = 80
        box_l = list(box_l)
    else:
        mpi_factor = min(2., float(np.amax(system.cell_system.node_grid)))
        box_l = (n_part * 4. / 3. * np.pi * (LJ_SIG / 2.)**3
                 / args.volume_fraction)**(1. / 3.)
        lb_grid = (n_part * args.lb_sites_per_particle)**(1. / 3.)
        lb_grid = int(mpi_factor * np.ceil(lb_grid / mpi_factor))
        agrid = box_l / lb_grid
        auto_steps = 40
        lb_grid = 3 * [lb_grid]
        box_l = 3 * [box_l]

    steps = measurement_steps if measurement_steps is not None else auto_steps

    if args.weak_scaling:
        box_l = list(np.array(box_l) * system.cell_system.node_grid)

    print(f"box length: {box_l}")
    print(f"LB shape: {lb_grid}")
    print(f"LB agrid: {agrid:.3f}")

    system.box_l = box_l
    system.cell_system.skin = INITIAL_SKIN if args.skin is None else args.skin

    if n_part:
        espressomd.assert_features(["LENNARD_JONES"])
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=LJ_EPS, sigma=LJ_SIG, cutoff=LJ_CUT, shift="auto")
        system.part.add(pos=np.random.random((n_part, 3)) * system.box_l)
        benchmarks.minimize(system, 1 / system.time_step)
        system.integrator.set_vv()
        system.thermostat.set_langevin(kT=KT, gamma=GAMMA, seed=SEED)
        print("Tune skin: {:.3f}".format(benchmarks.tune_skin_unless_fixed(
            system, args, min_skin, max_skin, tol=0.05, int_steps=100)))
        print("Equilibration")
        system.integrator.run(500)
        print("Tune skin: {:.3f}".format(benchmarks.tune_skin_unless_fixed(
            system, args, min_skin, max_skin, tol=0.05, int_steps=100)))
        print("Equilibration")
        system.integrator.run(500)
        system.thermostat.turn_off()

    meta = {
        "agrid": float(agrid), "tau": float(system.time_step),
        "kinematic_viscosity": 1., "density": 1.,
        "single_precision": bool(args.single_precision),
        "blocks_per_mpi_rank": [int(x) for x in args.blocks_per_mpi_rank],
        "gpu": bool(args.gpu), "multi_gpu": bool(args.multi_gpu),
    }
    lbf = make_lb(system, meta)
    system.lb = lbf
    if n_part:
        system.thermostat.set_lb(LB_fluid=lbf, gamma=1., seed=SEED)

    return {"n_part": n_part, "measurement_steps": steps,
            "lbf": lbf, "lb_meta": meta}


def save_lb_state(system, args, ctx):
    lbf = ctx["lbf"]
    meta = dict(ctx["lb_meta"])
    meta.update({
        "skin": float(system.cell_system.skin),
        "box_l": [float(x) for x in system.box_l],
        "n_part": int(ctx["n_part"]),
        "has_particles": bool(ctx["n_part"]),
        "time_step": float(system.time_step),
        "measurement_steps": int(ctx["measurement_steps"]),
        "n_iterations": n_iterations,
        "kT": KT, "gamma": GAMMA, "seed": SEED,
    })
    meta.update(benchmarks.topology_meta(system))
    arrays = {
        "lb_velocity": np.copy(lbf[:, :, :].velocity)
    }
    if ctx["n_part"]:
        p = system.part.all()
        arrays["pos"] = np.copy(p.pos)
        arrays["vel"] = np.copy(p.v)
    benchmarks.save_state(args.state_file, meta, **arrays)
    print(f"Saved state to {args.state_file}")


def run_from_state(system, args):
    meta, handle = benchmarks.load_state(args.state_file)
    benchmarks.verify_topology(system, meta)
    system.box_l = meta["box_l"]
    system.time_step = meta["time_step"]
    system.cell_system.skin = args.skin if args.skin is not None else meta["skin"]

    if meta["has_particles"]:
        espressomd.assert_features(["LENNARD_JONES"])
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=LJ_EPS, sigma=LJ_SIG, cutoff=LJ_CUT, shift="auto")
        system.part.add(pos=handle["pos"], v=handle["vel"])
        system.integrator.set_vv()

    lbf = make_lb(system, meta)
    system.lb = lbf
    lbf[:, :, :].velocity = handle["lb_velocity"]
    if meta["has_particles"]:
        system.thermostat.set_lb(
            LB_fluid=lbf, gamma=meta["gamma"], seed=int(meta["seed"]))
    return meta["measurement_steps"], meta["n_iterations"]


def time_and_report(system, args, measurement_steps, n_iterations):
    timings = benchmarks.get_timings(system, measurement_steps, n_iterations)
    avg, ci = benchmarks.get_average_time(timings)
    print(f"average: {1000 * avg:.2f} +/- {1000 * ci:.2f} ms (95% C.I.)")
    n_proc = system.cell_system.get_state()["n_nodes"]
    benchmarks.write_report(args.output, n_proc, timings, measurement_steps,
                            n_threads=benchmarks.get_omp_num_threads(system))


if args.mode == "run":
    measurement_steps, n_iterations = run_from_state(system, args)
    time_and_report(system, args, measurement_steps, n_iterations)
    sys.exit(0)

ctx = build_and_tune(system, args)

if args.state_file:
    save_lb_state(system, args, ctx)

if args.mode == "tune":
    sys.exit(0)

time_and_report(system, args, ctx["measurement_steps"], n_iterations)
