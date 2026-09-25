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
import espressomd.interactions
import benchmarks
import numpy as np
import argparse

LJ_EPS = 1.0  # LJ epsilon
LJ_SIG = 1.0  # particle diameter
LJ_CUT = LJ_SIG * 2**(1. / 6.)  # cutoff distance
DEFAULT_PARTICLES_PER_CORE = 1000
INITIAL_SKIN = 0.4
KT = 1.0
GAMMA = 1.0
SEED = 42
HARMONIC_K = 2.0

# Simulation parameters. The importlib-based benchmark tests (in
# testsuite/benchmarks) override these module-level globals to run a fast
# smoke test. ``measurement_steps = None`` means auto-compute from the
# particle count.
measurement_steps = None
n_iterations = 30
min_skin = 0.2
max_skin = 1.0

parser = argparse.ArgumentParser(description="Benchmark LJ simulations. "
                                 "Save the results to a CSV file.")
benchmarks.add_common_args(parser, DEFAULT_PARTICLES_PER_CORE)
parser.add_argument("--volume_fraction", metavar="FRAC", action="store",
                    type=float, default=0.50, required=False,
                    help="Fraction of the simulation box volume occupied by "
                    "particles (range: [0.01-0.74], default: 0.50)")
parser.add_argument("--bonds", action="store_true",
                    help="Add bonds between particle pairs, default: false")
parser.add_argument("--retune_skin_after", metavar="N", action="store",
                    type=int, default=None, required=False,
                    help="Retune skin every N timing iterations "
                    "(0 disables, default: 5). In 'run' mode this overrides "
                    "the value stored in the state file.")
group = parser.add_mutually_exclusive_group()
group.add_argument("--output", metavar="FILEPATH", action="store",
                   type=str, required=False, default="benchmarks.csv",
                   help="Output file (default: benchmarks.csv)")
group.add_argument("--visualizer", action="store_true",
                   help="Starts the visualizer (for debugging purposes)")

args = parser.parse_args()
benchmarks.validate_mode(args)

espressomd.assert_features(["LENNARD_JONES"])
np.random.seed(SEED)

system = espressomd.System(box_l=[1, 1, 1])
system.time_step = 0.01


def configure_lj(system):
    system.non_bonded_inter[0, 0].lennard_jones.set_params(
        epsilon=LJ_EPS, sigma=LJ_SIG, cutoff=LJ_CUT, shift="auto")


def resolve_retune(args, meta=None):
    """
    Resolve the skin-retune interval for the timing loop. An explicit
    --retune_skin_after always wins. Otherwise --skin (a fixed skin) disables
    retuning; failing that, fall back to the state-file value or the default 5.
    A value <= 0 disables retuning.
    """
    if args.retune_skin_after is not None:
        value = args.retune_skin_after
    elif args.skin is not None:
        value = 0
    elif meta is not None:
        value = int(meta["retune_skin_after"])
    else:
        value = 5
    return None if value <= 0 else value


def build_and_tune(system, args):
    """Build the LJ system, warm up, equilibrate and tune the skin."""
    assert args.volume_fraction > 0, "volume_fraction must be a positive number"
    assert args.volume_fraction < np.pi / (3 * np.sqrt(2)), \
        "volume_fraction exceeds the physical limit of sphere packing (~0.74)"
    assert not (args.bonds and args.volume_fraction > 0.5), \
        "volume_fraction too dense (>0.50) for a diatomic liquid"

    n_part = benchmarks.resolve_n_part(system, args)
    steps = measurement_steps
    if steps is None:
        steps = max(20, int(np.round(
            5e6 / n_part * system.cell_system.get_state()["n_nodes"], -2)))
    if not args.visualizer:
        assert steps >= 20, \
            f"{steps} steps per tick are too short"

    box_l = (n_part * 4. / 3. * np.pi * (LJ_SIG / 2.)**3
             / args.volume_fraction)**(1. / 3.)
    system.box_l = 3 * (box_l,)
    system.cell_system.skin = INITIAL_SKIN if args.skin is None else args.skin
    configure_lj(system)

    bond_pairs = []
    if not args.bonds:
        system.part.add(pos=np.random.random((n_part, 3)) * system.box_l)
    else:
        hb = espressomd.interactions.HarmonicBond(r_0=LJ_CUT, k=HARMONIC_K)
        system.bonded_inter.add(hb)
        for _ in range(0, n_part, 2):
            pos = np.random.random(3) * system.box_l
            vec = np.random.random(3)
            vec /= np.linalg.norm(vec)
            p1 = system.part.add(pos=pos)
            p2 = system.part.add(pos=pos + vec * hb.r_0)
            p1.add_bond((hb, p2))
            bond_pairs.append([p1.id, p2.id])

    benchmarks.minimize(system, 1 / system.time_step)
    system.integrator.set_vv()
    system.thermostat.set_langevin(kT=KT, gamma=GAMMA, seed=SEED)

    print("Equilibration")
    system.time_step /= 10.
    system.integrator.run(100)
    system.time_step *= 10.
    system.integrator.run(min(5 * steps, 60000))
    print("Tune skin: {:.3f}".format(benchmarks.tune_skin_unless_fixed(
        system, args, min_skin, max_skin, tol=0.025, int_steps=200)))
    print("Equilibration")
    system.integrator.run(min(10 * steps, 60000))

    return {"n_part": n_part, "measurement_steps": steps,
            "bond_pairs": bond_pairs}


def save_lj_state(system, args, ctx):
    resolved_retune = resolve_retune(args)
    p = system.part.all()
    meta = {
        "skin": float(system.cell_system.skin),
        "box_l": [float(x) for x in system.box_l],
        "n_part": int(ctx["n_part"]),
        "time_step": float(system.time_step),
        "measurement_steps": int(ctx["measurement_steps"]),
        "n_iterations": n_iterations,
        "volume_fraction": float(args.volume_fraction),
        "bonds": bool(args.bonds),
        "retune_skin_after": 0 if resolved_retune is None else resolved_retune,
        "harmonic_r_0": float(LJ_CUT),
        "harmonic_k": float(HARMONIC_K),
        "kT": KT, "gamma": GAMMA, "seed": SEED,
    }
    meta.update(benchmarks.topology_meta(system))
    benchmarks.save_state(
        args.state_file, meta,
        pos=np.copy(p.pos), vel=np.copy(p.v),
        bond_pairs=np.array(ctx["bond_pairs"], dtype=int))
    print(f"Saved state to {args.state_file}")


def run_from_state(system, args):
    meta, handle = benchmarks.load_state(args.state_file)
    benchmarks.verify_topology(system, meta)
    system.box_l = meta["box_l"]
    system.time_step = meta["time_step"]
    system.cell_system.skin = args.skin if args.skin is not None else meta["skin"]
    configure_lj(system)
    system.part.add(pos=handle["pos"], v=handle["vel"])
    if meta["bonds"]:
        hb = espressomd.interactions.HarmonicBond(
            r_0=meta["harmonic_r_0"], k=meta["harmonic_k"])
        system.bonded_inter.add(hb)
        for p1, p2 in handle["bond_pairs"]:
            system.part.by_id(int(p1)).add_bond((hb, int(p2)))
    system.integrator.set_vv()
    system.thermostat.set_langevin(
        kT=meta["kT"], gamma=meta["gamma"], seed=int(meta["seed"]))
    retune = resolve_retune(args, meta)
    return meta["measurement_steps"], meta["n_iterations"], retune


def time_and_report(system, args, measurement_steps, n_iterations, retune):
    timings = benchmarks.get_timings(
        system, measurement_steps, n_iterations,
        retune_skin_after_steps=retune)
    avg, ci = benchmarks.get_average_time(timings)
    print(f"average: {avg:.3e} +/- {ci:.3e} (95% C.I.)")
    n_proc = system.cell_system.get_state()["n_nodes"]
    benchmarks.write_report(args.output, n_proc, timings, measurement_steps,
                            n_threads=benchmarks.get_omp_num_threads(system))


if args.mode == "run":
    measurement_steps, n_iterations, retune = run_from_state(system, args)
    time_and_report(system, args, measurement_steps, n_iterations, retune)
    sys.exit(0)

ctx = build_and_tune(system, args)

if args.visualizer:
    import time
    import threading
    import espressomd.visualization
    visualizer = espressomd.visualization.openGLLive(system)

    def main_thread():
        while True:
            system.integrator.run(1)
            visualizer.update()
            time.sleep(1 / 60.)

    t = threading.Thread(target=main_thread)
    t.daemon = True
    t.start()
    visualizer.start()

if args.state_file:
    save_lj_state(system, args, ctx)

if args.mode == "tune":
    sys.exit(0)

time_and_report(system, args, ctx["measurement_steps"], n_iterations,
                resolve_retune(args))
