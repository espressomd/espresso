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
Benchmark a heavily bonded polymer system: bond creation, bond deletion,
repeated bond churn (remove/re-add random bonds with MD steps in between)
and the MD time per step.

The system consists of ``n_chains`` chains of ``chain_len`` beads, each
chain carrying harmonic, angle and dihedral bonds. The chains are placed
on a lattice without overlap and the box is scaled with the number of
chains, so that the particle density stays roughly constant and the number
of bonds can be increased systematically with ``--n_chains``.

The MD time per step is reported without the startup cost of
``integrator.run()`` (neighbor list construction etc.), which is reported
separately as ``startup``. ``system.min_global_cut`` is set explicitly:
without it, a bonded-only system (``--no-wca``) has an undefined
interaction range, a coarse cell grid and a neighbor list construction
that scales quadratically with the number of particles.

The resident memory of every MPI rank is reported before and after the
bonds are created (``rss_before_bonds_mb``, ``rss_after_bonds_mb``), read
from ``/proc``, which only works when all ranks run on the same Linux node.

The script works with both bond APIs (``system.bonds[bond].add()`` and the
older ``ParticleHandle.add_bond()``), so that it can be run unchanged on
older builds for comparison. The results are printed as a JSON line
prefixed with ``RESULT``.

Example::

    for n in 50 100 200 400 800 1600; do
        mpirun -n 4 ./pypresso ../maintainer/benchmarks/bonds.py --n_chains $n
    done
"""

import argparse
import json
import os
import time
import warnings
import numpy as np
import espressomd
import espressomd.interactions

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--n_chains", type=int, default=400,
                    help="number of chains (default: 400)")
parser.add_argument("--chain_len", type=int, default=60,
                    help="number of beads per chain (default: 60)")
parser.add_argument("--md_steps", type=int, default=200,
                    help="number of timed MD steps (default: 200)")
parser.add_argument("--churn_ops", type=int, default=5000,
                    help="number of bonds removed and re-added per churn "
                    "round, 0 disables the churn benchmark (default: 5000)")
parser.add_argument("--churn_rounds", type=int, default=5,
                    help="number of churn rounds (default: 5)")
parser.add_argument("--churn_steps", type=int, default=10,
                    help="MD steps after each removal/addition in a churn "
                    "round (default: 10)")
parser.add_argument("--no-wca", dest="wca", action="store_false",
                    help="do not add WCA interactions (bonded forces only)")
parser.add_argument("--bulk", action="store_true",
                    help="create the bonds with one bulk call per bond type, "
                    "if the bond API supports it")
args = parser.parse_args()

if args.wca:
    espressomd.assert_features(["WCA"])

# the deprecated per-particle bond API is used on older builds only
warnings.simplefilter("ignore", FutureWarning)


def rotate(vec, axis, phi):
    """Rodrigues rotation of ``vec`` around the unit vector ``axis``."""
    return (vec * np.cos(phi) + np.cross(axis, vec) * np.sin(phi) +
            axis * np.dot(axis, vec) * (1. - np.cos(phi)))


# all chains share the same helical conformation (non-degenerate bond and
# dihedral angles); they are placed on a lattice whose cells enclose one
# chain each, so that chains never overlap
theta, phi_step = np.deg2rad(109.5), np.deg2rad(50.)
template = []
pos = np.zeros(3)
direction, up = np.array([1., 0., 0.]), np.array([0., 1., 0.])
for _ in range(args.chain_len):
    template.append(pos.copy())
    direction = rotate(direction, up, np.pi - theta)
    direction /= np.linalg.norm(direction)
    up = rotate(up, direction, phi_step)
    up /= np.linalg.norm(up)
    pos = pos + direction
template = np.array(template) - np.min(template, axis=0)
cell = np.max(template, axis=0) + 1.5  # lattice cell with a safety margin
side = (args.n_chains * np.prod(cell))**(1. / 3.)
n_cells = np.maximum(1, np.round(side / cell)).astype(int)
while np.prod(n_cells) < args.n_chains:
    n_cells[np.argmin(n_cells * cell)] += 1
box_l = n_cells * cell
sites = np.array(list(np.ndindex(*n_cells))[:args.n_chains]) * cell
positions = (sites[:, np.newaxis, :] + template + 0.75).reshape(-1, 3)
chains = [list(range(c * args.chain_len, (c + 1) * args.chain_len))
          for c in range(args.n_chains)]

system = espressomd.System(box_l=box_l)
system.time_step = 0.001
system.cell_system.skin = 0.4
# the interaction range must cover the bond partners of angle and dihedral
# bonds (up to 3 bond lengths apart), otherwise partners can be missing
# from the ghost layer in MPI-parallel runs
system.min_global_cut = 3.0

harmonic = espressomd.interactions.HarmonicBond(k=30.0, r_0=1.0)
angle = espressomd.interactions.AngleHarmonic(bend=10.0, phi0=np.pi)
dihedral = espressomd.interactions.Dihedral(mult=1, bend=2.0, phase=0.0)
for bond in (harmonic, angle, dihedral):
    system.bonded_inter.add(bond)

system.part.add(id=np.arange(len(positions)), pos=positions)

NEW_API = hasattr(system, "bonds")


def add_bond(bond, ids):
    if NEW_API:
        system.bonds[bond].add(tuple(ids))
    else:
        system.part.by_id(ids[0]).add_bond((bond, *ids[1:]))


def remove_bond(bond, ids):
    if NEW_API:
        system.bonds[bond].remove(tuple(ids))
    else:
        system.part.by_id(ids[0]).delete_bond((bond, *ids[1:]))


def rss_per_rank_mb():
    """Resident memory of this process and its sibling MPI ranks."""
    def rss(pid):
        try:
            with open(f"/proc/{pid}/status") as f:
                for line in f:
                    if line.startswith("VmRSS:"):
                        return int(line.split()[1]) / 1024.
        except OSError:
            pass
        return None

    pids = [os.getpid()]
    if res_n_ranks > 1:
        ppid = os.getppid()
        exe = os.path.realpath(f"/proc/{os.getpid()}/exe")
        try:
            with open(f"/proc/{ppid}/task/{ppid}/children") as f:
                children = [int(x) for x in f.read().split()]
        except OSError:
            children = []
        pids = [pid for pid in children
                if os.path.realpath(f"/proc/{pid}/exe") == exe]
    return [x for x in map(rss, pids) if x is not None]


def timed(fn):
    tick = time.perf_counter()
    fn()
    return time.perf_counter() - tick


n = args.chain_len
all_bonds = []
for i in chains:
    all_bonds += [(harmonic, (i[m], i[m + 1])) for m in range(n - 1)]
for i in chains:
    all_bonds += [(angle, (i[m + 1], i[m], i[m + 2])) for m in range(n - 2)]
for i in chains:
    all_bonds += [(dihedral, tuple(i[m:m + 4])) for m in range(n - 3)]

res_n_ranks = system.cell_system.get_state()["n_nodes"]
res = {"bond_api": "system.bonds" if NEW_API else "ParticleHandle",
       "n_ranks": system.cell_system.get_state()["n_nodes"],
       "n_part": len(positions), "n_bonds": len(all_bonds),
       "wca": args.wca}

# bond creation
system.integrator.run(0)
res["rss_before_bonds_mb"] = rss_per_rank_mb()
if args.bulk and NEW_API and hasattr(system.bonds[harmonic], "_is_batch"):
    def add_all_bulk():
        for bond in (harmonic, angle, dihedral):
            rows = [ids for b, ids in all_bonds if b is bond]
            system.bonds[bond].add(np.array(rows))
    res["create_all"] = timed(add_all_bulk)
else:
    res["create_all"] = timed(
        lambda: [add_bond(b, ids) for b, ids in all_bonds])
system.integrator.run(0)
res["rss_after_bonds_mb"] = rss_per_rank_mb()

if args.wca:
    system.non_bonded_inter[0, 0].wca.set_params(epsilon=1., sigma=0.8)

# MD time per step, without the startup cost of integrator.run()
system.integrator.run(10)
res["startup"], res["md_step"] = [], []
for _ in range(2):
    t_startup = timed(lambda: system.integrator.run(0))
    t_run = timed(lambda: system.integrator.run(args.md_steps))
    res["startup"].append(t_startup)
    res["md_step"].append((t_run - t_startup) / args.md_steps)

# bond churn: remove random bonds of all types, run, re-add them, run
n_ops = min(args.churn_ops, len(all_bonds))
if n_ops > 0 and args.churn_rounds > 0:
    rng = np.random.default_rng(42)
    t_rm = t_add = t_md = 0.
    for _ in range(args.churn_rounds):
        sel = [all_bonds[k] for k in
               rng.choice(len(all_bonds), size=n_ops, replace=False)]
        t_rm += timed(lambda: [remove_bond(b, ids) for b, ids in sel])
        t_md += timed(lambda: system.integrator.run(args.churn_steps))
        t_add += timed(lambda: [add_bond(b, ids) for b, ids in sel])
        t_md += timed(lambda: system.integrator.run(args.churn_steps))
    n_total = args.churn_rounds * n_ops
    res["churn_remove_per_op"] = t_rm / n_total
    res["churn_add_per_op"] = t_add / n_total
    res["churn_md_per_step"] = t_md / (2 * args.churn_rounds *
                                       args.churn_steps)

# bond deletion
res["remove_all"] = timed(
    lambda: [remove_bond(b, ids) for b, ids in all_bonds])
assert sum(len(p.bonds) for p in system.part.all()) == 0

print("RESULT " + json.dumps(res))
