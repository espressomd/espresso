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
import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.interactions
import numpy as np
import itertools


@utx.skipIfMissingFeatures("BOND_CONSTRAINT")
class RigidBondTest(ut.TestCase):

    system = espressomd.System(box_l=[10., 10., 10.])

    def tearDown(self):
        self.system.part.clear()
        self.system.bonded_inter.clear()
        self.system.thermostat.turn_off()

    def test(self):
        target_acc = 1E-3
        tol = 1.2 * target_acc
        system = self.system
        system.cell_system.skin = 0.4
        system.time_step = 0.01
        system.thermostat.set_langevin(kT=1, gamma=1, seed=42)
        rigid_bond = espressomd.interactions.RigidBond(
            r=1.2, ptol=1E-3, vtol=target_acc)
        system.bonded_inter.add(rigid_bond)

        # create polymer
        last_p = None
        for i in range(5):
            p = system.part.add(pos=(i * 1.2, 0, 0))
            if last_p is not None:
                p.add_bond((rigid_bond, last_p))
            last_p = p

        system.integrator.run(5000)

        # check every bond
        p1_iter, p2_iter = itertools.tee(system.part)
        next(p2_iter, None)  # advance second iterator by 1 step
        for p1, p2 in zip(p1_iter, p2_iter):
            d = system.distance(p2, p1)
            v_d = system.distance_vec(p2, p1)
            self.assertAlmostEqual(d, 1.2, delta=tol)
            # Velocity projection on distance vector
            vel_proj = np.dot(p2.v - p1.v, v_d) / d
            self.assertLess(vel_proj, tol)

    def test_trajectory_independent_of_resorts(self):
        """
        RATTLE reads the positions at the start of the time step of ghost
        particles. They are part of RATTLE's own ghost update, since ghost
        particles are rebuilt when particles are resorted during that update.
        A free rigid dimer gas must follow the same trajectory with a skin
        that triggers resorts in most time steps and with a large skin.
        Only sensitive on more than one MPI rank: on one rank, bond partners
        are always resolved to the real particles, not to ghosts.
        """
        system = self.system
        system.time_step = 0.01
        rng = np.random.default_rng(seed=42)
        rigid_bond = espressomd.interactions.RigidBond(
            r=1., ptol=1e-10, vtol=1e-10)
        system.bonded_inter.add(rigid_bond)
        n_dimers = 200
        pos1 = rng.uniform(0., 10., (n_dimers, 3))
        direction = rng.normal(size=(n_dimers, 3))
        direction /= np.linalg.norm(direction, axis=1)[:, np.newaxis]
        vel = rng.normal(0., 2., (2 * n_dimers, 3))
        # remove the velocity components along the bonds
        for i in range(n_dimers):
            v_rel = vel[2 * i + 1] - vel[2 * i]
            v_par = np.dot(v_rel, direction[i]) * direction[i]
            vel[2 * i] += 0.5 * v_par
            vel[2 * i + 1] -= 0.5 * v_par

        def run(skin):
            system.part.clear()
            system.cell_system.skin = skin
            p1 = system.part.add(pos=pos1, v=vel[0::2])
            p2 = system.part.add(pos=pos1 + direction, v=vel[1::2])
            for a, b in zip(p1, p2):
                b.add_bond((rigid_bond, a))
            system.integrator.run(100)
            p = system.part.all()
            return np.copy(p.pos), np.copy(p.v)

        pos_ref, vel_ref = run(skin=1.5)
        pos, vel_ = run(skin=0.02)
        np.testing.assert_allclose(pos, pos_ref, rtol=0., atol=1e-8)
        np.testing.assert_allclose(vel_, vel_ref, rtol=0., atol=1e-8)


if __name__ == "__main__":
    ut.main()
