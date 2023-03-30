#
# Copyright (C) 2010-2022 The ESPResSo project
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
import numpy as np
import espressomd
import espressomd.interactions
import espressomd.lees_edwards
import espressomd.polymer


@utx.skipIfMissingFeatures("LENNARD_JONES")
class AnalyzeChain(ut.TestCase):
    system = espressomd.System(box_l=3 * [1.0])
    system.cell_system.set_n_square(use_verlet_lists=False)
    system.time_step = 0.01
    system.cell_system.skin = 0.1
    fene = espressomd.interactions.FeneBond(k=30., d_r_max=2.)
    system.bonded_inter.add(fene)

    np.random.seed(1234)

    def setUp(self):
        box_l = 20.0
        self.system.box_l = np.array(3 * [box_l])

    def tearDown(self):
        self.system.part.clear()

    def insert_polymers(self, num_poly, num_mono):
        positions = espressomd.polymer.linear_polymer_positions(
            n_polymers=num_poly, bond_length=0.9,
            beads_per_chain=num_mono, seed=42)
        for p in positions:
            for ndx, m in enumerate(p):
                partcl = self.system.part.add(pos=m)
                if ndx > 0:
                    partcl.add_bond((self.fene, partcl.id - 1))
        # bring two polymers to opposite corners:
        # far in cell centre, but mirror images are close
        pol_0_parts = self.system.part.by_ids(range(0, num_mono))
        cm = np.mean(pol_0_parts.pos, axis=0)
        pol_0_parts.pos += (- cm + self.system.box_l)
        head_id = num_mono
        tail_id = head_id + num_mono
        pol_1_parts = self.system.part.by_ids(range(head_id, tail_id))
        cm = np.mean(pol_1_parts.pos, axis=0)
        pol_1_parts.pos -= cm

    # python version of the espresso core function,
    # does not check mirror distances
    def calc_re(self, num_poly, num_mono):
        head_ids = np.arange(0, num_poly * num_mono, num_mono)
        tail_ids = head_ids + num_mono - 1
        dist = self.system.part.by_ids(
            head_ids).pos - self.system.part.by_ids(tail_ids).pos
        dist = np.sum(dist**2, axis=-1)
        return np.mean(np.sqrt(dist)), np.std(
            np.sqrt(dist)), np.mean(dist), np.std(dist)

    # python version of the espresso core function,
    # does not check mirror distances
    def calc_rg(self, num_poly, num_mono):
        head_ids = np.arange(0, num_poly * num_mono, num_mono)
        tail_ids = head_ids + num_mono
        rg2 = []
        for p in range(num_poly):
            ids = range(head_ids[p], tail_ids[p])
            r = np.copy(self.system.part.by_ids(ids).pos)
            rg2.append(np.var(r, axis=0))
        rg2 = np.array(rg2)
        rg2 = np.sum(rg2, axis=1)
        return np.mean(np.sqrt(rg2)), np.std(
            np.sqrt(rg2)), np.mean(rg2), np.std(rg2)

    # python version of the espresso core function,
    # does not check mirror distances
    def calc_rh(self, num_poly, num_mono):
        head_ids = np.arange(0, num_poly * num_mono, num_mono)
        tail_ids = head_ids + num_mono - 1
        rh = []
        for p in range(num_poly):
            ids = range(head_ids[p], tail_ids[p] + 1)
            r = np.copy(self.system.part.by_ids(ids).pos)
            # this generates indices for all i<j combinations
            ij = np.triu_indices(len(r), k=1)
            r_ij = r[ij[0]] - r[ij[1]]
            dist = np.linalg.norm(r_ij, axis=1)
            rh.append(1. / np.mean(1. / dist))
        rh = np.array(rh)
        return np.mean(rh), np.std(rh)

    # python version of the espresso core function,
    # does not check mirror distances
    # test core results versus python results (no PBC)
    def test_observables_no_pbc(self):
        num_poly = 2
        num_mono = 5
        self.insert_polymers(num_poly, num_mono)
        # increase PBC to remove interactions with periodic images
        all_partcls = self.system.part.all()
        old_pos = all_partcls.pos.copy()
        self.system.box_l = self.system.box_l * 2.
        all_partcls.pos = old_pos
        # compare calc_re()
        core_re = self.system.analysis.calc_re(chain_start=0,
                                               number_of_chains=num_poly,
                                               chain_length=num_mono)
        np.testing.assert_allclose(core_re, self.calc_re(num_poly, num_mono))
        # compare calc_rg()
        core_rg = self.system.analysis.calc_rg(chain_start=0,
                                               number_of_chains=num_poly,
                                               chain_length=num_mono)
        np.testing.assert_allclose(core_rg, self.calc_rg(num_poly, num_mono))
        # compare calc_rh()
        core_rh = self.system.analysis.calc_rh(chain_start=0,
                                               number_of_chains=num_poly,
                                               chain_length=num_mono)
        np.testing.assert_allclose(core_rh, self.calc_rh(num_poly, num_mono))

    def test_observables_lebc(self):
        lin_protocol = espressomd.lees_edwards.LinearShear(
            initial_pos_offset=-0.02, time_0=0., shear_velocity=0.)
        self.system.lees_edwards.set_boundary_conditions(
            shear_direction="y", shear_plane_normal="x", protocol=lin_protocol)

        # place particles on a line across a Lees-Edwards periodic boundary
        for x in np.arange(-0.2, 0.21, 0.01):
            x += 1e-7  # small offset due to Lees-Edwards edge case at x=0
            pos = [x, 0., 1.2 * x]
            vel = [0., 0., 0.]
            # for particles beyond the shear boundary, place them initially in
            # the central box with an x-velocity that will cause them to cross
            # the shear boundary in 1 time step; the shear y-offset is exactly
            # compensated to form a straight line in unfolded coordinates
            if x < 0.:
                pos[0] += 0.2
                vel[0] -= 0.2 / self.system.time_step
                pos[1] -= lin_protocol.initial_pos_offset
            self.system.part.add(pos=pos, v=vel)
        self.system.integrator.run(1)
        num_poly = 1
        num_mono = len(self.system.part)

        # compare calc_re()
        core_re = self.system.analysis.calc_re(chain_start=0,
                                               number_of_chains=num_poly,
                                               chain_length=num_mono)
        np.testing.assert_allclose(core_re, self.calc_re(num_poly, num_mono))
        # compare calc_rg()
        core_rg = self.system.analysis.calc_rg(chain_start=0,
                                               number_of_chains=num_poly,
                                               chain_length=num_mono)
        np.testing.assert_allclose(core_rg, self.calc_rg(num_poly, num_mono),
                                   atol=1e-8)
        # compare calc_rh()
        core_rh = self.system.analysis.calc_rh(chain_start=0,
                                               number_of_chains=num_poly,
                                               chain_length=num_mono)
        np.testing.assert_allclose(core_rh, self.calc_rh(num_poly, num_mono))

    def test_exceptions(self):
        num_poly = 2
        num_mono = 5
        self.insert_polymers(num_poly, num_mono)
        err_msg = ("Particle with id 10 does not exist; cannot perform "
                   "analysis on the range chain_start=0, number_of_chains=2, "
                   "chain_length=10. Please provide a contiguous range of "
                   "particle ids.")
        analysis = self.system.analysis
        for method in (analysis.calc_re, analysis.calc_rg, analysis.calc_rh):
            with self.assertRaisesRegex(RuntimeError, err_msg):
                method(chain_start=0, number_of_chains=num_poly,
                       chain_length=2 * num_mono)
        self.assertIsNone(analysis.call_method("unknown"))


if __name__ == "__main__":
    ut.main()
