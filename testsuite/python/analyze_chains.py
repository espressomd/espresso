# Copyright (C) 2010-2019 The ESPResSo project
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
import unittest as ut
import unittest_decorators as utx
import numpy as np
import espressomd
import espressomd.interactions
import espressomd.polymer


@utx.skipIfMissingFeatures("LENNARD_JONES")
class AnalyzeChain(ut.TestCase):
    system = espressomd.System(box_l=3 * [1.0])
    np.random.seed(1234)

    num_poly = 2
    num_mono = 5

    @classmethod
    def setUpClass(cls):
        box_l = 20.0
        # start with a small box
        cls.system.box_l = np.array(3 * [box_l])
        cls.system.cell_system.set_n_square(use_verlet_lists=False)
        fene = espressomd.interactions.FeneBond(k=30, d_r_max=2)
        cls.system.bonded_inter.add(fene)
        positions = espressomd.polymer.linear_polymer_positions(
            n_polymers=cls.num_poly, bond_length=0.9,
            beads_per_chain=cls.num_mono, seed=42)
        for p in positions:
            for ndx, m in enumerate(p):
                partcl = cls.system.part.add(pos=m)
                if ndx > 0:
                    partcl.add_bond((fene, partcl.id - 1))
        # bring two polymers to opposite corners:
        # far in cell centre, but mirror images are close
        pol_0_parts = cls.system.part.by_ids(range(0, cls.num_mono))
        cm = np.mean(pol_0_parts.pos, axis=0)
        pol_0_parts.pos += (- cm + cls.system.box_l)
        head_id = cls.num_mono
        tail_id = head_id + cls.num_mono
        pol_1_parts = cls.system.part.by_ids(range(head_id, tail_id))
        cm = np.mean(pol_1_parts.pos, axis=0)
        pol_1_parts.pos -= cm

    # python version of the espresso core function,
    # does not check mirror distances
    def calc_re(self):
        head_ids = np.arange(0, self.num_poly * self.num_mono, self.num_mono)
        tail_ids = head_ids + self.num_mono - 1
        dist = self.system.part.by_ids(
            head_ids).pos - self.system.part.by_ids(tail_ids).pos
        dist = np.sum(dist**2, axis=-1)
        return np.mean(np.sqrt(dist)), np.std(
            np.sqrt(dist)), np.mean(dist), np.std(dist)

    # python version of the espresso core function,
    # does not check mirror distances
    def calc_rg(self):
        head_ids = np.arange(0, self.num_poly * self.num_mono, self.num_mono)
        tail_ids = head_ids + self.num_mono
        rg2 = []
        for p in range(self.num_poly):
            rg2.append(
                np.var(self.system.part.by_ids(range(head_ids[p], tail_ids[p])).pos, axis=0))
        rg2 = np.array(rg2)
        rg2 = np.sum(rg2, axis=1)
        return np.mean(np.sqrt(rg2)), np.std(
            np.sqrt(rg2)), np.mean(rg2), np.std(rg2)

    # python version of the espresso core function,
    # does not check mirror distances
    def calc_rh(self):
        head_ids = np.arange(0, self.num_poly * self.num_mono, self.num_mono)
        tail_ids = head_ids + self.num_mono - 1
        rh = []
        for p in range(self.num_poly):
            r = np.array(
                self.system.part.by_ids(
                    range(
                        head_ids[p],
                        tail_ids[p] +
                        1)).pos)
            # this generates indices for all i<j combinations
            ij = np.triu_indices(len(r), k=1)
            r_ij = r[ij[0]] - r[ij[1]]
            dist = np.linalg.norm(r_ij, axis=1)
            # rh.append(self.num_mono*self.num_mono*0.5/(np.sum(1./dist)))
            # the other way do it, with the proper prefactor of N(N-1)
            rh.append(1. / np.mean(1. / dist))
        rh = np.array(rh)
        return np.mean(rh), np.std(rh)

    # python version of the espresso core function,
    # does not check mirror distances
    # test core results versus python variants (no PBC)
    def test_radii(self):
        # increase PBC for remove mirror images
        all_partcls = self.system.part.all()
        old_pos = all_partcls.pos.copy()
        self.system.box_l = self.system.box_l * 2.
        all_partcls.pos = old_pos
        # compare calc_re()
        core_re = self.system.analysis.calc_re(chain_start=0,
                                               number_of_chains=self.num_poly,
                                               chain_length=self.num_mono)
        np.testing.assert_allclose(core_re, self.calc_re())
        # compare calc_rg()
        core_rg = self.system.analysis.calc_rg(chain_start=0,
                                               number_of_chains=self.num_poly,
                                               chain_length=self.num_mono)
        np.testing.assert_allclose(core_rg, self.calc_rg())
        # compare calc_rh()
        core_rh = self.system.analysis.calc_rh(chain_start=0,
                                               number_of_chains=self.num_poly,
                                               chain_length=self.num_mono)
        np.testing.assert_allclose(core_rh, self.calc_rh())
        # restore PBC
        self.system.box_l = self.system.box_l / 2.
        all_partcls.pos = old_pos


if __name__ == "__main__":
    ut.main()
