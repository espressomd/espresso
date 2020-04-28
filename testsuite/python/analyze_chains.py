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
from espressomd.interactions import FeneBond
from espressomd import polymer


@utx.skipIfMissingFeatures("LENNARD_JONES")
class AnalyzeChain(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    np.random.seed(1234)

    num_poly = 2
    num_mono = 5

    @classmethod
    def setUpClass(cls):
        box_l = 20.0
        # start with a small box
        cls.system.box_l = np.array([box_l, box_l, box_l])
        cls.system.cell_system.set_n_square(use_verlet_lists=False)
        fene = FeneBond(k=30, d_r_max=2)
        cls.system.bonded_inter.add(fene)
        positions = polymer.linear_polymer_positions(n_polymers=cls.num_poly,
                                                     bond_length=0.9,
                                                     beads_per_chain=cls.num_mono,
                                                     seed=42)
        for p in positions:
            for ndx, m in enumerate(p):
                part_id = len(cls.system.part)
                cls.system.part.add(id=part_id, pos=m)
                if ndx > 0:
                    cls.system.part[part_id].add_bond((fene, part_id - 1))
        # bring two polymers to opposite corners:
        # far in cell centre, but mirror images are close
        head_id = 0
        tail_id = head_id + cls.num_mono
        cm = np.mean(cls.system.part[head_id:tail_id].pos, axis=0)
        cls.system.part[head_id:tail_id].pos = cls.system.part[
            head_id:tail_id].pos - cm + cls.system.box_l
        head_id = cls.num_mono + 1
        tail_id = head_id + cls.num_mono
        cm = np.mean(cls.system.part[head_id:tail_id].pos, axis=0)
        cls.system.part[head_id:tail_id].pos -= cm

    # python version of the espresso core function,
    # does not check mirror distances
    def calc_re(self):
        head_id = np.arange(0, self.num_poly * self.num_mono, self.num_mono)
        tail_id = head_id + self.num_mono - 1
        dist = self.system.part[head_id].pos - self.system.part[tail_id].pos
        dist = np.sum(dist**2, axis=-1)
        return np.mean(np.sqrt(dist)), np.std(
            np.sqrt(dist)), np.mean(dist), np.std(dist)

    # python version of the espresso core function,
    # does not check mirror distances
    def calc_rg(self):
        head_id = np.arange(0, self.num_poly * self.num_mono, self.num_mono)
        tail_id = head_id + self.num_mono - 1
        rg2 = []
        for p in range(self.num_poly):
            rg2.append(
                np.var(self.system.part[head_id[p]:tail_id[p] + 1].pos, axis=0))
        rg2 = np.array(rg2)
        rg2 = np.sum(rg2, axis=1)
        return np.mean(np.sqrt(rg2)), np.std(
            np.sqrt(rg2)), np.mean(rg2), np.std(rg2)

    # python version of the espresso core function,
    # does not check mirror distances
    def calc_rh(self):
        head_id = np.arange(0, self.num_poly * self.num_mono, self.num_mono)
        tail_id = head_id + self.num_mono - 1
        rh = []
        for p in range(self.num_poly):
            r = np.array(self.system.part[head_id[p]:tail_id[p] + 1].pos)
            # this generates indices for all i<j combinations
            ij = np.triu_indices(len(r), k=1)
            r_ij = r[ij[0]] - r[ij[1]]
            dist = np.sqrt(np.sum(r_ij**2, axis=1))
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
        old_pos = self.system.part[:].pos.copy()
        self.system.box_l = self.system.box_l * 2.
        self.system.part[:].pos = old_pos
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
        self.system.part[:].pos = old_pos


if __name__ == "__main__":
    ut.main()
