# Copyright (C) 2010-2018 The ESPResSo project
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
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
import espressomd.electrostatics
from espressomd import electrostatic_extensions


@ut.skipIf(not espressomd.has_features(["ELECTROSTATICS", "PARTIAL_PERIODIC"]),
           "Features not available, skipping test!")
class ELC_vs_MMM2D_neutral(ut.TestCase):
    # Handle to espresso system

    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    acc = 1e-6
    elc_gap = 5.0
    box_l = 10.0
    bl2 = box_l * 0.5
    system.time_step = 0.01
    system.cell_system.skin = 0.1

    def test_elc_vs_mmm2d(self):
        elc_param_sets = {
            "inert": {"gap_size": self.elc_gap, "maxPWerror": self.acc},
            "dielectric": {
                "gap_size": self.elc_gap,
                "maxPWerror": self.acc,
                "delta_mid_bot": 0.1,
                "delta_mid_top": 0.9},
            "const_pot_0": {
                "gap_size": self.elc_gap,
                "maxPWerror": self.acc,
                "const_pot": 1,
                "pot_diff": 0.0},
            "const_pot_1": {
                "gap_size": self.elc_gap,
                "maxPWerror": self.acc,
                "const_pot": 1,
                "pot_diff": 1.0},
            "const_pot_m1": {"gap_size": self.elc_gap, "maxPWerror": self.acc, "const_pot": 1, "pot_diff": -1.0}
        }

        mmm2d_param_sets = {
            "inert": {"prefactor": 1.0, "maxPWerror": self.acc},
            "dielectric": {
                "prefactor": 1.0,
                "maxPWerror": self.acc,
                "dielectric_contrast_on": 1,
                "delta_mid_bot": 0.1,
                "delta_mid_top": 0.9},
            "const_pot_0": {
                "prefactor": 1.0,
                "maxPWerror": self.acc,
                "const_pot": 1,
                "pot_diff": 0.0},
            "const_pot_1": {
                "prefactor": 1.0,
                "maxPWerror": self.acc,
                "const_pot": 1,
                "pot_diff": 1.0},
            "const_pot_m1": {"prefactor": 1.0, "maxPWerror": self.acc, "const_pot": 1, "pot_diff": -1.0}
        }

        self.system.box_l = [self.box_l, self.box_l, self.box_l]
        buf_node_grid = self.system.cell_system.node_grid
        self.system.cell_system.set_layered(
            n_layers=10, use_verlet_lists=False)

        self.system.periodicity = [1, 1, 0]

        q = 1.0
        self.system.part.add(id=0, pos=(5.0, 5.0, 5.0), q=-q)
        self.system.part.add(id=1, pos=(2.0, 2.0, 5.0), q=q / 3.0)
        self.system.part.add(id=2, pos=(2.0, 5.0, 2.0), q=q / 3.0)
        self.system.part.add(id=3, pos=(5.0, 2.0, 7.0), q=q / 3.0)

        #MMM2D
        mmm2d = espressomd.electrostatics.MMM2D(**mmm2d_param_sets["inert"])
        self.system.actors.add(mmm2d)
        mmm2d_res = {}
        mmm2d_res["inert"] = self.scan()

        mmm2d.set_params(**mmm2d_param_sets["dielectric"])
        mmm2d_res["dielectric"] = self.scan()

        mmm2d.set_params(**mmm2d_param_sets["const_pot_0"])
        mmm2d_res["const_pot_0"] = self.scan()

        mmm2d.set_params(**mmm2d_param_sets["const_pot_1"])
        mmm2d_res["const_pot_1"] = self.scan()

        mmm2d.set_params(**mmm2d_param_sets["const_pot_m1"])
        mmm2d_res["const_pot_m1"] = self.scan()

        self.system.actors.remove(mmm2d)

        #ELC
        self.system.box_l = [self.box_l, self.box_l, self.box_l + self.elc_gap]
        self.system.cell_system.set_domain_decomposition(
            use_verlet_lists=True)
        self.system.cell_system.node_grid = buf_node_grid
        self.system.periodicity = [1, 1, 1]
        p3m = espressomd.electrostatics.P3M(prefactor=1.0, accuracy=self.acc,
                                            mesh=[16, 16, 24], cao=6)
        self.system.actors.add(p3m)

        elc = electrostatic_extensions.ELC(**elc_param_sets["inert"])
        self.system.actors.add(elc)
        elc_res = {}
        elc_res["inert"] = self.scan()

        elc.set_params(**elc_param_sets["dielectric"])
        elc_res["dielectric"] = self.scan()

        elc.set_params(**elc_param_sets["const_pot_0"])
        elc_res["const_pot_0"] = self.scan()

        elc.set_params(**elc_param_sets["const_pot_1"])
        elc_res["const_pot_1"] = self.scan()

        elc.set_params(**elc_param_sets["const_pot_m1"])
        elc_res["const_pot_m1"] = self.scan()

        for run in elc_res:
            self.assertTrue(np.testing.assert_allclose(
                mmm2d_res[run], elc_res[run], rtol=0, atol=1e-4) is None)

    def scan(self):
        n = 10
        d = 0.5
        res = []
        for i in range(n + 1):
            z = self.box_l - d - 1.0 * i / n * (self.box_l - 2 * d)
            self.system.part[0].pos = [self.bl2, self.bl2, z]
            self.system.integrator.run(0)
            energy = self.system.analysis.energy()
            m = [z]
            m.extend(self.system.part[0].f)
            m.append(energy['coulomb'])
            res.append(m)

        return res

if __name__ == "__main__":
    ut.main()
