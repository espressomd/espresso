#
# Copyright (C) 2013-2019 The ESPResSo project
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
import espressomd.electrostatics
import espressomd.scafacos
import tests_common


@utx.skipIfMissingFeatures("ELECTROSTATICS")
class CoulombMixedPeriodicity(ut.TestCase):

    """Test mixed periodicity electrostatics"""

    system = espressomd.System(box_l=[10, 10, 10])
    data = np.genfromtxt(tests_common.abspath(
        "data/coulomb_mixed_periodicity_system.data"))

    tolerance_force = 5E-4
    tolerance_energy = 1.8E-3

    # Reference energy from MMM2D
    reference_energy = 216.640984711

    def setUp(self):
        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.

        # Add particles to system and store reference forces in hash
        # Input format: id pos q f
        self.system.part.add(pos=self.data[:, 1:4], q=self.data[:, 4])
        self.forces = self.data[:, 5:8]

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()

    def compare(self, method_name, energy=True):
        # Compare forces and energy now in the system to stored ones

        # Force
        force_diff = np.linalg.norm(
            self.system.part.all().f - self.forces, axis=1)
        self.assertLessEqual(
            np.mean(force_diff), self.tolerance_force,
            "Absolute force difference too large for method " + method_name)

        # Energy
        if energy:
            self.assertAlmostEqual(
                self.system.analysis.energy()["total"],
                self.reference_energy, delta=self.tolerance_energy,
                msg="Absolute energy difference too large for " + method_name)

    @utx.skipIfMissingFeatures(["P3M"])
    def test_elc(self):
        # Make sure, the data satisfies the gap
        for p in self.system.part:
            assert p.pos[2] >= 0. and p.pos[2] <= 9., f'particle {p.id} in gap'

        self.system.cell_system.set_domain_decomposition()
        self.system.cell_system.node_grid = sorted(
            self.system.cell_system.node_grid, key=lambda x: -x)
        self.system.periodicity = [1, 1, 1]

        p3m = espressomd.electrostatics.P3M(
            prefactor=1, accuracy=1e-6, mesh=(64, 64, 64))
        elc = espressomd.electrostatics.ELC(
            p3m_actor=p3m, maxPWerror=1E-6, gap_size=1)

        self.system.actors.add(elc)
        self.system.integrator.run(0)
        self.compare("elc", energy=True)

    @ut.skipIf(not espressomd.has_features("SCAFACOS")
               or 'p2nfft' not in espressomd.scafacos.available_methods(),
               'Skipping test: missing feature SCAFACOS or p2nfft method')
    def test_scafacos_p2nfft(self):
        self.system.periodicity = [1, 1, 0]
        self.system.cell_system.set_domain_decomposition()

        scafacos = espressomd.electrostatics.Scafacos(
            prefactor=1,
            method_name="p2nfft",
            method_params={
                "tolerance_field": 5E-5,
                "pnfft_n": "96,96,128",
                "pnfft_N": "96,96,128",
                "r_cut": 2.4,
                "pnfft_m": 3})
        self.system.actors.add(scafacos)
        self.system.integrator.run(0)
        self.compare("scafacos_p2nfft", energy=True)


if __name__ == "__main__":
    ut.main()
