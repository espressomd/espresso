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

"""Testmodule for the Widom Insertion.
"""
import unittest as ut
import unittest_decorators as utx
import numpy as np
import espressomd
import espressomd.reaction_ensemble
import tests_common


@utx.skipIfMissingFeatures(["LENNARD_JONES"])
class WidomInsertionTest(ut.TestCase):

    """Test the implementation of the widom insertion.

       The excess chemical potential is calculated for identical particles in
       a 20 cubed box with a single particle, interacting via a LJ-potential
       (cut-off at 5 sigma)."""

    N0 = 1
    TEMPERATURE = 0.5
    TYPE_HA = 0
    CHARGE_HA = 0
    LJ_EPS = 1.0
    LJ_SIG = 1.0
    LJ_CUT = 5
    BOX_L = 2 * LJ_CUT
    LJ_SHIFT = tests_common.lj_potential(
        LJ_CUT, LJ_EPS, LJ_SIG, LJ_CUT + 1.0, 0.0)

    radius = np.linspace(1e-10, LJ_CUT, 1000)
    # numerical integration for radii smaller than the cut-off in spherical
    # coordinates
    integrateUpToCutOff = 4 * np.pi * np.trapz(
        radius**2 * np.exp(-tests_common.lj_potential(radius,
                                                      LJ_EPS,
                                                      LJ_SIG,
                                                      LJ_CUT,
                                                      LJ_SHIFT) / TEMPERATURE),
        x=radius)
    # numerical solution for V_lj=0 => corresponds to the volume (as exp(0)=1)
    integreateRest = (BOX_L**3 - 4.0 / 3.0 * np.pi * LJ_CUT**3)

    # calculate excess chemical potential of the system, see Frenkel Smith,
    # p 174. Note: He uses scaled coordinates, which is why we need to divide
    # by the box volume
    target_mu_ex = -TEMPERATURE * \
        np.log((integrateUpToCutOff + integreateRest) / BOX_L**3)

    system = espressomd.System(box_l=np.ones(3) * BOX_L)
    system.cell_system.set_n_square()
    np.random.seed(69)  # make reaction code fully deterministic
    system.cell_system.skin = 0.4
    volume = system.volume()

    Widom = espressomd.reaction_ensemble.WidomInsertion(
        kT=TEMPERATURE, seed=1)

    # Set the hidden particle type to the lowest possible number to speed
    # up the simulation
    Widom.set_non_interacting_type(1)

    def setUp(self):
        self.system.part.add(pos=0.5 * self.system.box_l, type=self.TYPE_HA)

        self.system.non_bonded_inter[self.TYPE_HA, self.TYPE_HA].lennard_jones.set_params(
            epsilon=self.LJ_EPS, sigma=self.LJ_SIG, cutoff=self.LJ_CUT,
            shift="auto")

        self.Widom.add_reaction(
            reactant_types=[],
            reactant_coefficients=[],
            product_types=[self.TYPE_HA],
            product_coefficients=[1],
            default_charges={self.TYPE_HA: self.CHARGE_HA})

    def test_widom_insertion(self):

        num_samples = 10000
        particle_insertion_potential_energy_samples = []

        for _ in range(num_samples):
            # 0 for insertion reaction
            particle_insertion_potential_energy = self.Widom.calculate_particle_insertion_potential_energy(
                0)
            particle_insertion_potential_energy_samples.append(
                particle_insertion_potential_energy)

        mu_ex_mean, mu_ex_Delta = self.Widom.calculate_excess_chemical_potential(
            particle_insertion_potential_energy_samples=particle_insertion_potential_energy_samples)

        deviation_mu_ex = abs(np.mean(mu_ex_mean) - self.target_mu_ex)

        self.assertLess(
            deviation_mu_ex,
            1e-3,
            msg="\nExcess chemical potential for single LJ-particle computed via Widom insertion is wrong.\n"
            + f"  average mu_ex: {np.mean(mu_ex_mean):.4f}"
            + f"   mu_ex_std_err: {np.std(mu_ex_Delta):.5f}"
            + f"  target_mu_ex: {self.target_mu_ex:.4f}"
        )


if __name__ == "__main__":
    ut.main()
