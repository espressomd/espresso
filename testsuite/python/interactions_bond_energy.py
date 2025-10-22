#
# Copyright (C) 2013-2025 The ESPResSo project
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
import numpy as np

import espressomd
import espressomd.electrostatics
import espressomd.interactions

#
# Analytical expressions for interactions
#


def harmonic_potential(scalar_r, k, r_0):
    return 0.5 * k * (scalar_r - r_0)**2


def angle_harmonic_potential(phi, bend=1.0, phi0=np.pi):
    return 0.5 * bend * (phi - phi0)**2


def angle(a, b):
    """Get angle in radians."""
    dot_product = np.dot(a, b)
    norms = np.linalg.norm(a) * np.linalg.norm(b)
    # Clipping to avoid numerical errors
    cosine_angle = np.clip(dot_product / norms, -1.0, 1.0)
    return np.arccos(cosine_angle)


class InteractionsBondedTest(ut.TestCase):
    system = espressomd.System(box_l=[17.0, 9.0, 8.0])

    box_l = 20.

    def setUp(self):
        self.system.box_l = [self.box_l] * 3
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.2

    def test(self):
        system = self.system
        hb1 = espressomd.interactions.HarmonicBond(k=1, r_0=0.3, r_cut=3)
        hb2 = espressomd.interactions.HarmonicBond(k=3, r_0=3, r_cut=3)
        ab = espressomd.interactions.AngleHarmonic(bend=2, phi0=0)
        for b in [hb1, hb2, ab]:
            system.bonded_inter.add(b)

        p1 = system.part.add(pos=[0, 0, 0])
        p2 = system.part.add(pos=[1, 0, 0])
        p3 = system.part.add(pos=[-1, -1, -1])

        p1.bonds = ((hb1, p2), (ab, p3, p2))
        p3.bonds = ((hb2, p2),)

        np.testing.assert_allclose(  # hb1 p1 p2
            system.analysis.particle_bond_energy(p1, p1.bonds[0]),
            harmonic_potential(system.distance(p1, p2), hb1.k, hb1.r_0)
        )
        np.testing.assert_allclose(  # hb2 p3 p2
            system.analysis.particle_bond_energy(p3, p3.bonds[0]),
            harmonic_potential(system.distance(p3, p2), hb2.k, hb2.r_0)
        )

        r_13 = system.distance_vec(p1, p3)
        r_12 = system.distance_vec(p1, p2)
        np.testing.assert_allclose(  # ab p3-p1-p2
            system.analysis.particle_bond_energy(p1, p1.bonds[1]),
            angle_harmonic_potential(angle(r_13, r_12), ab.bend, ab.phi0)
        )


if __name__ == '__main__':
    ut.main()
