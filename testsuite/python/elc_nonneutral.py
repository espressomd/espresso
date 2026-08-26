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

import unittest as ut
import unittest_decorators as utx
import espressomd.electrostatics

import numpy as np


def image_series(z, h, q, delta_bot, delta_top, n_max=100000):
    """Energy and z-force of a charge between two dielectric interfaces.

    A charge `q` at height `z` in a slab `0 < z < h` generates the
    geometric image series of Tyagi et al. 2008, Eqs. (2.3)-(2.6)::

        q*db*D^m  at  -(2 m h + z)      m = 0,1,2,...
        q*D^m     at  -(2 m h - z)      m = 1,2,3,...
        q*dt*D^m  at  +(2(m+1) h - z)   m = 0,1,2,...
        q*D^m     at  +(2 m h + z)      m = 1,2,3,...

    with ``D = delta_bot*delta_top``. Returns the polarization energy
    ``1/2 sum q q_img / d`` and the force ``sum q q_img d / |d|^3``.
    """
    D = delta_bot * delta_top
    m = np.arange(n_max + 1, dtype=float)
    m1 = m[1:]
    q_img, z_img = [], []
    if delta_bot != 0.:
        q_img.append(q * delta_bot * D**m)
        z_img.append(-(2. * m * h + z))
    if delta_top != 0.:
        q_img.append(q * delta_top * D**m)
        z_img.append(2. * (m + 1.) * h - z)
    if D != 0.:
        q_img.append(q * D**m1)
        z_img.append(-(2. * m1 * h - z))
        q_img.append(q * D**m1)
        z_img.append(2. * m1 * h + z)
    q_img = np.concatenate(q_img)
    d = z - np.concatenate(z_img)
    return (0.5 * np.sum(q * q_img / np.abs(d)),
            float(np.sum(q * q_img * d / np.abs(d)**3)))


class ELCNonNeutralTest:
    BOX_L = 16.
    BOX_H = 4.
    GAP = 4.

    system = espressomd.System(box_l=[BOX_L, BOX_L, BOX_H + GAP],
                               time_step=0.01)
    system.cell_system.skin = 0.

    # a net-charged group used by several tests
    POSITIONS = np.array([[3., 3., 0.7], [9., 5., 1.9], [6., 11., 3.1]])
    CHARGES = np.array([1., 1., -1.])

    # walls: metallic (grounded electrodes) and a non-metallic pair
    WALLS = [(-1., -1., True), (0.6, 0.6, False)]

    def tearDown(self):
        self.system.part.clear()
        self.system.electrostatics.clear()

    def setup_elc(self, delta_mid_bot=-1., delta_mid_top=-1., const_pot=True,
                  pot_diff=0., mesh=56, cao=7, r_cut=3.57874842826277,
                  alpha=1.0174111158710029):
        p3m = espressomd.electrostatics.P3M(
            prefactor=1., accuracy=1e-7, mesh=mesh, cao=cao, r_cut=r_cut,
            alpha=alpha, tune=False, check_neutrality=False,
            **self.p3m_params)
        self.system.electrostatics.solver = espressomd.electrostatics.ELC(
            actor=p3m, gap_size=self.GAP, maxPWerror=1e-7,
            check_neutrality=False, delta_mid_bot=delta_mid_bot,
            delta_mid_top=delta_mid_top, const_pot=const_pot,
            pot_diff=pot_diff)

    def test_reflection_symmetry(self):
        """
        Two identical interfaces make the slab invariant under
        ``z -> box_h - z``. Reflecting every charge is therefore an isometry:
        the energy is unchanged, in-plane forces are unchanged and the normal
        forces flip sign.
        """
        system = self.system
        self.assertNotAlmostEqual(np.sum(self.CHARGES), 0.)

        for delta, const_pot in [(-1., True), (0.6, False)]:
            with self.subTest(delta_mid=delta):
                system.part.clear()
                system.electrostatics.clear()
                parts = system.part.add(pos=self.POSITIONS, q=self.CHARGES)
                self.setup_elc(delta_mid_bot=delta, delta_mid_top=delta,
                               const_pot=const_pot)

                energy = system.analysis.energy()['coulomb']
                system.integrator.run(0)
                forces = np.copy(parts.f)

                parts.pos = np.column_stack(
                    [self.POSITIONS[:, :2], self.BOX_H - self.POSITIONS[:, 2]])
                energy_mirror = system.analysis.energy()['coulomb']
                system.integrator.run(0)
                forces_mirror = np.copy(parts.f)

                np.testing.assert_allclose(energy_mirror, energy,
                                           atol=self.atol, rtol=0.)
                np.testing.assert_allclose(forces_mirror[:, :2],
                                           forces[:, :2],
                                           atol=self.atol, rtol=0.)
                np.testing.assert_allclose(forces_mirror[:, 2],
                                           -forces[:, 2],
                                           atol=self.atol, rtol=0.)

    def test_energy_force_consistency(self):
        system = self.system
        base = np.array([[8., 8., 1.7], [3., 5., 2.9]])
        eps = 1e-3
        neutral, charged = [1., -1.], [1., 1.]

        regimes = [
            ("neutral, metallic", neutral, -1., -1., True, 0.),
            ("charged, metallic", charged, -1., -1., True, 0.),
            ("charged, metallic, pot_diff", charged, -1., -1., True, -3.),
            ("charged, dielectric", charged, 0.6, 0.6, False, 0.),
        ]

        for name, charges, db, dt, const_pot, pot_diff in regimes:
            with self.subTest(regime=name):
                system.part.clear()
                system.electrostatics.clear()
                probe = system.part.add(pos=base[0], q=charges[0])
                system.part.add(pos=base[1], q=charges[1])
                self.setup_elc(delta_mid_bot=db, delta_mid_top=dt,
                               const_pot=const_pot, pot_diff=pot_diff,
                               r_cut=3.4805082753300667,
                               alpha=1.0386519879600566)

                energies = []
                for shift in (-eps, eps):
                    probe.pos = base[0] + [0., 0., shift]
                    energies.append(system.analysis.energy()['coulomb'])
                probe.pos = base[0]
                system.integrator.run(0)

                gradient = -(energies[1] - energies[0]) / (2. * eps)
                np.testing.assert_allclose(probe.f[2], gradient,
                                           rtol=1e-4, atol=1e-7)

    def test_grounded_plates_vs_image_series(self):
        """
        A single charge between two grounded metallic plates. Perfect
        conductors screen massively, so the 2D-periodic energy and force are
        almost exactly the image-charge results, with no PBC-correction.
        """
        system = self.system
        q = 1.
        p = system.part.add(pos=[self.BOX_L / 2., self.BOX_L / 2., 1.], q=q)
        self.setup_elc(r_cut=3.3050408844304116, alpha=1.080290324528693)

        for z in (0.8, 1.6, 2.4, 3.2):
            with self.subTest(z=z):
                p.pos = [self.BOX_L / 2., self.BOX_L / 2., z]
                energy = system.analysis.energy()['coulomb']
                system.integrator.run(0)
                ref_energy, ref_force = image_series(z, self.BOX_H, q,
                                                     -1., -1.)
                np.testing.assert_allclose(p.f[2], ref_force,
                                           atol=self.atol_image, rtol=0.)
                np.testing.assert_allclose(energy, ref_energy,
                                           atol=self.atol_image, rtol=0.)

    def test_energy_is_continuous(self):
        """
        ELC switches between an explicit image charge and an analytic image sum at
        a some layer height. Check that there is no jump between the layers.
        """
        system = self.system
        # pin r_cut so the internal switching height is reproducible
        r_cut = 2.0
        switch = min(self.GAP / 3., self.GAP - r_cut, self.BOX_H / 2.)
        offsets = np.array([-.04, -.03, -.02, -.01, .01, .02, .03, .04])

        for db, dt, const_pot in self.WALLS:
            with self.subTest(delta_mid_bot=db, delta_mid_top=dt):
                system.part.clear()
                system.electrostatics.clear()
                p = system.part.add(
                    pos=[self.BOX_L / 2., self.BOX_L / 2., 1.], q=1.)
                self.setup_elc(delta_mid_bot=db, delta_mid_top=dt,
                               const_pot=const_pot, mesh=98, r_cut=r_cut,
                               alpha=1.802701676542639)

                energies = []
                for off in offsets:
                    p.pos = [self.BOX_L / 2., self.BOX_L / 2., switch + off]
                    energies.append(system.analysis.energy()['coulomb'])
                energies = np.array(energies)

                below = np.polyval(np.polyfit(
                    offsets[:4], energies[:4], 1), 0.)
                above = np.polyval(np.polyfit(
                    offsets[4:], energies[4:], 1), 0.)
                np.testing.assert_allclose(above, below, atol=1e-5, rtol=0.)


@utx.skipIfMissingFeatures(["P3M"])
class ELCNonNeutralTestCPU(ELCNonNeutralTest, ut.TestCase):

    p3m_params = {"gpu": False}
    atol = 1e-9
    atol_image = 1e-4


if __name__ == "__main__":
    ut.main()
