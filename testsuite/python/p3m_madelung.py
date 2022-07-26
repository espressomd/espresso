#
# Copyright (C) 2013-2022 The ESPResSo project
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

import numpy as np
import scipy.special
import itertools
import unittest as ut
import unittest_decorators as utx
import espressomd.electrostatics
import espressomd.magnetostatics


class Test(ut.TestCase):
    """
    Check all P3M algorithms against the Madelung constant of 1D, 2D and 3D
    NaCl lattices. See user guide sections :ref:`Madelung electrostatics` and
    :ref:`Madelung magnetostatics` for more details.
    """

    system = espressomd.System(box_l=[1., 1., 1.])
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()

    def get_normalized_obs_per_ion(self, pressure=True):
        energy = self.system.analysis.energy()["coulomb"]
        if pressure:
            p_scalar = self.system.analysis.pressure()["coulomb"]
            p_tensor = self.system.analysis.pressure_tensor()["coulomb"]
        else:
            p_scalar = 0.
            p_tensor = np.zeros((3, 3))
        N = len(self.system.part)
        V = self.system.volume()
        return 2. * energy / N, 2. * p_scalar * V / N, 2. * p_tensor * V / N

    def get_reference_obs_per_ion(self, madelung, base_vector):
        base_tensor = base_vector * np.eye(3)
        ref_energy = madelung
        ref_pressure = madelung * base_tensor / np.trace(base_tensor)
        return ref_energy, ref_pressure

    def get_normalized_obs_per_dipole(self, dipm, spacing):
        energy = self.system.analysis.energy()["dipolar"]
        return energy / len(self.system.part) / dipm**2 * spacing**3

    @utx.skipIfMissingFeatures(["DP3M"])
    def test_infinite_magnetic_wire(self):
        n_part = 64
        spacing = 1.2
        dipm = 1.05
        self.system.box_l = 3 * [n_part * spacing]

        dipoles_parallel = []
        dipoles_antiparallel = []
        for j in range(n_part):
            dipoles_antiparallel.append([(-1)**j, 0., 0.])
            dipoles_parallel.append([1., 0., 0.])
            pos = np.array([j, 0, 0])
            self.system.part.add(pos=spacing * pos, dip=[1., 0., 0.])

        dp3m_params = dict(prefactor=1., accuracy=1e-9, mesh=48, r_cut=28.,
                           cao=6, tuning=False)
        mdlc_params = {'maxPWerror': 1e-9, 'gap_size': 16.}

        def check():
            # minimal energy configuration
            self.system.part.all().dip = dipm * np.array(dipoles_parallel)
            mc = -2. * scipy.special.zeta(3)
            energy_per_dip = self.get_normalized_obs_per_dipole(dipm, spacing)
            np.testing.assert_allclose(energy_per_dip, mc, atol=0., rtol=2e-9)

            # maximal energy configuration
            self.system.part.all().dip = dipm * np.array(dipoles_antiparallel)
            mc = 3. / 2. * scipy.special.zeta(3)
            energy_per_dip = self.get_normalized_obs_per_dipole(dipm, spacing)
            np.testing.assert_allclose(energy_per_dip, mc, atol=0., rtol=2e-9)

        dp3m = espressomd.magnetostatics.DipolarP3M(**dp3m_params)
        self.system.actors.add(dp3m)
        check()
        self.system.actors.remove(dp3m)
        mdlc = espressomd.magnetostatics.DLC(actor=dp3m, **mdlc_params)
        self.system.actors.add(mdlc)
        check()

    @utx.skipIfMissingFeatures(["DP3M"])
    def test_infinite_magnetic_sheet(self):
        n_part = 64
        spacing = 1.2
        dipm = 1.05
        self.system.box_l = 3 * [n_part * spacing]

        dipoles_parallel = []
        dipoles_antiparallel = []
        for j, k in itertools.product(range(n_part), repeat=2):
            if j % 2 == 0:
                theta = (-1)**(k + 0) * np.pi / 4. + np.pi
            else:
                theta = (-1)**(k + 1) * np.pi / 4.
            dipoles_antiparallel.append([np.cos(theta), np.sin(theta), 0.])
            dipoles_parallel.append([0., (-1)**j, 0.])
            pos = np.array([j, k, 0])
            self.system.part.add(pos=spacing * pos, dip=[0., 0., 1.])

        dp3m_params = dict(prefactor=1., accuracy=1e-9, mesh=48, r_cut=33.6,
                           cao=7, tuning=False)
        mdlc_params = {'maxPWerror': 1e-9, 'gap_size': 16.}

        def check():
            # minimal energy configuration
            self.system.part.all().dip = dipm * np.array(dipoles_parallel)
            mc = -2.54944
            energy_per_dip = self.get_normalized_obs_per_dipole(dipm, spacing)
            np.testing.assert_allclose(energy_per_dip, mc, atol=0., rtol=4e-6)

            # maximal energy configuration
            self.system.part.all().dip = dipm * np.array(dipoles_antiparallel)
            mc = +3.01716
            energy_per_dip = self.get_normalized_obs_per_dipole(dipm, spacing)
            np.testing.assert_allclose(energy_per_dip, mc, atol=0., rtol=4e-6)

        dp3m = espressomd.magnetostatics.DipolarP3M(**dp3m_params)
        self.system.actors.add(dp3m)
        check()
        self.system.actors.remove(dp3m)
        mdlc = espressomd.magnetostatics.DLC(actor=dp3m, **mdlc_params)
        self.system.actors.add(mdlc)
        check()

    @utx.skipIfMissingFeatures(["DP3M"])
    def test_infinite_magnetic_cube(self):
        n_part = 16
        spacing = 1.2
        dipm = 1.05
        self.system.box_l = 3 * [n_part * spacing]

        dipoles_parallel = []
        dipoles_antiparallel = []
        for i, j, k in itertools.product(range(n_part), repeat=3):
            if i % 2 == 0:
                theta = np.pi + (-1)**j * np.pi / 4.
            else:
                theta = (-1)**(j + 1) * np.pi / 4.
            dipoles_antiparallel.append([np.cos(theta), np.sin(theta), 0.])
            dipoles_parallel.append([0., 0., (-1)**(i + j)])
            pos = np.array([i, j, k])
            self.system.part.add(pos=spacing * pos, dip=[0., 0., 1.])

        dp3m_params = dict(
            prefactor=1., accuracy=5e-8, mesh=64, r_cut=8.2, cao=7, tuning=False)

        def check():
            # minimal energy configuration
            self.system.part.all().dip = dipm * np.array(dipoles_parallel)
            mc = -2.67679
            energy_per_dip = self.get_normalized_obs_per_dipole(dipm, spacing)
            np.testing.assert_allclose(energy_per_dip, mc, atol=0., rtol=1e-6)

            # maximal energy configuration
            self.system.part.all().dip = dipm * np.array(dipoles_antiparallel)
            mc = +4.84473  # imprecise
            energy_per_dip = self.get_normalized_obs_per_dipole(dipm, spacing)
            np.testing.assert_allclose(energy_per_dip, mc, atol=0., rtol=3e-4)

        dp3m = espressomd.magnetostatics.DipolarP3M(**dp3m_params)
        self.system.actors.add(dp3m)
        check()

    @utx.skipIfMissingFeatures(["P3M"])
    def test_infinite_ionic_wire(self):
        n_pairs = 64
        self.system.box_l = [n_pairs / 2, n_pairs / 2, 2 * n_pairs]
        for i in range(n_pairs):
            self.system.part.add(pos=[0., 0., 2. * i + 0.], q=+1.)
            self.system.part.add(pos=[0., 0., 2. * i + 1.], q=-1.)

        mc = -2. * np.log(2.)
        base = [0., 0., 1.]
        ref_energy, ref_pressure = self.get_reference_obs_per_ion(mc, base)
        p3m_params = dict(prefactor=1., accuracy=1e-8, mesh=32, r_cut=14.,
                          cao=7, tuning=False)

        def check():
            energy, p_scalar, p_tensor = self.get_normalized_obs_per_ion()
            np.testing.assert_allclose(energy, ref_energy, atol=0., rtol=5e-7)
            np.testing.assert_allclose(p_scalar, np.trace(ref_pressure) / 3.,
                                       atol=1e-12, rtol=1e-6)
            np.testing.assert_allclose(p_tensor, ref_pressure, atol=1e-12,
                                       rtol=1e-6)

        p3m = espressomd.electrostatics.P3M(**p3m_params)
        self.system.actors.add(p3m)
        check()
        if espressomd.has_features("CUDA"):
            self.system.actors.clear()
            p3m = espressomd.electrostatics.P3MGPU(**p3m_params)
            self.system.actors.add(p3m)
            check()

    @utx.skipIfMissingFeatures(["P3M"])
    def test_infinite_ionic_sheet(self):
        n_pairs = 32
        self.system.box_l = [2 * n_pairs, 2 * n_pairs, 4 * n_pairs]
        for j, k in itertools.product(range(2 * n_pairs), repeat=2):
            self.system.part.add(pos=[j, k, 0.], q=(-1)**(j + k))

        mc = -1.6155426267128247
        base = [1., 1., 0.]
        ref_energy, ref_pressure = self.get_reference_obs_per_ion(mc, base)
        p3m_params = dict(prefactor=1., accuracy=1e-8, mesh=48, r_cut=26.,
                          cao=7, tuning=False)
        elc_params = dict(maxPWerror=1e-6, gap_size=16)

        def check(pressure=True):
            energy, p_scalar, p_tensor = self.get_normalized_obs_per_ion(
                pressure)
            np.testing.assert_allclose(energy, ref_energy, atol=0., rtol=1e-7)
            if not pressure:
                return
            np.testing.assert_allclose(p_scalar, np.trace(ref_pressure) / 3.,
                                       atol=1e-12, rtol=5e-7)
            np.testing.assert_allclose(p_tensor, ref_pressure, atol=1e-12,
                                       rtol=5e-7)

        p3m = espressomd.electrostatics.P3M(**p3m_params)
        self.system.actors.add(p3m)
        check(pressure=True)
        self.system.actors.remove(p3m)
        elc = espressomd.electrostatics.ELC(actor=p3m, **elc_params)
        self.system.actors.add(elc)
        check(pressure=False)
        if espressomd.has_features("CUDA"):
            self.system.actors.clear()
            p3m = espressomd.electrostatics.P3MGPU(**p3m_params)
            self.system.actors.add(p3m)
            check(pressure=True)
            self.system.actors.remove(p3m)
            elc = espressomd.electrostatics.ELC(actor=p3m, **elc_params)
            self.system.actors.add(elc)
            check(pressure=False)

    @utx.skipIfMissingFeatures(["P3M"])
    def test_infinite_ionic_cube(self):
        n_pairs = 8
        self.system.box_l = 3 * [2 * n_pairs]
        for j, k, l in itertools.product(range(2 * n_pairs), repeat=3):
            self.system.part.add(pos=[j, k, l], q=(-1)**(j + k + l))

        mc = -1.74756459463318219
        base = [1., 1., 1.]
        ref_energy, ref_pressure = self.get_reference_obs_per_ion(mc, base)
        p3m_params = dict(prefactor=1., accuracy=3e-7, mesh=44, r_cut=7., cao=7,
                          tuning=False)

        def check():
            energy, p_scalar, p_tensor = self.get_normalized_obs_per_ion()
            np.testing.assert_allclose(energy, ref_energy, atol=0., rtol=1e-6)
            np.testing.assert_allclose(p_scalar, np.trace(ref_pressure) / 3.,
                                       atol=1e-12, rtol=5e-6)
            np.testing.assert_allclose(p_tensor, ref_pressure, atol=1e-12,
                                       rtol=5e-6)

        p3m = espressomd.electrostatics.P3M(**p3m_params)
        self.system.actors.add(p3m)
        check()
        if espressomd.has_features("CUDA"):
            self.system.actors.clear()
            p3m = espressomd.electrostatics.P3MGPU(**p3m_params)
            self.system.actors.add(p3m)
            check()


if __name__ == "__main__":
    ut.main()
