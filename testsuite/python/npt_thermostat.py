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
import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.propagation
import tests_common
import numpy as np


@utx.skipIfMissingFeatures("NPT")
class NPTThermostat:

    """Test NpT dynamics"""
    system = espressomd.System(box_l=[2.0, 2.0, 2.0])
    system.cell_system.skin = 0.
    system.periodicity = [True, True, True]

    def setUp(self):
        np.random.seed(42)
        self.system.box_l = [2., 2., 2.]
        self.system.time_step = 0.01

    def tearDown(self):
        self.system.non_bonded_inter.reset()
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.integrator.set_vv()
        if espressomd.has_features("ELECTROSTATICS"):
            self.system.electrostatics.clear()

    def test_01__rng(self):
        """Test for RNG consistency."""
        def reset_particle_and_box():
            self.system.part.clear()
            self.system.box_l = [4, 4, 4]
            p = self.system.part.add(pos=[0, 0, 0])
            return p

        system = self.system
        system.time_step = 0.01

        kT = 1.0
        gamma0 = 2.0
        gammav = 0.04
        p_ext = 2.0
        if self.barostat == "Andersen":
            piston = 0.01
        else:
            piston = 4.0
        vel2force = system.time_step / 2

        # No seed should throw exception
        if self.no_seed:
            with self.assertRaises(ValueError):
                system.thermostat.set_npt(kT=kT, gamma0=gamma0, gammav=gammav)

        # Negative seed should throw exception
        with self.assertRaises(ValueError):
            system.thermostat.set_npt(
                kT=kT, gamma0=gamma0, gammav=gammav, seed=-1)

        system.thermostat.set_npt(kT=kT, gamma0=gamma0, gammav=gammav, seed=41)
        system.integrator.set_isotropic_npt(
            ext_pressure=p_ext, piston=piston, barostat=self.barostat)

        # run(0) does not increase the philox counter and should give the same
        # force and box volume
        p = reset_particle_and_box()
        system.integrator.run(0, recalc_forces=True)
        force0 = np.copy(p.v) / vel2force
        system.integrator.run(0, recalc_forces=True)
        force1 = np.copy(p.v) / vel2force
        boxl1 = np.copy(system.box_l)
        np.testing.assert_almost_equal(force0, force1)
        np.testing.assert_almost_equal(boxl1, [4, 4, 4])

        # run(1) should give a different force and box volume
        p = reset_particle_and_box()
        system.integrator.run(1)
        force2 = np.copy(p.v) / vel2force
        boxl2 = np.copy(system.box_l)
        self.assertTrue(np.all(np.not_equal(force1, force2)))
        self.assertTrue(np.all(np.not_equal(boxl2, [4, 4, 4])))

        # Same seed should not give the same force and box volume with a
        # different counter state
        # force2: npt_iso.rng_counter() = 1, npt_iso.rng_seed() = 41
        # force3: npt_iso.rng_counter() = 2, npt_iso.rng_seed() = 41
        p = reset_particle_and_box()
        system.thermostat.set_npt(kT=kT, gamma0=gamma0, gammav=gammav, seed=41)
        system.integrator.run(1)
        force3 = np.copy(p.v) / vel2force
        boxl3 = np.copy(system.box_l)
        self.assertTrue(np.all(np.not_equal(force2, force3)))
        self.assertTrue(np.all(np.not_equal(boxl2, boxl3)))

        # Seed offset should not give the same force and box volume with a lag
        # force4: npt_iso.rng_counter() = 3, npt_iso.rng_seed() = 42
        # force5: npt_iso.rng_counter() = 4, npt_iso.rng_seed() = 41
        p = reset_particle_and_box()
        system.thermostat.set_npt(kT=kT, gamma0=gamma0, gammav=gammav, seed=42)
        system.integrator.run(1)
        force4 = np.copy(p.v) / vel2force
        boxl4 = np.copy(system.box_l)
        p = reset_particle_and_box()
        system.thermostat.set_npt(kT=kT, gamma0=gamma0, gammav=gammav, seed=41)
        system.integrator.run(1)
        force5 = np.copy(p.v) / vel2force
        boxl5 = np.copy(system.box_l)
        self.assertTrue(np.all(np.not_equal(force4, force5)))
        self.assertTrue(np.all(np.not_equal(boxl4, boxl5)))

    @utx.skipIfMissingFeatures("WCA")
    def test_02__direction(self):
        """Test for NpT constrained in one direction."""

        data = np.genfromtxt(tests_common.data_path("npt_lj_system.data"))
        ref_box_l = 1.01 * np.max(data[:, 0:3])

        system = self.system
        system.box_l = 3 * [ref_box_l]
        system.non_bonded_inter[2, 2].wca.set_params(epsilon=1., sigma=1.)
        system.time_step = 0.01
        if self.barostat == "Andersen":
            gamma0 = 2.0
            gammav = 0.004
            piston = 0.0001
            ext_pressure = 2.0
        else:
            gamma0 = 1.0
            gammav = 0.004
            piston = 4.0
            ext_pressure = 2.0

        for n in range(3):
            direction = np.roll([True, False, False], n)
            system.box_l = 3 * [ref_box_l]
            system.part.add(pos=data[:, 0:3], type=len(data) * [2])
            system.part.all().pos = data[:, 0:3]
            system.part.all().v = data[:, 3:6]
            system.thermostat.set_npt(
                kT=1.0, gamma0=gamma0, gammav=gammav, seed=42)
            system.integrator.set_isotropic_npt(ext_pressure=ext_pressure, piston=piston,
                                                direction=direction, barostat=self.barostat)
            system.integrator.run(20)
            box_l_rel = np.copy(system.box_l) / ref_box_l
            box_l_rel_ref = np.roll([np.max(box_l_rel), 1., 1.], n)
            np.testing.assert_allclose(box_l_rel, box_l_rel_ref, atol=1e-10)
            self.assertGreater(np.max(box_l_rel), 2)
            system.part.clear()

    def test_07__virtual(self):
        system = self.system
        self.system.time_step = 0.01
        if self.barostat == "Andersen":
            piston = 0.01
        else:
            piston = 4.0

        Propagation = espressomd.propagation.Propagation
        virtual = system.part.add(pos=[0, 0, 0], v=[1, 0, 0],
                                  propagation=Propagation.NONE)
        physical = system.part.add(pos=[0, 0, 0], v=[1, 0, 0])

        system.thermostat.set_npt(kT=1.0, gamma0=2.0, gammav=0.04, seed=42)
        system.integrator.set_isotropic_npt(
            ext_pressure=2.0, piston=piston, barostart=self.barostat)

        system.integrator.run(1)

        np.testing.assert_almost_equal(np.copy(virtual.v), [1, 0, 0])
        self.assertTrue(np.all(np.not_equal(np.copy(physical.v), [1, 0, 0])))

    def test_integrator_exceptions(self):
        system = self.system

        # invalid parameters should throw exceptions
        with self.assertRaises(Exception):
            system.integrator.set_isotropic_npt(
                ext_pressure=-1., piston=1., barostat=self.barostat)
        with self.assertRaises(Exception):
            system.integrator.set_isotropic_npt(
                ext_pressure=1., piston=-1., barostat=self.barostat)
        with self.assertRaises(Exception):
            system.integrator.set_isotropic_npt(
                ext_pressure=1., piston=0., barostat=self.barostat)
        with self.assertRaises(Exception):
            system.integrator.set_isotropic_npt(ext_pressure=1., piston=1.,
                                                direction=[0, 0, 0], barostat=self.barostat)

    @utx.skipIfMissingFeatures(["WCA", "P3M"])
    def test_pressure_with_p3m(self):
        """Test for NpT with P3M."""

        data = np.genfromtxt(tests_common.data_path("npt_lj_system.data"))
        ref_box_l = np.max(data[:, 0:3])
        p_ext = 1.0

        system = self.system
        system.box_l = 3 * [ref_box_l]
        system.time_step = 0.01
        system.non_bonded_inter[2, 2].wca.set_params(epsilon=1., sigma=1.)
        system.part.add(pos=data[:, 0:3], v=data[:, 3:6], type=len(data) * [2],
                        q=np.sign(np.arange(100) - 50 + 0.5))
        system.integrator.set_vv()
        system.electrostatics.solver = espressomd.electrostatics.P3M(
            prefactor=2.0, accuracy=1e-2, mesh=3 * [18], cao=5, tune=True)

        if self.barostat == "Andersen":
            system.thermostat.set_npt(kT=1.0, gamma0=0.2, gammav=0.01, seed=42)
            system.integrator.set_isotropic_npt(
                ext_pressure=p_ext, piston=0.0001)
        else:
            system.thermostat.set_npt(
                kT=1.0, gamma0=0.5, gammav=0.001, seed=42)
            system.integrator.set_isotropic_npt(
                ext_pressure=p_ext, piston=4.0, barostat=self.barostat)

        steps = int(0.1 / system.time_step)

        for _ in range(100):
            system.integrator.run(steps)
            p_sim = system.analysis.pressure()['total']
            p_kin = system.analysis.pressure()['kinetic']
            # virial of electrostatic force from system.analysis
            p_vir = p_sim - p_kin
            # virial of electrostatic force from instantaneous_pressure
            p_inst_vir = system.analysis.get_instantaneous_pressure_virial()

            np.testing.assert_allclose(p_vir, p_inst_vir, rtol=1e-2, atol=1e-7)


@utx.skipIfMissingFeatures("NPT")
class NPTThermostat_Andersen(NPTThermostat, ut.TestCase):
    barostat = "Andersen"
    no_seed = True


@utx.skipIfMissingFeatures("NPT")
class NPTThermostat_MTK(NPTThermostat, ut.TestCase):
    barostat = "MTK"
    no_seed = False


if __name__ == "__main__":
    ut.main()
