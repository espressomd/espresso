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
import numpy as np
import espressomd
import tests_common


@utx.skipIfMissingFeatures(["NPT", "LENNARD_JONES"])
class IntegratorNPT:

    """
    Compare the pressure and compressibility of a WCA fluid
    against expected values. For p_ext = 1.0 and kT = 1.0 in LJ units,
    compressibility should be around 0.5.
    It can be estimated by the Carnahan-Starling equation of state
    (J. Chem. Phys. 51, 635-636 (1969))
    because this equation describes well the pressure of WCA fluid
    at kT=1.0 (J. Chem. Phys. 124, 164507 (2006)).
    """
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    def setUp(self):
        self.system.box_l = [5] * 3
        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.25

    def tearDown(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.integrator.set_vv()

    def test_compressibility_and_pressure(self):
        system = self.system
        system.box_l = [5.86326165] * 3

        data = np.genfromtxt(tests_common.data_path("npt_lj_system.data"))
        p_ext = 1.0

        system.part.add(pos=data[:, :3], v=data[:, 3:])
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1, sigma=1, cutoff=1.12246, shift=0.25)

        if self.barostat == "Andersen":
            system.thermostat.set_npt(kT=1.0, gamma0=1, gammav=0.004, seed=42)
            system.integrator.set_isotropic_npt(
                ext_pressure=p_ext, piston=0.0001)
        else:
            system.thermostat.set_npt(
                kT=1.0, gamma0=0.5, gammav=0.001, seed=42)
            system.integrator.set_isotropic_npt(
                ext_pressure=p_ext, piston=4.0, barostat=self.barostat)

        system.integrator.run(800)
        # averaged pressure by system.analysis.pressure()
        avp = 0
        # averaged p_{kin}*V by system.analysis.pressure(),
        # where p_{kin} stands for kinetic pressure and V is system volume
        # p_{kin}*V should be same as N*kT=100
        avpV_sim = 0
        # averaged virial pressure by system.analysis.pressure()
        avp_sim_vir = 0
        # averaged p_{kin}*V by get_instantaneous_pressure()
        # it also should be N*kT = 100
        avpV_inst = 0
        # averaged virial pressure by get_instantaneous_pressure_virial()
        avp_inst_vir = 0
        # large sample size needed due to precision loss in fast-math builds
        n = 55000
        skip_p = 8
        ls = np.zeros(n)
        for t in range(n):
            system.integrator.run(2)
            if t % skip_p == 0:
                volume = float(np.prod(system.box_l))
                pressure = system.analysis.pressure()
                avp += pressure['total']
                avpV_sim += pressure['kinetic'] * volume
                avp_sim_vir += pressure['non_bonded']

                p_inst_vir = system.analysis.get_instantaneous_pressure_virial()
                avp_inst_vir += p_inst_vir
                p_inst = system.analysis.get_instantaneous_pressure()
                p_inst_kin = p_inst - p_inst_vir
                avpV_inst += p_inst_kin * volume
            ls[t] = system.box_l[0]

        avp /= (n / skip_p)
        avpV_sim /= (n / skip_p)
        avp_sim_vir /= (n / skip_p)
        avpV_inst /= (n / skip_p)
        avp_inst_vir /= (n / skip_p)
        Vs = np.array(ls)**3
        compressibility = np.var(Vs) / np.average(Vs)

        self.assertAlmostEqual(avp, p_ext, delta=0.02)
        self.assertAlmostEqual(compressibility, 0.5, delta=0.05)
        np.testing.assert_allclose(avp_sim_vir, avp_inst_vir, atol=1e-10)
        self.assertAlmostEqual(avpV_sim, 100., delta=1.)
        self.assertAlmostEqual(avpV_inst, 100., delta=1.)

    def test_negative_volume(self):
        """Test for NpT with bad parameters."""

        data = np.genfromtxt(tests_common.data_path("npt_lj_system.data"))
        ref_box_l = np.max(data[:, 0:3])

        system = self.system
        system.box_l = 3 * [ref_box_l]
        dt = 0.01
        system.time_step = dt
        if self.barostat == "Andersen":
            piston = 0.0001
        else:
            piston = 4.0

        direction = [True] * 3
        ext_pressure = 100.0  # Too large external pressure
        system.box_l = 3 * [ref_box_l]
        system.part.add(pos=data[:, 0:3], type=len(data) * [2])
        system.part.all().pos = data[:, 0:3]
        system.part.all().v = data[:, 3:6]
        self.system.integrator.set_vv()

        system.thermostat.set_npt(kT=1.0, gamma0=0.1, gammav=0.001, seed=42)
        system.integrator.set_isotropic_npt(ext_pressure=ext_pressure,
                                            piston=piston,
                                            direction=direction,
                                            barostat=self.barostat)

        if self.barostat == "Andersen":
            with self.assertRaises(Exception):
                system.integrator.run(10)
        elif self.barostat == "MTK":
            with self.assertRaises(Exception):
                system.integrator.run(10)
            # Volume cannot be negative within NPT ensemble based on MTK equation
            self.assertTrue(float(np.prod(system.box_l)) > 0.)


@utx.skipIfMissingFeatures("NPT")
class IntegratorNPT_Andersen(IntegratorNPT, ut.TestCase):
    barostat = "Andersen"


@utx.skipIfMissingFeatures("NPT")
class IntegratorNPT_MTK(IntegratorNPT, ut.TestCase):
    barostat = "MTK"


if __name__ == "__main__":
    ut.main()
