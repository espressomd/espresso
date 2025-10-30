#
# Copyright (C) 2023 The ESPResSo project
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

import espressomd
import espressomd.observables
import numpy as np
import unittest as ut
import unittest_decorators as utx
import espressomd.polymer
import tests_common
import espressomd.propagation
Propagation = espressomd.propagation.Propagation
from espressomd.observables import MagneticDipoleMoment


def generate_random_unit_vectors(N_PART):
    z = np.random.uniform(-1, 1, N_PART)
    r = np.sqrt(1 - z*z)
    phi = np.random.uniform(0, 2*np.pi, N_PART)
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return np.column_stack((x, y, z))

@utx.skipIfMissingFeatures(["THERMAL_STONER_WOHLFARTH"])
class Test(ut.TestCase):
    """
    Check the total dipole field for a magnetic LJ fluid (500 particles,
    density approx 0.002, mu^2=1, no PBC).
    """
    res_dict_fluid={3.4283694213261087:0.9, 1.1427898071087026:0.6, 0.28569745177717565:0.2}
    res_dict_solid={3.4283694213261087:0.8, 1.1427898071087026:0.54, 0.28569745177717565:0.2}
    skin = 0.4
    time_step = 0.001
    temperature = 1
    kT_KVm_inv=5
    dt_incr= 0.001*3.437060795580368e-08
    HK_inv = 0.17501031139401407
    dip_reduced = 1.7501031139401464
    gamma_T = 74.86576383782938
    gamma_R = 24.955254612609792
    tau0_inv=735412234.8230474
    SNAPSHOT_SEPARATION = 12477
    n_part = 100
    system = espressomd.System(box_l=(29.69314567, 29.69314567, 29.69314567))
    system.cell_system.skin = 0.4
    system.time_step = 0.001
    system.periodicity = [True, True, True]
    orientor_list = generate_random_unit_vectors(N_PART=n_part)
    dip_mom_list = dip_reduced*orientor_list
    positions = espressomd.polymer.linear_polymer_positions(n_polymers=n_part, beads_per_chain=1, min_distance=1., bond_length=1., seed=42)
    positions = np.reshape(positions, (-1, 3))

    def tearDown(self):
        self.system.part.clear()

    def setUp(self):
        system = self.system
        particles = system.part.add(pos=self.positions, director=self.orientor_list)
        particles.sw_real = True
        particles.rotation = (True, True, True)
        particles.kT_KVm_inv = self.kT_KVm_inv
        particles.dt_incr = self.dt_incr
        particles.tau0_inv = self.tau0_inv
        for p1,dipm_el in zip(list(particles),self.dip_mom_list):
            p2=system.part.add(
            pos=p1.pos, dip=dipm_el, rotation=[False, False, False],
            sw_virt=True, Hkinv=self.HK_inv,
            sat_mag=self.dip_reduced)
            p2.vs_auto_relate_to(p1)
            p2.propagation = Propagation.TRANS_VS_RELATIVE | Propagation.ROT_VS_INDEPENDENT

    def test_tSW_fluid(self):
        system = self.system
        for h_reduced,res in self.res_dict_fluid.items():
            print(h_reduced)
            for x in system.constraints:
                system.constraints.remove(x)
            ExtH = espressomd.constraints.HomogeneousMagneticField(H=(0, 0, h_reduced))
            system.constraints.add(ExtH)
            dipm_tot = MagneticDipoleMoment(ids=system.part.select(lambda p: p.sw_virt == True).id)
            norm = 1/(self.dip_reduced*self.n_part)
            system.integrator.run(self.SNAPSHOT_SEPARATION)
            mag_el = dipm_tot.calculate()*norm
            self.assertAlmostEqual(mag_el[-1],res,delta=0.05)
    
    def test_tSW_solid(self):
        system = self.system
        part_slice=system.part.select(lambda p: p.sw_real== True)
        part_slice.rotation = [False, False, False]
        part_slice.fix = [True, True, True]

        for h_reduced,res in self.res_dict_solid.items():
            print(h_reduced)
            for x in system.constraints:
                system.constraints.remove(x)
            ExtH = espressomd.constraints.HomogeneousMagneticField(H=(0, 0, h_reduced))
            system.constraints.add(ExtH)
            dipm_tot = MagneticDipoleMoment(ids=system.part.select(lambda p: p.sw_virt == True).id)
            norm = 1/(self.dip_reduced*self.n_part)
            system.integrator.run(self.SNAPSHOT_SEPARATION)
            mag_el = dipm_tot.calculate()*norm
            self.assertAlmostEqual(mag_el[-1],res,delta=0.05) 
        
if __name__ == "__main__":
    ut.main()
