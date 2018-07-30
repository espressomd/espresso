

#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
# Tests particle property setters/getters
from __future__ import print_function
import os
import unittest as ut
import numpy as np
import espressomd
import espressomd.electrostatics as el
import espressomd.electrostatic_extensions as el_ext
from espressomd import scafacos
import tests_common


@ut.skipIf(not espressomd.has_features(["ELECTROSTATICS"]),
           "Features not available, skipping test!")
class CoulombMixedPeriodicity(ut.TestCase):
    """"Test mixed periodicity electrostatics"""
     
    S = espressomd.System(box_l=[1.0, 1.0, 1.0])
    buf_node_grid = S.cell_system.node_grid
    S.thermostat.turn_off()
    forces = {}
    tolerance_force = 5E-4
    tolerance_energy = 1.8E-3
    generate_data=False

    # Reference energy from MMM2D
    reference_energy = 216.640984711 

    def setUp(self):
        self.S.box_l = (10, 10, 10)
        self.S.time_step = 0.01
        self.S.cell_system.skin = 0.
        while len(self.S.actors)>0:
            del self.S.actors[0]

        #  Clear actors that might be left from prev tests
        self.S.part.clear()
        data = np.genfromtxt(tests_common.abspath(
                "data/coulomb_mixed_periodicity_system.data"))

        # Add particles to system and store reference forces in hash
        # Input format: id pos q f
        for particle in data:
            id = particle[0]
            pos = particle[1:4]
            q = particle[4]
            f = particle[5:]
            self.S.part.add(id=int(id), pos=pos, q=q)
            self.forces[id] = f

    def compare(self, method_name, energy=True):
        # Compare forces and energy now in the system to stored ones

        # Force
        rms_force_diff = 0.
        for p in self.S.part:
            rms_force_diff += np.sum((np.copy(p.f) - self.forces[p.id])**2)
        rms_force_diff=np.sqrt(rms_force_diff/len(self.S.part))

        print(method_name, "rms force difference", rms_force_diff)

        # Energy
        if energy:
            energy_abs_diff = abs(self.S.analysis.energy()["total"] - self.reference_energy)
            print(method_name, "energy", self.S.analysis.energy()["total"])
            print(method_name, "energy difference", energy_abs_diff)
            self.assertLessEqual(
                energy_abs_diff,
                self.tolerance_energy,
                "Absolte energy difference " +
                str(energy_abs_diff) +
                " too large for " +
                method_name)
        self.assertLessEqual(
            rms_force_diff,
            self.tolerance_force,
            "Asbolute force difference " +
            str(rms_force_diff) +
            " too large for method " +
            method_name)

    # Tests for individual methods

    if espressomd.has_features(["P3M"]):
        def test_zz_p3mElc(self):
            # Make sure, the data satisfies the gap 
            for p in self.S.part:
                if p.pos[2]<0 or p.pos[2]>9.:
                    raise Exception("Particle z pos invalid")
            
            self.S.cell_system.set_domain_decomposition()
            self.S.cell_system.node_grid=sorted(self.S.cell_system.node_grid,key=lambda x:-x)
            self.S.periodicity=1,1,1
            self.S.box_l = (10, 10, 10)

            p3m=el.P3M(prefactor=1, accuracy = 1e-6,mesh=(64,64,64))
            
            self.S.actors.add(p3m)
            elc=el_ext.ELC(maxPWerror=1E-6,gap_size=1)
            self.S.actors.add(elc)
            self.S.integrator.run(0)
            self.compare("elc", energy=True)
            self.S.actors.remove(p3m)

    def test_MMM2D(self):
        self.S.box_l = (10, 10, 10)
        self.S.cell_system.set_layered(n_layers=10,use_verlet_lists=False)
        self.S.periodicity=1,1,0
        mmm2d=(el.MMM2D(prefactor=1,maxPWerror=1E-7))
            
        self.S.actors.add(mmm2d)
        self.S.integrator.run(0)
        if self.generate_data:
            n=len(self.S.part)
            data=np.hstack((self.S.part[:].id.reshape((n,1)),self.S.part[:].pos_folded,self.S.part[:].q.reshape((n,1)),self.S.part[:].f))
            np.savetxt(tests_common.abspath(
                "data/coulomb_mixed_periodicity_system.data"),data)
        self.compare("mmm2d (compared to stored data)", energy=True)
        self.S.actors.remove(mmm2d)

    if espressomd.has_features("SCAFACOS"):
        if "p2nfft" in scafacos.available_methods():
            def test_scafacos_p2nfft(self):
                self.S.periodicity=1,1,0
                self.S.cell_system.set_domain_decomposition()
                self.S.box_l=10,10,10
                
                scafacos=el.Scafacos(
                        prefactor=1,
                        method_name="p2nfft",
                        method_params={
                            "tolerance_field": 5E-5,
                            "pnfft_n":"96,96,128",
                            "pnfft_N":"96,96,128",
                            "r_cut":2.4,
                            "pnfft_m":3})
                self.S.actors.add(scafacos)
                self.S.integrator.run(0)
                self.compare("scafacos_p2nfft", energy=True)
                self.S.actors.remove(scafacos)


if __name__ == "__main__":
    ut.main()
