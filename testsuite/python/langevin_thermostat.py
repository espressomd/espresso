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
import unittest as ut
import espressomd
import numpy as np
from espressomd.interactions import FeneBond


class LangevinThermostat(ut.TestCase):
    """Tests the velocity distribution created by the Langevin thermostat agaisnt 
       the single component Maxwell distribution."""

    s = espressomd.System()

    def single_component_maxwell(self,x1,x2,kT):
        """Integrate the probability density from x1 to x2 using the trapez rule"""
        x=np.linspace(x1,x2,1000)
        return np.trapz(np.exp(-x**2/(2.*kT)),x)/np.sqrt(2.*np.pi*kT)
    
    def check_velocity_distribution(self,vel,minmax,n_bins,error_tol,kT):
        """check the recorded particle distributions in vel againsta histogram with n_bins bins. Drop velocities outside minmax. Check individual histogram bins up to an accuracy of error_tol agaisnt the analytical result for kT."""
        for i in range(3):
            hist=np.histogram(vel[:,i],range=(-minmax,minmax),bins=n_bins,normed=False)
            data=hist[0]/float(vel.shape[0])
            bins=hist[1]
            for j in range(n_bins):
                found=data[j]
                expected=self.single_component_maxwell(bins[j],bins[j+1],kT)
                self.assertLessEqual(abs(found-expected),error_tol)
            



    
    
    def test_aa_verify_single_component_maxwell(self):
        """Verifies the normalization of the analytical expression."""
        self.assertLessEqual(abs(self.single_component_maxwell(-10,10,4.)-1.),1E-4)   
    
    def test_global_langevin(self):
        """Test for global Langevin parameters."""
        N=200
        s=self.s
        s.part.clear()
        s.cell_system.skin=0
        s.cell_system.set_domain_decomposition(use_verlet_lists=False)
        s.min_global_cut=0.2
        s.time_step=0.1
        s.part.add(pos=np.random.random((N,3)))
        kT=2.3
        gamma=1.5
        s.thermostat.set_langevin(kT=kT,gamma=gamma)
        s.integrator.run(50)
        loops=4000
        v_stored=np.zeros((N*loops,3))
        for i in range(loops):
            s.integrator.run(5)
            v_stored[i*N:(i+1)*N,:]=s.part[:].v
        v_minmax=5
        bins=5
        error_tol=0.01
        self.check_velocity_distribution(v_stored,v_minmax,bins,error_tol,kT)


    @ut.skipIf(not espressomd.has_features("LANGEVIN_PER_PARTICLE"),"Test requires LANGEVIN_PER_PARTICLE")
    def test_langevin_per_particle(self):
        """Test for Langevin particle. Covers all combinations of
           particle specific gamma and temp set or not set.
        """
        N=200
        s=self.s
        s.part.clear()
        s.cell_system.skin=0
        s.cell_system.set_domain_decomposition(use_verlet_lists=False)
        s.min_global_cut=0.2
        s.time_step=0.1
        s.part.add(pos=np.random.random((N,3)))
        kT=2.3
        gamma=1.5
        gamma2=2.3
        kT2=1.5
        s.thermostat.set_langevin(kT=kT,gamma=gamma)
        # Set different kT on 2nd half of particles
        s.part[N/2:].temp =kT2
        # Set different gamma on half of the partiles (overlap over both kTs)
        if espressomd.has_features("PARTICLE_ANISOTROPY"):
            s.part[N/4:3*N/4].gamma=gamma2,gamma2,gamma2
        else:    
            s.part[N/4:3*N/4].gamma=gamma2


        s.integrator.run(50)
        loops=4000
        
        v_kT=np.zeros((N/2*loops,3))
        v_kT2=np.zeros((N/2*loops,3))

        for i in range(loops):
            s.integrator.run(5)
            v_kT[i*N/2:(i+1)*N/2,:]=s.part[:N/2].v
            v_kT2[i*N/2:(i+1)*N/2,:]=s.part[N/2:].v
        v_minmax=5
        bins=5
        error_tol=0.014
        self.check_velocity_distribution(v_kT,v_minmax,bins,error_tol,kT)
        self.check_velocity_distribution(v_kT2,v_minmax,bins,error_tol,kT2)





            

            


        
if __name__ == "__main__":
    #print("Features: ", espressomd.features())
    ut.main()
