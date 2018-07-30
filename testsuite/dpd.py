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
from time import time
from itertools import product

@ut.skipIf(not espressomd.has_features("DPD"),"Skipped because feature is disabled")
class DPDThermostat(ut.TestCase):
    """Tests the velocity distribution created by the dpd thermostat against 
       the single component Maxwell distribution."""

    s = espressomd.System(box_l=[1.0, 1.0, 1.0])
    s.box_l = 3 * [10]
    s.time_step = 0.01
    s.cell_system.skin=0.4

    def setUp(self):
        self.s.seed = range(self.s.cell_system.get_state()["n_nodes"])
        np.random.seed(16)

    def tearDown(self):
        s = self.s
        s.part.clear()

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

    def check_total_zero(self):
        v_total = np.sum(self.s.part[:].v, axis=0)
        self.assertTrue(v_total[0] < 1e-11)
        self.assertTrue(v_total[1] < 1e-11)
        self.assertTrue(v_total[2] < 1e-11)

    def test_single(self):
        """Test velocity distribution of a dpd fluid with a single type."""
        N=200
        s=self.s
        s.part.add(pos=s.box_l * np.random.random((N,3)))
        kT=2.3
        gamma=1.5
        s.thermostat.set_dpd(kT=kT)
        s.non_bonded_inter[0,0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.5,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.5)
        s.integrator.run(100)
        loops=250
        v_stored=np.zeros((N*loops,3))
        for i in range(loops):
            s.integrator.run(10)
            v_stored[i*N:(i+1)*N,:]=s.part[:].v
        v_minmax=5
        bins=5
        error_tol=0.01
        self.check_velocity_distribution(v_stored,v_minmax,bins,error_tol,kT)
        self.check_total_zero()

    def test_binary(self):
        """Test velocity distribution of binary dpd fluid"""
        N=200
        s=self.s
        s.part.add(pos=s.box_l * np.random.random((N//2,3)), type=N//2*[0])
        s.part.add(pos=s.box_l * np.random.random((N//2,3)), type=N//2*[1])
        kT=2.3
        gamma=1.5
        s.thermostat.set_dpd(kT=kT)
        s.non_bonded_inter[0,0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.0,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.0)
        s.non_bonded_inter[1,1].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.0,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.0)
        s.non_bonded_inter[0,1].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.5,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.5)
        s.integrator.run(100)
        loops=300
        v_stored=np.zeros((N*loops,3))
        for i in range(loops):
            s.integrator.run(10)
            v_stored[i*N:(i+1)*N,:]=s.part[:].v
        v_minmax=5
        bins=5
        error_tol=0.01
        self.check_velocity_distribution(v_stored,v_minmax,bins,error_tol,kT)
        self.check_total_zero()

    def test_disable(self):
        N=200
        s=self.s
        s.time_step=0.01
        s.part.add(pos=np.random.random((N,3)))
        kT=2.3
        gamma=1.5
        s.thermostat.set_dpd(kT=kT)
        s.non_bonded_inter[0,0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.5,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.5)
        s.integrator.run(10)

        s.thermostat.turn_off()

        # Reset velocities
        s.part[:].v = [1.,2.,3.]

        s.integrator.run(10)

        # Check that there was neither noise nor friction
        for v in s.part[:].v:
            for i in range(3):
                self.assertTrue(v[i] == float(i+1))

        # Turn back on
        s.thermostat.set_dpd(kT=kT)

        # Reset velocities for faster convergence
        s.part[:].v = [0.,0.,0.]

        # Equilibrate
        s.integrator.run(250)

        loops=250
        v_stored=np.zeros((N*loops,3))
        for i in range(loops):
            s.integrator.run(10)
            v_stored[i*N:(i+1)*N,:]=s.part[:].v
        v_minmax=5
        bins=5
        error_tol=0.01
        self.check_velocity_distribution(v_stored,v_minmax,bins,error_tol,kT)

    def test_const_weight_function(self):
        s=self.s
        kT=0.
        gamma=1.42
        s.thermostat.set_dpd(kT=kT)
        s.non_bonded_inter[0,0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.2,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.4)

        s.part.add(id=0, pos=[5,5,5], type = 0, v=[0,0,0])
        v=[.5,.8,.3]
        s.part.add(id=1, pos=[3,5,5], type = 0, v = v)

        s.integrator.run(0)

        # Outside of both cutoffs, forces should be 0
        for f in s.part[:].f:
            self.assertTrue(f[0] == 0.)
            self.assertTrue(f[1] == 0.)
            self.assertTrue(f[2] == 0.)

        # Only trans
        s.part[1].pos = [5. - 1.3, 5, 5]

        s.integrator.run(0)

        # Only trans, so x component should be zero
        self.assertTrue(s.part[0].f[0] == 0.)
        # f = gamma * v_ij
        self.assertTrue(abs(s.part[0].f[1] - gamma*v[1]) < 1e-11)
        self.assertTrue(abs(s.part[0].f[2] - gamma*v[2]) < 1e-11)
        # Momentum conservation
        self.assertTrue(s.part[1].f[0] == 0.)
        self.assertTrue(abs(s.part[1].f[1] + gamma*v[1]) < 1e-11)
        self.assertTrue(abs(s.part[1].f[2] + gamma*v[2]) < 1e-11)

        # Trans and parallel
        s.part[1].pos = [5. - 1.1, 5, 5]

        s.integrator.run(0)

        self.assertTrue(abs(s.part[0].f[0] - gamma*v[0]) < 1e-11)
        self.assertTrue(abs(s.part[0].f[1] - gamma*v[1]) < 1e-11)
        self.assertTrue(abs(s.part[0].f[2] - gamma*v[2]) < 1e-11)
        self.assertTrue(abs(s.part[1].f[0] + gamma*v[0]) < 1e-11)
        self.assertTrue(abs(s.part[1].f[1] + gamma*v[1]) < 1e-11)
        self.assertTrue(abs(s.part[1].f[2] + gamma*v[2]) < 1e-11)

    def test_linear_weight_function(self):
        s=self.s
        kT=0.
        gamma=1.42
        s.thermostat.set_dpd(kT=kT)
        s.non_bonded_inter[0,0].dpd.set_params(
            weight_function=1, gamma=gamma, r_cut=1.2,
            trans_weight_function=1, trans_gamma=gamma, trans_r_cut=1.4)

        def omega(dist, r_cut):
            return (1. - dist / r_cut)

        s.part.add(id=0, pos=[5,5,5], type = 0, v=[0,0,0])
        v=[.5,.8,.3]
        s.part.add(id=1, pos=[3,5,5], type = 0, v = v)

        s.integrator.run(0)

        # Outside of both cutoffs, forces should be 0
        for f in s.part[:].f:
            self.assertTrue(f[0] == 0.)
            self.assertTrue(f[1] == 0.)
            self.assertTrue(f[2] == 0.)

        # Only trans
        s.part[1].pos = [5. - 1.3, 5, 5]

        s.integrator.run(0)

        # Only trans, so x component should be zero
        self.assertTrue(s.part[0].f[0] == 0.)
        # f = gamma * v_ij
        self.assertTrue(abs(s.part[0].f[1] - omega(1.3, 1.4)**2*gamma*v[1]) < 1e-11)
        self.assertTrue(abs(s.part[0].f[2] - omega(1.3, 1.4)**2*gamma*v[2]) < 1e-11)
        # Momentum conservation
        self.assertTrue(s.part[1].f[0] == 0.)
        self.assertTrue(abs(s.part[1].f[1] + omega(1.3, 1.4)**2*gamma*v[1]) < 1e-11)
        self.assertTrue(abs(s.part[1].f[2] + omega(1.3, 1.4)**2*gamma*v[2]) < 1e-11)

        # Trans and parallel
        s.part[1].pos = [5. - 1.1, 5, 5]

        s.integrator.run(0)

        self.assertTrue(abs(s.part[0].f[0] - omega(1.1, 1.2)**2*gamma*v[0]) < 1e-11)
        self.assertTrue(abs(s.part[0].f[1] - omega(1.1, 1.4)**2*gamma*v[1]) < 1e-11)
        self.assertTrue(abs(s.part[0].f[2] - omega(1.1, 1.4)**2*gamma*v[2]) < 1e-11)
        self.assertTrue(abs(s.part[1].f[0] + omega(1.1, 1.2)**2*gamma*v[0]) < 1e-11)
        self.assertTrue(abs(s.part[1].f[1] + omega(1.1, 1.4)**2*gamma*v[1]) < 1e-11)
        self.assertTrue(abs(s.part[1].f[2] + omega(1.1, 1.4)**2*gamma*v[2]) < 1e-11)

        # Trans and parallel 2nd point
        s.part[1].pos = [5. - 0.5, 5, 5]

        s.integrator.run(0)

        self.assertTrue(abs(s.part[0].f[0] - omega(0.5, 1.2)**2*gamma*v[0]) < 1e-11)
        self.assertTrue(abs(s.part[0].f[1] - omega(0.5, 1.4)**2*gamma*v[1]) < 1e-11)
        self.assertTrue(abs(s.part[0].f[2] - omega(0.5, 1.4)**2*gamma*v[2]) < 1e-11)
        self.assertTrue(abs(s.part[1].f[0] + omega(0.5, 1.2)**2*gamma*v[0]) < 1e-11)
        self.assertTrue(abs(s.part[1].f[1] + omega(0.5, 1.4)**2*gamma*v[1]) < 1e-11)
        self.assertTrue(abs(s.part[1].f[2] + omega(0.5, 1.4)**2*gamma*v[2]) < 1e-11)

    def test_ghosts_have_v(self):
        s=self.s

        s.box_l = 3 * [10.]

        r_cut=1.5
        dx = 0.25 * r_cut

        def f(i):
            if i == 0:
                return dx
            else:
                return 10. - dx

        # Put a particle in every corner
        for ind in product([0, 1], [0, 1], [0, 1]):
            pos = [f(x) for x in ind]
            v = ind
            s.part.add(pos=pos, v=v)

        gamma=1.0
        s.thermostat.set_dpd(kT=0.0)
        s.non_bonded_inter[0,0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=r_cut,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=r_cut)

        s.integrator.run(0)

        id = 0
        for ind in product([0, 1], [0, 1], [0, 1]):
            for i in ind:
                if ind[i] == 0:
                    sgn = 1
                else:
                    sgn = -1

                self.assertAlmostEqual(sgn*4.0, s.part[id].f[i])
            id += 1

    def test_constraint(self):
        import espressomd.shapes

        s = self.s

        s.constraints.add(shape=espressomd.shapes.Wall(dist=0, normal=[1,0,0]), particle_type=0, particle_velocity=[1,2,3])

        s.thermostat.set_dpd(kT=0.0)
        s.non_bonded_inter[0,0].dpd.set_params(
            weight_function=0, gamma=1., r_cut=1.0,
            trans_weight_function=0, trans_gamma=1., trans_r_cut=1.0)

        p = s.part.add(pos=[0.5,0,0], type=0, v=[0,0,0])

        s.integrator.run(0)

        self.assertAlmostEqual(p.f[0], 1.)
        self.assertAlmostEqual(p.f[1], 2.)
        self.assertAlmostEqual(p.f[2], 3.)

        for c in s.constraints:
            s.constraints.remove(c)

if __name__ == "__main__":
    ut.main()
