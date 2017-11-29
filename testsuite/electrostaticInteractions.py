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
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
from espressomd import electrostatics
from tests_common import *
import matplotlib.pyplot as plt



@ut.skipIf(not espressomd.has_features(["ELECTROSTATICS", "EXTERNAL_FORCES"]),
           "Features not available, skipping test!")
class ElectrostaticInteractionsTests(ut.TestCase):
    # Handle to espresso system
    system = espressomd.System()
    system.seed = system.cell_system.get_state()['n_nodes'] * [1234]

    def setUp(self):
        self.system.box_l = [20, 20, 20]
        self.system.time_step = 0.01

        if not self.system.part.exists(0):
            self.system.part.add(id=0, pos=(1.0, 2.0, 2.0), q=1, fix=[1, 1, 1])
        if not self.system.part.exists(1):
            self.system.part.add(
                id=1, pos=(3.0, 2.0, 2.0), q=-1, fix=[1, 1, 1])
        print("ut.TestCase setUp")
    def calc_cdh_potential(self, r, cdf_params):
        kT=1.0
        q1=self.system.part[0].q
        q2=self.system.part[1].q
        u = np.zeros_like(r)

        # r>r_cut
        i = np.where(r>cdf_params['r_cut'])[0]
        u[i]=0

        # r_cut>r>r1
        i = np.where( (cdf_params['r_cut']>r) & (r>cdf_params['r1']))
        u[i]=cdf_params['bjerrum_length']*kT*q1*q2*np.exp(-cdf_params['kappa']*r[i])/(cdf_params['eps_ext']*r[i])

        # r0<r<r1
        i = np.where( (cdf_params['r0']<r) & (r<cdf_params['r1']))
        u[i]=cdf_params['bjerrum_length']*kT*q1*q2*np.exp(-cdf_params['alpha']*(r[i]-cdf_params['r0']))/(cdf_params['eps_int']*r[i])

        # r<r0
        i = np.where(r<cdf_params['r0'])[0]
        u[i]=cdf_params['bjerrum_length']*kT*q1*q2/(cdf_params['eps_int']*r[i])
        
        return u
       

    @ut.skipIf(not espressomd.has_features(["P3M"]),
               "Features not available, skipping test!")
    def test_p3m(self):
        self.system.part[0].pos = [1.0, 2.0, 2.0]
        self.system.part[1].pos =[3.0, 2.0, 2.0]
        # results,
        p3m_energy = -0.501062398379
        p3m_force = 2.48921612e-01
        test_P3M = generate_test_for_class(
            self.system,
            electrostatics.P3M,
            dict(
                accuracy=9.910945054074526e-08,
                 mesh=[22, 22, 22],
                 cao=7,
                 r_cut=8.906249999999998,
                 alpha=0.387611049779351,
                 tune=False))
        p3m = espressomd.electrostatics.P3M(bjerrum_length=1.0,
                                            accuracy=9.910945054074526e-08,
                                            mesh=[22, 22, 22],
                                            cao=7,
                                            r_cut=8.906249999999998,
                                            alpha=0.387611049779351,
                                            tune=False)
        self.system.actors.add(p3m)
        self.assertAlmostEqual(self.system.analysis.energy()['coulomb'],
                               p3m_energy)
        # need to update forces
        self.system.integrator.run(0)
        self.assertTrue(np.allclose(self.system.part[0].f,
                                    [p3m_force, 0, 0]))
        self.assertTrue(np.allclose(self.system.part[1].f,
                                    [-p3m_force, 0, 0]))
        self.system.actors.remove(p3m)

    @ut.skipIf(not espressomd.has_features(["COULOMB_DEBYE_HUECKEL"]),
               "Features not available, skipping test!")
    def test_dh(self):
        cdh_params = dict(bjerrum_length=1.0,
                          kappa=2.3,
                          r_cut=5.0,
                          r0=1.0,
                          r1=1.9,
                          eps_int=0.8,
                          eps_ext=1.0,
                          alpha=2.0)
        test_CDH = generate_test_for_class(
            self.system,
            electrostatics.CDH,
            cdh_params)
        cdh = espressomd.electrostatics.CDH(
            bjerrum_length=cdh_params['bjerrum_length'],
                                            kappa=cdh_params['kappa'],
                                            r_cut=cdh_params['r_cut'],
                                            r0=cdh_params['r0'],
                                            r1=cdh_params['r1'],
                                            eps_int=cdh_params['eps_int'],
                                            eps_ext=cdh_params['eps_ext'],
                                            alpha=cdh_params['alpha'])
        self.system.actors.add(cdh)
        print("yo, CDH:")
        dr=0.001
        r = np.arange(.5, 1.01*cdh_params['r_cut'], dr)
        u_cdh = self.calc_cdh_potential(r, cdh_params)
        f_cdh = -np.gradient(u_cdh, dr)
        #print(u_cdh)
        #print(self.system.analysis.energy()['coulomb'])
        # self.assertAlmostEqual(self.system.analysis.energy()['coulomb'],
                               # self.p3m_energy)

        u_cdh_core=np.zeros_like(r)
        f_cdh_core=np.zeros_like(r)
        # need to update forces
        for i,ri in enumerate(r):
            self.system.part[1].pos=self.system.part[0].pos+[ri, 0, 0]
            self.system.integrator.run(0)
#            print(self.system.part[0].f[0], f_cdh[i], ri)
            u_cdh_core[i]=self.system.analysis.energy()['coulomb']
            f_cdh_core[i]=self.system.part[0].f[0]
            #print (i, r[i], u_cdh_core[i], u_cdh[i], f_cdh_core[i],  f_cdh[i])
            #self.assertAlmostEqual(f_cdh_core[i], f_cdh[i], delta=1e-2)

        #self.assertTrue( np.allclose(u_cdh_core,
        #                             u_cdh))
        
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(r, u_cdh_core, label='core')
        ax.plot(r, u_cdh, ls=':', label='python')
        ax.legend(loc='best')
        ax.set_xlabel(r'distance $r$')
        ax.set_ylabel(r'potential $U_\mathrm{CDH}$')
        #ax.axis(ymin=0)
        #ax.axis(xmin=0)
        fig.savefig('cdh_u.png')
        plt.close(fig)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(r, f_cdh_core, label='core')
        ax.plot(r, f_cdh, ls=':', label='python')
        ax.legend(loc='best')
        ax.set_xlabel(r'distance $r$')
        ax.set_ylabel(r'force $f_\mathrm{CDH}$')
        ax.axis(ymin=-10)
        ax.axis(ymax=1)
        #ax.axis(xmin=0)
        fig.savefig('cdh_f.png')
        plt.close(fig)

        # self.assertTrue( np.allclose(self.system.part[0].f,
                                     #[self.p3m_force, 0, 0]))
        # self.assertTrue( np.allclose(self.system.part[1].f,
                                     #[-self.p3m_force, 0, 0]))
        self.system.actors.remove(cdh)


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
