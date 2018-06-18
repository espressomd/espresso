from __future__ import print_function
import sys
import unittest as ut
import numpy as np
import espressomd

class AnalyzeGyration(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
    np.random.seed(1234)
    type_mono=0

    @classmethod
    def setUpClass(self):
        box_l = 20.0
        self.system.box_l = np.array([box_l, box_l, box_l])
        self.system.cell_system.set_n_square(use_verlet_lists=False)
        #4x4 cube
        for x,y,z  in np.ndindex((4,4,4)):
            self.system.part.add(pos=[x,y,z], type=self.type_mono)
        # long stick in z, force z as principal axis
        for x,y,z  in np.ndindex((1,1,10)):
           self.system.part.add(pos=[x+1.5, y+1.5, z+4], type=self.type_mono)
        # two small nubs in x, force y as secondary axis
        self.system.part.add(pos=[1.5, 4, 1.5], type=self.type_mono)
        self.system.part.add(pos=[1.5, -1, 1.5], type=self.type_mono)
    def test_gyration_tensor(self):
        # get results
        res=self.system.analysis.gyration_tensor(p_type=self.type_mono)
        rg=self.system.analysis.calc_rg(chain_start=0, number_of_chains=1, chain_length=len(self.system.part[:]))[0]
        #test if principal and secondary  axis is [0,0,1] and [0,1,0]
        self.assertTrue( np.allclose(np.abs(res['eva0'][1]), [0.,0.,1.], atol=1e-6))
        self.assertTrue( np.allclose(np.abs(res['eva1'][1]), [0.,1.,0.], atol=1e-6))
        self.assertTrue( np.allclose(np.abs(res['eva2'][1]), [1.,0.,0.], atol=1e-6))
        self.assertTrue( np.allclose(rg**2, res['Rg^2'], atol=1e-6))


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
