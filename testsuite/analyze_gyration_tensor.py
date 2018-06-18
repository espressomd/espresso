from __future__ import print_function
import sys
import unittest as ut
import numpy as np
import espressomd

class AnalyzeGyration(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
    np.random.seed(1234)
    cube_len=4
    type_cube=0
    type_stick=1

    @classmethod
    def setUpClass(self):
        box_l = 20.0
        cube_centre=0.5*(self.cube_len-1)
        self.system.box_l = np.array([box_l, box_l, box_l])
        self.system.cell_system.set_n_square(use_verlet_lists=False)
        #4x4 cube
        for x,y,z  in np.ndindex((self.cube_len, self.cube_len, self.cube_len)):
            self.system.part.add(pos=[x,y,z], type=self.type_cube)
        # long stick in z, force z as principal axis
        for x,y,z  in np.ndindex((1,1,10)):
           self.system.part.add(pos=[x+cube_centre, y+cube_centre, z+self.cube_len], type=self.type_stick)
        # two small nubs in y, force y as secondary axis
        self.system.part.add(pos=[cube_centre, self.cube_len, cube_centre], type=self.type_stick)
        self.system.part.add(pos=[cube_centre, -1, cube_centre], type=self.type_stick)

    def test_gyration_tensor_cube(self):
        # get results
        res=self.system.analysis.gyration_tensor(p_type=self.type_cube)
        rg=self.system.analysis.calc_rg(chain_start=0, number_of_chains=1, chain_length=self.cube_len**3)[0]
        # make sure all eigenvalues (for the cube) are identical
        self.assertTrue( np.allclose(np.abs(res['eva0'][0]), np.abs(res['eva1'][0]), np.abs(res['eva2'][0]), atol=1e-6))
        self.assertTrue( np.allclose(rg**2, res['Rg^2'], atol=1e-6))

    def test_gyration_tensor(self):
        # get results
        res=self.system.analysis.gyration_tensor(p_type=[self.type_stick, self.type_cube])
        rg=self.system.analysis.calc_rg(chain_start=0, number_of_chains=1, chain_length=len(self.system.part[:]))[0]
        #test if principal and secondary  axis is [0,0,1] and [0,1,0]
        self.assertTrue( np.allclose(np.abs(res['eva0'][1]), [0.,0.,1.], atol=1e-6))
        self.assertTrue( np.allclose(np.abs(res['eva1'][1]), [0.,1.,0.], atol=1e-6))
        self.assertTrue( np.allclose(np.abs(res['eva2'][1]), [1.,0.,0.], atol=1e-6))
        self.assertTrue( np.allclose(rg**2, res['Rg^2'], atol=1e-6))

    def test_mom_intertia(self):
        
        sqr_dist = np.sum((self.system.analysis.center_of_mass(p_type=0)-self.system.part.select(type=0).pos)**2, axis=0)
        mom_I = self.system.analysis.moment_of_inertia_matrix(p_type=0)
        # the cube case should have zero as off- diagonal components
        self.assertTrue( np.allclose([mom_I[0,1], mom_I[0,2], mom_I[1,2], mom_I[1,0], mom_I[2,0], mom_I[2,1]] , np.zeros(6), atol=1e-6))
        self.assertTrue( np.allclose([mom_I[0,0], mom_I[1,1], mom_I[2,2]], [sqr_dist[1]+sqr_dist[2], sqr_dist[0]+sqr_dist[2], sqr_dist[1]+sqr_dist[2]], atol=1e-6))

        
        

if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
