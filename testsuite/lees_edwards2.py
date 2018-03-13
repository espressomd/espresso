
#!/usr/bin/env python


from __future__ import print_function

import espressomd
from espressomd.interactions import HarmonicBond
from espressomd.virtual_sites import VirtualSitesRelative
import unittest as ut
import numpy as np


@ut.skipIf(not espressomd.has_features(['LEES_EDWARDS']),
  'Feature not available, skipping test!')

class LeesEdwardsTest(ut.TestCase):

    system = espressomd.System(box_l=[1,1,1])
    bond = HarmonicBond(k=1,r_0=0)
    system.bonded_inter.add(bond)
    system.time_step=0.01
    system.cell_system.skin =0.15
    system.min_global_cut =0.15

    system.virtual_sites=VirtualSitesRelative()



    def test_distance_calc(self):
        print("normal particles")
        s=self.system
        bond=self.bond

        s.part.clear()
        s.part.add(id=0,pos=(0.5,0.9,0.5))
        s.part.add(id=1,pos=(0.5,0.1,0.5),bonds=(bond,0))
        s.integrator.run(0)
        np.testing.assert_allclose(s.part[0].f,[0,0.2,0])

        s.lees_edwards_offset=0.1
        s.integrator.run(0)
        np.testing.assert_allclose(s.part[0].f,[0.1,0.2,0])
        print(s.part[1].pos,s.part[1].pos_folded,s.part[1].image_box)
    
    def test_vs(self):
        print("VS")
        s=self.system
        bond=self.bond

        s.lees_edwards_offset=0
        s.part.clear()
        s.part.add(id=0,pos=(0.5,0.9,0.5))
        s.part.add(id=1,pos=(0.5,0.8,0.5))
        s.part.add(id=2,pos=(0.5,0.9,0.5),bonds=(bond,0))
        s.part[2].vs_auto_relate_to(1)
        s.integrator.run(0)
        np.testing.assert_allclose(s.part[0].f,[0,0.0,0])
        
        s.part[1].pos=0.5,0.95,0.5
        s.integrator.run(0)
        np.testing.assert_allclose(s.part[0].f,[0.0,0.15,0])


        s.lees_edwards_offset=0.1
        s.integrator.run(0)
        print(s.part[2].pos,s.part[2].bare_position,s.part[2].image_box)
        np.testing.assert_allclose(s.part[0].f,[0.1,0.15,0])



    
if __name__ == "__main__":
  ut.main()
