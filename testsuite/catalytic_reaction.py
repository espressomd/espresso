from espressomd import has_features, System
from espressomd.reaction import Reaction

import unittest as ut
import numpy as np

@ut.skipIf(not has_features(["ROTATION","CATALYTIC_REACTIONS"]),
           "Features missing")
class ReactionTest(ut.TestCase):

    def test_reaction(self):
        system = System()
        system.box_l = [10,10,10]
        system.cell_system.skin = 0.1
        system.time_step = 0.01

        system.part.add(pos=0.5*system.box_l,type=0)
        reactant = system.part.add(pos=0.5*system.box_l-[0,0,1],type=1)
        product = system.part.add(pos=0.5*system.box_l+[0,0,1],type=2)

        # Set up the reaction, rate is infinity to always trigger a reaction
        Reaction(reactant_type=1,catalyzer_type=0,product_type=2,ct_range=3.0,ct_rate=np.inf)

        system.integrator.run(1)
        
        # Check that particles have switched type
        self.assertEqual(reactant.type,2)
        self.assertEqual(product.type,1)
        
if __name__ == '__main__':
    ut.main()

