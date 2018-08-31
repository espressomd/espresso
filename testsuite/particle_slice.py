from __future__ import print_function
import unittest as ut
import espressomd
from espressomd import has_features
import numpy as np


class ParticleSliceTest(ut.TestCase):

    state = [[0, 0, 0], [0, 0, 1]]
    system = espressomd.System(box_l=[10, 10, 10])

    def __init__(self, *args, **kwargs):
        super(ParticleSliceTest, self).__init__(*args, **kwargs)
        self.system.part.clear()
        self.system.part.add(pos=[0, 0, 0])
        self.system.part.add(pos=[0, 0, 1])
        self.system.part.add(pos=[0, 0, 2])
        self.system.part.add(pos=[0, 0, 3])

        if has_features(["EXTERNAL_FORCES"]):
            self.system.part[1].fix = self.state[1]
            self.assertTrue(np.array_equal(
                self.system.part[0].fix, self.state[0]))
            self.assertTrue(np.array_equal(
                self.system.part[1].fix, self.state[1]))
            self.assertTrue(np.array_equal(
                self.system.part[:2].fix, self.state))
        xs = self.system.part[:].pos
        for i in range(len(xs)):
            self.assertTrue(np.array_equal(xs[i], self.system.part[i].pos))

    @ut.skipIf(
        not has_features(
            ["EXTERNAL_FORCES"]),
        "Features not available, skipping test!")
    def test_1_set_different_values(self):
        self.state[0] = [1, 0, 0]
        self.state[1] = [1, 0, 0]
        self.system.part[:2].fix = self.state
        self.assertTrue(np.array_equal(self.system.part[:2].fix, self.state))

    @ut.skipIf(
        not has_features(
            ["EXTERNAL_FORCES"]),
        "Features not available, skipping test!")
    def test_2_set_same_value(self):
        self.state[0] = [0, 1, 0]
        self.state[1] = [0, 1, 0]
        self.system.part[:2].fix = self.state[1]
        self.assertTrue(np.array_equal(self.system.part[:2].fix, self.state))

    @ut.skipIf(
        not has_features(
            ["EXTERNAL_FORCES"]),
        "Features not available, skipping test!")
    def test_3_set_one_value(self):
        self.state[1] = [0, 0, 1]
        self.system.part[1:2].fix = self.state[1]
        self.assertTrue(np.array_equal(self.system.part[:2].fix, self.state))

    @ut.skipIf(
        not has_features(
            ["EXTERNAL_FORCES"]),
        "Features not available, skipping test!")
    def test_4_str(self):
        self.assertEqual(
            repr(self.system.part[0].fix), repr(np.array([0, 1, 0])))
        self.assertEqual(repr(self.system.part[:2].fix), repr(
            np.array([[0, 1, 0], [0, 0, 1]])))

    def test_pos_str(self):
        self.system.part[0].pos = [0,0,0]
        self.system.part[1].pos = [0,0,1]
        self.assertEqual(repr(self.system.part[0].pos), repr(
                np.array([0.0, 0.0, 0.0])))
        self.assertEqual(repr(self.system.part[:2].pos), repr(
                np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])))

    @ut.skipIf(
        not has_features(
            ["ELECTROSTATICS"]),
        "Features not available, skipping test!")
    def test_scalar(self):
        self.system.part[:1].q = 1.3
        self.assertEqual(self.system.part[0].q, 1.3)
        self.system.part[:2].q = 2.0
        self.assertEqual(self.system.part[0].q, 2)
        self.assertEqual(self.system.part[1].q, 2)
        self.system.part[:2].q = 3
        self.assertEqual(self.system.part[0].q, 3)
        self.assertEqual(self.system.part[1].q, 3)
        self.system.part[:2].q = [-1, 1.0]
        self.assertEqual(self.system.part[0].q, -1)
        self.assertEqual(self.system.part[1].q, 1)
        qs = self.system.part[:2].q
        self.assertEqual(qs[0], -1)
        self.assertEqual(qs[1], 1)
   
    def test_bonds(self):

        fene = espressomd.interactions.FeneBond(k=1, d_r_max=1, r_0=1)
        self.system.bonded_inter.add(fene)

        # Setter

        # tuple
        b = fene,0
        self.system.part[2].bonds = b
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),), () ] )

        # list
        self.system.part[:].bonds = []
        b = [fene,0]
        self.system.part[2].bonds = b
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),), () ] )

        # nested list single
        self.system.part[:].bonds = []
        b = [[fene,0]]
        self.system.part[2].bonds = b
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),), () ] )

        # nested list multi
        self.system.part[:].bonds = []
        b = [[fene,0], [fene,1]]
        self.system.part[2].bonds = b
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),(fene,1)), () ] )

        # nested tuple single
        self.system.part[:].bonds = []
        b = ((fene,0),)
        self.system.part[2].bonds = b
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),), () ] )

        # nested tuple multi
        self.system.part[:].bonds = []
        b = ((fene,0), (fene,1))
        self.system.part[2].bonds = b
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),(fene,1)), () ] )

        #Add/Del bonds
        self.system.part[:].bonds = []
        self.system.part[2].add_bond((fene,0))
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),), () ] )
        self.system.part[2].add_bond((fene,1))
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),(fene,1)), () ] )
        self.system.part[2].delete_bond((fene,1))
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),), () ] )
        self.system.part[2].delete_bond((fene,0))
        self.assertTrue( self.system.part[:].bonds == [ (), (), (), () ] )
        
        self.system.part[:].bonds = []
        self.system.part[2].add_bond((fene,0))
        self.system.part[2].add_bond((fene,1))
        self.system.part[2].delete_all_bonds()
        self.assertTrue( self.system.part[:].bonds == [ (), (), (), () ] )


        # Slices

        # tuple for all
        self.system.part[:].bonds = []
        b = fene,0
        self.system.part[2:].bonds = b
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),), ((fene,0),) ] )

        # list for all
        self.system.part[:].bonds = []
        b = [fene,0]
        self.system.part[2:].bonds = b
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),), ((fene,0),) ] )

        # nested list single for all
        self.system.part[:].bonds = []
        b = [[fene,0]]
        self.system.part[2:].bonds = b
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),), ((fene,0),) ] )

        # nested list multi for all
        self.system.part[:].bonds = []
        b = [[fene,0], [fene,1]]
        self.system.part[2:].bonds = b
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),(fene,1)), ((fene,0),(fene,1)) ] )

        # tuples for each
        self.system.part[:].bonds = []
        b = (((fene,0),),((fene,1),))
        self.system.part[2:].bonds = b
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),), ((fene,1),) ] )

        # lists for each
        self.system.part[:].bonds = []
        b = [[[fene,0]],[[fene,1]]]
        self.system.part[2:].bonds = b
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),), ((fene,1),) ] )

        # multi tuples for each
        self.system.part[:].bonds = []
        b = (((fene,0),(fene,1)),((fene,0),(fene,1)))
        self.system.part[2:].bonds = b
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),(fene,1)), ((fene,0),(fene,1)) ] )

        # multi lists for each
        self.system.part[:].bonds = []
        b = [[[fene,0],[fene,1]],[[fene,0],[fene,1]]]
        self.system.part[2:].bonds = b
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),(fene,1)), ((fene,0),(fene,1)) ] )

        #Add/Del bonds

        self.system.part[:].bonds = []
        self.system.part[2:].add_bond((fene,0))
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),), ((fene,0),) ] )
        self.system.part[2:].add_bond((fene,1))
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),(fene,1)), ((fene,0),(fene,1)) ] )
        self.system.part[2:].delete_bond((fene,1))
        self.assertTrue( self.system.part[:].bonds == [ (), (), ((fene,0),), ((fene,0),) ] )
        self.system.part[2:].delete_bond((fene,0))
        self.assertTrue( self.system.part[:].bonds == [ (), (), (), () ] )
        
        b = [[[fene,0],[fene,1]],[[fene,0],[fene,1]]]
        self.system.part[2:].bonds = b
        self.system.part[:].delete_all_bonds()
        self.assertTrue( self.system.part[:].bonds == [ (), (), (), () ] )
 

    def cmp_array_like(self, A,B):
        return all(a.tolist() == b for a,b in zip(A,B))

    @ut.skipIf(
        not has_features(
            ["EXCLUSIONS"]),
        "Features not available, skipping test!")
    def test_exclusions(self):

        # Setter
        # int
        self.system.part[:].exclusions = []
        b = 1
        self.system.part[2].exclusions = b
        self.assertTrue( self.cmp_array_like(self.system.part[:].exclusions, [ [], [2], [1], [] ] ) )

        # single list
        self.system.part[:].exclusions = []
        b = [1]
        self.system.part[2].exclusions = b
        self.assertTrue( self.cmp_array_like(self.system.part[:].exclusions, [ [], [2], [1], [] ] ) )

        # tuple
        self.system.part[:].exclusions = []
        b = (0,1)
        self.system.part[2].exclusions = b
        self.assertTrue( self.cmp_array_like(self.system.part[:].exclusions, [ [2], [2], [0,1], [] ] ) )

        # list
        self.system.part[:].exclusions = []
        b = [0,1]
        self.system.part[2].exclusions = b
        self.assertTrue( self.cmp_array_like(self.system.part[:].exclusions, [ [2], [2], [0,1], [] ] ) )

        # Add/Del exclusions 
        self.system.part[:].exclusions = []
        self.system.part[2].add_exclusion(1)
        self.assertTrue( self.cmp_array_like(self.system.part[:].exclusions, [ [], [2], [1], [] ] ) )
        self.system.part[2].add_exclusion(0)
        self.assertTrue( self.cmp_array_like(self.system.part[:].exclusions, [ [2], [2], [1,0], [] ] ) )
        self.system.part[2].delete_exclusion(0)
        self.assertTrue( self.cmp_array_like(self.system.part[:].exclusions, [ [], [2], [1], [] ] ) )
        self.system.part[2].delete_exclusion(1)
        self.assertTrue( self.cmp_array_like(self.system.part[:].exclusions, [ [], [], [], [] ] ) )

        # Slices

        # single list for all
        self.system.part[:].exclusions = []
        b = [1]
        self.system.part[2:].exclusions = b
        self.assertTrue( self.cmp_array_like(self.system.part[:].exclusions, [ [], [2,3], [1], [1] ] ) )

        # list for all
        self.system.part[:].exclusions = []
        b = [0,1]
        self.system.part[2:].exclusions = b
        self.assertTrue( self.cmp_array_like(self.system.part[:].exclusions, [ [2,3], [2,3], [0,1], [0,1] ] ) )

        # single list for each
        self.system.part[:].exclusions = []
        b = [[0],[0]]
        self.system.part[2:].exclusions = b
        self.assertTrue( self.cmp_array_like(self.system.part[:].exclusions, [ [2,3], [], [0], [0] ] ) )

        # multi list for each
        self.system.part[:].exclusions = []
        b = [[0,1],[0,1]]
        self.system.part[2:].exclusions = b
        self.assertTrue( self.cmp_array_like(self.system.part[:].exclusions, [ [2,3], [2,3], [0,1], [0,1] ] ) )

        #Add/Del exclusions
        self.system.part[:].exclusions = []
        self.system.part[2:].add_exclusion(1)
        self.assertTrue( self.cmp_array_like(self.system.part[:].exclusions, [ [], [2,3], [1], [1] ] ) )
        self.system.part[2:].add_exclusion(0)
        self.assertTrue( self.cmp_array_like(self.system.part[:].exclusions, [ [2,3], [2,3], [1,0], [1,0] ] ) )
        self.system.part[2:].delete_exclusion(0)
        self.assertTrue( self.cmp_array_like(self.system.part[:].exclusions, [ [], [2,3], [1], [1] ] ) )
        self.system.part[2:].delete_exclusion(1)
        self.assertTrue( self.cmp_array_like(self.system.part[:].exclusions, [ [], [], [], [] ] ) )

    @ut.skipIf(not espressomd.has_features("VIRTUAL_SITES_RELATIVE"),
        "Features not available, skipping test!")
    def test_vs_relative(self):
        
        self.system.part.clear()
        self.system.part.add(pos=[0,0,0])
        self.system.part.add(pos=[0,0,0])
        self.system.part.add(pos=[0,0,0])
        self.system.part.add(pos=[0,0,0])
        self.system.part[0].vs_relative = [1,1.0,(1.0,1.0,1.0,1.0)]

        self.assertTrue( repr(self.system.part[:].vs_relative) == repr([(1, 1.0, np.array([ 1.,  1.,  1.,  1.])), (0, 0.0, np.array([ 0.,  0.,  0.,  0.])), (0, 0.0, np.array([ 0.,  0.,  0.,  0.])), (0, 0.0, np.array([ 0.,  0.,  0.,  0.]))]) )

        self.system.part[:].vs_relative = [1,1.0,(1.0,1.0,1.0,1.0)]
        
        self.assertTrue( repr(self.system.part[:].vs_relative) == repr([(1, 1.0, np.array([ 1.,  1.,  1.,  1.])), (1, 1.0, np.array([ 1.,  1.,  1.,  1.])), (1, 1.0, np.array([ 1.,  1.,  1.,  1.])), (1, 1.0, np.array([ 1.,  1.,  1.,  1.]))]) )

        self.system.part[:].vs_relative = [ [1,1.0,(1.0,1.0,1.0,1.0)], [1,2.0,(1.0,1.0,1.0,1.0)], [1,3.0,(1.0,1.0,1.0,1.0)], [1,4.0,(1.0,1.0,1.0,1.0)] ]

        self.assertTrue( repr(self.system.part[:].vs_relative) == repr([(1, 1.0, np.array([ 1.,  1.,  1.,  1.])), (1, 2.0, np.array([ 1.,  1.,  1.,  1.])), (1, 3.0, np.array([ 1.,  1.,  1.,  1.])), (1, 4.0, np.array([ 1.,  1.,  1.,  1.]))]) )

    def test_multiadd(self):
        self.system.part.clear()
        positions = [[1,2,3],[4,5,6],[7,8,9]]
        self.system.part.add(pos=positions)
        for p in positions:
            self.system.part.add(pos=p)
        np.testing.assert_allclose(np.copy(self.system.part[:3].pos), np.copy(self.system.part[3:6].pos))
        
        with self.assertRaises(ValueError):
            self.system.part.add(pos=([1,1,1],[2,2,2]),type = 0)
        with self.assertRaises(ValueError):
            self.system.part.add(pos=([1,1,1],[2,2,2]),type = (0,1,2))

        self.system.part.clear()
        self.system.part.add(pos=([1,1,1],[2,2,2]),type = (0,1))
        self.assertEqual(self.system.part[0].type,0)
        self.assertEqual(self.system.part[1].type,1)

    def test_empty(self):
        self.assertTrue(np.array_equal(self.system.part[0:0].pos, np.empty(0)))

    def test_len(self):
        self.assertEqual(len(self.system.part[0:0]), 0)
        self.assertEqual(len(self.system.part[0:1]), 1)
        self.assertEqual(len(self.system.part[0:2]), 2)

    def test_non_existing_property(self):
        with self.assertRaises(AttributeError):
            self.system.part[:].thispropertydoesnotexist = 1.0

if __name__ == "__main__":
    ut.main()
