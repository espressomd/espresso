from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np

class ParticleSlice(ut.TestCase):

	state = [[0,0,0], [0,0,1]]
	system = espressomd.System()

	system.part.add(pos=[0,0,0])
	system.part.add(pos=[0,0,1], fix=state[1])
	assert np.array_equal(system.part[0].fix, state[0])
	assert np.array_equal(system.part[1].fix, state[1])
	assert np.array_equal(system.part[:].fix, state)

	def test_1_set_different_values(self):
		self.state[0] = [1,0,0]
		self.state[1] = [1,0,0]
		self.system.part[:].fix = self.state
		assert np.array_equal(self.system.part[:].fix, self.state)

	def test_2_set_same_value(self):
		self.state[0] = [0,1,0]
		self.state[1] = [0,1,0]
		self.system.part[:].fix = self.state[1]
		assert np.array_equal(self.system.part[:].fix, self.state)

	def test_3_set_one_value(self):
		self.state[1] = [0,0,1]
		self.system.part[1:].fix = self.state[1]
		assert np.array_equal(self.system.part[:].fix, self.state)

	def test_4_str(self):
		self.assertEqual( repr(self.system.part[0].pos), repr(np.array([0.0,0.0,0.0])) )
		self.assertEqual( repr(self.system.part[:].pos), repr(np.array([[0.0,0.0,0.0], [0.0,0.0,1.0]])) )
		self.assertEqual( repr(self.system.part[0].fix), repr(np.array([0,1,0])) )
		self.assertEqual( repr(self.system.part[:].fix), repr([np.array([0,1,0]), np.array([0,0,1])]) )


if __name__ == "__main__":
    ut.main()
