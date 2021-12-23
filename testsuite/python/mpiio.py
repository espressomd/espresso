#
# Copyright (C) 2013-2019 The ESPResSo project
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

"""
Testmodule for MPI-IO.
"""

import espressomd
import espressomd.io
import espressomd.interactions
import numpy as np
import unittest as ut
import random
import os
import dataclasses

# Number of particles
npart = 10
# Number of different bond types
nbonds = 100

filename = "testdata.mpiio"
exts = ["head", "pref", "id", "type", "pos", "vel", "boff", "bond"]
filenames = [filename + "." + ext for ext in exts]


@dataclasses.dataclass
class MockParticle:
    id: int
    type: int
    pos: np.ndarray
    v: np.ndarray
    bonds: list


def clean_files():
    for f in filenames:
        if os.path.isfile(f):
            os.remove(f)


def randint_different_from(a, b, n):
    """Returns a random integer in [a, b) that is not n."""
    r = n
    while r == n:
        r = random.randint(a, b - 1)
    return r


def get_random_mock_particles():
    """Returns a list of random particle descriptions."""
    parts = []
    for i in range(npart):
        p = MockParticle(id=i,
                         type=random.randint(0, 100),
                         pos=np.random.rand(3),
                         v=np.random.rand(3),
                         bonds=[])
        # Up to 20 bonds; otherwise this test will take ages
        for _ in range(random.randint(0, 20)):
            btype = random.randint(0, nbonds - 1)
            # Don't create loops, i.e. exclude "i" itself
            p1 = randint_different_from(0, npart, i)
            p2 = randint_different_from(0, npart, i)
            # Don't add the same bond twice
            if (btype, p1, p2) not in p.bonds:
                p.bonds.append((btype, p1, p2))
        parts.append(p)
    return parts


class MPIIOTest(ut.TestCase):

    """
    Test class for the MPI-IO core functionality.
    Generates random particles, dumps them, reads them in,
    again and then checks the input against the initially created random
    particles.
    """
    system = espressomd.system.System(box_l=[1, 1, 1])
    # Just a bunch of random interactions such that add_bond does not throw
    for i in range(nbonds):
        system.bonded_inter[i] = espressomd.interactions.AngleHarmonic(
            bend=i, phi0=i)
    test_mock_particles = get_random_mock_particles()

    def setUp(self):
        """Sets up a system from test_mock_particles and prepares environment
        for the tests."""
        clean_files()  # Prior call might not have completed successfully
        for p in self.test_mock_particles:
            self.system.part.add(id=p.id, type=p.type, pos=p.pos, v=p.v)
            for b in p.bonds:
                self.system.part.by_id(p.id).add_bond(b)

    def tearDown(self):
        clean_files()

    def check_files_exist(self):
        """Checks if all necessary files have been written."""
        for fn in filenames:
            self.assertTrue(os.path.isfile(fn))

    def check_sample_system(self):
        """Checks the particles in the ESPResSo system "self.s" against the
        true values in "self.test_particles"."""
        for p, q in zip(self.system.part, self.test_mock_particles):
            self.assertEqual(p.id, q.id)
            self.assertEqual(p.type, q.type)
            np.testing.assert_array_equal(np.copy(p.pos), q.pos)
            np.testing.assert_array_equal(np.copy(p.v), q.v)
            self.assertEqual(len(p.bonds), len(q.bonds))
            # Check all bonds
            for bp, bq in zip(p.bonds, q.bonds):
                # Bond type - "bend" stores the index of the bond
                self.assertEqual(bp[0].params["bend"], bq[0])
                # Bond partners
                np.testing.assert_array_equal(bp[1:], bq[1:])

    def test_mpiio(self):
        espressomd.io.mpiio.mpiio.write(
            filename, types=True, positions=True, velocities=True, bonds=True)

        self.check_files_exist()

        self.system.part.clear()  # Clear to be on the safe side
        espressomd.io.mpiio.mpiio.read(
            filename, types=True, positions=True, velocities=True, bonds=True)

        self.check_sample_system()


if __name__ == '__main__':
    ut.main()
