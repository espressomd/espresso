#
# Copyright (C) 2013-2022 The ESPResSo project
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
import tempfile

# Number of particles
npart = 10
# Number of different bond types
nbonds = 100


@dataclasses.dataclass
class MockParticle:
    id: int
    type: int
    pos: np.ndarray
    v: np.ndarray
    bonds: list


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

    @classmethod
    def setUpClass(cls):
        cls.temp_dir = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        cls.temp_dir.cleanup()

    def add_particles(self):
        """Sets up a system from test_mock_particles and prepares environment
        for the tests."""
        for p in self.test_mock_particles:
            self.system.part.add(id=p.id, type=p.type, pos=p.pos, v=p.v)
            for b in p.bonds:
                self.system.part.by_id(p.id).add_bond(b)

    def tearDown(self):
        self.system.part.clear()

    def generate_prefix(self, test_id):
        return os.path.join(self.temp_dir.name, test_id.rsplit('.')[-1])

    def build_list_of_expected_files(self, prefix, **fields):
        exts = {'head', 'pref', 'id'}
        if fields.get('types', False):
            exts.add('type')
        if fields.get('positions', False):
            exts.add('pos')
        if fields.get('velocities', False):
            exts.add('vel')
        if fields.get('bonds', False):
            exts.add('boff')
            exts.add('bond')
        return {f'{prefix}.{ext}' for ext in exts}

    def check_files_exist(self, prefix, **fields):
        """Checks if all necessary files have been written."""
        filepaths = self.build_list_of_expected_files(prefix, **fields)
        for filepath in filepaths:
            self.assertTrue(os.path.isfile(filepath))

    def check_sample_system(self, **fields):
        """Checks particles in the system against the reference values."""
        self.assertEqual(len(self.system.part), len(self.test_mock_particles))
        for p, q in zip(self.system.part, self.test_mock_particles):
            ref_t = q.type if fields.get('types', False) else 0
            ref_p = q.pos if fields.get('positions', False) else [0., 0., 0.]
            ref_v = q.v if fields.get('velocities', False) else [0., 0., 0.]
            ref_b = q.bonds if fields.get('bonds', False) else []
            self.assertEqual(p.id, q.id)
            self.assertEqual(p.type, ref_t)
            np.testing.assert_array_equal(np.copy(p.pos), ref_p)
            np.testing.assert_array_equal(np.copy(p.v), ref_v)
            self.assertEqual(len(p.bonds), len(ref_b))
            # Check all bonds
            for bp, bq in zip(p.bonds, ref_b):
                # Bond type - "bend" stores the index of the bond
                self.assertEqual(bp[0].params["bend"], bq[0])
                # Bond partners
                np.testing.assert_array_equal(bp[1:], bq[1:])

    def test_mpiio(self):
        fields = {
            'types': True,
            'positions': True,
            'velocities': True,
            'bonds': True}
        prefix = self.generate_prefix(self.id())
        mpiio = espressomd.io.mpiio.Mpiio(system=self.system)

        self.add_particles()
        mpiio.write(prefix, **fields)
        self.check_files_exist(prefix, **fields)

        self.system.part.clear()
        mpiio.read(prefix, **fields)
        self.check_sample_system(**fields)

    def test_mpiio_without_positions(self):
        prefix = self.generate_prefix(self.id())
        mpiio = espressomd.io.mpiio.Mpiio(system=self.system)
        self.add_particles()
        mpiio.write(prefix, types=True, positions=False)
        self.system.part.clear()
        mpiio.read(prefix, types=True, positions=False)
        self.check_sample_system(types=True, positions=False)

    def test_mpiio_without_types(self):
        prefix = self.generate_prefix(self.id())
        mpiio = espressomd.io.mpiio.Mpiio(system=self.system)
        self.add_particles()
        mpiio.write(prefix, types=False, positions=True)
        self.system.part.clear()
        mpiio.read(prefix, types=False, positions=True)
        self.check_sample_system(types=False, positions=True)

    def test_mpiio_multiple_instances(self):
        fields1 = {
            'types': True,
            'positions': True,
            'velocities': True,
            'bonds': True}
        fields2 = {
            'types': True,
            'positions': True,
            'velocities': False,
            'bonds': False}
        prefix1 = self.generate_prefix(self.id()) + '.1'
        prefix2 = self.generate_prefix(self.id()) + '.2'
        mpiio1 = espressomd.io.mpiio.Mpiio(system=self.system)
        mpiio2 = espressomd.io.mpiio.Mpiio(system=self.system)

        self.add_particles()
        mpiio1.write(prefix1, **fields1)
        mpiio2.write(prefix2, **fields2)
        self.check_files_exist(prefix1, **fields1)
        self.check_files_exist(prefix2, **fields2)

        self.system.part.clear()
        mpiio1.read(prefix1, **fields1)
        self.check_sample_system(**fields1)

        self.system.part.clear()
        mpiio2.read(prefix2, **fields2)
        self.check_sample_system(**fields2)

    def test_mpiio_exceptions(self):
        mpiio = espressomd.io.mpiio.Mpiio(system=self.system)
        prefix = self.generate_prefix(self.id())
        msg_prefix = "Need to supply output prefix via the 'prefix' argument."
        with self.assertRaisesRegex(ValueError, msg_prefix):
            mpiio.write(None, positions=True)
        with self.assertRaisesRegex(ValueError, msg_prefix):
            mpiio.read(None, positions=True)
        with self.assertRaisesRegex(ValueError, "No output fields chosen."):
            mpiio.write(prefix)
        with self.assertRaisesRegex(ValueError, "No output fields chosen."):
            mpiio.read(prefix)


if __name__ == '__main__':
    ut.main()
