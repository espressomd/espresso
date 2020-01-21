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
Testmodule for the H5MD interface.
"""
import os
import sys
import unittest as ut
import unittest_decorators as utx
import numpy as np
import espressomd
import h5py  # h5py has to be imported *after* espressomd (MPI)
from espressomd.interactions import Virtual

npart = 26


class CommonTests(ut.TestCase):

    """
    Class that holds common test methods.
    """
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    # avoid particles to be set outside of the main box, otherwise particle
    # positions are folded in the core when writing out and we cannot directly
    # compare positions in the dataset and where particles were set. One would
    # need to unfold the positions of the hdf5 file.
    box_l = npart / 2.0
    system.box_l = [box_l, box_l, box_l]
    system.cell_system.skin = 0.4
    system.time_step = 0.01

    for i in range(npart):
        system.part.add(id=i, pos=np.array(3 * [i], dtype=float),
                        v=np.array([1.0, 2.0, 3.0]), type=23)
        if espressomd.has_features(['MASS']):
            system.part[i].mass = 2.3
        if espressomd.has_features(['EXTERNAL_FORCES']):
            system.part[i].ext_force = [0.1, 0.2, 0.3]

    vb = Virtual()
    system.bonded_inter.add(vb)

    for i in range(npart - 1):
        system.part[i].add_bond((vb, i + 1))

    system.integrator.run(steps=0)

    @classmethod
    def setUpClass(cls):
        if os.path.isfile('test.h5'):
            os.remove('test.h5')
        cls.py_file = cls.py_pos = cls.py_vel = cls.py_f = cls.py_id = cls.py_img = None

    def test_metadata(self):
        """Test if the H5MD metadata has been written properly."""
        self.assertEqual(self.py_file['h5md'].attrs['version'][0], 1)
        self.assertEqual(self.py_file['h5md'].attrs['version'][1], 1)
        self.assertIn('creator', self.py_file['h5md'])
        self.assertIn('name', self.py_file['h5md/creator'].attrs)
        self.assertIn('version', self.py_file['h5md/creator'].attrs)
        self.assertEqual(
            self.py_file['h5md/creator'].attrs['name'][:], b'ESPResSo')
        self.assertIn('author', self.py_file['h5md'])
        self.assertIn('name', self.py_file['h5md/author'].attrs)

    def test_pos(self):
        """Test if positions have been written properly."""
        np.testing.assert_allclose(
            np.array([3 * [float(i) % self.box_l] for i in range(npart)]),
            np.array([x for (_, x) in sorted(zip(self.py_id, self.py_pos))]))

    def test_img(self):
        """Test if images have been written properly."""
        images = np.append(np.zeros((int(npart / 2), 3)),
                           np.ones((int(npart / 2), 3)))
        images = images.reshape(npart, 3)
        np.testing.assert_allclose(
            [x for (_, x) in sorted(zip(self.py_id, self.py_img))], images)

    def test_vel(self):
        """Test if velocities have been written properly."""
        np.testing.assert_allclose(
            np.array([[1.0, 2.0, 3.0] for _ in range(npart)]),
            np.array([x for (_, x) in sorted(zip(self.py_id, self.py_vel))]),
            err_msg="Velocities not written correctly by H5md!")

    @utx.skipIfMissingFeatures(['EXTERNAL_FORCES'])
    def test_f(self):
        """Test if forces have been written properly."""
        np.testing.assert_allclose(
            np.array([[0.1, 0.2, 0.3] for _ in range(npart)]),
            np.array([x for (_, x) in sorted(zip(self.py_id, self.py_f))]),
            err_msg="Forces not written correctly by H5md!")

    def test_bonds(self):
        """Test if bonds have been written properly."""
        self.assertEqual(len(self.py_bonds), npart - 1)

        for i in range(npart - 1):
            bond = [x for x in self.py_bonds if x[0] == i][0]
            self.assertEqual(bond[0], i + 0)
            self.assertEqual(bond[1], i + 1)


@utx.skipIfMissingFeatures(['H5MD'])
class H5mdTestOrdered(CommonTests):

    """
    Test the core implementation of writing hdf5 files if written ordered.
    """

    @classmethod
    def setUpClass(cls):
        write_ordered = True
        from espressomd.io.writer import h5md
        h5 = h5md.H5md(
            filename="test.h5",
            write_pos=True,
            write_vel=True,
            write_force=True,
            write_species=True,
            write_mass=True,
            write_ordered=write_ordered)
        h5.write()
        h5.flush()
        h5.close()
        cls.py_file = h5py.File("test.h5", 'r')
        cls.py_pos = cls.py_file['particles/atoms/position/value'][0]
        cls.py_img = cls.py_file['particles/atoms/image/value'][0]
        cls.py_vel = cls.py_file['particles/atoms/velocity/value'][0]
        cls.py_f = cls.py_file['particles/atoms/force/value'][0]
        cls.py_id = cls.py_file['particles/atoms/id/value'][0]
        cls.py_bonds = cls.py_file['connectivity/atoms']

    @classmethod
    def tearDownClass(cls):
        os.remove("test.h5")

    def test_ids(self):
        """Test if ids have been written properly."""
        np.testing.assert_allclose(np.array(range(npart)), self.py_id,
                                   err_msg="ids incorrectly ordered and written by H5md!")


@utx.skipIfMissingFeatures(['H5MD'])
class H5mdTestUnordered(CommonTests):

    """
    Test the core implementation of writing hdf5 files if written un-ordered.
    """

    @classmethod
    def setUpClass(cls):
        write_ordered = False
        from espressomd.io.writer import h5md
        h5 = h5md.H5md(
            filename="test.h5",
            write_pos=True,
            write_vel=True,
            write_force=True,
            write_species=True,
            write_mass=True,
            write_ordered=write_ordered)
        h5.write()
        h5.flush()
        h5.close()
        cls.py_file = h5py.File("test.h5", 'r')
        cls.py_pos = cls.py_file['particles/atoms/position/value'][0]
        cls.py_img = cls.py_file['particles/atoms/image/value'][0]
        cls.py_vel = cls.py_file['particles/atoms/velocity/value'][0]
        cls.py_f = cls.py_file['particles/atoms/force/value'][0]
        cls.py_id = cls.py_file['particles/atoms/id/value'][0]
        cls.py_bonds = cls.py_file['connectivity/atoms']

    @classmethod
    def tearDownClass(cls):
        os.remove("test.h5")


if __name__ == "__main__":
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(H5mdTestUnordered))
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(H5mdTestOrdered))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
