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
import espressomd.interactions
import espressomd.io.writer
try:
    import h5py  # h5py has to be imported *after* espressomd (MPI)
    skipIfMissingPythonPackage = utx.no_skip
except ImportError:
    skipIfMissingPythonPackage = ut.skip(
        "Python module h5py not available, skipping test!")


N_PART = 26


@utx.skipIfMissingFeatures(['H5MD'])
@skipIfMissingPythonPackage
class H5mdTests(ut.TestCase):
    """
    Test the core implementation of writing hdf5 files.

    """
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    # avoid particles to be set outside of the main box, otherwise particle
    # positions are folded in the core when writing out and we cannot directly
    # compare positions in the dataset and where particles were set. One would
    # need to unfold the positions of the hdf5 file.
    box_l = N_PART / 2.0
    system.box_l = [box_l, box_l, box_l]
    system.cell_system.skin = 0.4
    system.time_step = 0.01

    for i in range(N_PART):
        p = system.part.add(id=i, pos=np.array(3 * [i], dtype=float),
                            v=np.array([1.0, 2.0, 3.0]), type=23)
        if espressomd.has_features(['MASS']):
            p.mass = 2.3
        if espressomd.has_features(['EXTERNAL_FORCES']):
            p.ext_force = [0.1, 0.2, 0.3]
        if espressomd.has_features(['ELECTROSTATICS']):
            p.q = i

    vb = espressomd.interactions.Virtual()
    system.bonded_inter.add(vb)

    for i in range(N_PART - 1):
        system.part.by_id(i).add_bond((vb, i + 1))

    system.integrator.run(steps=0)
    system.time = 12.3

    @classmethod
    def setUpClass(cls):
        if os.path.isfile('test.h5'):
            os.remove('test.h5')
        h5_units = espressomd.io.writer.h5md.UnitSystem(
            time='ps', mass='u', length='m', charge='e')
        h5 = espressomd.io.writer.h5md.H5md(
            file_path="test.h5", unit_system=h5_units)
        h5.write()
        h5.write()
        h5.flush()
        h5.close()
        cls.py_file = h5py.File("test.h5", 'r')
        cls.py_pos = cls.py_file['particles/atoms/position/value'][1]
        cls.py_img = cls.py_file['particles/atoms/image/value'][1]
        cls.py_mass = cls.py_file['particles/atoms/mass/value'][1]
        cls.py_vel = cls.py_file['particles/atoms/velocity/value'][1]
        cls.py_charge = cls.py_file['particles/atoms/charge/value'][1]
        cls.py_f = cls.py_file['particles/atoms/force/value'][1]
        cls.py_id = cls.py_file['particles/atoms/id/value'][1]
        cls.py_id_time = cls.py_file['particles/atoms/id/time'][1]
        cls.py_id_step = cls.py_file['particles/atoms/id/step'][1]
        cls.py_bonds = cls.py_file['connectivity/atoms/value'][1]
        cls.py_box = cls.py_file['particles/atoms/box/edges/value'][1]

    @classmethod
    def tearDownClass(cls):
        os.remove("test.h5")

    def test_opening(self):
        h5 = espressomd.io.writer.h5md.H5md(file_path="test.h5")
        h5.close()

    def test_box(self):
        np.testing.assert_allclose(self.py_box, self.box_l)

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
            np.array([3 * [float(i) % self.box_l] for i in range(N_PART)]),
            np.array([x for (_, x) in sorted(zip(self.py_id, self.py_pos))]))

    def test_time(self):
        """Test for time dataset."""
        self.assertEqual(self.py_id_time, 12.3)

    def test_img(self):
        """Test if images have been written properly."""
        images = np.append(np.zeros((int(N_PART / 2), 3)),
                           np.ones((int(N_PART / 2), 3)))
        images = images.reshape(N_PART, 3)
        np.testing.assert_allclose(
            [x for (_, x) in sorted(zip(self.py_id, self.py_img))], images)

    @utx.skipIfMissingFeatures("MASS")
    def test_mass(self):
        """Test if masses have been written correct."""
        np.testing.assert_allclose(self.py_mass, 2.3)

    @utx.skipIfMissingFeatures(['ELECTROSTATICS'])
    def test_charge(self):
        """Test if charges have been written properly."""
        charges = np.arange(N_PART)
        np.testing.assert_allclose(
            [x for (_, x) in sorted(zip(self.py_id, self.py_charge))], charges)

    def test_vel(self):
        """Test if velocities have been written properly."""
        np.testing.assert_allclose(
            np.array([[1.0, 2.0, 3.0] for _ in range(N_PART)]),
            np.array([x for (_, x) in sorted(zip(self.py_id, self.py_vel))]),
            err_msg="Velocities not written correctly by H5md!")

    @utx.skipIfMissingFeatures(['EXTERNAL_FORCES'])
    def test_f(self):
        """Test if forces have been written properly."""
        np.testing.assert_allclose(
            np.array([[0.1, 0.2, 0.3] for _ in range(N_PART)]),
            np.array([x for (_, x) in sorted(zip(self.py_id, self.py_f))]),
            err_msg="Forces not written correctly by H5md!")

    def test_bonds(self):
        """Test if bonds have been written properly."""
        self.assertEqual(len(self.py_bonds), N_PART - 1)

        for i in range(N_PART - 1):
            bond = [x for x in self.py_bonds if x[0] == i][0]
            self.assertEqual(bond[0], i + 0)
            self.assertEqual(bond[1], i + 1)

    def test_script(self):
        with open(sys.argv[0], 'r') as f:
            ref = f.read()
        data = self.py_file['parameters/files'].attrs['script'].decode('utf-8')
        self.assertEqual(data, ref)

    def test_units(self):
        self.assertEqual(
            self.py_file['particles/atoms/id/time'].attrs['unit'], b'ps')
        self.assertEqual(
            self.py_file['particles/atoms/position/value'].attrs['unit'], b'm')
        if espressomd.has_features(['ELECTROSTATICS']):
            self.assertEqual(
                self.py_file['particles/atoms/charge/value'].attrs['unit'], b'e')
        if espressomd.has_features(['MASS']):
            self.assertEqual(
                self.py_file['particles/atoms/mass/value'].attrs['unit'], b'u')
        self.assertEqual(
            self.py_file['particles/atoms/force/value'].attrs['unit'],
            b'm u ps-2')
        self.assertEqual(
            self.py_file['particles/atoms/velocity/value'].attrs['unit'],
            b'm ps-1')

    def test_links(self):
        time_ref = self.py_id_time
        step_ref = self.py_id_step
        for group in "position", "velocity", "force", "charge", "mass", "image":
            time = self.py_file['particles/atoms/' + group + '/time'][1]
            step = self.py_file['particles/atoms/' + group + '/step'][1]
            self.assertEqual(time, time_ref)
            self.assertEqual(step, step_ref)

        bond_time = self.py_file['connectivity/atoms/time'][1]
        self.assertEqual(bond_time, time_ref)
        bond_step = self.py_file['connectivity/atoms/step'][1]
        self.assertEqual(bond_step, step_ref)
        box_time = self.py_file['particles/atoms/box/edges/time'][1]
        self.assertEqual(box_time, time_ref)
        box_step = self.py_file['particles/atoms/box/edges/step'][1]
        self.assertEqual(box_step, step_ref)


if __name__ == "__main__":
    ut.main()
