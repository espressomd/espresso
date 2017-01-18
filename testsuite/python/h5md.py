#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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

"""Testmodule for the H5MD interface.
"""
import os
import unittest as ut
import numpy as np
import espressomd  # pylint: disable=import-error
import h5py  # h5py has to be imported *after* espressomd (MPI)


@ut.skipIf('H5MD' not in espressomd.code_info.features(),
           "H5MD not compiled in, can not check functionality.")
class H5mdTest(ut.TestCase):
    """Test the core implementation of writing hdf5 files."""
    @classmethod
    def setUpClass(self):
        """Prepare a testsystem."""
        from espressomd.io.writer import h5md  # pylint: disable=import-error
        self.system = espressomd.System()
        self.system.box_l = [20.0, 20.0, 20.0]
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.01
        for i in range(10):
            self.system.part.add(id=i, pos=np.array([float(i),
                                                     float(i),
                                                     float(i)]),
                                 v=np.array([1.0, 2.0, 3.0]), type=23)
            if 'MASS' in espressomd.code_info.features():
                self.system.part[i].mass = 2.3
            if 'EXTERNAL_FORCE' in espressomd.code_info.features():
                self.system.part[i].ext_force = [0.1, 0.2, 0.3]
        self.system.integrator.run(steps=0)
        self.h5 = h5md.H5md(filename="test.h5", write_pos=True, write_vel=True,
                            write_force=True, write_type=True, write_mass=True, write_ordered=True)
        self.h5.write()
        self.h5.close()
        del self.h5
        self.py_file = h5py.File("test.h5", 'r')
        self.py_pos = self.py_file['particles/atoms/position/value']
        self.py_vel = self.py_file['particles/atoms/velocity/value']
        self.py_f = self.py_file['particles/atoms/force/value']
        self.py_type = self.py_file['particles/atoms/type/value']
        self.py_mass = self.py_file['particles/atoms/mass/value']
        self.py_id = self.py_file['particles/atoms/id/value']

    @classmethod
    def tearDownClass(self):
        os.remove("test.h5")

    def test_pos(self):
        """Test if positions have been written properly."""
        self.assertTrue(np.allclose(
            np.array([(float(i), float(i), float(i)) for i in range(10)]),
            np.array([x for (y, x) in sorted(zip(self.py_id, self.py_pos))])),
                        msg="Positions not written correctly by H5md!")

    def test_vel(self):
        """Test if velocities have been written properly."""
        self.assertTrue(np.allclose(
            np.array([[1.0, 2.0, 3.0] for i in range(10)]),
            np.array([x for (y, x) in sorted(zip(self.py_id, self.py_vel))])),
                        msg="Velocities not written correctly by H5md!")

    @ut.skipIf('EXTERNAL_FORCE' not in espressomd.code_info.features(),
               "EXTERNAL_FORCE not compiled in, can not check writing forces.")
    def test_f(self):
        """Test if forces have been written properly."""
        self.assertTrue(np.allclose(
            np.array([[0.1, 0.2, 0.3] for i in range(10)]),
            np.array([x for (y, x) in sorted(zip(self.py_id, self.py_f))])),
                        msg="Forces not written correctly by H5md!")


if __name__ == "__main__":
    suite = ut.TestLoader().loadTestsFromTestCase(H5mdTest)
    ut.TextTestRunner(verbosity=2).run(suite)
