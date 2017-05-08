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

"""
Testmodule for the VTF file writing.
"""
import os
import sys
import unittest as ut
import numpy as np
import espressomd  # pylint: disable=import-error
from espressomd.io.writer import vtf

npart = 50


class CommonTests(ut.TestCase):
    """
    Class that holds common test methods.
    """
    system = espressomd.System()
    # avoid particles to be set outside of the main box, otherwise particle
    # positions are folded in the core when writing out and we cannot directly
    # compare positions in the dataset and where particles were set. One would
    # need to unfold the positions of the hdf5 file.
    system.box_l = [npart, npart, npart]
    system.cell_system.skin = 0.4
    system.time_step = 0.01
    written_pos = None
    types_to_write = None
    for i in range(npart):
        system.part.add(id=i, pos=np.array([float(i),
                                            float(i),
                                            float(i)]),
                        v=np.array([1.0, 2.0, 3.0]), type=1+(-1)**i)
    system.integrator.run(steps=0)


    def test_pos(self):
        """Test if positions have been written properly."""
        if self.types_to_write=='all': simulation_pos=np.array([((i), float(i), float(i), float(i)) for i in range(npart)]),
        elif (2 in self.types_to_write): simulation_pos=np.array([((i*2), float(i*2), float(i*2), float(i*2)) for i in range(npart//2)]),
        self.assertTrue(np.allclose(
            simulation_pos, self.written_pos),
            msg="Positions not written correctly by writevcf!")


class VCFTestAll(CommonTests):
    """
    Test the writing VTF files.
    """
    @classmethod
    def tearDownClass(cls):
        os.remove("test.vtf")

    @classmethod
    def setUpClass(cls):
        """Prepare a testsystem."""
        cls.types_to_write='all'
        with open('test.vtf','w') as fp:
            vtf.writevcf(cls.system, fp, types=cls.types_to_write)
        cls.written_pos=np.loadtxt("test.vtf",comments="t")

class VCFTestType(CommonTests):
    """
    Test the writing VTF files.
    """
    @classmethod
    def tearDownClass(cls):
        os.remove("test.vtf")

    @classmethod
    def setUpClass(cls):
        """Prepare a testsystem."""
        cls.types_to_write=[2, 23]
        with open('test.vtf','w') as fp:
            vtf.writevcf(cls.system, fp, types=cls.types_to_write)
        cls.written_pos=np.loadtxt("test.vtf",comments="t")



if __name__ == "__main__":
    suite = ut.TestLoader().loadTestsFromTestCase(VCFTestAll)
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(VCFTestType))

    result = ut.TextTestRunner(verbosity=4).run(suite)
    #if os.path.isfile("test.vtf"):
    #    os.remove("test.vtf")
    sys.exit(not result.wasSuccessful())

