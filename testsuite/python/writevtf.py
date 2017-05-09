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
    system.set_random_state_PRNG()
    # make sure the box is small enough for some to be printed outside
    system.box_l = [npart*.5, npart*.2, npart*.2]
    system.cell_system.skin = 0.4
    system.time_step = 0.01
    written_pos = None
    written_bonds = None
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


    def test_bonds(self):
        """Test if bonds have been written properly."""
        
        if self.types_to_write=='all': 
            types=[0,2]
        elif (2 in self.types_to_write): 
            types=[2]
            sim_bonds=[]
            for p in system.part:
                if p.type in types:
                    for b in p.bonds:
                      if (b[1].type in types)
                          sim_bonds.append([p.id, b[1].type])
            sim_bonds=np.array(sim_bonds)
            print ("sim bonds")
            print (sim_bonds)
            print ("Written bonds")
            print (written_bonds)

            ##########
            # WRITE CODE HERE!
            #########
        #self.assertTrue(np.allclose(
        #    sim_bonds, written_bonds),
        #    msg="Bonds not written correctly by writevsf!")


class VCFTestAll(CommonTests):
    """
    Test the writing VTF files.
    """
    @classmethod
    def tearDownClass(cls):
        os.remove("test.vcf")
#        os.remove("test.vsf")

    @classmethod
    def setUpClass(cls):
        """Prepare a testsystem."""
        cls.types_to_write='all'
        with open('test.vcf','w') as fp:
            vtf.writevcf(cls.system, fp, types=cls.types_to_write)
        cls.written_pos=np.loadtxt("test.vcf",comments="t")

        with open('test.vsf','w') as fp:
            vtf.writevsf(cls.system, fp, types=cls.types_to_write)
        
        ###########
        ## CODE HERE!
        ###########
        cls.written_bonds=np.loadtxt("test.vsf", skiplines=1, comments="a", delimiter=':')

class VCFTestType(CommonTests):
    """
    Test the writing VTF files.
    """
    @classmethod
    def tearDownClass(cls):
        os.remove("test.vcf")

    @classmethod
    def setUpClass(cls):
        """Prepare a testsystem."""
        cls.types_to_write=[2, 23]
        with open('test.vcf','w') as fp:
            vtf.writevcf(cls.system, fp, types=cls.types_to_write)
        cls.written_pos=np.loadtxt("test.vcf",comments="t")



if __name__ == "__main__":
    suite = ut.TestLoader().loadTestsFromTestCase(VCFTestAll)
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(VCFTestType))

    result = ut.TextTestRunner(verbosity=4).run(suite)
    if os.path.isfile("test.vcf"):
        os.remove("test.vcf")
    #if os.path.isfile("test.vsf"):
    #    os.remove("test.vsf")
    sys.exit(not result.wasSuccessful())

