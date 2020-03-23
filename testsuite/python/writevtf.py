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
Testmodule for the VTF file writing.
"""
import sys
import unittest as ut
import numpy as np
import espressomd
from espressomd import interactions
from espressomd.io.writer import vtf
import tempfile

npart = 50


class CommonTests(ut.TestCase):

    """
    Class that holds common test methods.
    """
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    # avoid particles to be set outside of the main box, otherwise particle
    # positions are folded in the core when writing out and we cannot directly
    # compare positions in the dataset and where particles were set. One would
    # need to unfold the positions of the hdf5 file.
    system.box_l = 3 * [npart]
    system.cell_system.skin = 0.4
    system.time_step = 0.01
    written_pos = None
    written_bonds = None
    written_atoms = None

    types_to_write = None
    for i in range(npart):
        system.part.add(id=i, pos=np.array(3 * [i], dtype=float),
                        v=np.array([1.0, 2.0, 3.0]), type=1 + (-1)**i)

    system.bonded_inter.add(interactions.FeneBond(k=1., d_r_max=10.0))
    system.part[0].add_bond((0, 1))
    system.part[0].add_bond((0, 2))
    system.part[0].add_bond((0, 3))

    system.integrator.run(steps=0)

    def test_pos(self):
        """Test if positions have been written properly."""
        if self.types_to_write == 'all':
            simulation_pos = np.array(
                [((i), float(i), float(i), float(i)) for i in range(npart)])
        elif 2 in self.types_to_write:
            simulation_pos = np.array(
                [((i * 2), float(i * 2), float(i * 2), float(i * 2)) for i in range(npart // 2)])

        np.testing.assert_allclose(
            simulation_pos[:, 1:], self.written_pos[:, 1:],
            err_msg="Positions not written correctly by writevcf!")

    def test_bonds(self):
        """Test if bonds have been written properly: just look at number of bonds"""
        if self.types_to_write == 'all':
            simulation_bonds = np.array([1, 2, 3])  # the two bonded particles
        elif 2 in self.types_to_write:
            simulation_bonds = np.array(2)  # only this one is type 2

        np.testing.assert_allclose(
            np.shape(simulation_bonds), np.shape(self.written_bonds),
            err_msg="Bonds not written correctly by writevsf!")

    def test_atoms(self):
        """Test if atom declarations have been written properly."""
        if self.types_to_write == 'all':
            simulation_atoms = np.array(
                [((i), (1 + (-1)**i)) for i in range(npart)])
        elif 2 in self.types_to_write:
            simulation_atoms = np.array([((i * 2), 2)
                                         for i in range(npart // 2)])

        np.testing.assert_allclose(
            simulation_atoms[:, 1], self.written_atoms[:, 1],
            err_msg="Atoms not written correctly by writevsf!")


class VCFTestAll(CommonTests):

    """
    Test the writing VTF files: all particle types.
    """

    @classmethod
    def setUpClass(cls):
        """Prepare a testsystem."""
        cls.types_to_write = 'all'

        with tempfile.TemporaryFile(mode='w+t') as fp:
            vtf.writevcf(cls.system, fp, types=cls.types_to_write)
            fp.flush()
            fp.seek(0)
            cls.written_pos = np.loadtxt(fp, comments="t")

        with tempfile.TemporaryFile(mode='w+t') as fp:
            vtf.writevsf(cls.system, fp, types=cls.types_to_write)
            fp.flush()
            fp.seek(0)
            cls.written_bonds = np.loadtxt(
                fp,
                skiprows=1,
                comments="a",
                delimiter=":",
                usecols=[1])  # just the second bonded member
            fp.seek(0)
            cls.written_atoms = np.loadtxt(
                fp, skiprows=1, comments="b", usecols=[
                    1, 7])  # just the part_ID and type_ID


class VCFTestType(CommonTests):

    """
    Test the writing VTF files: only particle types 2 and 23.
    """

    @classmethod
    def setUpClass(cls):
        """Prepare a testsystem."""
        cls.types_to_write = [2, 23]
        with tempfile.TemporaryFile(mode='w+') as fp:

            vtf.writevcf(cls.system, fp, types=cls.types_to_write)
            fp.flush()
            fp.seek(0)
            cls.written_pos = np.loadtxt(fp, comments="t")

        with tempfile.TemporaryFile(mode='w+') as fp:
            vtf.writevsf(cls.system, fp, types=cls.types_to_write)
            fp.flush()
            fp.seek(0)
            cls.written_bonds = np.loadtxt(
                fp,
                skiprows=1,
                comments="a",
                delimiter=":",
                usecols=[1])  # just the second bonded member
            fp.seek(0)
            cls.written_atoms = np.loadtxt(
                fp, skiprows=1, comments="b",
                usecols=[1, 7])  # just the part_ID and type_ID


if __name__ == "__main__":
    suite = ut.TestLoader().loadTestsFromTestCase(VCFTestAll)
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(VCFTestType))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
