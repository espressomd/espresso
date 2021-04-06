#
# Copyright (C) 2017-2021 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published byss
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

import unittest as ut
import espressomd
import numpy as np
import itertools


class StructureFactorTest(ut.TestCase):
    '''
    Test structure factor analysis against rectangular lattices.
    We do not check the wavevectors directly, but rather the
    quantity (wavevectors * 2 / pi)^2, which is more readable.
    '''

    box_l = 16
    part_ty = 0
    sf_order = 16
    system = espressomd.System(box_l=[box_l, box_l, box_l])

    def tearDown(self):
        self.system.part.clear()

    def generate_peaks(self, a, b, c, condition, cutoff=box_l):
        '''
        Generate the main diffraction peaks for crystal structures.

        Parameters
        ----------
        a: :obj:`float`
            Length of the unit cell on the x-axis.
        b: :obj:`float`
            Length of the unit cell on the y-axis.
        c: :obj:`float`
            Length of the unit cell on the z-axis.
        condition: :obj:`function`
            Reflection conditions for the crystal lattice.
        cutoff: :obj:`float` (optional)
            Cutoff value for the generation of reflections.
        '''
        reflections = [np.linalg.norm([h * a, k * b, l * c])
                       for (h, k, l) in itertools.product(np.arange(0, 10), repeat=3)
                       if condition(h, k, l)]
        return (np.array([x for x in sorted(set(reflections))
                          if 0. < x < 1.001 * cutoff]) / 4) ** 2

    def test_tetragonal(self):
        """Check tetragonal lattice."""
        xen = range(0, self.box_l, 2)
        yen = range(0, self.box_l, 4)
        zen = range(0, self.box_l, 8)
        for i, j, k in itertools.product(xen, yen, zen):
            self.system.part.add(type=self.part_ty, pos=(i, j, k))
        wavevectors, intensities = self.system.analysis.structure_factor(
            sf_types=[self.part_ty], sf_order=self.sf_order)
        intensities_int = np.around(intensities).astype(int)
        # check all reflections
        # no condition on (h,k,l)
        peaks_ref = self.generate_peaks(2, 4, 8, lambda h, k, l: True)
        peaks = (wavevectors[np.nonzero(intensities_int)] * 2 / np.pi)**2
        np.testing.assert_array_almost_equal(peaks, peaks_ref)

    def test_sc(self):
        """Check simple cubic lattice."""
        xen = range(0, self.box_l, 4)
        for i, j, k in itertools.product(xen, repeat=3):
            self.system.part.add(type=self.part_ty, pos=(i, j, k))
        wavevectors, intensities = self.system.analysis.structure_factor(
            sf_types=[self.part_ty], sf_order=self.sf_order)
        intensities_int = np.around(intensities).astype(int)
        np.testing.assert_array_almost_equal(intensities, intensities_int)
        intensities_ref = np.zeros(intensities.shape)
        intensities_ref[np.nonzero(intensities_int)] = len(self.system.part)
        np.testing.assert_array_equal(intensities_int, intensities_ref)
        # check all reflections
        # no condition on (h,k,l)
        peaks_ref = self.generate_peaks(4, 4, 4, lambda h, k, l: True)
        peaks = (wavevectors[np.nonzero(intensities_int)] * 2 / np.pi)**2
        np.testing.assert_array_almost_equal(peaks, peaks_ref)

    def test_bcc(self):
        """Check body-centered cubic lattice."""
        xen = range(0, self.box_l, 4)
        for i, j, k in itertools.product(xen, repeat=3):
            self.system.part.add(type=self.part_ty, pos=(i, j, k))
            self.system.part.add(type=self.part_ty, pos=(i + 2, j + 2, k + 2))
        wavevectors, intensities = self.system.analysis.structure_factor(
            sf_types=[self.part_ty], sf_order=self.sf_order)
        intensities_int = np.around(intensities).astype(int)
        np.testing.assert_array_almost_equal(intensities, intensities_int)
        intensities_ref = np.zeros(intensities.shape)
        intensities_ref[np.nonzero(intensities_int)] = len(self.system.part)
        np.testing.assert_array_equal(intensities_int, intensities_ref)
        # check all reflections
        # (h+k+l) even => F = 2f, otherwise F = 0
        peaks_ref = self.generate_peaks(
            4, 4, 4, lambda h, k, l: (h + k + l) % 2 == 0)
        peaks = (wavevectors[np.nonzero(intensities_int)] * 2 / np.pi)**2
        np.testing.assert_array_almost_equal(peaks, peaks_ref)

    def test_fcc(self):
        """Check face-centered cubic lattice."""
        xen = range(0, self.box_l, 4)
        for i, j, k in itertools.product(xen, repeat=3):
            self.system.part.add(type=self.part_ty, pos=(i, j, k))
            self.system.part.add(type=self.part_ty, pos=(i + 2, j + 2, k))
            self.system.part.add(type=self.part_ty, pos=(i + 2, j, k + 2))
            self.system.part.add(type=self.part_ty, pos=(i, j + 2, k + 2))
        wavevectors, intensities = self.system.analysis.structure_factor(
            sf_types=[self.part_ty], sf_order=self.sf_order)
        intensities_int = np.around(intensities).astype(int)
        np.testing.assert_array_almost_equal(intensities, intensities_int)
        intensities_ref = np.zeros(intensities.shape)
        intensities_ref[np.nonzero(intensities_int)] = len(self.system.part)
        np.testing.assert_array_equal(intensities_int, intensities_ref)
        # check all reflections
        # (h,k,l) all even or odd => F = 4f, otherwise F = 0
        peaks_ref = self.generate_peaks(
            4, 4, 4, lambda h, k, l:
            h % 2 == 0 and k % 2 == 0 and l % 2 == 0 or
            h % 2 == 1 and k % 2 == 1 and l % 2 == 1)
        peaks = (wavevectors[np.nonzero(intensities_int)] * 2 / np.pi)**2
        np.testing.assert_array_almost_equal(peaks, peaks_ref)

    def test_cco(self):
        """Check c-centered orthorhombic lattice."""
        xen = range(0, self.box_l, 4)
        for i, j, k in itertools.product(xen, repeat=3):
            self.system.part.add(type=self.part_ty, pos=(i, j, k))
            self.system.part.add(type=self.part_ty, pos=(i + 2, j + 2, k))
        wavevectors, intensities = self.system.analysis.structure_factor(
            sf_types=[self.part_ty], sf_order=self.sf_order)
        intensities_int = np.around(intensities).astype(int)
        # check all reflections
        peaks_ref = self.generate_peaks(4, 4, 4, lambda h, k, l: True)
        peaks = (wavevectors[np.nonzero(intensities_int)] * 2 / np.pi)**2
        np.testing.assert_array_almost_equal(peaks, peaks_ref)
        # check main reflections
        # (h+k) even => F = 2f, otherwise F = 0
        main_peaks_ref = 4 * self.generate_peaks(
            4, 4, 4, lambda h, k, l: (h + k) % 2 == 0, self.box_l / 2)
        peaks = (wavevectors * 2 / np.pi)**2
        for main_peak in main_peaks_ref:
            idx = (np.abs(peaks - main_peak)).argmin()
            self.assertAlmostEqual(peaks[idx], main_peak)
            self.assertAlmostEqual(intensities[idx], len(self.system.part))

    def test_exceptions(self):
        with self.assertRaisesRegex(ValueError, 'order has to be a strictly positive number'):
            self.system.analysis.structure_factor(sf_types=[0], sf_order=0)


if __name__ == "__main__":
    ut.main()
