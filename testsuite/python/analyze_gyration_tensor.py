# Copyright (C) 2010-2019 The ESPResSo project
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
import unittest as ut
import numpy as np
import espressomd


class AnalyzeGyration(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    np.random.seed(1234)
    cube_len = 4
    type_cube = 0
    type_stick = 1

    @classmethod
    def setUpClass(cls):
        box_l = 20.0
        cube_centre = 0.5 * (cls.cube_len - 1)
        cls.system.box_l = np.array([box_l, box_l, box_l])
        cls.system.cell_system.set_n_square(use_verlet_lists=False)
        # 4x4 cube
        for x, y, z in np.ndindex((cls.cube_len, cls.cube_len, cls.cube_len)):
            cls.system.part.add(pos=[x, y, z], type=cls.type_cube)
        # long stick in z, force z as principal axis
        for x, y, z in np.ndindex((1, 1, 10)):
            cls.system.part.add(
                pos=[x + cube_centre, y + cube_centre, z + cls.cube_len], type=cls.type_stick)
        # two small nubs in y, force y as secondary axis
        cls.system.part.add(
            pos=[cube_centre, cls.cube_len, cube_centre], type=cls.type_stick)
        cls.system.part.add(
            pos=[cube_centre, -1, cube_centre], type=cls.type_stick)

    def test_gyration_tensor_cube(self):
        # get results
        res = self.system.analysis.gyration_tensor(p_type=self.type_cube)
        rg = self.system.analysis.calc_rg(
            chain_start=0, number_of_chains=1, chain_length=self.cube_len**3)[0]
        # check eigenvalues are identical
        np.testing.assert_allclose(
            [np.abs(res['eva' + x][0]) for x in '012'], 3 * [1.25], atol=1e-6)
        np.testing.assert_allclose(rg**2, res['Rg^2'], atol=1e-6)

    def test_gyration_tensor(self):
        # get results
        res = self.system.analysis.gyration_tensor(
            p_type=[self.type_stick, self.type_cube])
        rg = self.system.analysis.calc_rg(
            chain_start=0, number_of_chains=1, chain_length=len(self.system.part[:]))[0]
        # check eigenvectors
        np.testing.assert_allclose(
            np.abs(res['eva0'][1]), [0., 0., 1.], atol=1e-6)
        np.testing.assert_allclose(
            np.abs(res['eva1'][1]), [0., 1., 0.], atol=1e-6)
        np.testing.assert_allclose(
            np.abs(res['eva2'][1]), [1., 0., 0.], atol=1e-6)
        np.testing.assert_allclose(rg**2, res['Rg^2'], atol=1e-6)

    def test_mom_intertia(self):
        sqr_dist = np.sum(
            (self.system.analysis.center_of_mass(p_type=0) - self.system.part.select(type=0).pos)**2, axis=0)
        mom_I = self.system.analysis.moment_of_inertia_matrix(p_type=0)
        # the cube case should have zero as off-diagonal components
        np.testing.assert_allclose(
            mom_I, np.diag(np.diag(mom_I)), rtol=0, atol=1e-6)
        np.testing.assert_allclose(
            np.diag(mom_I), sqr_dist[(1, 0, 1), ] + sqr_dist[2], atol=1e-6)


if __name__ == "__main__":
    ut.main()
