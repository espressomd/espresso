#
# Copyright (C) 2026 The ESPResSo project
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

import espressomd
import espressomd.interactions
import numpy as np
import unittest as ut
import unittest_decorators as utx
import contextlib
import itertools

with contextlib.suppress(ImportError):
    import object_in_fluid as oif


@utx.skipIfMissingModules("object_in_fluid")
class Test(ut.TestCase):

    system = espressomd.System(box_l=[10., 10., 10.])
    system.time_step = 0.01
    system.cell_system.skin = 0.1

    def tearDown(self):
        self.system.thermostat.turn_off()
        self.system.part.clear()

    def test_oif_local_forces(self):
        # create a square of side 1 made of two right triangles
        p1, p2, p3, p4 = self.system.part.add(
            pos=[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.]])
        p2.v = [-1., +3., +16.]
        p3.v = [+2., -6., -16.]
        r0 = np.sqrt(2.)
        phi0 = np.pi
        A01 = 0.5
        A02 = 0.5

        def setup_bond(**kwargs):
            defs = dict(
                r0=r0, ks=0., kslin=1., kb=1., kal=1., kvisc=0.,
                phi0=phi0, A01=A01, A02=A02)
            bond = espressomd.interactions.OifLocalForces(**{**defs, **kwargs})
            self.system.bonded_inter.add(bond)
            p2.bonds = ((bond, p1, p3, p4),)
            return bond

        # relax surface to a bent shape
        for kb in [0.5, 1., 6.]:
            bond = setup_bond(phi0=np.pi / 2., kb=kb)
            self.system.integrator.run(recalc_forces=True, steps=0)
            forces_sim = np.copy(self.system.part.all().f)
            forces_ref = kb * np.outer([-1, 1, 1, -1],
                                       [0., 0., np.pi / np.sqrt(2.)])
            np.testing.assert_almost_equal(forces_sim, forces_ref)
            forces_utils = np.array(oif.oif_calc_bending_force(
                kb, np.copy(p1.pos), np.copy(p2.pos), np.copy(p3.pos),
                np.copy(p4.pos), bond.phi0, phi0)).reshape((-1, 3)) * (np.pi / np.sqrt(2.))
            np.testing.assert_almost_equal(forces_utils, forces_ref[(0, 3), :])

        # stretch the square diagonal to make a lozenge
        for scale, kvisc in itertools.product([0.5, 2.], [0., 3.]):
            vij = np.copy(p3.v - p2.v)
            eij = np.copy(p3.pos - p2.pos) / r0
            forces_stretch = np.outer([0, -1, 1, 0], [1, -1, 0])
            forces_visc = np.outer([0, -1, 1, 0], np.dot(vij, eij) * eij)
            # linear elastic
            for kslin in [0.5, 1., 6.]:
                bond = setup_bond(r0=r0 * scale, kslin=kslin, kvisc=kvisc)
                self.system.integrator.run(recalc_forces=True, steps=0)
                forces_sim = np.copy(self.system.part.all().f)
                forces_ref = kslin * forces_stretch * (1. - scale)
                forces_ref -= kvisc * forces_visc
                np.testing.assert_almost_equal(forces_sim, forces_ref)
                if kvisc == 0.:
                    force_utils = oif.oif_calc_linear_stretching_force(
                        kslin, np.copy(p2.pos), np.copy(p3.pos), r0 * scale, r0)
                    np.testing.assert_almost_equal(force_utils, forces_ref[1])
            # neo-Hookean
            for ks in [0.5, 1., 6.]:
                bond = setup_bond(r0=r0 * scale, kslin=0., ks=ks, kvisc=kvisc)
                self.system.integrator.run(recalc_forces=True, steps=0)
                s = r0 / bond.r0
                fac = (s**0.5 + s**-2.5) / (s + s**-3.)
                forces_sim = np.copy(self.system.part.all().f)
                forces_ref = ks * fac * forces_stretch * (1. - scale)
                forces_ref -= kvisc * forces_visc
                np.testing.assert_almost_equal(forces_sim, forces_ref)
                if kvisc == 0.:
                    fac_utils = oif.oif_neo_hookean_nonlin(s)
                    np.testing.assert_almost_equal(fac_utils, fac)
                    force_utils = oif.oif_calc_stretching_force(
                        ks, np.copy(p2.pos), np.copy(p3.pos), r0 * scale, r0)
                    np.testing.assert_almost_equal(force_utils, forces_ref[1])

        # relax each triangle to a different surface area
        for scale in list(1 / n for n in range(2, 7)) + list(range(2, 7)):
            # first triangle
            bond = setup_bond(A01=A01 * scale, kal=48.)
            self.system.integrator.run(recalc_forces=True, steps=0)
            forces_sim = np.copy(self.system.part.all().f)
            forces_ref = np.zeros((4, 3))
            forces_ref[0:3, 0:2] = [[1, 1], [-2, 1], [1, -2]]
            if scale > 1.:
                forces_ref *= -2. * (scale - 1.)
            else:
                forces_ref *= 2. * (1. / scale - 1.) * scale
            forces_utils = np.array(oif.oif_calc_local_area_force(
                48., np.copy(p1.pos), np.copy(p2.pos), np.copy(p3.pos),
                bond.A01, A01)).reshape((-1, 3)) / 3.
            np.testing.assert_almost_equal(forces_sim, forces_ref)
            np.testing.assert_almost_equal(forces_utils, forces_ref[:3, :])
            # second triangle
            bond = setup_bond(A02=A02 * scale, kal=48.)
            self.system.integrator.run(recalc_forces=True, steps=0)
            forces_sim = np.copy(self.system.part.all().f)
            forces_ref = np.zeros((4, 3))
            forces_ref[1:4, 0:2] = [[-1, 2], [2, -1], [-1, -1]]
            if scale > 1.:
                forces_ref *= -2. * (scale - 1.)
            else:
                forces_ref *= 2. * (1. / scale - 1.) * scale
            forces_utils = np.array(oif.oif_calc_local_area_force(
                48., np.copy(p2.pos), np.copy(p3.pos), np.copy(p4.pos),
                bond.A02, A02)).reshape((-1, 3)) / 3.
            np.testing.assert_almost_equal(forces_sim, forces_ref)
            np.testing.assert_almost_equal(forces_utils, forces_ref[1:, :])

    def test_utils(self):
        a = np.array([-1.1, 2.3, 4.6])
        b = np.array([-3.2, 5.3, 2.8])
        c = np.array([-5.4, 1.7, 9.6])
        np.testing.assert_almost_equal(oif.norm(a), np.linalg.norm(a))
        np.testing.assert_almost_equal(oif.norm(b), np.linalg.norm(b))
        np.testing.assert_almost_equal(oif.norm(c), np.linalg.norm(c))
        ref_n = np.cross(b - a, c - a)
        np.testing.assert_almost_equal(oif.get_triangle_normal(a, b, c), ref_n)
        ref_area = 0.5 * np.linalg.norm(np.cross(b - a, c - a))
        np.testing.assert_almost_equal(oif.area_triangle(a, b, c), ref_area)
        ref_vec = np.linalg.norm(b - a)
        np.testing.assert_almost_equal(oif.vec_distance(a, b), ref_vec)


if __name__ == "__main__":
    ut.main()
