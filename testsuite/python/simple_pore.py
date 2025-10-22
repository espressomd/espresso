#
# Copyright (C) 2010-2022 The ESPResSo project
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

import unittest as ut
import unittest_decorators as utx
import numpy as np
import espressomd
import espressomd.shapes


class SimplePoreConstraint(ut.TestCase):

    box_yz = 15.
    box_x = 20.
    system = espressomd.System(box_l=[box_x, box_yz, box_yz])
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.constraints.clear()
        self.system.part.clear()
        self.system.non_bonded_inter.reset()

    def test_orientation(self):
        pore = espressomd.shapes.SimplePore(axis=[1., 0., 0.], radius=2., smoothing_radius=.1,
                                            length=2., center=[5., 5., 5.])

        d, _ = pore.calc_distance(position=[.0, .0, .0])
        self.assertGreater(d, 0.)

        d, _ = pore.calc_distance(position=[5., 5., .0])
        self.assertLess(d, 0.)

    def test_distance_calculation(self):
        """
        Check the distance calculation function for all 3 sections of
        a simple pore. The test works in cylindrical coordinates.
        """
        twopi = 2. * np.pi
        shape = espressomd.shapes.SimplePore(
            axis=[0., 0., 1.], radius=3., smoothing_radius=.1, length=6.,
            center=self.system.box_l / 2.)

        def calc_cartesian_coord(body, pos):
            """
            Convert cylindrical coordinates in the body frame
            to Cartesian coordinates in the lab frame.
            """
            origin = body.center
            r, phi, z = pos
            x = r * np.cos(phi) + origin[0]
            y = r * np.sin(phi) + origin[1]
            z = z + origin[2]
            return [x, y, z]

        def calc_distance(shape, pos):
            """Compute the distance vector from a shape surface."""
            return shape.calc_distance(
                position=calc_cartesian_coord(shape, pos))

        def angledist(val, ref):
            """Compute the absolute distance between two angles."""
            twopi = 2. * np.pi
            d1 = np.abs((ref - val) % twopi)
            d2 = np.abs((ref - val) % twopi - twopi)
            return min(d1, d2)

        tols = {"atol": 1e-10, "rtol": 1e-7}

        # scan points along the cylinder main axis
        for h in np.linspace(-1., 1., 21):
            z = h * (shape.length / 2. - shape.smoothing_radius)
            for phi in np.linspace(0., twopi, 31):
                dist, vec = calc_distance(shape, (0., phi, z))
                np.testing.assert_allclose(dist, shape.radius, **tols)
                np.testing.assert_allclose(vec[0], -shape.radius, **tols)
                np.testing.assert_allclose(vec[2], 0., **tols)

        # scan cylinder section
        for ref_dist in np.linspace(-0.99, 0.99, 21) * shape.smoothing_radius:
            r = shape.radius - ref_dist
            for h in np.linspace(-0.99, 0.99, 21):
                z = h * (shape.length / 2. - shape.smoothing_radius)
                for phi in np.linspace(0., twopi, 31):
                    pos = calc_cartesian_coord(shape, (r, phi, z))
                    dist, vec = shape.calc_distance(position=pos)
                    is_inside = shape.is_inside(position=pos)
                    cur_r = np.linalg.norm(vec[:2])
                    np.testing.assert_allclose(dist, ref_dist, **tols)
                    np.testing.assert_allclose(cur_r, np.abs(ref_dist), **tols)
                    np.testing.assert_allclose(vec[2], 0., **tols)
                    if np.abs(ref_dist) > 1e-4:
                        ref_phi = phi if ref_dist < 0. else phi + np.pi
                        cur_phi = np.arctan2(vec[1], vec[0])
                        angle_diff = angledist(cur_phi, ref_phi)
                        np.testing.assert_allclose(angle_diff, 0., **tols)
                        self.assertEqual(is_inside, ref_dist < 0.)

        # scan wall section
        for r in np.linspace(2.01, 3.01, 21) + shape.radius:
            for ref_dist in np.linspace(-2.0, 2.0, 21):
                for sgn in [+1, -1]:
                    z = ref_dist + shape.length / 2.
                    for phi in np.linspace(0., twopi, 31):
                        pos = calc_cartesian_coord(shape, (r, phi, sgn * z))
                        dist, vec = shape.calc_distance(position=pos)
                        is_inside = shape.is_inside(position=pos)
                        ref_z = sgn * ref_dist
                        np.testing.assert_allclose(dist, ref_dist, **tols)
                        np.testing.assert_allclose(vec[0], 0., **tols)
                        np.testing.assert_allclose(vec[1], 0., **tols)
                        np.testing.assert_allclose(vec[2], ref_z, **tols)
                        if np.abs(ref_dist) > 1e-4:
                            self.assertEqual(is_inside, ref_dist < 0.)

        # scan torus section
        for phi in np.linspace(0., twopi, 21):
            r = shape.radius + shape.smoothing_radius
            z = shape.length / 2. - shape.smoothing_radius
            origin = np.array(calc_cartesian_coord(shape, (r, phi, z)))
            for rho in np.linspace(0.01, 2.01, 21) * shape.smoothing_radius:
                ref_dist = rho - shape.smoothing_radius
                for theta in np.linspace(np.pi / 2. + 0.01, np.pi - 0.01, 31):
                    dr = rho * np.cos(theta)
                    dz = rho * np.sin(theta)
                    rdr = ref_dist * np.cos(theta)
                    rdz = ref_dist * np.sin(theta)
                    pos_cyl = (r + dr, phi, z + dz)
                    normal = (r + rdr, phi, z + rdz)
                    pos = calc_cartesian_coord(shape, pos_cyl)
                    dist, vec = shape.calc_distance(position=pos)
                    is_inside = shape.is_inside(position=pos)
                    ref_vec = calc_cartesian_coord(shape, normal) - origin
                    np.testing.assert_allclose(dist, ref_dist, **tols)
                    np.testing.assert_allclose(np.copy(vec), ref_vec, **tols)
                    if np.abs(ref_dist) > 1e-4:
                        self.assertEqual(is_inside, ref_dist < 0.)

    @utx.skipIfMissingFeatures(["LENNARD_JONES"])
    def test_scattering(self):
        """
        Create a line of particles along the pore main axis with a velocity
        vector pointing up. Particles inside the pore and outside the pore
        will bounce up and down (force vector is perpendicular to main axis).
        Particles at the pore mouth will scatter along the main axis due
        to reflection by the torus. Two peaks will appear in the trajectory
        of the particles x-position, which roughly follow a sine wave
        (we'll use a polynomial fit). Since we don't use a thermostat,
        the two peaks are perfectly symmetric.
        """
        system = self.system
        system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=1., sigma=1., cutoff=2**(1. / 6.), shift="auto")
        xpos_init = np.arange(800) * (self.box_x / 800.)

        def get_diffusion():
            return np.copy(self.system.part.all().pos[:, 0]) - xpos_init

        def get_series():
            ydata = get_diffusion()
            xdata = np.arange(250, 340)
            xdata = xdata[np.nonzero(ydata[250:340])[0]]
            ydata = ydata[xdata[0]:xdata[-1] + 1]
            xdata = np.array(xdata, dtype=float) / system.box_l[0]
            return (xdata, ydata)

        system.constraints.add(
            particle_type=0, penetrable=False, only_positive=False,
            shape=espressomd.shapes.SimplePore(
                axis=[1., 0., 0.], radius=3., smoothing_radius=.1, length=5.,
                center=system.box_l / 2.))

        for x0 in xpos_init:
            rpos = [x0, 0.5 * self.box_yz, 0.5 * self.box_yz]
            system.part.add(pos=rpos, type=1, v=[0., 1., 0.])

        # integrate until particles interact with the simple pore
        system.integrator.run(300)
        system.time = 0.
        xdata, ydata = get_series()
        coefs_init = np.polyfit(xdata, ydata, 4)[::-1]
        coefs_growth = np.array([5398.954, -1483.69746, 152.9315,
                                 -7.0115546, 0.1207058])

        # verify particles trajectory
        for _ in range(10):
            system.integrator.run(20)
            # check symmetry
            ydata = get_diffusion()[1:]
            np.testing.assert_allclose(ydata, -ydata[::-1], rtol=0., atol=1e-9)
            # check one peak
            coefs = coefs_init + system.time * coefs_growth
            xdata, ydata = get_series()
            ydata_hat = np.sum([coefs[j] * xdata**j for j in range(5)], axis=0)
            magnitude = -np.min(ydata_hat)
            deviation = (ydata - ydata_hat)[2:-2]  # remove outliers
            rmsd = np.sqrt(np.mean(deviation**2))
            self.assertLess(rmsd / magnitude, 0.01)

    @utx.skipIfMissingFeatures(["LENNARD_JONES"])
    def test_stability(self):
        """
        Stability test for a simple pore tilted along the box diagonal.

        The rationale is to hit the pore everywhere with particles
        and check that no singularity occurs. The cylinder is needed
        because the pore is tilted with respect to the box, without
        it particles could enter the constraint over the periodic
        boundaries, leading to force jumps.
        """
        box_yz = self.box_yz
        box_x = self.box_x
        system = self.system
        lj_eps = 1.0
        lj_sig = 1.0
        lj_cut = lj_sig * 2**(1. / 6.)

        system.constraints.add(
            particle_type=0, penetrable=False, only_positive=False,
            shape=espressomd.shapes.SimplePore(
                axis=[1., 0.5, 0.5], radius=3., smoothing_radius=.1,
                length=5, center=[.5 * box_x, .5 * box_yz, .5 * box_yz]))
        system.constraints.add(
            particle_type=0, penetrable=False, only_positive=False,
            shape=espressomd.shapes.Cylinder(
                axis=[1., 0, 0], radius=0.5 * box_yz, length=4 * lj_cut + box_x,
                center=[.5 * box_x, .5 * box_yz, .5 * box_yz], direction=-1))

        system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")

        for i in range(200):
            rpos = [i * (box_x / 200.), 0.5 * box_yz, 0.5 * box_yz]
            system.part.add(pos=rpos, type=1, v=[1., 1., 1.])

        start_energy = system.analysis.energy()['total']
        system.integrator.run(1000)
        end_energy = system.analysis.energy()['total']
        rel_diff = abs(end_energy - start_energy) / start_energy

        self.assertLess(rel_diff, 1e-3)


if __name__ == "__main__":
    ut.main()
