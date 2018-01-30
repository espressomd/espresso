from __future__ import division, print_function

import unittest as ut

import espressomd

from espressomd import interactions
from espressomd.shapes import Ellipsoid

import numpy


@ut.skipIf(not espressomd.has_features(["CONSTRAINTS", "LENNARD_JONES_GENERIC"]),
           "Features not available, skipping test!")
class EllipsoidTest(ut.TestCase):

    box_l = 30.

    def prepare(self, system):
        system.box_l = [self.box_l, self.box_l, self.box_l]
        system.time_step = 0.01
        system.cell_system.skin = 0.4
        system.part.add(pos=[0., 0., 0.], type=0)

        # abuse generic LJ to measure distance via the potential V(r) = r
        system.non_bonded_inter[0, 1].generic_lennard_jones.set_params(
            epsilon=1., sigma=1., cutoff=7., shift=0., offset=0., e1=-1, e2=0, b1=1., b2=0.)

    def pos_on_surface(self, theta, v, semiaxis0, semiaxis1, semiaxis2, center=numpy.array([15, 15, 15])):
        """Return postion on ellipsoid surface."""
        pos = numpy.array([semiaxis0 * numpy.sqrt(1. - v * v) * numpy.cos(theta),
                           semiaxis1 *
                           numpy.sqrt(1. - v * v) * numpy.sin(theta),
                           semiaxis2 * v])
        return pos + center

    def test_distance(self):
        system = espressomd.System()
        self.prepare(system)

        N = 10

        # check oblate ellipsoid

        semiaxes = [2.18, 5.45]
        e = Ellipsoid(a=semiaxes[0], b=semiaxes[1],
                      center=[self.box_l / 2., self.box_l / 2., self.box_l / 2.], direction=+1)

        constraint_e = espressomd.constraints.ShapeBasedConstraint(
            shape=e, particle_type=1, penetrable=1)
        const1 = system.constraints.add(constraint_e)

        for i in range(N):
            for j in range(N):
                theta = 2. * i / float(N) * numpy.pi
                v = j / float(N - 1) * 2. - 1
                pos = self.pos_on_surface(
                    theta, v, semiaxes[0], semiaxes[1], semiaxes[1])
                system.part[0].pos = pos
                system.integrator.run(recalc_forces=True, steps=0)
                energy = system.analysis.energy()
                self.assertAlmostEqual(energy["total"], 0., places=6)

        system.constraints.remove(const1)

        # check prolate ellipsoid

        semiaxes = [3.61, 2.23]
        e = Ellipsoid(a=semiaxes[0], b=semiaxes[1],
                      center=[self.box_l / 2., self.box_l / 2., self.box_l / 2.], direction=+1)

        constraint_e = espressomd.constraints.ShapeBasedConstraint(
            shape=e, particle_type=1, penetrable=1)
        const1 = system.constraints.add(constraint_e)

        for i in range(N):
            for j in range(N):
                theta = 2. * i / float(N) * numpy.pi
                v = j / float(N - 1) * 2. - 1
                pos = self.pos_on_surface(
                    theta, v, semiaxes[0], semiaxes[1], semiaxes[1])
                system.part[0].pos = pos
                system.integrator.run(recalc_forces=True, steps=0)
                energy = system.analysis.energy()
                self.assertAlmostEqual(energy["total"], 0., places=6)

        # check sphere (multiple distances from surface)

        # change ellipsoid paramters instead of creating a new constraint
        e.a = 1.
        e.b = 1.

        radii = numpy.linspace(1., 6.5, 7)

        for i in range(N):
            for j in range(N):
                theta = 2. * i / float(N) * numpy.pi
                v = j / float(N - 1) * 2. - 1
                for r in radii:
                    pos = self.pos_on_surface(theta, v, r, r, r)
                    system.part[0].pos = pos
                    system.integrator.run(recalc_forces=True, steps=0)
                    energy = system.analysis.energy()
                    self.assertAlmostEqual(energy["total"], r - 1.)


if __name__ == "__main__":
    ut.main()
