from __future__ import division, print_function

import unittest as ut

import espressomd

from espressomd import interactions
from espressomd.shapes import Ellipsoid

import numpy


@ut.skipIf(not espressomd.has_features(["CONSTRAINTS", "LENNARD_JONES_GENERIC"]),
           "Features not available, skipping test!")
class EllipsoidTest(ut.TestCase):

    def prepare(self, S):
        S.box_l = [30., 30., 30.]
        S.time_step = 0.01
        S.cell_system.skin = 0.4
        S.part.add(pos=[0., 0., 0.], type=0)

        # abuse generic LJ to measure distance via the potential V(r) = r -- (b1 = m)
        S.non_bonded_inter[0, 1].generic_lennard_jones.set_params(
            epsilon=1., sigma=1., cutoff=7., shift=0., offset=0., e1=-1, e2=0, b1=1., b2=0.)

    def pos_on_surface(self, theta, v, semiaxis0, semiaxis1, semiaxis2, center=numpy.array([15., 15., 15.])):
        """Return postion on ellipsoid surface."""
        pos = numpy.array([semiaxis0 * numpy.sqrt(1. - v * v) * numpy.cos(theta),
               semiaxis1 * numpy.sqrt(1. - v * v) * numpy.sin(theta),
               semiaxis2 * v])
        return pos + center

    def test_distance(self):
        S = espressomd.System()
        self.prepare(S)

        N = 10

        # check oblate ellipsoid

        semiaxes = [2.18, 5.45]
        e = Ellipsoid(a=semiaxes[0], b=semiaxes[1], center = [15,15,15], direction=+1)

        constraint_e = espressomd.constraints.ShapeBasedConstraint(
            shape=e, particle_type=1, penetrable=1)
        const1 = S.constraints.add(constraint_e)

        for i in range(N):
            for j in range(N):
                theta = 2. * i / float(N) * numpy.pi
                v = j / float(N-1) * 2. - 1
                pos = self.pos_on_surface(theta, v, semiaxes[0], semiaxes[1], semiaxes[1])
                S.part[0].pos = pos
                S.integrator.run(recalc_forces=True, steps=0)
                energy = S.analysis.energy()
                self.assertAlmostEqual(energy["total"], 0., places=6)

        S.constraints.remove(const1)

        # check prolate ellipsoid

        semiaxes = [3.61, 2.23]
        e = Ellipsoid(a=semiaxes[0], b=semiaxes[1], center = [15,15,15], direction=+1)

        constraint_e = espressomd.constraints.ShapeBasedConstraint(
            shape=e, particle_type=1, penetrable=1)
        const1 = S.constraints.add(constraint_e)

        for i in range(N):
            for j in range(N):
                theta = 2. * i / float(N) * numpy.pi
                v = j / float(N-1) * 2. - 1
                pos = self.pos_on_surface(theta, v, semiaxes[0], semiaxes[1], semiaxes[1])
                S.part[0].pos = pos
                S.integrator.run(recalc_forces=True, steps=0)
                energy = S.analysis.energy()
                self.assertAlmostEqual(energy["total"], 0., places=6)

        # check sphere (multiple distances from surface)

        # change ellipsoid paramters instead of creating a new constraint
        e.a = 1.
        e.b = 1.

        radii = numpy.linspace(1., 6.5, 7)

        for i in range(N):
            for j in range(N):
                theta = 2. * i / float(N) * numpy.pi
                v = j / float(N-1) * 2. - 1
                for r in radii:
                    pos = self.pos_on_surface(theta, v, r, r, r)
                    S.part[0].pos = pos
                    S.integrator.run(recalc_forces=True, steps=0)
                    energy = S.analysis.energy()
                    self.assertAlmostEqual(energy["total"], r-1.)#, places=7)


if __name__ == "__main__":
    ut.main()
