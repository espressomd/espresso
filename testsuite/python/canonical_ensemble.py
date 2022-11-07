#
# Copyright (C) 2022 The ESPResSo project
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
import espressomd.constraints
import espressomd.observables
import espressomd.accumulators
import espressomd.reaction_methods

import numpy as np
import scipy.optimize
import unittest as ut
import unittest_decorators as utx


class Test(ut.TestCase):

    system = espressomd.System(box_l=[1., 1., 1.])
    system.cell_system.skin = 0.4
    system.setup_type_map(type_list=[0])

    def tearDown(self):
        self.system.part.clear()
        self.system.constraints.clear()

    @utx.skipIfMissingFeatures("ELECTROSTATICS")
    def test_linear_potential(self):
        """
        Set a particle in a box with a linear potential along the x-axis.
        The particle distribution resulting from accepted Monte Carlo moves
        should follow a Maxwell-Boltzmann distribution.
        """

        method = espressomd.reaction_methods.ReactionEnsemble(
            kT=0.2, seed=42, exclusion_range=0., search_algorithm="order_n")
        method.set_non_interacting_type(type=1)

        p = self.system.part.add(pos=[0., 0., 0.], q=1, type=0)
        obs_pos = espressomd.observables.ParticlePositions(ids=(p.id,))
        acc_pos = espressomd.accumulators.TimeSeries(obs=obs_pos)

        E = np.array([-1., 0., 0.])
        field = espressomd.constraints.LinearElectricPotential(E=E, phi0=0.)
        self.system.constraints.add(field)

        for _ in range(100000):
            method.do_particle_displacement_MC_move()
            acc_pos.update()

        # the x-position should follow an exponential distribution
        # -> mean = kT, median = kT x ln(2), variance = kT^2
        series = acc_pos.time_series()[:, p.id, 0]
        ydata, xbins = np.histogram(series, bins=15, range=[0., 1.])
        xdata = (xbins[1:] + xbins[:-1]) / 2.
        ydata = ydata / float(ydata[0])
        (a, b, c), _ = scipy.optimize.curve_fit(
            lambda x, a, b, c: a * np.exp(-b * x) + c, xdata, ydata)
        # check histogram profile is roughly exponential

        self.assertAlmostEqual(a, 1., delta=0.2)
        self.assertAlmostEqual(b, 1. / method.kT, delta=0.3)
        self.assertAlmostEqual(c, 0., delta=0.01)
        # check distribution parameters with high accuracy
        ln2 = np.log(2)
        self.assertAlmostEqual(np.mean(series), method.kT, delta=0.02)
        self.assertAlmostEqual(np.median(series) / ln2, method.kT, delta=0.02)
        self.assertAlmostEqual(np.sqrt(np.var(series)), method.kT, delta=0.02)

        # the y- and z-position should follow a uniform distribution
        for axis in (1, 2):
            series = acc_pos.time_series()[:, p.id, axis]
            ydata, _ = np.histogram(series, bins=10, range=[0., 1.])
            ydata = ydata / np.mean(ydata)
            np.testing.assert_allclose(ydata, 1., atol=0.25)


if __name__ == "__main__":
    ut.main()
