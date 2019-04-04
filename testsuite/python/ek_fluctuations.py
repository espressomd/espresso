#Copyright(C) 2011, 2012, 2013, 2014, 2015, 2016 The ESPResSo project
#
#This file is part of ESPResSo.
#
#ESPResSo is free software : you can redistribute it and / or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#ESPResSo is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
#GNU General Public License for more details.
#

from __future__ import print_function
import unittest as ut
import espressomd
import espressomd.electrokinetics
from espressomd import *
import numpy as np
import sys
import math

##########################################################################
#Set up the System #
##########################################################################
#Build a fluctuating ek species.


@ut.skipIf(not espressomd.gpu_available() or 
           not espressomd.has_features(
    ["ELECTROKINETICS"]),
    "Features or gpu not available, skipping test!")
class ek_fluctuations(ut.TestCase):

    es = espressomd.System(box_l=[1.0, 1.0, 1.0])

    def test(self):
        system = self.es

#Set parameters

        box_x = 16
        box_y = 16
        box_z = 16
        time_step = 0.001
        rho0 = 27.0
        diff = 1.0
        agrid = 1.0

        system.box_l = box_l = [box_x, box_y, box_z]
        system.time_step = time_step
        system.cell_system.skin = 0.2
        system.thermostat.turn_off()

#Setup the Fluid

        ek = electrokinetics.Electrokinetics(agrid=agrid,
                                             lb_density=1.0,
                                             viscosity=1.0,
                                             friction=0.0,
                                             T=1.0,
                                             prefactor=1.0,
                                             stencil='linkcentered',
                                             advection=False,
                                             fluctuations=True,
                                             fluctuation_amplitude=1.0)

        species = electrokinetics.Species(density=rho0,
                                          D=diff,
                                          valency=0.0)

        ek.add_species(species)
        system.actors.add(ek)

#Warm - Up

        system.integrator.run(1000)

#Set integration and binning parameters

        n_min = 10.0
        n_max = 44.0
        bin_size = 0.25
        bins = np.zeros(int((n_max - n_min) / bin_size))
        x_range = np.linspace(n_min, n_max, int((n_max - n_min) / bin_size))
        sample_steps = 100
        integration_steps = 200
        count = 0

#Integrate

        for t in range(sample_steps):
            system.integrator.run(integration_steps)
            for i in range(box_x):
                for j in range(box_y):
                    for k in range(box_z):
                        dens = species[i, j, k].density
                        if(dens < n_max and dens > n_min):
                            x = int((dens - n_min) / bin_size)
                            bins[x] += 1
                            count += 1

        bins = bins / count / bin_size

#Analysis

        p = []
        for i in x_range:
            p.append(1.0 / (math.sqrt(2.0 * math.pi * i))
                     * math.pow(rho0 / i, i) * math.exp(i - rho0))

        max_diff = 0.0
        for i in range(len(x_range)):
            max_diff = max(math.fabs(p[i] - bins[i]), max_diff)

        self.assertLess(max_diff, 5.0e-03,
                        "Density distribution accuracy not achieved, allowed deviation: 5.0e-03, measured: {}".format(max_diff))

if __name__ == "__main__":
    ut.main()
