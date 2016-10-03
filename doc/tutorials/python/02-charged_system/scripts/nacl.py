#
# Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
#   Max-Planck-Institute for Polymer Research, Theory Group
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
from espressomd import thermostat
from espressomd import electrostatics
from espressomd import code_info
from espressomd import integrate
import numpy

# Setup system 
n_part = 200
density = 0.7
box_l = (n_part / density)**(1. / 3.)

system = espressomd.System()
system.box_l = [box_l, box_l, box_l]
system.periodicity = [1, 1, 1]
system.time_step = 0.01
system.skin = 0.3

# Place particles
q = 1.
ptype = 0
for i in xrange(n_part):
    q *= -1
    ptype = 1 - ptype
    system.part.add(id=i, 
                    type=ptype, 
                    pos=numpy.random.random(3) * system.box_l, 
                    q=q)


# Thermostat
temp = 1.
gamma = 1.
thermostat.Thermostat().set_langevin(kT=temp, gamma=gamma)


# Lennard-Jones interactions
lj_sig = 1.
lj_eps = 1.
lj_cut = 2.**(1. / 6.)
system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")
system.non_bonded_inter[1, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")
system.non_bonded_inter[1, 1].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")


# Electrostatic interaction
p3m = electrostatics.P3M(bjerrum_length=10.0, accuracy=1e-3)
system.actors.add(p3m)
print("P3M parameter:\n")
p3m_params = p3m.get_params()
for key in p3m_params.keys():
    print("{} = {}".format(key, p3m_params[key]))


# Warmup integration loop
print("\nWarmup")
for cap in xrange(20, 200, 20):
    system.non_bonded_inter.set_force_cap(cap)
    integrate.integrate(100)
system.non_bonded_inter.set_force_cap(0)


# Main integration loop
integ_steps = 200
int_n_times = 20

print("\nIntegration")
for i in xrange(int_n_times):
    temp = system.analysis.energy()['ideal'] / ((3.0 / 2.0) * n_part)
    print("t={0:.1f}, E_total={1:.2f}, E_coulomb={2:.2f}, T={3:.4f}".format(system.time,
                                       system.analysis.energy()['total'],
                                       system.analysis.energy()['coulomb'],
                                       temp))
    integrate.integrate(integ_steps)

    # Interally append particle configuration
    system.analysis.append()


# Calculate the averaged rdfs
rdf_bins = 100
r_min  = 0.0
r_max  = system.box_l[0]/2.0
r,rdf_00 = system.analysis.rdf(rdf_type='<rdf>', 
                            type_list_a=[0],
                            type_list_b=[0], 
                            r_min=r_min,
                            r_max=r_max, 
                            r_bins=rdf_bins)

r,rdf_01 = system.analysis.rdf(rdf_type='<rdf>', 
                            type_list_a=[0],
                            type_list_b=[1], 
                            r_min=r_min,
                            r_max=r_max, 
                            r_bins=rdf_bins)

# Write out the data
rdf_fp = open('rdf.data', 'w')
for i in range(rdf_bins):
    rdf_fp.write("%1.5e %1.5e %1.5e\n" % (r[i], rdf_00[i], rdf_01[i]))
rdf_fp.close()
