#
# Copyright (C) 2013,2014 The ESPResSo project
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
from espressomd import code_info
from espressomd import thermostat
from espressomd import analyze
from espressomd import integrate
from espressomd import interactions
from espressomd import electrostatics
import sys
import numpy as np
import cPickle as pickle
import os

print(code_info.features())

# Seed
#############################################################
np.random.seed(42)

# System parameters
#############################################################

system = espressomd.System()
system.time_step = 0.01
system.skin = 0.4
system.box_l = [100, 100, 100]
system.periodic = [1,1,1]
system.thermostat.set_langevin(kT=1.0, gamma=1.0)
# system.cell_system.set_n_square(use_verlet_lists=False)
system.max_num_cells = 2744

# Non-bonded interactions
###############################################################
# WCA between monomers
system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=1, sigma=1,
    cutoff=2**(1. / 6), shift="auto")

# WCA counterions - polymer
system.non_bonded_inter[0, 1].lennard_jones.set_params(
    epsilon=1, sigma=1,
    cutoff=2**(1. / 6), shift="auto")

# WCA coions - polymer
system.non_bonded_inter[0, 2].lennard_jones.set_params(
    epsilon=1, sigma=1,
    cutoff=2**(1. / 6), shift="auto")

# WCA between ions
system.non_bonded_inter[1, 2].lennard_jones.set_params(
    epsilon=1, sigma=1,
    cutoff=2**(1. / 6), shift="auto")


# Bonded interactions
################################################################
# fene = interactions.FeneBond(k=10, d_r_max=2)
# system.bonded_inter.add(fene)
harmonic = interactions.HarmonicBond(k=10, r_0=2)
harmonicangle = interactions.Angle_Harmonic(bend=10, phi0=np.pi)
system.bonded_inter.add(harmonic)
system.bonded_inter.add(harmonicangle)


# Create Monomer beads and bonds
#########################################################################################
n_monomers = 20

init_polymer_pos=np.dstack((np.arange(n_monomers),np.zeros(n_monomers),np.zeros(n_monomers)))[0]+np.array([system.box_l[0]/2-n_monomers/2, system.box_l[1]/2, system.box_l[2]/2])

system.part.add(id=np.arange(n_monomers), pos=init_polymer_pos)


system.part[:-1].add_bond((harmonic, np.arange(n_monomers)[1:]))
system.part[1:-1].add_bond((harmonicangle, np.arange(n_monomers)[:-2], np.arange(n_monomers)[2:]))

# Particle creation with loops:
# for i in range(n_monomers):
#     if i > 0:
#         system.part[i].add_bond((harmonic, i - 1))

# for i in range(1,n_monomers-1):
#     system.part[i].add_bond((harmonicangle,i - 1, i + 1))

system.part[:n_monomers].q = -np.ones(n_monomers)

# Create counterions
###################################################################
system.part.add(pos=np.random.random((n_monomers,3)) * system.box_l,
                q=np.ones(n_monomers),
                type=np.ones(n_monomers))

# Create ions
###############################################################
n_ions = 100

system.part.add(pos=np.random.random((n_ions,3)) * system.box_l, 
                q=np.hstack((np.ones(n_ions/2),-np.ones(n_ions/2))),
                type=np.hstack((np.ones(n_ions/2),2*np.ones(n_ions/2))))


# Sign charges to particles after the particle creation:
# system.part[2*n_monomers:2*n_monomers+n_ions/2] = np.ones(n_ions/2)
# system.part[2*n_monomers+n_ions/2:] = -np.ones(n_ions/2)

print("types:", system.part[:].type)
print("")
print("Q_tot:", np.sum(system.part[:].q))



#############################################################
#      Warmup                                               #
#############################################################

system.non_bonded_inter.set_force_cap(10)

for i in range(1000):
    sys.stdout.write("\rWarmup: %03i"%i)
    sys.stdout.flush()
    integrate.integrate(1)
    system.non_bonded_inter.set_force_cap(10*i)

system.non_bonded_inter.set_force_cap(0)

print("\nWarmup finished!\n")

#############################################################
#      Sampling                                             #
#############################################################
#
# Activate electostatic with checkpoint example
#############################################################
read_checkpoint = False
# Load checkpointed p3m class
if os.path.isfile("p3m_checkpoint") and read_checkpoint == True:
    print("reading p3m from file")
    p3m = pickle.load(open("p3m_checkpoint","r"))
else:
    p3m = electrostatics.P3M(bjerrum_length=1.0, accuracy=1e-2)
    print("Tuning P3M")
    
system.actors.add(p3m)

# Checkpoint AFTER tuning (adding method to actors)
pickle.dump(p3m,open("p3m_checkpoint","w"),-1)

print("P3M parameter:\n")
p3m_params = p3m.get_params()
for key in p3m_params.keys():
    print("{} = {}".format(key, p3m_params[key]))

print(system.actors)

# Apply external Force
#############################################################
n_part = system.n_part
system.part[:].ext_force = np.dstack((system.part[:].q * np.ones(n_part), np.zeros(n_part), np.zeros(n_part)))[0]

# print(system.part[:].ext_force)


# Activate LB
############################################################
# lbf = lb.LBF(dens=1, tau=0.01, visc=1, fric=1, agrid=1)
# system.actors.add(lbf)

# Data arrays
v_list = []
pos_list = []

# Sampling Loop
for i in range(4000):
    sys.stdout.write("\rSampling: %04i"%i)
    sys.stdout.flush()
    integrate.integrate(1)

    v_list.append(system.part[:n_monomers].v)
    pos_list.append(system.part[:n_monomers].pos)
    # other observales:
    # energies = analyze.energy(system=system)
    # linear_momentum = analyze.analyze_linear_momentum(system=system)

print("\nSampling finished!\n")

# Data evaluation
############################################################
# Convert data to numpy arrays
# shape = [time_step, monomer, coordinate]!
v_list = np.array(v_list)
pos_list = np.array(pos_list)

# Calculate COM and COM velocity
COM = pos_list.sum(axis=1)/n_monomers
COM_v = (COM[1:] - COM[:-1])/system.time_step

# Calculate the Mobility mu = v/E
##################################
mu = COM_v.mean()/1.0
print("MOBILITY", mu)

# Calculate the Persistence length
# fits better for longer sampling
##################################
# this calculation method requires
# numpy 1.10 or higher
if float(np.version.version.split(".")[1]) >= 10:
    from scipy.optimize import curve_fit
    from numpy.linalg import norm

    # First get bond vectors
    bond_vec = pos_list[:,1:,:] - pos_list[:,:-1,:]
    bond_abs = norm(bond_vec, axis=2, keepdims=True)
    bond_abs_avg = bond_abs.mean(axis=0)[:,0]

    c_length = bond_abs_avg
    for i in range(1,len(bond_abs_avg)):
        c_length[i] += c_length[i-1]

    bv_norm = bond_vec / bond_abs


    bv_zero = np.empty_like(bv_norm)
    for i in range(bv_zero.shape[1]):
        bv_zero[:,i,:] = bv_norm[:,0,:]
    
    # Calculate <cos(theta)>
    cos_theta = (bv_zero*bv_norm).sum(axis=2).mean(axis=0)

    def decay(x,lp):
        return np.exp(-x/lp)

    fit,_ = curve_fit(decay, c_length, cos_theta)

    print c_length.shape, cos_theta.shape
    print "PERSISTENCE LENGTH", fit[0]

# Plot Results
############################################################
import matplotlib.pyplot as pp

direction = ["x", "y", "z"]
fig1=pp.figure()
ax=fig1.add_subplot(111)
for i in range(3):
    ax.plot(COM[:-500,i], label="COM pos %s" %direction[i])
ax.legend(loc="best")
ax.set_xlabel("time_step")
ax.set_ylabel("r")

fig2=pp.figure()
ax=fig2.add_subplot(111)
for i in range(3):
    ax.plot(COM_v[:-500,i], label="COM v %s" %direction[i])
ax.legend(loc="best")
ax.set_xlabel("time_step")
ax.set_ylabel("v")

if float(np.version.version.split(".")[1]) >= 10:
    fig3=pp.figure()
    ax=fig3.add_subplot(111)
    ax.plot(c_length, cos_theta, label="sim data")
    ax.plot(c_length, decay(c_length, fit[0]), label="fit")
    ax.legend(loc="best")
    ax.set_xlabel("contour length")
    ax.set_ylabel("<cos(theta)>")

pp.show()


print("\nJob finished!\n")
