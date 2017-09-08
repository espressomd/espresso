#!/usr/bin/env python

from __future__ import print_function

import espressomd as md
from espressomd.interactions import HarmonicBond
from espressomd import analyze

import numpy as np
import scipy as sp

print("Program Information:")
print(md.features())

system = md.System()

#############################################################
# Calculate the shearstress_bonded

def calculatestress_bonded(n_parts, pos, box_l, wall_displacement):
  r = pos[1] - pos[0]
  f = -1.0*r
  return np.einsum('i,j', f, r) / box_l**3

#############################################################
# Calculate the shearstress_ideal

def calculatestress_ideal(n_parts, pos, vel, box_l, wall_displacement):
  return np.einsum('ij,ik->jk', vel, vel) / box_l**3

#############################################################
# System parameters

box_l        = 1.0
system.box_l =(box_l, box_l, box_l)

skin        = 0.4
time_step   = 0.01

#############################################################
# Thermostat and cell system

system.thermostat.set_langevin(kT=0.0, gamma=1.0)
system.cell_system.skin = skin
system.time_step = time_step
system.periodicity=1,1,1

#############################################################
# Placement of particles

system.part.add(id=0,pos=(0.5,0.75,0.5))
system.part.add(id=1,pos=(0.5,1.25,0.5))
n_parts = len(system.part)

print("Number of particles")
print(n_parts)

#############################################################
# Bonding of particles

harmonic=HarmonicBond(k=1, r_0=0)
system.bonded_inter.add(harmonic)
system.part[0].add_bond((harmonic, 1))

#############################################################
# Calculate forces and stress

system.integrator.run(steps=0)

#############################################################
# Without Lees-Edwards
wall_displacement = 0.0

total_pressure = system.analysis.pressure()['total']
total_stress_tensor = system.analysis.stress_tensor()['total']
ideal_pressure = system.analysis.pressure()['ideal']
ideal_stress_tensor = system.analysis.stress_tensor()['ideal']
bonded_pressure = system.analysis.pressure()['bonded']
bonded_stress_tensor = system.analysis.stress_tensor()['bonded']

pos = system.part[:].pos
vel = system.part[:].v

print("Particle positions")
print(pos)
print("Particle velocities")
print(vel)

print("Total pressure")
print(total_pressure)
print("Trace of stress-tensor")
print(np.matrix.trace(total_stress_tensor)/3)
print("Total stress tensor")
print(total_stress_tensor)
print("Ideal stress tensor")
print(ideal_stress_tensor)
print("Bonded stress tensor")
print(bonded_stress_tensor)

print("Manual stress total case:\n", calculatestress_bonded(n_parts, pos, box_l, wall_displacement) + calculatestress_ideal(n_parts, pos, vel, box_l,wall_displacement))
print("Manual ideal stress: \n", calculatestress_ideal(n_parts, pos, vel, box_l, wall_displacement))
print("Manual bonded stress: \n", calculatestress_bonded(n_parts, pos, box_l, wall_displacement))

#############################################################
# With Lees-Edwards

wall_displacement = 0.3
system.lees_edwards_offset = wall_displacement

print("Lees Edwards offset")

system.integrator.run(steps=1)

system.part[0].v = [1.0,1.0,1.0]
system.part[1].v = [1.0,1.0,1.0]

print(system.lees_edwards_offset)
print(wall_displacement)

total_pressure = system.analysis.pressure()['total']
total_stress_tensor = system.analysis.stress_tensor()['total']
ideal_pressure = system.analysis.pressure()['ideal']
ideal_stress_tensor = system.analysis.stress_tensor()['ideal']
bonded_pressure = system.analysis.pressure()['bonded']                                                                                                     
bonded_stress_tensor = system.analysis.stress_tensor()['bonded']

pos = system.part[:].pos
vel = system.part[:].v

print("Particle positions")
print(pos)
print("Particle velocities")
print(vel)

print("Total pressure:")
print(total_pressure)
print("Trace of stress-tensor")
print(np.matrix.trace(total_stress_tensor)/3)
print("Total stress tensor")
print(total_stress_tensor)
print("Ideal stress tensor")
print(ideal_stress_tensor)
print("Bonded stress tensor")
print(bonded_stress_tensor)

print("Manual stress total case:\n", calculatestress_bonded(n_parts, pos, box_l, wall_displacement) + calculatestress_ideal(n_parts, pos, vel, box_l,wall_displacement))
print("Manual stress ideal case: \n",calculatestress_ideal(n_parts, pos, vel, box_l, wall_displacement))
print("Manual stress bonded case: \n", calculatestress_bonded(n_parts, pos, box_l, wall_displacement))

#Compare manual and automatic calculation

#############################################################
# End of test
exit
