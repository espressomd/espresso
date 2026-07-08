#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 17:00:50 2026

@author: abhinav
"""

import sys
# sys.path.insert(0, "/home/abhinav/espresso_github/build/src/python")

from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]   # /home/abku051h/espresso_int if script is in codes/
sys.path.insert(0, str(ROOT / "build/src/python"))

import numpy as np
import espressomd

# ------------------------------------------------------------
# Two-particle test for anisotropic Gaussian interaction
# ------------------------------------------------------------

system = espressomd.System(box_l=[10.0, 10.0, 10.0])
system.time_step = 0.01
system.cell_system.skin = 0.4

# Clear particles/interactions if needed
system.part.clear()
system.non_bonded_inter.reset()

# Parameters
eps = -2.0
sig_x = 0.5
sig_y = 1.0
sig_z = 2.0
cutoff = 1.0

# Particle positions
pos0 = np.array([1.0, 1.0, 1.0])
pos1 = np.array([1.2, 1.0, 1.0])

d = pos1 - pos0
dx, dy, dz = d
r = np.linalg.norm(d)

# Add particles
p0 = system.part.add(pos=pos0, type=0)
p1 = system.part.add(pos=pos1, type=0)

# Set interaction
system.non_bonded_inter[0, 0].gaussian_aniso.set_params(
    eps=eps,
    sig_x=sig_x,
    sig_y=sig_y,
    sig_z=sig_z,
    cutoff=cutoff,
)

    
print("features:", espressomd.features())
print("params:", system.non_bonded_inter[0, 0].gaussian_aniso.get_params())
print("p0 type:", p0.type, "p1 type:", p1.type)

# Compute forces and energy without integrating
# system.integrator.run(0)
system.integrator.run(0, recalc_forces=True)

# Analytic energy
expected_energy = eps * np.exp(
    -0.5 * (
        dx**2 / sig_x**2 +
        dy**2 / sig_y**2 +
        dz**2 / sig_z**2
    )
)

expected_force_on_0 = -expected_energy * np.array([
    dx / sig_x**2,
    dy / sig_y**2,
    dz / sig_z**2,
])

energy_dict = system.analysis.energy()
measured_energy = energy_dict["non_bonded"]

f0 = np.array(p0.f)
f1 = np.array(p1.f)

print("------------------------------------------------------------")
print("Two-particle anisotropic Gaussian test")
print("------------------------------------------------------------")
print(f"d = {d}")
print(f"r = {r}")
print()
print("Energy:")
print(f"  expected = {expected_energy:.16e}")
print(f"  measured = {measured_energy:.16e}")
print(f"  diff     = {measured_energy - expected_energy:.16e}")
print()
print("Force on particle 0:")
print(f"  expected = {expected_force_on_0}")
print(f"  measured = {f0}")
print(f"  diff     = {f0 - expected_force_on_0}")
print()
print("Force on particle 1:")
print(f"  expected = {-expected_force_on_0}")
print(f"  measured = {f1}")
print(f"  diff     = {f1 + expected_force_on_0}")
print()
print("Newton's third law check:")
print(f"  f0 + f1 = {f0 + f1}")
print()

# Assertions
np.testing.assert_allclose(measured_energy, expected_energy, rtol=1e-12, atol=1e-12)
np.testing.assert_allclose(f0, expected_force_on_0, rtol=1e-12, atol=1e-12)
np.testing.assert_allclose(f1, -expected_force_on_0, rtol=1e-12, atol=1e-12)
np.testing.assert_allclose(f0 + f1, np.zeros(3), rtol=1e-12, atol=1e-12)

print("PASS: energy and force match analytic anisotropic Gaussian.")