#
# Copyright (C) 2023 The ESPResSo project
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

import enum


class Propagation(enum.IntFlag):
    """
    Flags for propagation modes. Use the pipe operator to combine multiple
    propagation modes. Not all combinations of propagation modes are allowed.
    Flags for virtual sites are special and instruct the propagator to skip
    integration; virtual sites can still be coupled to thermostats (Langevin,
    lattice-Boltzmann) to apply friction and noise to their forces and torques.
    """
    NONE = 0
    """No propagation."""
    SYSTEM_DEFAULT = 2**0
    """Use the system default propagator for motion and rotation."""
    TRANS_NEWTON = 2**1
    """Velocity-Verlet algorithm that integrates Newton's equations of motion."""
    TRANS_LANGEVIN = 2**2
    """Velocity-Verlet algorithm that integrates Langevin's equations of motion."""
    TRANS_LANGEVIN_NPT = 2**3
    """Velocity-Verlet algorithm that integrates Langevin's equations of motion coupled to a piston."""
    TRANS_VS_RELATIVE = 2**4
    """Algorithm for virtual sites relative motion."""
    TRANS_LB_MOMENTUM_EXCHANGE = 2**5
    """Algorithm for momentum exchange between particles and LB fluid cells."""
    TRANS_LB_TRACER = 2**6
    """Algorithm for LB inertialess tracers advection."""
    TRANS_BROWNIAN = 2**7
    """Euler algorithm that integrates Brownian's equations of motion."""
    TRANS_STOKESIAN = 2**8
    """Euler algorithm that integrates Stoke's equations of motion."""
    ROT_EULER = 2**10
    """Velocity-Verlet algorithm that integrates Euler's equations of rotation."""
    ROT_LANGEVIN = 2**11
    """Velocity-Verlet algorithm that integrates Langevin's equations of rotation."""
    ROT_VS_RELATIVE = 2**12
    """Algorithm for virtual sites relative rotation."""
    ROT_BROWNIAN = 2**13
    """Euler algorithm that integrates Brownian's equations of rotation."""
    ROT_STOKESIAN = 2**14
    """Euler algorithm that integrates Stokes' equations of rotation."""
