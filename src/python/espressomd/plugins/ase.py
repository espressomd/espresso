#
# Copyright (C) 2024-2025 The ESPResSo project
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

import typing
import ase
import numpy as np
import espressomd
if typing.TYPE_CHECKING:
    from espressomd.system import System
    from espressomd.particle_data import ParticleSlice


class ASEInterface:
    """
    Interface for ASE :cite:`hjorthlarsen17a` with calculator support.

    Parameters
    ----------
    system : :obj:`espressomd.system.System`
        The ESPResSo system object.
    particle_slice : :obj:`espressomd.particle_data.ParticleSlice`
        The particle slice to work on.
    export_charges : :obj:`bool`, optional
        Whether to make particle charges available to ASE.
    export_masses : :obj:`bool`, optional
        Whether to make particle masses available to ASE.
    export_momenta : :obj:`bool`, optional
        Whether to make particle momenta available to ASE.
    export_forces : :obj:`bool`, optional
        Whether to make particle forces available to ASE.
    assume_constant_charges : :obj:`bool`, optional
        Assume that the particles' charges won't change while this instance
        is valid (faster update).
    assume_constant_masses : :obj:`bool`, optional
        Assume that the particles' masses won't change while this instance
        is valid (faster update).
    assume_constant_types : :obj:`bool`, optional
        Assume that the particles' types won't change while this instance
        is valid (faster update).
    use_folded_positions : :obj:`bool`, optional
        If True, use folded positions (particles.pos_folded) which are always
        within the simulation box. If False, use unfolded positions (particles.pos).

    """

    def __init__(self, system: "System", particle_slice: "ParticleSlice",
                 export_charges: bool = False, export_masses: bool = False,
                 export_momenta: bool = False, export_forces: bool = False,
                 assume_constant_charges: bool = False,
                 assume_constant_masses: bool = False,
                 assume_constant_types: bool = False,
                 use_folded_positions: bool = True):

        espressomd.assert_features(["EXTERNAL_FORCES"])
        self._system = system
        self.particle_slice = particle_slice
        self.export_charges = export_charges
        self.export_masses = export_masses
        self.export_momenta = export_momenta
        self.export_forces = export_forces
        self.assume_constant_charges = assume_constant_charges
        self.assume_constant_masses = assume_constant_masses
        self.assume_constant_types = assume_constant_types
        self.use_folded_positions = use_folded_positions
        self.atoms = None

        self.reset()

    def __getstate__(self):
        raise NotImplementedError("ASE plugin doesn't support checkpointing")

    def __setstate__(self, state):
        raise NotImplementedError("ASE plugin doesn't support checkpointing")

    def reset(self):
        """
        Re-create the ASE atoms object using ESPResSo system properties.

        New ASE atom objects are created from the current particle slice
        with positions, types, periodicity and box dimensions.
        If export flags are set, charges, masses, momenta and/or forces 
        are also initialized.
        """
        pos_attr = "pos_folded" if self.use_folded_positions else "pos"
        particles = self.particle_slice
        positions = np.copy(getattr(particles, pos_attr))
        types = np.copy(particles.type)

        # Create new atoms object
        self.atoms = ase.Atoms(
            positions=positions,
            numbers=types,
            pbc=np.copy(self._system.periodicity),
            cell=np.copy(self._system.box_l),
        )

        # Initialize charges if requested (always set on reset)
        if self.export_charges:
            charges = np.copy(particles.q)
            self.atoms.set_initial_charges(charges)

        # Initialize masses if requested (always set on reset)
        if self.export_masses:
            masses = np.copy(particles.mass)
            self.atoms.set_masses(masses)

        # Initialize momenta if requested (always set on reset)
        if self.export_momenta:
            momenta = np.copy(particles.v) * \
                np.copy(particles.mass)[:, np.newaxis]
            self.atoms.set_momenta(momenta)

        # Initialize forces if requested (always set on reset)
        if self.export_forces:
            self.atoms.arrays["forces"] = np.copy(particles.f)

        if not self.assume_constant_types:
            self.atoms.numbers = np.copy(particles.type)

    def set_slice(self, particles):
        """
        Set the slice of particles to work on.
        This results in the re-creation of the ASE atom objects.
        """
        self.particles = particles
        self.reset()

    def update_ase(self):
        """
        Update the arrays in the atom objects based on the desired properties.

        The ``assume_constant_*`` flags determine which properties to update.
        """
        if self.atoms is None:
            raise RuntimeError(
                "atoms object not initialized, call reset() first")

        pos_attr = "pos_folded" if self.use_folded_positions else "pos"
        particles = self.particle_slice

        # Always update positions (pos -> positions)
        self.atoms.positions[:] = np.copy(getattr(particles, pos_attr))

        # Update charges if requested and not assumed constant
        if self.export_charges and not self.assume_constant_charges:
            charges = np.copy(particles.q)
            self.atoms.set_initial_charges(charges)

        # Update masses if requested and not assumed constant
        if self.export_masses and not self.assume_constant_masses:
            masses = np.copy(particles.mass)
            self.atoms.set_masses(masses)

        # Update momenta if requested
        if self.export_momenta:
            momenta = np.copy(particles.v) * \
                np.copy(particles.mass)[:, np.newaxis]
            self.atoms.set_momenta(momenta)

        if self.export_forces:
            forces = np.copy(particles.f)
            self.atoms.arrays["forces"] = forces

        if not self.assume_constant_types:
            self.atoms.numbers = np.copy(particles.type)

    def integrate(self, steps: int, calculator) -> int:
        """
        Integrate the system for the specified number of steps.

        For each step:

        1. Update ASE with current particle data
        2. Get forces from calculator
        3. Set external forces on particles
        4. Run one integration step

        Since the particle property :attr:`~espressomd.particle_data.ParticleHandle.ext_force`
        is updated at every time step, cache invalidation is inevitable,
        and the symplectic Euler method is needed to conserve momentum
        (:class:`espressomd.integrate.SymplecticEuler`).

        Parameters
        ----------
        steps : :obj:`int`
            Number of integration steps to perform.
        calculator : :obj:`ase.calculators.calculator.Calculator`
            ASE calculator to use for force calculations.

        Returns
        -------
        int
            Number of steps actually performed
        """
        if calculator is None:
            raise RuntimeError(
                "No calculator assigned. Set self.calculator before integrating.")
        if self.atoms is None:
            raise RuntimeError(
                "atoms object not initialized, call reset() first")
        if not np.allclose(self._system.box_l, self.atoms.cell.lengths()):
            raise RuntimeError(
                "Box size has changed. Call the reset() method on the ase interface")
        self.atoms.calc = calculator
        for _ in range(steps):
            # Update ASE with current particle data
            self.update_ase()

            # Get forces from ASE calculator
            forces = self.atoms.get_forces()

            # Set external forces on particle slice
            self.particle_slice.ext_force = forces

            # Run one integration step
            self._system.integrator.run(1)

        return steps
