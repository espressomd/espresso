#
# Copyright (C) 2024 The ESPResSo project
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

import dataclasses
import typing
import ase
from ase.calculators.singlepoint import SinglePointCalculator
import numpy as np
if typing.TYPE_CHECKING:
    from espressomd.system import System


@dataclasses.dataclass
class ASEInterface:
    """
    ASE interface for ESPResSo.
    """

    _system: typing.Union["System", None] = None

    def register_system(self, system):
        """Register the system."""
        self._system = system

    def get(self, folded=False) -> ase.Atoms:
        """Export the ESPResSo system particle data to an ASE atoms object."""
        particles = self._system.part.all()
        positions = np.copy(particles.pos_folded if folded else particles.pos)
        types = np.copy(particles.type)
        forces = np.copy(particles.f)
        if any(p.is_virtual() for p in particles):
            raise RuntimeError("ASE doesn't support virtual sites")

        atoms = ase.Atoms(
            positions=positions,
            numbers=types,
            pbc=np.copy(self._system.periodicity),
            cell=np.copy(self._system.box_l),
        )
        atoms.calc = SinglePointCalculator(atoms, forces=forces)
        return atoms
