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

    type_mapping: dict
    """
    Mapping of ESPResSo particle types to ASE symbols. E.g. ``{0: "H", 1: "O"}``.
    """
    _system: typing.Union["System", None] = None

    def register_system(self, system):
        """Register the system."""
        self._system = system

    def __getstate__(self):
        return {"type_mapping": self.type_mapping}

    def get(self) -> ase.Atoms:
        """Export the ESPResSo system particle data to an ASE atoms object."""
        particles = self._system.part.all()
        positions = np.copy(particles.pos)
        types = np.copy(particles.type)
        forces = np.copy(particles.f)
        unknown_types = set(types) - set(self.type_mapping)
        if unknown_types:
            raise RuntimeError(
                f"Particle types '{unknown_types}' haven't been registered in the ASE type map"  # nopep8
            )
        if any(p.is_virtual() for p in particles):
            raise RuntimeError("ASE doesn't support virtual sites")

        atoms = ase.Atoms(
            positions=positions,
            symbols=[self.type_mapping[t] for t in types],
            pbc=np.copy(self._system.periodicity),
            cell=np.copy(self._system.box_l),
        )
        atoms.calc = SinglePointCalculator(atoms, forces=forces)
        return atoms
