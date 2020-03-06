# Copyright (C) 2010-2019 The ESPResSo project
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
from .script_interface import ScriptInterfaceHelper, script_interface_register
from .particle_data import ParticleHandle


class _PairCriterion(ScriptInterfaceHelper):

    def decide(self, p1, p2):
        """Makes a decision based on the two particles specified.

        Parameters
        ----------
        p1, p2 : :obj:`espressomd.particle_data.ParticleHandle` or :obj:`int` containing the particle id.
            Particle pair.
        """
        id1 = None
        id2 = None
        if isinstance(p1, ParticleHandle) and isinstance(p2, ParticleHandle):
            id1 = p1.id
            id2 = p2.id
        elif isinstance(p1, int) and isinstance(p2, int):
            id1 = p1
            id2 = p2
        else:
            raise ValueError(
                "arguments must be instances of int or ParticleHandle")
        return self.call_method("decide", id1=id1, id2=id2)


@script_interface_register
class DistanceCriterion(_PairCriterion):

    """Pair criterion returning true, if particles are closer than a cutoff.
    Periodic boundaries are treated via minimum image convention.

    The following parameters can be passed to the constructor, changed via
    ``set_params()`` and retrieved via ``get_params()``.

    Parameters
    ----------
    cut_off : :obj:`float`
        distance cutoff for the criterion
    """
    _so_name = "PairCriteria::DistanceCriterion"
    _so_creation_policy = "LOCAL"


@script_interface_register
class EnergyCriterion(_PairCriterion):

    """Pair criterion returning true, if the short range energy between the
    particles is superior or equal to the cutoff.

    Be aware that the short range energy contains the short range part of
    dipolar and electrostatic interactions, but not the long range part.

    The following parameters can be passed to the constructor, changed via
    ``set_params()`` and retrieved via ``get_params()``.

    Parameters
    ----------
    cut_off : :obj:`float`
        energy cutoff for the criterion
    """
    _so_name = "PairCriteria::EnergyCriterion"
    _so_creation_policy = "LOCAL"


@script_interface_register
class BondCriterion(_PairCriterion):

    """Pair criterion returning true, if a pair bond of given type exists between them

    The following parameters can be passed to the constructor, changed via
    ``set_params()`` and retrieved via ``get_params()``.

    Parameters
    ----------
    bond_type : :obj:`int`
        numeric type of the bond
    """
    _so_name = "PairCriteria::BondCriterion"
    _so_creation_policy = "LOCAL"
