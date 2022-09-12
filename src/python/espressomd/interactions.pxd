#
# Copyright (C) 2013-2022 The ESPResSo project
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

# Handling of interactions

from libc cimport stdint

from .thermostat cimport thermalized_bond

cdef extern from "script_interface/interactions/bonded.hpp":
    int bonded_ia_params_zero_based_type(int bond_id) except +

# Map the boost::variant type indices to python type identifiers. These enum
# values must be in the same order as in the definition of the boost::variant.
cdef enum enum_bonded_interaction:
    BONDED_IA_NONE = 0,
    BONDED_IA_FENE,
    BONDED_IA_HARMONIC,
    BONDED_IA_QUARTIC,
    BONDED_IA_BONDED_COULOMB,
    BONDED_IA_BONDED_COULOMB_SR,
    BONDED_IA_ANGLE_HARMONIC,
    BONDED_IA_ANGLE_COSINE,
    BONDED_IA_ANGLE_COSSQUARE,
    BONDED_IA_DIHEDRAL,
    BONDED_IA_TABULATED_DISTANCE,
    BONDED_IA_TABULATED_ANGLE,
    BONDED_IA_TABULATED_DIHEDRAL,
    BONDED_IA_THERMALIZED_DIST,
    BONDED_IA_RIGID_BOND,
    BONDED_IA_IBM_TRIEL,
    BONDED_IA_IBM_VOLUME_CONSERVATION,
    BONDED_IA_IBM_TRIBEND,
    BONDED_IA_OIF_GLOBAL_FORCES,
    BONDED_IA_OIF_LOCAL_FORCES,
    BONDED_IA_VIRTUAL_BOND

cdef extern from "thermostat.hpp":
    void thermalized_bond_set_rng_counter(stdint.uint64_t counter)

cdef extern from "immersed_boundary/ImmersedBoundaries.hpp":
    cppclass ImmersedBoundaries:
        double get_current_volume(int softID)

cdef extern from "immersed_boundaries.hpp":
    extern ImmersedBoundaries immersed_boundaries
