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
from libcpp cimport bool

from .utils cimport Vector3i, Vector3d

cdef extern from "grid.hpp":
    Vector3i node_grid

    cppclass BoxGeometry:
        void set_periodic(unsigned coord, bool value)
        bool periodic(unsigned coord)
        const Vector3d & length()
        void set_length(Vector3d)

    BoxGeometry box_geo

    Vector3d get_mi_vector(Vector3d, Vector3d, const BoxGeometry & )
    Vector3d folded_position(Vector3d, const BoxGeometry &)
    Vector3d unfolded_position(Vector3d, Vector3i, const Vector3d & )
