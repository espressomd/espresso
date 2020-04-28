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
# Interface to the scafacos library. These are the methods shared between
# dipolar and electrostatics methods

include "myconfig.pxi"

from libcpp.string cimport string
from libcpp cimport bool
from libcpp.list cimport list
IF SCAFACOS:
    cdef extern from "electrostatics_magnetostatics/scafacos.hpp" namespace "Scafacos":
        void set_parameters(string & method_name, string & params, bool dipolar)
        string get_method_and_parameters()
        list[string] available_methods_core "Scafacos::available_methods" ()
        void free_handle()
