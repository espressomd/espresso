#
# Copyright (C) 2019 The ESPResSo project
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

from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "boost/utility/string_ref.hpp" namespace "boost":
    cdef cppclass string_ref:
        string_ref()
        string_ref(const char * )
        string_ref(const string & )

        bool operator == (const string_ref & )
        bool operator != (const string_ref & )

        string to_string()
        const char * data()
