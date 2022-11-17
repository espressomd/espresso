#
# Copyright (C) 2016-2022 The ESPResSo project
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

include(unit_test)
unit_test(NAME Wall_test SRC Wall_test.cpp DEPENDS Espresso::shapes
          Espresso::utils)
unit_test(NAME HollowConicalFrustum_test SRC HollowConicalFrustum_test.cpp
          DEPENDS Espresso::shapes Espresso::utils)
unit_test(NAME Union_test SRC Union_test.cpp DEPENDS Espresso::shapes
          Espresso::utils)
unit_test(NAME Ellipsoid_test SRC Ellipsoid_test.cpp DEPENDS Espresso::shapes
          Espresso::utils)
unit_test(NAME Sphere_test SRC Sphere_test.cpp DEPENDS Espresso::shapes
          Espresso::utils)
unit_test(NAME NoWhere_test SRC NoWhere_test.cpp DEPENDS Espresso::shapes
          Espresso::utils)
