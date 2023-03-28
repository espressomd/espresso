#
# Copyright (C) 2015-2022 The ESPResSo project
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

set(logo_FILES
    logo.svg logo.pdf logo_500x500.png logo-small.svg logo-small.pdf
    logo_32x32.png logo_48x48.png logo_100x100.png logo-animated-100.gif
    logo-animated-200.gif logo-animated-500.gif)

foreach(lf logo_FILES)
  add_custom_command(
    OUTPUT ${lf}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/${lf} ${CMAKE_CURRENT_BINARY_DIR})
endforeach()
