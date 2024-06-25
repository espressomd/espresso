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

include(espresso_resource_files)

add_custom_target(object_in_fluid)
file(GLOB python_SRC RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.py)

espresso_target_resources(object_in_fluid ${python_SRC})

# Install resource files (Python files, text files, etc.)
get_property(OIF_RESOURCE_FILES TARGET object_in_fluid
             PROPERTY EspressoResourceFiles)
foreach(RESOURCE_ABSPATH ${OIF_RESOURCE_FILES})
  cmake_path(RELATIVE_PATH RESOURCE_ABSPATH BASE_DIRECTORY
             ${CMAKE_CURRENT_BINARY_DIR} OUTPUT_VARIABLE RESOURCE_RELPATH)
  cmake_path(GET RESOURCE_RELPATH PARENT_PATH RESOURCE_RELPARENT)
  install(
    FILES "${RESOURCE_ABSPATH}"
    DESTINATION
      "${ESPRESSO_INSTALL_PYTHON}/object_in_fluid/${RESOURCE_RELPARENT}")
endforeach()
