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

define_property(TARGET PROPERTY EspressoResourceFiles
                BRIEF_DOCS "List of resource files to be deployed with target")

# Register resource files (Python files, text files, etc.) that need to be
# deployed alongside a target. If the file exists in the project source
# directory, it is configured with COPYONLY. If not, it is assumed to be a
# generated file.
function(espresso_target_resources)
  list(POP_FRONT ARGV TARGET_NAME)
  foreach(RESOURCE_RELPATH ${ARGV})
    if(IS_ABSOLUTE ${RESOURCE_RELPATH})
      message(
        FATAL_ERROR
          "function espresso_target_resources() only supports relative paths, could not process \"${RESOURCE_RELPATH}\""
      )
    endif()
    set(RESOURCE_SOURCE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${RESOURCE_RELPATH}")
    set(RESOURCE_BINARY_PATH "${CMAKE_CURRENT_BINARY_DIR}/${RESOURCE_RELPATH}")
    if(EXISTS ${RESOURCE_SOURCE_PATH})
      configure_file(${RESOURCE_SOURCE_PATH} ${RESOURCE_BINARY_PATH} COPYONLY)
    endif()
    set_property(TARGET ${TARGET_NAME} APPEND
                 PROPERTY EspressoResourceFiles "${RESOURCE_BINARY_PATH}")
  endforeach()
endfunction()
