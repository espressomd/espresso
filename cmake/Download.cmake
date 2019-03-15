###########################################################################
#   Copyright 2017 Florian Reiterer
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
###########################################################################

set(download_template [=[
cmake_minimum_required(VERSION 2.8.2)
project(@name@ NONE)
include(ExternalProject)
ExternalProject_Add(@name@
    BINARY_DIR "@download_dir@"
    SOURCE_DIR "@source_dir@"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    TEST_COMMAND ""
    @ARGN@
)
]=])

function(download directory)
    get_filename_component(name "${directory}" NAME)
    message(STATUS "Downloading/updating ${name}")
    get_filename_component(directory "${directory}" ABSOLUTE)

    set(download_dir "${directory}/download")
    set(source_dir "${directory}/src")

    string(CONFIGURE "${download_template}" cmake_lists @ONLY)
    file(WRITE "${download_dir}/CMakeLists.txt" "${cmake_lists}")
    execute_process(COMMAND ${CMAKE_COMMAND} .
                    RESULT_VARIABLE result
                    WORKING_DIRECTORY "${download_dir}"
    )
    if(result)
        message(FATAL_ERROR "CMake step for ${name} failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
                    RESULT_VARIABLE result
                    WORKING_DIRECTORY "${download_dir}"
    )
    if(result)
        message(FATAL_ERROR "Build step for ${name} failed: ${result}")
    endif()
endfunction()