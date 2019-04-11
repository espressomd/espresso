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

set(collect_snippet [=[
set(_all_targets "" CACHE INTERNAL "")
macro(add_library name)
    set(_all_targets "${_all_targets};${name}" CACHE INTERNAL "")
    _add_library(${name} ${ARGN})
endmacro()

macro(add_executable name)
    set(_all_targets "${_all_targets};${name}" CACHE INTERNAL "")
    _add_executable(${name} ${ARGN})
endmacro()
]=])

set(export_snippet [=[
set(_targets_file "@binary_dir@/_external_targets.cmake")
file(WRITE "${_targets_file}" "")
foreach(target ${_all_targets}) 
    if (TARGET ${target})
        get_target_property(imported ${target} IMPORTED)
        get_target_property(aliased ${target} ALIASED_TARGET)
        get_target_property(target_type ${target} TYPE)
        if ("${target_type}" MATCHES "^(STATIC_LIBRARY|SHARED_LIBRARY|EXECUTABLE)$")
            if(imported)
                get_target_property(location ${target} LOCATION)
                if(location)
                    file(APPEND "${_targets_file}" "add_library(${target} STATIC IMPORTED)\nset_property(TARGET ${target} PROPERTY IMPORTED_LOCATION \"${location}\")\n")
                endif()
            elseif(aliased)
                get_target_property(location ${target} LOCATION)
                if(location)
                    file(APPEND "${_targets_file}" "add_library(@namespace@::${target} STATIC IMPORTED)\nset_property(TARGET @namespace@::${target} PROPERTY IMPORTED_LOCATION \"${location}\")\n")
                endif()
            else()
                list(APPEND export_targets ${target}) 
            endif() 
        endif()
    endif()
endforeach() 
file(APPEND "${_targets_file}" "set(_@namespace@_targets ${export_targets})\n\n")
export(TARGETS ${export_targets} NAMESPACE @namespace@:: APPEND FILE "${_targets_file}")
]=])

function(add_external namespace source_dir)
    get_filename_component(source_dir "${source_dir}" ABSOLUTE)
    if(NOT EXISTS "${source_dir}")
        message(FATAL_ERROR "Source directory \"${source_dir}\" not found")
    endif()
    set(binary_dir "${CMAKE_BINARY_DIR}/external/${namespace}")
    file(MAKE_DIRECTORY "${binary_dir}")

    # patch CMakeLists.txt if not patched or updated
    if(EXISTS "${source_dir}/CMakeLists.txt.hash")
        file(READ "${source_dir}/CMakeLists.txt.hash" should_hash)
        file(SHA1 "${source_dir}/CMakeLists.txt" is_hash)
    endif()
    if(NOT is_hash OR NOT(is_hash STREQUAL should_hash))
        file(READ "${source_dir}/CMakeLists.txt" cmake_list)
        file(WRITE "${source_dir}/CMakeLists.txt" "${collect_snippet}\n${cmake_list}\n\ninclude(_target_export.cmake)\n")
        file(SHA1 "${source_dir}/CMakeLists.txt" hash)
        file(WRITE "${source_dir}/CMakeLists.txt.hash" "${hash}")
    endif()
    string(CONFIGURE "${export_snippet}" target_export @ONLY)
    file(WRITE "${source_dir}/_target_export.cmake" "${target_export}")

    foreach(param IN LISTS ARGN)
        if(param MATCHES "^[A-Za-z_]+=.+$")
            list(APPEND cmake_params "-D${param}")
        elseif(param MATCHES "^-.+$")
            list(APPEND cmake_params "${param}")
        else()
            message(FATAL_ERROR "Invalid parameter \"${param}\"")
        endif()
    endforeach()

    execute_process(COMMAND "${CMAKE_COMMAND}" -G${CMAKE_GENERATOR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -Wno-dev ${cmake_params} "${source_dir}"
                    RESULT_VARIABLE result
                    WORKING_DIRECTORY "${binary_dir}"
    )
    if(result)
        message(FATAL_ERROR "CMake step for ${namespace} failed: ${result}")
    endif()

    include("${binary_dir}/_external_targets.cmake")
    foreach(target IN LISTS _${namespace}_targets)
        if(TARGET ${namespace}::${target})
            get_target_property(location ${namespace}::${target} LOCATION)
            add_custom_command(
                OUTPUT "${location}"
                COMMAND cmake --build . --config ${CMAKE_BUILD_TYPE} --target ${target}
                WORKING_DIRECTORY "${binary_dir}"
                VERBATIM USES_TERMINAL
                COMMENT "Generating ${namespace}::${target}"
            )
            add_custom_target(_${namespace}_${target}_compile
                DEPENDS ${location}
            )
            add_dependencies(${namespace}::${target} _${namespace}_${target}_compile)
        endif()
    endforeach()
endfunction()