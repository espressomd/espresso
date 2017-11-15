  # Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  # Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
  #   Max-Planck-Institute for Polymer Research, Theory Group
  
  # This file is part of ESPResSo.
  
  # ESPResSo is free software: you can redistribute it and/or modify
  # it under the terms of the GNU General Public License as published by
  # the Free Software Foundation, either version 3 of the License, or
  # (at your option) any later version.
  
  # ESPResSo is distributed in the hope that it will be useful,
  # but WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  # GNU General Public License for more details.
  
  # You should have received a copy of the GNU General Public License
  # along with this program.  If not, see <http://www.gnu.org/licenses/>. 

# Build tests only when testing
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL ON)

# Target for the unit tests
add_custom_target(check_unit_tests COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure)

# Run unit tests on check
add_dependencies(check check_unit_tests)

# unit_test function
function(UNIT_TEST)
  cmake_parse_arguments(TEST "" "NAME;NUM_PROC" "SRC" ${ARGN} )
  add_executable(${TEST_NAME} ${TEST_SRC})
  target_link_libraries(${TEST_NAME} ${LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  add_dependencies(${TEST_NAME} EspressoConfig)

  # If NUM_PROC is given, set up MPI parallel test case
  if( TEST_NUM_PROC )
    if(DEFINED TEST_NP)
      if(${TEST_NUM_PROC} GREATER ${TEST_NP})
        set(TEST_NUM_PROC ${TEST_NP})
      endif()
    endif()

    add_test(${TEST_NAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${TEST_NUM_PROC} ${TEST_NAME})
  else( )
    add_test(${TEST_NAME} ${TEST_NAME})
  endif( )

  add_dependencies(check_unit_tests ${TEST_NAME})
endfunction(UNIT_TEST)

# Add tests here
unit_test(NAME Vector_test SRC Vector_test.cpp)
unit_test(NAME RuntimeError_test SRC RuntimeError_test.cpp)
unit_test(NAME RunningAverage_test SRC RunningAverage_test.cpp)
unit_test(NAME RuntimeErrorCollector_test SRC RuntimeErrorCollector_test.cpp ../RuntimeErrorCollector.cpp ../RuntimeError.cpp)

unit_test(NAME ScriptInterface_test SRC ScriptInterface_test.cpp)
unit_test(NAME Wall_test SRC Wall_test.cpp ../shapes/Wall.cpp)
unit_test(NAME Factory_test SRC Factory_test.cpp)
unit_test(NAME NumeratedContainer_test SRC NumeratedContainer_test.cpp)

unit_test(NAME MpiCallbacks_test SRC MpiCallbacks_test.cpp ../MpiCallbacks.cpp NUM_PROC 2)

unit_test(NAME ParallelScriptInterface_test SRC ParallelScriptInterface_test.cpp
  ../MpiCallbacks.cpp ../../script_interface/ScriptInterfaceBase.cpp ../../script_interface/ParallelScriptInterface.cpp
  ../../script_interface/ParallelScriptInterfaceSlave.cpp ../../script_interface/Variant.cpp NUM_PROC 2)

unit_test(NAME gather_buffer_test SRC gather_buffer_test.cpp NUM_PROC 4)
unit_test(NAME scatter_buffer_test SRC scatter_buffer_test.cpp NUM_PROC 4)
unit_test(NAME AutoParameters_test SRC AutoParameters_test.cpp)
unit_test(NAME make_function_test SRC make_function_test.cpp)
unit_test(NAME AutoParameter_test SRC AutoParameter_test.cpp ../../script_interface/Variant.cpp)
unit_test(NAME Variant_test SRC Variant_test.cpp ../../script_interface/Variant.cpp)
unit_test(NAME ParticleIterator_test SRC ParticleIterator_test.cpp)
unit_test(NAME Range_test SRC Range_test.cpp)
unit_test(NAME link_cell_test SRC link_cell_test.cpp)
unit_test(NAME verlet_ia_test SRC verlet_ia_test.cpp)
unit_test(NAME Batch_test SRC Batch_test.cpp)
unit_test(NAME ParticleCache_test SRC ParticleCache_test.cpp ../MpiCallbacks.cpp NUM_PROC 2)
unit_test(NAME SkipIterator_test SRC SkipIterator_test.cpp)
unit_test(NAME List_test SRC List_test.cpp)
unit_test(NAME Particle_test SRC Particle_test.cpp)
unit_test(NAME for_each_pair_test SRC for_each_pair_test.cpp)
unit_test(NAME all_compare_test SRC all_compare_test.cpp NUM_PROC 3)
unit_test(NAME None_test SRC None_test.cpp)
unit_test(NAME keys_test SRC keys_test.cpp)
