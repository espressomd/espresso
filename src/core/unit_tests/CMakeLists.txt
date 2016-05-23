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
add_custom_target(check_unit_tests COMMAND ${CMAKE_CTEST_COMMAND})

# Run unit tests on check
add_dependencies(check check_unit_tests)

# unit_test function
function(UNIT_TEST TEST_NAME TEST_SRC)
  add_executable(${TEST_NAME} ${TEST_SRC})
  target_link_libraries(${TEST_NAME} ${LIBRARIES})
  add_dependencies(${TEST_NAME} EspressoConfig)
  
  add_test(${TEST_NAME} ${TEST_NAME})
  add_dependencies(check_unit_tests ${TEST_NAME})  
endfunction(UNIT_TEST)

unit_test(Vector_test Vector_test.cpp)
unit_test(RuntimeError_test RuntimeError_test.cpp)
unit_test(RunningAverage_test RunningAverage_test.cpp)

set(RuntimeErrorCollector_test_SRC RuntimeErrorCollector_test.cpp ../RuntimeErrorCollector.cpp ../RuntimeError.cpp)
unit_test(RuntimeErrorCollector_test "${RuntimeErrorCollector_test_SRC}")

# ScriptInterface
add_executable(ScriptInterface_test ScriptInterface_test.cpp)
target_link_libraries(ScriptInterface_test ${LIBRARIES})
add_test(script_interface ScriptInterface_test)
add_dependencies(check_unit_tests ScriptInterface_test)

# Wall
add_executable(Wall_test Wall_test.cpp ../shapes/Wall.cpp)
target_link_libraries(Wall_test ${LIBRARIES})
add_test(wall Wall_test)
add_dependencies(check_unit_tests Wall_test)

# Factory
add_executable(Factory_test Factory_test.cpp)
target_link_libraries(Factory_test ${LIBRARIES})
add_test(factory Factory_test)
add_dependencies(check_unit_tests Factory_test)

