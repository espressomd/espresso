enable_testing()

file(GLOB TESTS *.tcl)
foreach(TEST ${TESTS})
  add_test(${TEST} runtest.sh ${TEST})
endforeach(TEST)
