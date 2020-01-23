#!/usr/bin/env bash
# Copyright (C) 2018-2019 The ESPResSo project
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.

## @file    BashUnitTests.sh
## @brief   Bash unit testing
##
## How to use this library:
##
## 1. Create a file called `test_cmake.sh` containing:
##    @code{.sh}
##    #!/bin/bash
##    # load Unit Testing library
##    source /path/to/espresso/testsuite/cmake/BashUnitTests.sh
##
##    # setup
##    function set_up() {
##      cat >> myconfig.hpp << EOF
##      /* global features */
##      #define ELECTROSTATICS
##      #define EXTERNAL_FORCES
##    EOF
##      cmake .. -DCMAKE_INSTALL_PREFIX=/tmp/espresso
##      make
##      make install
##    }
##
##    # cleanup
##    function tear_down() {
##      rm -rf /tmp/espresso
##      make dist-clean
##    }
##
##    # write tests
##    function test_install() {
##      assert_file_exists bin/main
##      assert_file_exists lib/python3.6/site-packages/espressomd/__init__.py
##    }
##    function test_run() {
##      assert_return_code python3 -c "import espressomd"
##    }
##
##    # run tests
##    run_test_suite
##    @endcode
##
## 2. Add execution permissions: `chmod +x test_cmake.sh`
##
## 3. Run it in a subshell with `(./test_cmake.sh)`. The library always exits
##    from the current shell with a non-zero error code in case of failure,
##    or zero in case of success. This is meant to be used in Makefiles:
##    @code{.sh}
##    test_cmake:
##    	./test_cmake.sh
##    @endcode
##
## 4. Output in case of success:
##    @code{.txt}
##    install ..
##    run .
##    Successfully ran 3 tests
##    @endcode
##
##    Output in case of failure:
##    @code{.txt}
##    install .x
##    run x
##    Found 2 errors in 3 tests
##    @endcode
##
## Safety mechanism: the library will not run if the current environment
## already defines functions named `test_*`, because the @ref run_test_suite
## function loops over all functions starting with prefix `test_`.
##
## Tests run in alphabetical order. If a specific order is required, simply
## name the functions `test_01_python`, `test_02_cpp`, ...
##
## Users may declare new assertions using the following template:
##
## @code{.sh}
## function assert_<type>() {
##   local result=$1
##   local reference=$2
##   local message=$3
##   if [ -z "${message}" ]
##   then
##     message="<generic message>"
##   fi
##   if [ <condition> ]
##   then
##     logSuccess
##   else
##     logFailure "${message}"
##   fi
## }
## @endcode

## @defgroup Private Private functions
##
## Functions that collect and run setup, teardown and unit tests.
## Not meant to be used directly.
## @{

## @brief Print message to stderr
## @param $1 Message to display
stderr() {
  local message=$1
  echo "${message}" 1>&2
}

## @brief Check for possible side effects
##
## If functions named `test_*` are already declared in the current environment,
## stop unit testing now. We don't want to run them and create side effects.
check_namespace() {
  if [ ! -z "$(declare -F | sed -r 's/declare +-f +//' | grep -P '^test_')" ]; then
    stderr 'Functions named test_* already exist:'
    for test in $(declare -F | sed -r 's/declare +-f +//' | grep -P "^test_"); do
      stderr "  ${test}"
    done
    exit 1
  fi
}

# exit now if test_* functions already exist
check_namespace

## @defgroup TryCatch Bash try/catch statement
##
## Bash implementation of a try/catch statement: run a command in a subshell,
## optionally with modifiers, exiting the subshell on first error.
## Returns the subshell error code.
##
## When dispatching unit test scripts with Open MPI, there is a risk that a
## script freezes indefinitely if two calls are made to a try/catch statement.
## For that reason, try/catch functions will exit gracefully if Open MPI
## environmental variables are detected.
## @{


## @brief Detect if the bash script is run by an Open MPI program
##
## Open MPI programs (orterun, mpirun, mpiexec) declare [environmental
## variables](https://www.open-mpi.org/faq/?category=running#mpi-environmental-variables)
## that can easily be checked for.
## @returns 1 if MPI is used, 0 otherwise
detect_open_mpi() {
  if [ -z "${OMPI_COMM_WORLD_SIZE}" ] && [ -z "${OMPI_COMM_WORLD_RANK}" ] \
  && [ -z "${OMPI_COMM_WORLD_LOCAL_RANK}" ] && [ -z "${OMPI_UNIVERSE_SIZE}" ] \
  && [ -z "${OMPI_COMM_WORLD_LOCAL_SIZE}" ] && [ -z "${OMPI_COMM_WORLD_NODE_RANK}" ]; then
    return 0
  else
    error_log+=("Runtime error: cannot run this job from MPI, there is a conflict with")
    error_log+=("${FUNCNAME[1]}() that can cause this job to freeze indefinitely.")
    return 1
  fi
}

## @brief Try/Catch statement in Bash
##
## @param $@ Command to run, possibly with modifiers
## @returns Error code returned by the command
try_catch() {
  detect_open_mpi || return 1
  (
    set -e  # exit from current subshell on first error
    "${@}"
  )
  return ${?}
}

## @brief Try/Catch statement in Bash without output to the terminal
##
## Run a command in a subshell, exiting the subshell on first error.
## @param $@ Command to run, possibly with modifiers
## @returns Error code returned by the command
try_catch_silent() {
  detect_open_mpi || return 1
  (
    set -e  # exit from current subshell on first error
    "${@}" 1>/dev/null 2>/dev/null
  )
  return ${?}
}

## @brief Try/Catch statement in Bash while logging stdout/stderr
##
## Run a command in a subshell, exiting the subshell on first error,
## logging stdout/stderr to the temporary file at #TMPNAME.
## @param $@ Command to run, possibly with modifiers
## @returns Error code returned by the command
try_catch_capture_output() {
  detect_open_mpi || return 1
  rm -f "${TMPNAME}"
  (
    set -e  # exit from current subshell on first error
    "${@}" 2>>"${TMPNAME}" 1>>"${TMPNAME}"
  )
  return ${?}
}

## @}

## @brief Run the set_up() function if it exists
## @returns Error code returned by setUp()
run_set_up() {
  if [ "$(type -t set_up)" = "function" ]; then
    try_catch set_up
    local -r retcode=${?}
    if [ "${retcode}" -ne "0" ]
    then
      stderr "Failed to run set_up()"
      exit ${retcode}
    fi
  fi
  return 0
}

## @brief Run the tear_down() function if it exists
## @returns Error code returned by tear_down()
run_tear_down() {
  rm -f "${TMPNAME}"
  if [ "$(type -t tear_down)" = "function" ]; then
    try_catch tear_down
    local retcode=${?}
    if [ "${retcode}" -ne "0" ]; then
      stderr "Failed to run tear_down()"
      exit ${retcode}
    fi
  fi
  return 0
}

## @brief Run the tests
run_tests() {
  for test in $(declare -F | sed -r 's/declare +-f +//' | grep -P "^test_"); do
    start_test_block "${test#test_}"
    ${test}
    end_test_block
  done
}

## @var error_log
## @brief Array of assertion error messages
error_log=()

## @var total_tests
## @brief Number of tests
total_tests=0

## @var error_counter
## @brief Number of errors
error_counter=0

## @var TMPNAME
## @brief Name of a temporary file where to log Try/Catch output
readonly TMPNAME=$(mktemp -u)

## @brief Start a test
##
## Print the test name and clear @ref error_log
## @param $1 Test name
function start_test_block() {
  local label="${1}"
  echo -n "${label} "
  error_log=()
}

## @brief End a test
##
## Print a newline character and print all error messages in @ref error_log
end_test_block() {
  echo ""
  for (( i=0; i<${#error_log[@]}; i++ )); do
    stderr "    ${error_log[i]}"
  done
}

## @brief Log a successful assertion
log_success() {
  echo -n '.'
  total_tests=$((total_tests + 1))
}

## @brief Log a failed assertion
## @param $* Description of the failure
log_failure() {
  local message="${*}"
  echo -n 'x'
  error_log+=("${message}")
  error_counter=$((error_counter + 1))
  total_tests=$((total_tests + 1))
}

## @}

## @defgroup Public Unit Testing
##
## Functions available to the user to write and run unit tests in bash script.
##
## @copydetails BashUnitTests.sh
##
## @{

## @brief Run the setup, tests, teardown
run_test_suite() {
  run_set_up
  run_tests
  run_tear_down
  if [ "${error_counter}" -ne 0 ]; then
    stderr "Found ${error_counter} errors in ${total_tests} tests"
    exit 1
  else
    echo "Successfully ran ${total_tests} tests"
    exit 0
  fi
}

## @defgroup Assertions Assertion statements
##
## Assert functions used to write unit tests. Checks return values, exit codes,
## files, etc.
## @{

## @brief Check if a file exists
## @param $1 Filepath
## @param $2 Message on failure (optional)
assert_file_exists() {
  local -r filepath="${1}"
  local message="${2}"
  if [ -z "${message}" ]; then
    message="file not found: ${filepath}"
  fi
  if [ -f "${filepath}" ]; then
    log_success
  else
    log_failure "${message}"
  fi
}

## @brief Check if two strings are identical
## @param $1 Obtained result
## @param $2 Expected result
## @param $3 Message on failure (optional)
assert_string_equal() {
  local -r result="${1}"
  local -r expected="${2}"
  local message="${3}"
  if [ -z "${message}" ]; then
    message="${result} != ${expected}"
  fi
  if [ "${result}" = "${expected}" ]; then
    log_success
  else
    log_failure "${message}"
  fi
}

## @brief Check if two integers are identical
## @param $1 Obtained result
## @param $2 Expected result
## @param $3 Message on failure (optional)
assert_equal() {
  local -r result="${1}"
  local -r expected="${2}"
  local message="${3}"
  if [ -z "${message}" ]; then
    message="${result} != ${expected}"
  fi
  if [ "${result}" -eq "${expected}" ]; then
    log_success
  else
    log_failure "${message}"
  fi
}

## @brief Check if a variable is non-zero
## @param $1 Obtained result
## @param $2 Message on failure (optional)
assert_non_zero() {
  local -r result="${1}"
  local message="${2}"
  if [ -z "${message}" ]; then
    message="${result} == 0"
  fi
  if [ "${result}" -ne "0" ]; then
    log_success
  else
    log_failure "${message}"
  fi
}

## @brief Check if a variable is zero
## @param $1 Obtained result
## @param $2 Message on failure (optional)
assert_zero() {
  local -r result="${1}"
  local message="${2}"
  if [ -z "${message}" ]; then
    message="${result} != 0"
  fi
  if [ "${result}" -eq "0" ]; then
    log_success
  else
    log_failure "${message}"
  fi
}

## @brief Check if a return code is zero
##
## Cannot be used in a script run by Open MPI (see \ref TryCatch).
## @param $@ Command to run, possibly with modifiers
assert_return_code() {
  try_catch_capture_output "${@}"
  local -r retcode=${?}
  if [ "${retcode}" -eq "0" ]; then
    log_success
  else
    local -r message="non-zero return code (${retcode}) for command \`$*\`"
    local -r logfile=$(sed 's/^/        /' < "${TMPNAME}")
    log_failure "${message}"$'\n'"${logfile}"
  fi
}

## @}

## @}

