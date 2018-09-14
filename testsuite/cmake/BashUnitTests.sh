#!/bin/bash
## @file    BashUnitTests.sh
## @brief   Bash unit testing
##
## How to use this library:
##
## 1. Create a file called `test_cmake.sh` containing:
##    @code{.sh}
##    #!/bin/bash
##    # load Unit Testing library
##    source /work/jgrad/espresso-fork/benchmarks/BashUnitTests.sh
##
##    # setup
##    function set_up() {
##      cat >> myconfig.hpp << EOF
##      /* global features */
##      #define PARTIAL_PERIODIC
##      #define ELECTROSTATICS
##      #define EXTERNAL_FORCES
##    EOF
##      make
##      make install DESTDIR="install"
##    }
##
##    # cleanup
##    function tear_down() {
##      rm -rf install
##      make dist-clean
##    }
##
##    # write tests
##    function test_install() {
##      assert_file_exists bin/main
##      assert_file_exists lib/python2.7/site-packages/espressomd/__init__.py
##    }
##    function test_run() {
##      assert_return_code python2 -c "import espressomd"
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
function stderr() {
  local message=$1
  echo "${message}" 1>&2
}

## @brief Check for possible side effects
##
## If functions named `test_*` are already declared in the current environment,
## stop unit testing now. We don't want to run them and create side effects.
function check_namespace() {
  if [ ! -z "$(declare -F | sed -r 's/declare +-f +//' | grep -P '^test_')" ]
  then
    stderr 'Functions named test_* already exist:'
    for test in $(declare -F | sed -r 's/declare +-f +//' | grep -P "^test_")
    do
      stderr "  ${test}"
    done
    exit 1
  fi
}

# exit now if test_* functions already exist
check_namespace

## @brief Try/Catch statement in Bash
##
## Run a command in a subshell, exiting the subshell on first error.
## @param $@ Command to run, possibly with modifiers
## @returns Error code returned by the command
function try_catch() {
  (
    set -e  # exit from current subshell on first error
    "$@"
  )
  return $?
}

## @brief Try/Catch statement in Bash without output to the terminal
##
## Run a command in a subshell, exiting the subshell on first error.
## @param $@ Command to run, possibly with modifiers
## @returns Error code returned by the command
function try_catch_silent() {
  (
    set -e  # exit from current subshell on first error
    "$@" 1>/dev/null 2>/dev/null
  )
  return $?
}

## @brief Run the set_up() function if it exists
## @returns Error code returned by setUp()
function run_set_up() {
  if [ "$(type -t set_up)" = "function" ]
  then
    try_catch set_up
    local retcode=$?
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
function run_tear_down() {
  if [ "$(type -t tear_down)" = "function" ]
  then
    try_catch tear_down
    local retcode=$?
    if [ "${retcode}" -ne "0" ]
    then
      stderr "Failed to run tear_down()"
      exit ${retcode}
    fi
  fi
  return 0
}

## @brief Run the tests
function run_tests() {
  for test in $(declare -F | sed -r 's/declare +-f +//' | grep -P "^test_")
  do
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

## @brief Start a test
##
## Print the test name and clear @ref error_log
## @param $1 Test name
function start_test_block() {
  local label=$1
  echo -n "${label} "
  error_log=()
}

## @brief End a test
##
## Print a newline character and print all error messages in @ref error_log
function end_test_block() {
  echo ""
  for (( i=0; i<${#error_log[@]}; i++ ));
  do
    stderr "    ${error_log[i]}"
  done
}

## @brief Log a successful assertion
function log_success() {
  echo -n '.'
  total_tests=$((total_tests + 1))
}

## @brief Log a failed assertion
## @param $1 Description of the failure
function log_failure() {
  local message=$1
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
function run_test_suite() {
  run_set_up
  run_tests
  run_tear_down
  if [ "${error_counter}" -ne 0 ]
  then
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
function assert_file_exists() {
  local filepath=$1
  local message=$2
  if [ -z "${message}" ]
  then
    message="file not found: ${filepath}"
  fi
  if [ -f "${filepath}" ]
  then
    log_success
  else
    log_failure "${message}"
  fi
}

## @brief Check if two variables are identical
## @param $1 Obtained result
## @param $2 Expected result
## @param $3 Message on failure (optional)
function assert_equal() {
  local result=$1
  local expected=$2
  local message=$3
  if [ -z "${message}" ]
  then
    message="${result} != ${expected}"
  fi
  if [ "${result}" -eq "${expected}" ]
  then
    log_success
  else
    log_failure "${message}"
  fi
}

## @brief Check if a variable is non-zero
## @param $1 Obtained result
## @param $2 Message on failure (optional)
function assert_non_zero() {
  local result=$1
  local message=$2
  if [ -z "${message}" ]
  then
    message="${result} == 0"
  fi
  if [ "${result}" -ne "0" ]
  then
    log_success
  else
    log_failure "${message}"
  fi
}

## @brief Check if a variable is zero
## @param $1 Obtained result
## @param $2 Message on failure (optional)
function assert_zero() {
  local result=$1
  local message=$2
  if [ -z "${message}" ]
  then
    message="${result} != 0"
  fi
  if [ "${result}" -eq "0" ]
  then
    log_success
  else
    log_failure "${message}"
  fi
}

## @brief Check if a return code is zero
## @param $@ Command to run, possibly with modifiers
function assert_return_code() {
  try_catch_silent "$@"
  local retcode=$?
  local message="non-zero return code (${retcode}) for command \`$@\`"
  if [ "${retcode}" -eq "0" ]
  then
    log_success
  else
    log_failure "${message}"
  fi
}

## @}

## @}

