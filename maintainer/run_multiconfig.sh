#!/bin/bash

# parallelism levels to check
parallelisms="34core nompi"
# default build dir
builddir=build_multiconfig
# continue after encountering an error
contaftererror=n
# number of processors to build on
build_procs=8

# parse options
########################################################
configure_params=""
while test -n "$1"; do
    case "$1" in
        --myconfigs) myconfigs="$2"; shift ;;
        --continue) contaftererror=y ;;
        --build-dir) builddir="$2"; shift ;;
        --build-procs) build_procs="$2"; shift ;;
        --help) echo "usage: $0 [--myconfigs \"config1 config2...\"] [--continue] [--build-dir <dir>] [--build-procs <number of processors for building>] <configure arguments>" 1>&2 ; exit 1 ;;
        *) configure_params="$configure_params $1" ;;
    esac
    shift
done

# get the topdir
########################################################
topdir="`dirname $0`/.."

if ! test -f "$topdir/src/initialize.cpp"; then
    echo "Cannot determine the source directory." 1>&2
    echo "I am >$0< and think the source directory is >$topdir<." 1>&2
    echo "However, >$topdir/src/initialize.cpp< is missing." 1>&2
    exit 1
fi

# check the build directory
########################################################

if test -d "$builddir"; then
    if ! test -f "$builddir/.multiconfig_stamp"; then
        echo "Refusing to build in existing directory $builddir, which was not created by $0." 1>&2
        echo "Remove the directory first." 1>&2
        exit 1
    fi
else
    mkdir -p "$builddir"
fi
touch "$builddir/.multiconfig_stamp"

# make builddir absolute path
builddir="`cd "$builddir"; pwd -P`"

# get the configs we have
########################################################

if test -z "$myconfigs"; then
    myconfigs=""
    pushd "$topdir/maintainer/jenkins/configs"
    for file in *.hpp; do
        myconfigs="$myconfigs ${file%.hpp}"
    done
    popd
fi

# run the tests
########################################################

logdir="$builddir/log"
mkdir -p "$logdir"

export builddir
export insource=false
export build_procs
export configure_params

log="$logdir/configure"
set -o pipefail
bash "$topdir"/maintainer/jenkins/configure.sh 2>&1 | tee "$log"
# don't put that into if, pipefail doesn't work then, and we don't really want to continue
if test "$?" != 0 ; then
    echo "Could not configure Espresso, see $log and $builddir/config.log for details" 1>&2
    echo "Exiting..." 1>&2
    exit 1
fi

for myconfig in $myconfigs; do
    for parallelism in $parallelisms; do
        log="$logdir/$myconfig-$parallelism"
        export myconfig
        export parallelism
        # cleanup first to force rebuild
        ( cd "$builddir" && make clean )

        bash "$topdir"/maintainer/jenkins/build.sh 2>&1  | tee "$log"
        if test "$?" != 0 ; then
            buildfails="$buildfails $myconfig-$parallelism"
        else
            bash "$topdir"/maintainer/jenkins/check.sh 2>&1  | tee -a "$log"
            if test "$?" != 0 ; then
                checkfails="$checkfails $myconfig-$parallelism"
            else
                # failed tests are reported differently, figure out from the log
                if grep -q "Not all test cases were successful!" "$log"; then
                    checkfails="$checkfails $myconfig-$parallelism"
                fi
            fi
        fi
        # save runtest.log
        mv "$builddir/testsuite/runtest.log" "$log-runtest"
    done
done

if test -n "$buildfails"; then
    echo ""
    echo "The following test configurations could not be compiled, see $logdir/<configuration> for details:"
    for conf in $buildfails; do echo $conf; done
fi
if test -n "$checkfails"; then
    echo ""
    echo "The following test configurations could be compiled, but tests failed, see $logdir/<configuration> and $logdir/<configuration>-runtest for details:"
    for conf in $checkfails; do echo $conf; done
fi