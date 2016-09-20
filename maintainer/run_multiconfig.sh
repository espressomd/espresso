#!/bin/bash
# Copyright (C) 2016 The ESPResSo project
# Copyright (C) 2013,2014 Axel Arnold
# Copyright (C) 2015 Florian Weik
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

# parallelism levels to check
parallelisms="34core nompi"
# default build dir
builddir=build_multiconfig
# continue after encountering an error
contaftererror=n
# number of processors to build on
build_procs=8
# compilers to run
compilers="auto"

# parse options
########################################################
my_configure_params=""
while test -n "$1"; do
    case "$1" in
        --myconfigs) myconfigs="$2"; shift ;;
        --continue) contaftererror=y ;;
        --build-dir) builddir="$2"; shift ;;
        --build-procs) build_procs="$2"; shift ;;
	--compilers) compilers="$2"; shift;;
	--parallelisms) parallelisms= "$2"; shift ;;
        --help) echo "usage: $0 [--myconfigs \"config1 config2...\"] [--continue] [--build-dir <dir>] [--build-procs <number of processors for building>] <configure arguments>" 1>&2 ; exit 1 ;;
        *) my_configure_params="$my_configure_params $1" ;;
    esac
    shift
done

# get the topdir
########################################################
topdir="`dirname $0`/.."

sourcefile="$topdir/src/features.def"
if ! test -f $sourcefile; then
    echo "Cannot determine the source directory." 1>&2
    echo "I am >$0< and think the source directory is >$topdir<." 1>&2
    echo "However, >$sourcefile< is missing." 1>&2
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

for compiler in $compilers; do    
    if [ $compiler == "auto" ]; then
	confgiure_params="$my_configure_params"
    else	
	configure_params="CXX=$compiler $my_configure_params"
    fi
    export configure_params

    log="$logdir/configure-$compiler"
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
            log="$logdir/$myconfig-$compiler-$parallelism"
            export myconfig
            export parallelism
            # cleanup first to force rebuild
            ( cd "$builddir" && make clean )

            bash "$topdir"/maintainer/jenkins/build.sh 2>&1  | tee "$log"
            if test "$?" != 0 ; then
		buildfails="$buildfails $myconfig-$compiler-$parallelism"
            else
		bash "$topdir"/maintainer/jenkins/check.sh 2>&1  | tee -a "$log"
		if test "$?" != 0 ; then
                    checkfails="$checkfails $myconfig-$compiler-$parallelism"
		else
                    # failed tests are reported differently, figure out from the log
                    if grep -q "Not all test cases were successful!" "$log"; then
			checkfails="$checkfails $myconfig-$compiler-$parallelism"
                    fi
		fi
            fi
            # save runtest.log
            mv "$builddir/testsuite/runtest.log" "$log-runtest"
	done
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
