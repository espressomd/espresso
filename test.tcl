#!/bin/sh
# tricking...\
    PLATFORM=`uname -s`;
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np 8 $PLATFORM/tcl_md $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $PLATFORM/tcl_md $0 $* -procs 8
# Linux \
    else lamboot; exec /usr/lib/lam/bin/mpirun -np 8 -nsigs $PLATFORM/tcl_md $0 $*;

# \
    fi;

##################################################
# settings
##################################################
set tcl_precision 5

# data initialization
##################################################
puts "nproc = [setmd nproc]"

setmd box_l 10.0 8.0 16.0
puts "box =\{[setmd box]\}"
if {[setmd nproc] == 8} {
    setmd procgrid 2 2 2
}
puts "grid = \{[setmd proc]\}"

if { ![file exists "config"]} {
    error "please generate a configuration file with setup.tcl"
    exit
}

# read particles
##################################################

set f [open "|gzip -cd config.gz" r]
readmd $f

# setup interactions
##################################################
inter 0 0 lennard-jones 1 1 1.2 0 0
inter 1 1 lennard-jones 1 1 1.2 0 0
inter 2 2 lennard-jones 1 1 1.2 0 0
inter 0 1 lennard-jones 1 1 1.2 0 0
inter 0 2 lennard-jones 3 1 1.2 0 0
inter 1 2 lennard-jones 2 1 1.2 0 0
puts "nptypes = [setmd nptypes]"

# integration
##################################################
integrate init
integrate 2
integrate exit

# exit
##################################################
puts "finished"
exit
