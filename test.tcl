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

# external tcl files 
##################################################
set write yes
source polywr.tcl

# data initialization
##################################################
puts "n_node = [setmd n_node]"

setmd box_l 10.0 10.0 10.0
puts "box =\{[setmd box]\}"
if {[setmd n_node] == 8} {
    setmd node_grid 2 2 2
}
puts "grid = \{[setmd node_grid]\}"

if { ![file exists "config.gz"]} {
    error "please generate a configuration file with setup.tcl"
    exit
}

# read particles
##################################################

set f [open "|gzip -cd config.gz" r]
readmd $f

puts "read [expr [setmd maxpart] + 1] particles from config.gz"

# setup interactions
##################################################
inter 0 0 lennard-jones 1 1 2.0 0 0
inter 1 1 lennard-jones 1 1 2.0 0 0
inter 2 2 lennard-jones 1 1 2.0 0 0
inter 0 1 lennard-jones 1 1 2.0 0 0
inter 0 2 lennard-jones 3 1 2.0 0 0
inter 1 2 lennard-jones 2 1 1.2 0 0
puts "nptypes = [setmd nptypes]"

# integration
##################################################
integrate init
set write_steps 10
set configs 50
for {set i 0} { $i < $configs } { incr i } {
    puts "step [expr $i*$write_steps]"
    integrate $write_steps
    if {"$write" == "yes" } {
	polywrite [format "configs/t%04d.poly" $i]
    }
}

integrate exit

# write
set f [open "|gzip -c - >tconfig.gz" w]
writemd $f posx posy posz type q
close $f

if {"$write" == "yes" } {
    eval exec poly2pdb -poly [glob "configs/t*"] \
	-per -psf "movie/t" >/dev/null 2>&1
}


# exit
##################################################
puts "finished"
exit
