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

setmd bjerrum 1.54
setmd p3m_alpha 0.27
setmd p3m_r_cut 3.0
setmd p3m_mesh 20 20 20
setmd p3m_cao 5
setmd p3m_epsilon 0.1
setmd p3m_mesh_off 1.0 1.0 1.0


# integration
##################################################

puts "imd: [imd connect 12345]"
puts "imd: [imd stall 1]"

integrate init
set write_steps 2
set configs 2
for {set i 0} { $i < $configs } { incr i } {
    puts "step [expr $i*$write_steps]"
    integrate $write_steps
    if {"$write" == "yes" } {
	polywrite [format "configs/t%04d.poly" $i]
    }
    #puts "imd: [imd pos]"
}

integrate exit

# puts "imd: [imd disconnect]"

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
