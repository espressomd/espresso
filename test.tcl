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
set use_imd y

set write_steps 10
set configs 5000

setmd periodic 1 0 0
setmd bjerrum 0
setmd box_l 7.0 7.0 7.0

# external tcl files 
##################################################
set write no
source polywr.tcl

# data initialization
##################################################
puts "n_node = [setmd n_node]"

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

puts "read particles"

set f [open "|gzip -cd config.gz" r]

puts "opened file"

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

setmd skin 0.200001

setmd p3m_alpha 0.27
setmd p3m_r_cut 3.0
setmd p3m_mesh 8 8 8
setmd p3m_cao 3 5000
setmd p3m_epsilon 0.1
setmd p3m_mesh_off 0.5 0.5 0.5


# integration
##################################################

if { $use_imd == "y" } {
	for {set port 10000} { $port < 65000 } { incr port } {
		catch {imd connect $port} res
		if {$res == ""} break
		puts "imd port $port: result $res"
	}
	puts "opened port $port"


	while { [imd listen 10] != "connected" } {
		puts "wating for vmd to connect..."
	}

	while { [setmd transfer_rate] == 0 } {
		imd listen 10
	}
}

integrate init

set write_steps 2

for {set i 0} { $i < $configs } { incr i } {
    puts "step [expr $i*$write_steps]"
    integrate $write_steps
    if {"$write" == "yes" } {
	polywrite [format "configs/t%04d.poly" $i]
    }
    imd pos
}

integrate exit

puts "imd: [imd disconnect]"

# write
set f [open "|gzip -c - >tconfig.gz" w]
writemd $f posx posy posz type q
close $f

if {"$write" == "yes" } {
    eval exec poly2pdb -poly [glob -nocomplain "configs/t*"] \
	-per -psf "movie/t" >/dev/null 2>&1
}


# exit
##################################################
puts "finished"
exit
