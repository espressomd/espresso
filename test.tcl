#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ $# = 1 ]; then NP=$1; else NP=8; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $PLATFORM/tcl_md $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $PLATFORM/tcl_md $0 $* -procs $NP
# Linux \
    else lamboot; exec /usr/lib/lam/bin/mpirun -np $NP -nsigs $PLATFORM/tcl_md $0 $*;

# \
    fi;

##################################################
# settings
##################################################
set tcl_precision 5
set use_imd n

set write_steps 10
set configs 400

setmd periodic 1 1 1
setmd bjerrum 1.0
setmd box_l 10.0 10.0 10.0
setmd max_num_cells 512
setmd skin 0.4

# external tcl files 
##################################################
set write no
source polywr.tcl

# data initialization
##################################################
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

puts "max_seen_part = [setmd maxpart] (read from config.gz)"
puts "n_part = [part number]"
puts "grid   = \{[setmd node_grid]\}"
puts "n_node = [setmd n_node]"
puts "box    =\{[setmd box]\}"
# setup interactions
##################################################
inter 0 angle 2
inter 1 fene 2 0
inter 0 0 lennard-jones 1 1 2.0 0 0
inter 1 1 lennard-jones 1 1 2.0 0 0
inter 2 2 lennard-jones 1 1 2.0 0 0
inter 0 1 lennard-jones 1 1 2.0 0 0
inter 0 2 lennard-jones 3 1 2.0 0 0
inter 1 2 lennard-jones 2 1 1.2 0 0
puts "nptypes = [setmd nptypes]"
puts "interactions [inter]"
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
