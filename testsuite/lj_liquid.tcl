#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ "$1" != "" ]; then NP=$1; else NP=2; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $* -procs $NP
# Linux \
    else export EF_ALLOW_MALLOC_0=1; lamboot; exec mpirun -np $NP -nsigs $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $*;
# \
    fi;

#############################################################
#                                                           #
#  Test System 1: Lennard Jones Liquid                      #
#                                                           #
#                                                           #
#  Created:       17.03.2003 by HL                          #
#  Last modified: 18.03.2003 by HL                          #
#                                                           #
#############################################################

puts " "
puts "==================================================="
puts "=                lj_liquid.tcl                    ="
puts "==================================================="
puts " "

puts "Program Information: \n[code_info]\n"

#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "lj_liquid"
set ident "_t1"

# 
set vmd_output "yes"

# System parameters
#############################################################

set box_l   10.0
set density 0.6

# Interaction parameters (repulsive Lennard Jones)
#############################################################

set lj1_eps     1.0
set lj1_sig     1.0
set lj1_cut     1.12246
set lj1_shift   $lj1_eps

# Integration parameters
#############################################################

setmd time_step 0.01
setmd skin      0.4
setmd gamma     1.0
setmd temp      1.0

# warmup integration (with capped LJ potential)
set warm_steps   200
set warm_n_times 30
# do the warmup until the particles have at least the distance min__dist
set min_dist     0.9

# integration
set int_steps    200
set int_n_times  100

# Other parameters
#############################################################
set tcl_precision 6

#############################################################
#  Setup System                                             #
#############################################################

# Interaction setup
#############################################################

setmd box_l $box_l $box_l $box_l

inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift 0

# Particle setup
#############################################################

set volume [expr $box_l*$box_l*$box_l]
set n_part [expr floor($volume*$density)]

for {set i 0} { $i < $n_part } {incr i} {
    set posx [expr $box_l*[t_random]]
    set posy [expr $box_l*[t_random]]
    set posz [expr $box_l*[t_random]]
 
    part $i pos $posx $posy $posz type 0
}

puts "simulate $n_part particles in a cubic simulation box "
puts "[setmd box_l] at density $density"
puts "Interactions:\n[inter]"
set act_min_dist [analyze mindist]
puts "Start with minimal distance $act_min_dist"

#############################################################
#  Warmup Integration                                       #
#############################################################

# prepare vmd connection
if { $vmd_output=="yes" } {
    puts "\nWrite psf and pdb for VMD connection"
    writepsf "$name$ident.psf"
    writepdb "$name$ident.pdb"

    for {set port 10000} { $port < 65000 } { incr port } {
	catch {imd connect $port} res
	if {$res == ""} break
    }
    puts "What you have to do for an VMD connection:"
    puts "1: Start vmd in the current directory (best before you start the script)."
    puts "2: vmd command line: mol load psf $name$ident.psf pdb $name$ident.pdb"
    set HOSTNAME [exec hostname]
    puts "3: vmd command line: imd connect $HOSTNAME $port"
    imd listen 0
}

#open Observable file
set obs_file [open "$name$ident.obs" "w"]

puts "\nStart warmup integration:"
puts "At maximum $warm_n_times times $warm_steps steps"
puts "Stop if minimal distance is larger than $min_dist"

# set LJ cap
set cap 20
inter ljforcecap $cap

set i 0
while { $i < $warm_n_times && $act_min_dist < $min_dist } {
    puts -nonewline "run $i at time=[setmd time] (LJ cap=$cap) "

    integrate $warm_steps
    if { $vmd_output=="yes" } { imd positions }

    set act_min_dist [analyze mindist]
    puts -nonewline "minimal distance = $act_min_dist\r"
    flush stdout
#   write observables
    puts $obs_file "{ time [setmd time] } [analyze energy]"
#   Increase LJ cap
    set cap [expr $cap+10]
    inter ljforcecap $cap
    incr i
}

puts "\nMinimal distance [analyze mindist]"

# Just to see what else we may get from the c code
puts "\nro variables:"
puts "cell_grid     [setmd cell_grid]" 
puts "cell_size     [setmd cell_size]" 
puts "local_box_l   [setmd local_box_l]" 
puts "max_cut       [setmd max_cut]" 
puts "max_part      [setmd max_part]" 
puts "max_range     [setmd max_range]" 
puts "max_skin      [setmd max_skin]" 
puts "n_nodes       [setmd n_nodes]" 
puts "n_part        [setmd n_part]" 
puts "n_part_types  [setmd n_part_types]" 
puts "periodicity   [setmd periodicity]" 
puts "transfer_rate [setmd transfer_rate]" 
puts "verlet_reuse  [setmd verlet_reuse]" 

# write parameter file
polyBlockWrite "$name$ident.set" {box_l time_step skin temp gamma } "" 

#############################################################
#      Integration                                          #
#############################################################

inter ljforcecap 0

puts "\nStart integration: run $int_n_times times $int_steps steps"

# write start configuration
polyBlockWrite "$name$ident.start" {time box_l} {id pos type}

set j 0
for {set i 0} { $i < $int_n_times } { incr i} {
    puts -nonewline "run $i at time=[setmd time] Temp= [expr [analyze energy kin]/$n_part/1.5]\r"
    flush stdout

    integrate $int_steps
    if { $vmd_output=="yes" } { imd positions }

#   write observables
    puts $obs_file "{ time [setmd time] } [analyze energy]"
#   write intermediate configuration
    if { $i%5==0 } {
	polyBlockWrite "$name$ident.[format %04d $j]" {time box_l} {id pos type}
	incr j
    }
}

# write end configuration
polyBlockWrite "$name$ident.end" {time box_l} {id pos type}

close $obs_file

puts "\n\nFinished"