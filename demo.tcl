#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ "$1" != "" ]; then NP=$1; else NP=2; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $PLATFORM/tcl_md $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $PLATFORM/tcl_md $0 $* -procs $NP
# Linux \
    else export EF_ALLOW_MALLOC_0=1; lamboot; exec mpirun -np $NP -nsigs $PLATFORM/tcl_md $0 $*;
# \
    fi;

#############################################################
#                                                           #
#  Setup a Polyelectrolyte solution                         #
#                                                           #
#                                                           #
#  Created:       23.01.2003 by HL                          #
#  Last modified: 23.01.2003 by HL                          #
#                                                           #
#############################################################

#############################################################
#  Setup GUI                                                #
#############################################################

# tagging for demonstration
proc codetag {tag} {
    toplevel .tag
    button .tag.btn -text "$tag reached!" -command "destroy .tag"
    pack .tag.btn
    tkwait window .tag
    update
}

# simulation start/stop for demonstration
set simulation_start 0
proc simstart {} {
    global simulation_start
    .sim conf -text stop -command simstop
    set simulation_start 1
}
proc simstop {} {
    global simulation_start
    .sim conf -text integrate -command simstart
    set simulation_start 0
}
button .sim -text "integrate" -command simstart
pack .sim

# buttons
proc mdparam {var min max init text cmd} {
    global $var
    eval "set $var $init"
    label .lab_$var -text $text
    scale .sc_$var -orient h -from $min -to $max -resolution [expr ($max - $min)/10000] -command $cmd
    .sc_$var set $init
    pack .lab_$var .sc_$var
}

proc set_time_step {ts} {
    setmd time_step [expr pow(10, $ts)]
}

mdparam time_step -6 -2 -6 "log time step" set_time_step

#############################################################
#  Parameters                                               #
#############################################################

set vmd_output "yes"

# System parameters
#############################################################

# Number of polymers
set n_polymers 4
# length of polymers
set l_polymers 49
# Bond length
set bond_length 1.0
# distance between charged monomers
set cm_distance 2
# valency of charged monomers
set cm_valency 1.0

# valency of counterions
set ci_valency -1.0

# system density
set density 0.001

# Interaction parameters
#############################################################

set n_part_types 3

set fene1_cut   2.0
set fene1_k     7.0

setmd skin      0.3

# Integration parameters
#############################################################

set intsteps 100
set maxtime 1000

setmd gamma 1.0
setmd temp 1.0

# Other parameters
#############################################################
set tcl_precision 5
set pi 3.141592653589793

#############################################################
#  Setup System                                             #
#############################################################

puts "number of Polymers       : $n_polymers"
set n_monomers [expr ($n_polymers*$l_polymers)]
puts "number of monomers       : $n_monomers"
set polymer_charge [expr (($l_polymers+$cm_distance-1)/$cm_distance)]
puts "polymer_charge           : $polymer_charge"
set n_counterions [expr -int($polymer_charge*$n_polymers/$ci_valency)]
puts "n_counterions            : $n_counterions"
set npart [expr $n_counterions+$n_monomers]
puts "total number of particles: $npart"
set box_length [expr pow(($npart/$density),(1.0/3.0))]
puts "cubic box length         : $box_length"
setmd box_l $box_length $box_length $box_length 

# Interaction setup
#############################################################
codetag "Interaction setup"

# fene
inter 1 fene $fene1_k $fene1_cut

# pairwise lennard_jones for all particles
proc setup_ia {what val} {
    global lj1_epsilon lj1_cut n_part_types
    eval "set $what $val"
    for {set ia1 0} { $ia1 <= $n_part_types } { incr ia1 } {
	for {set ia2 0} { $ia2 <= $n_part_types } { incr ia2 } {
	    inter $ia1 $ia2 lennard-jones $lj1_epsilon 1 $lj1_cut 0 0
	}
    }
}
mdparam lj1_cut 1 2.5 1.12246 "epsilon for LJ" {setup_ia lj1_cut}
mdparam lj1_epsilon 0.5 2 1.0 "epsilon for LJ" {setup_ia lj1_epsilon}

# LJ force capping
proc setup_force_cap {lfc} {
    if {$lfc == 200} { set lfc 0 }
    setmd lj_force_cap $lfc
}
mdparam lj_force_cap 20 200 20 "LJ force cap" setup_force_cap

# electrostatics
proc setup_p3m {bj} {
    if { $bj == 0} {
	setmd bjerrum 0
    } {
	setmd bjerrum $bj
	setmd p3m_alpha 0.1138
	setmd p3m_r_cut 15.0
	setmd p3m_mesh 8 8 8
	setmd p3m_cao 5 5000
	setmd p3m_epsilon 0.1
	setmd p3m_mesh_off 0.5 0.5 0.5
    }
}
mdparam bjerrum 0 1 0 "Bjerrum length" setup_p3m

# Particle setup
#############################################################

# check existing warmup configuration
if { [file exists "pe_initial.gz"]} {
    # read configuration file
    codetag "Reading particles from file"
    set f [open "|gzip -cd pe_initial.gz" r]
    while {[blockfile $f read auto] != "eof" } {}
    close $f

    puts "Found file pe_initial.gz with:"
    puts "n_part = [part number]"
    puts "grid   = \{[setmd node_grid]\}"
    puts "box    = \{[setmd box]\}" 
    if { $npart != [part number] } {
	error "Configuration differs from this setup!"
    }
} {
    # create new start configuration
    codetag "Creating new configuration"
    # Polymer setup:
    set part_id 0
    for {set i 0} { $i < $n_polymers } {incr i} {
	# polymer start ion
	set posx [expr $box_length*[tcl_rand]]
	set posy [expr $box_length*[tcl_rand]]
	set posz [expr $box_length*[tcl_rand]]
	part $part_id pos $posx $posy $posz q $cm_valency type 1
	incr part_id 
	for {set n 1} { $n < $l_polymers } {incr n} {
	    # place next ion within bond radius
	    set phi [expr 2.0*$pi*[tcl_rand]]
	    set theta [expr $pi*[tcl_rand]]
	    set vecx [expr $bond_length*sin($phi)*sin($theta)]
	    set vecy [expr $bond_length*cos($phi)*sin($theta)]
	    set vecz [expr $bond_length*cos($theta)]
	    set posx [expr $posx+$vecx]
	    set posy [expr $posy+$vecy]
	    set posz [expr $posz+$vecz]
	    part $part_id pos $posx $posy $posz \
		bond 1 [expr $part_id-1]

	    if { $n % $cm_distance == 0 } {
		part $part_id q $cm_valency type  1
	    } {
		part $part_id q 0 type 0
	    }
	    incr part_id 
	}
    }

    # Counterion setup
    for {set i 0} { $i < $n_counterions } {incr i} {
	set posx [expr $box_length*[tcl_rand]]
	set posy [expr $box_length*[tcl_rand]]
	set posz [expr $box_length*[tcl_rand]]
	part $part_id pos $posx $posy $posz q $ci_valency type 2
	incr part_id 
    }

    # write
    set f [open "|gzip -c - >pe_initial.gz" w]
    blockfile $f write variable box_l
    blockfile $f write particles "id pos type q" all
    blockfile $f write bonds all
    close $f
    puts "wrote [part number] particles to pe_initial.gz"
}

#############################################################
#  Prepare IMD Connection                                   #
#############################################################
codetag "Setup of IMD"

if { $vmd_output=="yes" } {
    writepsf pe_sol.psf
    writepdb pe_sol.pdb

    for {set port 10000} { $port < 65000 } { incr port } {
	catch {imd connect $port} res
	if {$res == ""} break
    }
    label .imd -text "IMD listening on port $port"
    pack .imd
}

#############################################################
#      Integration                                          #
#############################################################
codetag "Integration"

#analyze
set f [open "pe_obs.dat" w]
puts $f "\#    RE            RG"

label .lab_md -text "mindist: *"
label .lab_RE -text "RE: *"
label .lab_RG -text "RG: *"
pack .lab_md .lab_RE .lab_RG

while { 1 } {
    if {$simulation_start == 1} {
	integrate $intsteps
	if { $vmd_output=="yes" } { puts "imd: [imd positions]" }
	set md [mindist]
	set RE [analyze 0 $n_polymers $l_polymers $n_counterions]
	set RG [analyze 1 $n_polymers $l_polymers $n_counterions]
	.lab_md conf -text "mindist: $md"
	.lab_RE conf -text "RE: $RE"
	.lab_RG conf -text "RG: $RG"
	puts $f "[setmd time]\t$md $RE $RG"
	sdfsfd time
	flush $f
    } {
	imd listen 1
	after 100
    }
    update
}
close $f

if { $vmd_output=="yes" } { puts "imd: [imd disconnect]" }