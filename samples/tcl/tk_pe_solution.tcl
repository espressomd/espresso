#############################################################
#                                                           #
#  Setup a Polyelectrolyte solution                         #
#                                                           #
#############################################################
#
# Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
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


#############################################################
#  Setup GUI                                                #
#############################################################

# tagging for demonstration
proc codetag {tag} {
    toplevel .tag
    button .tag.btn -text "$tag reached!" -command "destroy .tag"
    pack .tag.btn
    wm geometry .tag +60+0
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

    frame .fr_$var -borderwidth 4 -relief raised

    label .fr_$var.lab -text $text

    scale .fr_$var.sc -orient h -from $min -to $max \
	-resolution [expr ($max - $min)/10000] -command $cmd
    .fr_$var.sc set $init

    pack .fr_$var.lab .fr_$var.sc -in .fr_$var

    pack .fr_$var -fill x
}

proc set_time_step {ts} {
    setmd time_step [expr pow(10, $ts)]
}

mdparam time_step -4 -1.7 -4 "log time step" set_time_step

#############################################################
#  Parameters                                               #
#############################################################

set vmd_output "yes"

# System parameters
#############################################################

# Number of polymers
set n_polymers 2
# length of polymers
set l_polymers 100

# Bond length
set bond_length 1.0
# distance between charged monomers
set cm_distance 3
# valency of charged monomers
set cm_valency 1.0

# valency of counterions
set ci_valency -1.0

# system density
set density 0.001

# counterion lj cutoff
set fixed_lj_cut 1.12246
set fixed_lj_eps 1.0

# Interaction parameters
#############################################################

set n_part_types 3

set fene1_cut   2.0
set fene1_k     7.0

setmd skin      0.3

# Integration parameters
#############################################################

set intsteps 20
set maxtime 1000

thermostat langevin 1.0 1.0

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
    global fixed_lj_eps fixed_lj_cut lj1_epsilon lj1_cut
    eval "set $what $val"
    inter 0 0 lennard-jones $lj1_epsilon 1 $lj1_cut 0 0
    inter 0 1 lennard-jones $lj1_epsilon 1 $lj1_cut 0 0
    inter 1 1 lennard-jones $lj1_epsilon 1 $lj1_cut 0 0

    inter 0 2 lennard-jones $fixed_lj_eps 1 $fixed_lj_cut 0 0
    inter 1 2 lennard-jones $fixed_lj_eps 1 $fixed_lj_cut 0 0
    inter 2 2 lennard-jones $fixed_lj_eps 1 $fixed_lj_cut 0 0
}
mdparam lj1_cut 1 2.5 1.12246 "cutoff for LJ" {setup_ia lj1_cut}
mdparam lj1_epsilon 0.5 2 1.0 "epsilon for LJ" {setup_ia lj1_epsilon}

# LJ force capping
proc setup_force_cap {lfc} {
    if {$lfc == 200} { set lfc 0 }
    inter forcecap $lfc
}
mdparam lj_force_cap 20 200 20 "LJ force cap" setup_force_cap

# electrostatics
proc setup_p3m {bj} {
    if { $bj == 0} {
	inter coulomb 0
    } {
	inter coulomb $bj p3m 15.0 8 5 0.118 
    }
}
mdparam bjerrum 0 4 0 "Bjerrum length" setup_p3m

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
    puts "n_part = [setmd n_part]"
    puts "grid   = \{[setmd node_grid]\}"
    puts "box    = \{[setmd box]\}" 
    if { $npart != [setmd n_part] } {
	error "Configuration differs from this setup!"
    }
} {
    # create new start configuration
    codetag "Creating new configuration"
    # Polymer setup:
    set part_id 0
    for {set i 0} { $i < $n_polymers } {incr i} {
	# polymer start ion
	set posz [expr $box_length*[t_random]]
	set posy [expr $box_length*[t_random]]
	set posx [expr $box_length*[t_random]]
	part $part_id pos $posx $posy $posz q $cm_valency type 1
	incr part_id 
	for {set n 1} { $n < $l_polymers } {incr n} {
	    # place next ion within bond radius
	    set phi [expr 2.0*$pi*[t_random]]
	    set theta [expr $pi*[t_random]]
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
	set posx [expr $box_length*[t_random]]
	set posy [expr $box_length*[t_random]]
	set posz [expr $box_length*[t_random]]
	part $part_id pos $posx $posy $posz q $ci_valency type 2
	incr part_id 
    }

    # write
    set f [open "|gzip -c - >pe_initial.gz" w]
    blockfile $f write variable box_l
    blockfile $f write particles "id pos type q" all
    blockfile $f write bonds all
    close $f
    puts "wrote [setmd n_part] particles to pe_initial.gz"
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
    analyze set chains 0 $n_polymers $l_polymers
    if {$simulation_start == 1} {
	integrate $intsteps
	if { $vmd_output=="yes" } { puts "imd: [imd positions]" }
	set md [analyze mindist]
	set RE [lindex [analyze re] 0]
	set RG [lindex [analyze rg] 0]
	.lab_md conf -text "mindist: $md"
	.lab_RE conf -text "RE: $RE"
	.lab_RG conf -text "RG: $RG"
	puts $f "[setmd time]\t$md $RE $RG"
	flush $f
    } {
	imd listen 1
	after 100
    }
    update
}
close $f

if { $vmd_output=="yes" } { puts "imd: [imd disconnect]" }
