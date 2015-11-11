# Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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
# This script proves wether the bond angle potentials are working
# properly by checking the energy and forces for a simple configuration.
#
########################################
#                                      #
# Test of dihedral bond:               #
#                                      #
# Check equilibration angle            #
#                                      #
########################################

source "tests_common.tcl"

puts "----------------------------------------------"
puts "- Testcase dihedral.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------------"

require_feature THERMOSTAT_IGNORE_NON_VIRTUAL off

# Integration parameters
##########################################
set epsilon 1e-4; #allowed error
thermostat langevin 0.0 1.0
setmd time 1
setmd time_step 0.01
setmd skin 0.4

setmd box 10.0 10.0 10.0
setmd periodic 1 1 1

# create 4 particles
#########################################
part 0 pos 0.0 0.0 0.0	
part 1 pos 1.0 1.0 0.0
part 2 pos 1.0 1.0 2.0		
part 3 pos 2.0 1.0 2.0

#define bond
#########################################
inter 0 dihedral 1 1.0 [expr [PI]/4]
part 1 bond 0 0 2 3

inter 1 harmonic 1.0 1.0
part 0 bond 1 1
part 1 bond 1 2
part 2 bond 1 3

puts [bond_dihedral 0 1 2 3]



# load file for comparison and process data
##########################################
# if { [file exists "dihedral.data"] } {
#     set infile [open "dihedral.data"]
# } else {
#     error_exit "could not find file \"dihedral.data\"" 
# }
# set infile_data [read $infile]

# # initialize i/o and temporary variables
# set in1 0.0
# set in2 0.0
# set in3 0.0
# set i 0
# set j 0
# set k 0

# set data [split $infile_data "\n"]
# foreach line $data {
#     scan $line "%s%f%f%f" name in1 in2 in3
#     if { ![ string compare $name "angle:"]  } {
# 	set angles_ref($i) $in1
# 	incr i 1
#     } 
#     if { ![ string compare $name "energy:"] } {
# 	set energies_ref($j) $in1
# 	incr j 1 
#     } 
#     if { ![ string compare $name "force:"]  } {
# 	set fx_ref($k) $in1; set fy_ref($k) $in2; set fz_ref($k) $in3
# 	incr k 1 
#     } 
# }  


#integrate for testing potential
##########################################
set asteps 10.0

for { set i 0 } { $i <= $asteps } { incr i } {

    integrate 1000

    # analyze
    set angles($i)   [ bond_dihedral 0 1 2 3 ]
    set energies($i) [ analyze energy bonded 0 ]
    set forces($i)   [ part 0 print f ]

    scan $forces($i) "%f%f%f" fx($i) fy($i) fz($i)

    # set delta_a($i) [expr $angles($i) - $angles_ref($i)]
    # set delta_e($i) [expr $energies($i) - $energies_ref($i)]
    # set delta_fx($i) [expr $fx($i) - $fx_ref($i)]
    # set delta_fy($i) [expr $fy($i) - $fy_ref($i)]
    # set delta_fz($i) [expr $fz($i) - $fz_ref($i)]

    #output
    puts "angle: [expr [PI]/4-$angles($i)]"
    puts "energy: $energies($i)"
    puts "force: $forces($i)"
}
    
# Check if all deltas are smaller than the allowed error
##########################################################
# for { set i 0 } { $i <= $asteps } {incr i} {
#     if { [expr abs($delta_e($i))] > $epsilon } { error_exit "energy deviation too large" }
#     if { [expr abs($delta_fx($i))] > $epsilon } { error_exit "force deviation too large" }
#     if { [expr abs($delta_fy($i))] > $epsilon } { error_exit "force deviation too large" }
#     if { [expr abs($delta_fz($i))] > $epsilon } { error_exit "force deviation too large" }
# }

set dPhi [expr abs([PI]/4- [ bond_dihedral 0 1 2 3 ])]

if { $dPhi > $epsilon } {
    error_exit "angle deviation $dPhi to big."
}

puts ""
puts "Test for dihedral successful!"
puts ""

# close $infile
exit 0
