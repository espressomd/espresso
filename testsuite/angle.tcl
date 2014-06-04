# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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


source "tests_common.tcl"

require_feature "BOND_ANGLE"

puts "----------------------------------------"
puts "- Testcase angle.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------"


foreach type { harmonic cosine cossquare }  {

  # system and integration parameters
  set epsilon 1e-4
  thermostat off 
  setmd time 1
  setmd time_step 1
  setmd skin 0

  setmd box_l 10.0 10.0 10.0

  # remove previous configuration
  part deleteall 

  # create 3 particles
  part 0 pos 0.0 0.0 0.0	
  part 1 pos 1.0 0.0 0.0
  part 2 pos 1.0 1.0 0.0		

  # define bond
  inter 0 angle_$type 3.0 [expr [PI]]
  part 1 bond 0 0 2

  # integrate

  set asteps 10.0

  # load files for comparison
  if { [file exists "angle_$type.data"] } { set infile [open "angle_$type.data"]
  } else { puts "Error: Could not find file \"angle_$type.data\".Aborting..."; exit }
  set infile_data [read $infile]

  # initialize i/o and temporary variables
  set in1 0.0
  set in2 0.0
  set in3 0.0
  set i 0
  set j 0
  set k 0

  # process data file
  set data [split $infile_data "\n"]
  foreach line $data {
    scan $line "%s%f%f%f" name in1 in2 in3
  #  puts "read \" $name $in1 $in2 $in3 \" from file"
    if { ![ string compare $name "angle:"]  } { set angles_ref($i) $in1;   incr i 1 } 
    if { ![ string compare $name "energy:"] } { set energies_ref($j) $in1; incr j 1 } 
    if { ![ string compare $name "force:"]  } { set fx_ref($k) $in1; set fy_ref($k) $in2; set fz_ref($k) $in3;  incr k 1 } 
  }  

  # testing potential
  for { set i 0 } { $i <= $asteps } { incr i } {
    part 2 pos [expr {1+sin($i/$asteps*[PI])} ] [expr {cos($i/$asteps*[PI])} ] 0.0
    integrate 0

    # analyze
    set angles($i)   [ bond_angle 0 1 2 ]
    set energies($i) [ analyze energy bonded 0 ]
    set forces($i)   [ part 2 print f ]

    scan $forces($i) "%f%f%f" fx($i) fy($i) fz($i)

    set delta_a($i) [expr $angles($i) - $angles_ref($i)]
    set delta_e($i) [expr $energies($i) - $energies_ref($i)]
    set delta_fx($i) [expr $fx($i) - $fx_ref($i)]
    set delta_fy($i) [expr $fy($i) - $fy_ref($i)]
    set delta_fz($i) [expr $fz($i) - $fz_ref($i)]

#    #output
#    puts "angle: [bond_angle 0 1 2]"
#    puts "energy: $energies($i)"
#    puts "force: $forces($i)"

    # output
#    puts "angle: $delta_a($i)	energies: $delta_e($i)	forces: $delta_fx($i) $delta_fy($i) $delta_fz($i)" 
  }

  # Check if all deltas are smaller than the allowed error
  for { set i 0 } { $i <= $asteps } {incr i} {
    if { [expr abs($delta_e($i))] > $epsilon } { error "Error: energy deviation too large" ; exit }
    if { [expr abs($delta_fx($i))] > $epsilon } { puts "Error: force deviation too large"  ; exit}
    if { [expr abs($delta_fy($i))] > $epsilon } { puts "Error: force deviation too large"  ; exit}
    if { [expr abs($delta_fz($i))] > $epsilon } { puts "Error: force deviation too large"  ; exit}
  }

  puts ""
  puts "Test for angle_$type successful!"
  puts ""
}

close $infile
exit 0
