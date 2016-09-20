# Copyright (C) 2012,2013,2014,2015,2016 The ESPResSo project
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

#    TEST DESCRIPTION
#
#    This test case loads the model of a blood cell, 
#    puts it into a fluid, sets the nonzero fluid velocity 
#    on the left inflow and lets the blood cell flows for 100 timesteps.
#    
#    In the beginning, the configuration is loaded 
#    from object_in_fluid_system.data.init
#    The model consists of a triangular mesh with 400 nodes.
#
#    After 100 timesteps, the positions, velocities and forces of the 400 particles 
#    are stored in arrays POS, VEL and FOR.
# 
#    Then, the reference configuration is loaded 
#    from object_in_fluid_system.data.final
#    and a check is performed whether computed configuration
#    stored in FOR, VEL and POS corresponds to the reference configuration. 

source "tests_common.tcl"

require_feature "ELECTROKINETICS"
require_feature "EK_BOUNDARIES"
require_feature "EK_REACTION"

set nn [format %02d [setmd n_nodes]]

puts "###############################################################"
puts "#              Testcase ek_reaction.tcl running on            #"
puts "#                           $nn nodes                          #"
puts "###############################################################\n"

# Macro to compare two vtk files

proc calculate_vtk_max_pointwise_difference {file1 file2} {

  set fp1 [open $file1 r]
  set fp2 [open $file2 r]
  
  set file1 [read $fp1]
  set file2 [read $fp2]
  
  close $fp1
  close $fp2
  
  set file1_data [split $file1 "\n"]
  set file2_data [split $file2 "\n"]
  
  set file1_data [lreplace $file1_data 0 11]
  set file2_data [lreplace $file2_data 0 11]
  
  if {[llength $file1_data] != [llength $file2_data]} {
  
    puts "different #lines in vtk files"
    return -1;
  }
  
  set max_diff 0.0
  
  for {set i 0} {$i < [llength $file1_data] && [lindex $file1_data $i] != ""} {incr i} {
  
    set file1_line [split [lindex $file1_data $i] " "]
    set file2_line [split [lindex $file2_data $i] " "]
    
    if {[llength $file1_line] != [llength $file2_line]} {
    
      puts "different #elements in lines of vtk files"
      return -1;
    }
    
    for {set k 0} {$k < [llength $file1_line]} {incr k} {
    
      set x1 [lindex $file1_line $k]
      set x2 [lindex $file2_line $k]
      
      if {[string is double -strict $x1] && [string is double -strict $x2]} {

        if { $x1 != 0.0 || $x2 != 0.0 } {

          set current_diff [expr abs( $x2 - $x1 )/( abs($x2) + abs($x1) ) ]
        
          if { $current_diff > $max_diff} {
          
            set max_diff $current_diff
          }

        }
      } elseif {$x1 != $x2} {
      
        puts "at least one element is not double, and they are different in vtk files ($x1 ::: $x2)"
        return -1
      }
    }
  }
  
  return $max_diff
}

# Set to 1 if you need a new 
# comparison configuration

set new_configuration 0

# Set the precision with which 
# the Tcl values are output

set tcl_precision 15

set length 12.0
set box_x $length
set box_y $length
set box_z $length
setmd box_l $box_x $box_y $box_z

# Set the cell system

cellsystem domain_decomposition

# Set periodicity

setmd periodic 1 1 1

# Set the skin size

setmd skin 0.1

# The time step

set dt 0.2
setmd time_step $dt

# Set fluid species properties

set density_reactant_eq 0.1
set density_product0_eq 0.9
set density_product1_eq 0.0

set diff_val 0.3
set diff_reactant_eq $diff_val
set diff_product0_eq $diff_val
set diff_product1_eq $diff_val

set val_val 0
set val_reactant_eq $val_val
set val_product0_eq $val_val
set val_product1_eq $val_val

# Set the EK-LB properties

set agrid 1.0
set viscosity 1.0
set frict 1.0
set temp 1.0
set bjer_length 1.0
set rho_fluid [expr $density_reactant_eq + $density_product0_eq + $density_product1_eq]
set viscosity_kinematic [expr $viscosity/$rho_fluid]

# Inititalize thermostat

thermostat off

# Initialize fluid 

electrokinetics agrid $agrid lb_density 1.0 viscosity $viscosity_kinematic \
                friction $frict T $temp bjerrum_length $bjer_length

# Set diffusion properties

electrokinetics 0 density $density_reactant_eq \
                D $diff_reactant_eq valency $val_reactant_eq 
electrokinetics 1 density $density_product0_eq \
                D $diff_product0_eq valency $val_product0_eq
electrokinetics 2 density $density_product1_eq \
                D $diff_product1_eq valency $val_product1_eq

# Set up reaction parameters

set re_index 0
set p0_index 1
set p1_index 2
set rate 0.5
set fraction_0 1.0
set fraction_1 0.5
set mass_react 34.0
set mass_prod0 18.0
set mass_prod1 32.0

# Initialize reaction

electrokinetics reaction reactant_index $re_index \
                         product0_index $p0_index \
                         product1_index $p1_index \
                         reactant_resrv_density $density_reactant_eq \
                         product0_resrv_density $density_product0_eq \
                         product1_resrv_density $density_product1_eq \
                         reaction_rate $rate \
                         mass_reactant $mass_react \
                         mass_product0 $mass_prod0 \
                         mass_product1 $mass_prod1 \
                         reaction_fraction_pr_0 $fraction_0 \
                         reaction_fraction_pr_1 $fraction_1

# Set up the boundary

electrokinetics boundary charge_density 0.0 \
                         sphere center [expr $box_x/2.0] \
                                       [expr $box_y/2.0] \
                                       [expr $box_z/2.0] \
                         radius 3 direction outside

# Tag the reactive regions

electrokinetics reaction region 0 box

electrokinetics reaction region 1 \
                         sphere center [expr $box_x/2.0] \
                                       [expr $box_y/2.0] \
                                       [expr $box_z/2.0] \
                         radius 4 direction outside

electrokinetics reaction region 1 \
                         sphere center [expr $box_x/2.0] \
                                       [expr $box_y/2.0] \
                                       [expr $box_z/2.0] \
                         radius 3 direction outside

electrokinetics reaction region 0 wall normal 0.0 0.0 1.0 dist [expr $box_z/2.0]

electrokinetics reaction region 2 wall normal 1.0 0.0 0.0 dist 1.0
electrokinetics reaction region 2 wall normal 0.0 1.0 0.0 dist 1.0
electrokinetics reaction region 2 wall normal 0.0 0.0 1.0 dist 1.0
electrokinetics reaction region 2 wall normal -1.0  0.0  0.0 dist [expr 1.0 - $box_x]
electrokinetics reaction region 2 wall normal  0.0 -1.0  0.0 dist [expr 1.0 - $box_y]
electrokinetics reaction region 2 wall normal  0.0  0.0 -1.0 dist [expr 1.0 - $box_z]

integrate 1000

# Make a new configuration file, otherwise do the comparison

if { $new_configuration != 0 } {
  electrokinetics 2 print density vtk "ek_reaction_density.vtk"
} else {
  electrokinetics 2 print density vtk "ek_reaction_density_tmp.vtk"
  set difference [calculate_vtk_max_pointwise_difference "./ek_reaction_density.vtk" "./ek_reaction_density_tmp.vtk"]
  file delete "./ek_reaction_density_tmp.vtk"

  puts "Maximum deviation to the reference point is: $difference\n"
  puts "Minor deviations are to be expected due to the sensitive"
  puts "dependence of the value of the mass flux for this system"
  puts "on the numerical precision that is attainable on the GPU."

  if { abs($difference) > 1.0e-02 } {
    error_exit "There is a significant difference with a previous result.\nPlease verify if this is correct."
  }
}

exit 0
