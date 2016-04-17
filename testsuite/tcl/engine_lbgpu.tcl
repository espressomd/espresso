# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
source "tests_common.tcl"

puts "------------------------------------------------"
puts "- Testcase engine_lbgpu.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------------"

require_feature "ENGINE"
require_feature "LB_GPU"

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


set boxl 12
set sampsteps 2000
set tstep 0.01
set temp 0

setmd box_l $boxl $boxl $boxl 
setmd periodic 1 1 1
setmd skin 0.1
setmd time_step $tstep

part 0 pos 6.0 3.0 2.0 q 0 swimming pusher v_swim 0.10 dipole_length 1.0 rotational_friction  2.0 quat [expr sqrt(.5)] [expr sqrt(.5)]               0               0
part 1 pos 2.0 3.0 6.0 q 0 swimming pusher f_swim 0.03 dipole_length 2.0 rotational_friction 20.0 quat [expr sqrt(.5)]               0 [expr sqrt(.5)]               0
part 2 pos 3.0 2.0 6.0 q 0 swimming puller v_swim 0.15 dipole_length 0.5 rotational_friction 15.0 quat [expr sqrt(.5)]               0               0 [expr sqrt(.5)]
part 3 pos 3.0 6.0 2.0 q 0 swimming puller f_swim 0.05 dipole_length 1.5 rotational_friction  6.0 quat               0               0 [expr sqrt(.5)] [expr sqrt(.5)]

lbfluid gpu agrid 1.0 tau $tstep friction 0.5 viscosity 1.0 density 1.0 couple 2pt
thermostat lb $temp

integrate $sampsteps

# Make a new configuration file, otherwise do the comparison

if { $new_configuration != 0 } {
  lbfluid print vtk velocity "engine_lbgpu_2pt.vtk"
} else {
  lbfluid print vtk velocity "engine_lbgpu_tmp.vtk"
  set difference [calculate_vtk_max_pointwise_difference "./engine_lbgpu_2pt.vtk" "./engine_lbgpu_tmp.vtk"]
  file delete "./engine_lbgpu_tmp.vtk"

  puts "Maximum deviation to the reference point is: $difference\n"

  if { abs($difference) > 1.0e-07 } {
    error_exit "There is a significant difference with a previous result.\nPlease verify if this is correct."
  }
}

part 0 pos 6.0 3.0 2.0 q 0 swimming pusher v_swim 0.10 dipole_length 1.0 rotational_friction  2.0 quat [expr sqrt(.5)] [expr sqrt(.5)]               0               0
part 1 pos 2.0 3.0 6.0 q 0 swimming pusher f_swim 0.03 dipole_length 2.0 rotational_friction 20.0 quat [expr sqrt(.5)]               0 [expr sqrt(.5)]               0
part 2 pos 3.0 2.0 6.0 q 0 swimming puller v_swim 0.15 dipole_length 0.5 rotational_friction 15.0 quat [expr sqrt(.5)]               0               0 [expr sqrt(.5)]
part 3 pos 3.0 6.0 2.0 q 0 swimming puller f_swim 0.05 dipole_length 1.5 rotational_friction  6.0 quat               0               0 [expr sqrt(.5)] [expr sqrt(.5)]

lbfluid gpu agrid 1.0 tau $tstep friction 0.5 viscosity 1.0 density 1.0 couple 3pt
thermostat lb $temp

integrate $sampsteps

# Make a new configuration file, otherwise do the comparison

if { $new_configuration != 0 } {
  lbfluid print vtk velocity "engine_lbgpu_3pt.vtk"
} else {
  lbfluid print vtk velocity "engine_lbgpu_tmp.vtk"
  set difference [calculate_vtk_max_pointwise_difference "./engine_lbgpu_3pt.vtk" "./engine_lbgpu_tmp.vtk"]
  file delete "./engine_lbgpu_tmp.vtk"

  puts "Maximum deviation to the reference point is: $difference\n"

  if { abs($difference) > 2.0e-07 } {
    error_exit "There is a significant difference with a previous result.\nPlease verify if this is correct."
  }
}

ok_exit
