# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#    Max-Planck-Institute for Polymer Research, Theory Group
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
#=======================================================================================================
# test-1:  compute trajectories for a dipolar system of two particles using the dipolar P3M-algorithm 
#          and compare against well-know trajectories:
#          aim: discover errors in DIPOLES, ROTATION, INTEGRATION and
#               CONSTRAINT for MAGNETIC FIELD
#
# NB:      In this first test, we also test the ability of changing epsilon to metallic BC, which
#          has no sense for true magnetic simulations but is a nice feature for future electic dipoles
#=======================================================================================================
source "tests_common.tcl"

require_feature "DIPOLES"
require_feature "FFTW"
require_feature "CONSTRAINTS"
require_feature "ROTATION"

#require_max_nodes_per_side 1

puts "--------------------------------------------------------------------"
puts "- Testcase p3m_magnetostatics.tcl for magnetic dipoles running on [format %02d [setmd n_nodes]] nodes: -"
puts "--------------------------------------------------------------------"
puts "CAREFUL: tests do not check PRESSURES NOR ENERGIES"

if { [catch {

    set tcl_precision 15
    set accuracy 1e-8

    puts "Creating the magnetic system"
 
    setmd time_step 0.003
    setmd skin 0.05
    setmd max_num_cells 2744
    thermostat off
    setmd box_l 20.0 20.0 20.0

    part 0 pos 4 4 4 
    part 1 pos 6 6 6
    part 0 dip 1 0 0 
    part 1 dip 0 0 1 
    part 0  v  0 0 0
    part 1  v  0 0 0
    part 0 torque_lab 0 0 0
    part 1 torque_lab 0 0 0

    puts "Imposing magnetostatics (epsilon changed on purpose, disregard warning message about epsilon)"

 
    inter magnetic 1 p3m 5 32 7  0.3
    inter magnetic n_interpol  0
    inter magnetic mesh_off 0 0 0
    inter magnetic epsilon metallic    ;#see reasons above
 
    constraint ext_magn_field 0 0 15
   

    #------ create the trajectory of the two particles --------   
    puts "Start integration of the 2-particle system ..." 
   set outp [open "dummy.dat" w]
    for {set i_step 0} { $i_step <= 5 } {incr i_step} {   
      integrate 100
      puts "dipolar P3M test-1 Step: $i_step of 5 "
      for {set i 0} { $i < 2} {incr i} {
         puts $outp  [part $i print id pos quat v  ]
      }
    }
    close $outp
    
    #----- validate the trajectories against the known correct solution ---------------
    set outp  [open "p3m_magnetostatics.data" r]  
    set outp1 [open "dummy.dat" r]  
    for {set i_step 0} { $i_step <= 5 } {incr i_step} {   
      for {set i 0} { $i < 2} {incr i} {
          gets $outp chain
          gets $outp1 chain1
          set k [split $chain {" "}]
          set k1 [split $chain1 {" "}]
          for {set j 0} { $j < 11} {incr j} {
	      set a  [lindex $k  $j]
	      set a1 [lindex $k1 $j]
	      set diff [expr abs($a-$a1)]
	       if {$diff > $accuracy} {
	           puts "mismatch in the trajectories $i_step $i $j >>> $diff"
	           error "mismatch in the trajectories $i_step $i $j"
	       }
	  }
      }
    }
    close $outp
    close $outp1

     #end this part of the p3m-checks by cleaning the system ...................... 
  
    part deleteall
    inter magnetic 0.0
  
    exec rm -f "dummy.dat"

} res ] } {
    error_exit $res
} else {
   puts " P3M magnetostatics: test-1  oK"
}

exit 0
