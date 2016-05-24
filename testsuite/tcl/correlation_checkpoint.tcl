# Copyright (C) 2012,2014,2015,2016 The ESPResSo project
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

puts "---------------------------------------------------------------"
puts "- Testcase correlation_checkpoint.tcl running on [format %02d [setmd n_nodes]] nodes"
puts "---------------------------------------------------------------"

setmd box_l 6. 6. 6.
thermostat langevin 1.0 1.0
setmd skin 1.0

set dt 0.05
setmd time_step $dt

part 0 pos 4.100  3.000  3.000
part 1 pos 1.000  3.000  2.000

set pp  [observable new particle_positions all]
set rdf [observable new rdf 0 0 ]
set orgc1 [correlation new obs1 $pp corr_operation square_distance_componentwise dt $dt tau_max 1000 tau_lin 16]
set orgc2 [correlation new obs1 $pp obs2 $rdf corr_operation tensor_product dt $dt tau_max 1000 tau_lin 16]

correlation $orgc1 autoupdate start
correlation $orgc2 autoupdate start

for { set i 0 } { $i < 100 } { incr i } {
  integrate 200
}

correlation $orgc1 write_checkpoint_binary "correlation_checkpoint_c1.bin"
correlation $orgc1 write_checkpoint_ascii  "correlation_checkpoint_c1.txt"
correlation $orgc2 write_checkpoint_binary "correlation_checkpoint_c2.bin"
correlation $orgc2 write_checkpoint_ascii  "correlation_checkpoint_c2.txt"

set pp2  [observable new particle_positions all]
set rdf2 [observable new rdf 0 0 ]

set binc1 [correlation new obs1 $pp corr_operation square_distance_componentwise dt $dt tau_max 1000 tau_lin 16]
set binc2 [correlation new obs1 $pp obs2 $rdf corr_operation tensor_product dt $dt tau_max 1000 tau_lin 16]

set ascc1 [correlation new obs1 $pp corr_operation square_distance_componentwise dt $dt tau_max 1000 tau_lin 16]
set ascc2 [correlation new obs1 $pp obs2 $rdf corr_operation tensor_product dt $dt tau_max 1000 tau_lin 16]

correlation $binc1 autoupdate start
correlation $ascc1 autoupdate start

correlation $binc1 read_checkpoint_binary "correlation_checkpoint_c1.bin"
correlation $ascc1 read_checkpoint_ascii  "correlation_checkpoint_c1.txt"
correlation $binc2 read_checkpoint_binary "correlation_checkpoint_c2.bin"
correlation $ascc2 read_checkpoint_ascii  "correlation_checkpoint_c2.txt"

correlation $binc2 autoupdate start
correlation $ascc2 autoupdate start

for { set i 0 } { $i < 100 } { incr i } {
  integrate 200
}

correlation $orgc1 finalize
correlation $orgc2 finalize
correlation $binc1 finalize
correlation $ascc1 finalize
correlation $binc2 finalize
correlation $ascc2 finalize

correlation $orgc1 write_to_file "correlation_checkpoint_c1.txt"
correlation $orgc2 write_to_file "correlation_checkpoint_c2.txt"
correlation $binc1 write_to_file "correlation_checkpoint_c1.bin"
correlation $binc2 write_to_file "correlation_checkpoint_c2.bin"

if {[catch {exec diff  "correlation_checkpoint_c1.txt" "correlation_checkpoint_c1.bin"}]} {
    puts "Error in the correlation checkpointing."
    puts "First check failed"
    error_exit
}

if {[catch {exec diff  "correlation_checkpoint_c2.txt" "correlation_checkpoint_c2.bin"}]} {
    puts "Error in the correlation checkpointing."
    puts "Second check failed"
    error_exit
}

correlation $ascc1 write_to_file "correlation_checkpoint_c1.bin"
correlation $ascc2 write_to_file "correlation_checkpoint_c2.bin"

set f1 [open  "correlation_checkpoint_c1.bin" r]
set f2 [open  "correlation_checkpoint_c1.txt" r]

while {[gets $f1 line1] != -1 } {
    gets $f2 line2
    set length [llength $line1]
    for {set i 0} {$i < $length} {incr i} {
	if { [expr abs(([lindex $line1 $i] - [lindex $line2 $i])/ ([lindex $line2 $i] +0.0001) )] > 0.01 } {
	    puts $line1 
	    puts $line2
	    puts $i
	    puts "Error in the correlation checkpointing."
	    puts "Third check failed"
	    error_exit
	}
    }
}

set f1 [open  "correlation_checkpoint_c2.bin" r]
set f2 [open  "correlation_checkpoint_c2.txt" r]

while {[gets $f1 line1] != -1 } {
    gets $f2 line2
    set length [llength $line1]
    for {set i 0} {$i < $length} {incr i} {
	if { [expr abs(([lindex $line1 $i] - [lindex $line2 $i])/ ([lindex $line2 $i] +0.0001) )] > 0.01 } {
	    puts "Error in the correlation checkpointing."
	    puts "Forth check failed"
	    error_exit
	}
    }
}

puts "Test successful"
ok_exit
