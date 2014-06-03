# Copyright (C) 2010,2012,2013 The ESPResSo project
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
#
# This file can be used to generate the data for the testcase fene.tcl.
#

puts "------------------------------------------"
puts "- Generating the testcase for fene.tcl"
puts "------------------------------------------"

source "tests_common.tcl"
require_feature "ROTATION" off

##################################################
# system setup
set file "fene_system.data.gz"

set L 10
set temperature 1.0
set gamma 0.5

set time_step 0.001
set num_steps 10000
set skin 0.0

set length_polymer 300

set fene_k 30.0
set fene_dr 0.5
set fene_r0 1.0

set bl_histo_bins 30

# set up global parameters
setmd box_l $L $L $L
setmd time_step $time_step
setmd skin $skin
thermostat langevin $temperature $gamma

# set up fene interaction
inter 0 fene $fene_k $fene_dr $fene_r0

# set up polymer
polymer 1 $length_polymer $fene_r0 mode RW bond 0

##################################################
# simulation

# run a number of integration steps
puts "Running simulation..."
integrate $num_steps

##################################################
# compute the bond length distribution at the end
set blmin [expr $fene_r0 - $fene_dr]
set blmax [expr $fene_r0 + $fene_dr]
set blstep [expr ($blmax-$blmin)/$bl_histo_bins]

# reset the histogram
for { set i 0 } { $i < $bl_histo_bins } { incr i } { set blh($i) 0 }

# compute the histogram
for { set i 0 } { $i < [expr $length_polymer - 1] } { incr i } {
    set bl [bond_length $i [expr $i+1]]
    incr blh([expr int(($bl-$blmin)/$blstep)])
}

# print the histogram
puts "Bond length distribution..."
for { set i 0 } { $i < $bl_histo_bins } { incr i } {
    puts -nonewline [format "%5.3f: " [expr $blmin+$i*$blstep]]

    for { set j 0 } { $j < $blh($i) } { incr j } {
	puts -nonewline "#"
    }
    puts ""
}

# set f [open "gen_fene.vtf" w]
# writevsf $f
# writevcf $f
# close $f

##################################################
# compute and write the testcase data
thermostat off
velocities 0
integrate 0

# write the data
set f [open "|gzip -c - >$file" w]
set energy [analyze energy total]
set pressure [analyze pressure total]
blockfile $f write tclvariable {energy pressure}
blockfile $f write variable box_l
blockfile $f write interactions
blockfile $f write particles {id pos f}
blockfile $f write bonds all
close $f

exit 0
