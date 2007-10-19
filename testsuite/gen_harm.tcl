#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    exec $ESPRESSO_SOURCE/Espresso $0 $*
# 
#  This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
# 
#
# This file can be used to generate the data for the testcase harm.tcl.
#

puts "------------------------------------------"
puts "- Generating the testcase for harm.tcl"
puts "------------------------------------------"

##################################################
# system setup
set file "harm_system.data.gz"

set L 10
set temperature 1.0
set gamma 0.5

set time_step 0.001
set num_steps 10000
set skin 0.0

set length_polymer 300

set harm_k 30.0
set harm_r0 1.0

set blbins 30
set blmin 0.5
set blmax 1.5

# set up global parameters
setmd box_l $L $L $L
setmd time_step $time_step
setmd skin $skin
thermostat langevin $temperature $gamma

# set up harm interaction
inter 0 harm $harm_k $harm_r0

# set up polymer
polymer 1 $length_polymer $harm_r0 mode RW bond 0

##################################################
# simulation

# run a number of integration steps
puts "Running simulation..."
integrate $num_steps

##################################################
# compute the bond length distribution at the end
set blstep [expr ($blmax-$blmin)/$blbins]

# reset the histogram
for { set i 0 } { $i < $blbins } { incr i } { set blh($i) 0 }

# compute the histogram
for { set i 0 } { $i < [expr $length_polymer - 1] } { incr i } {
    set bl [bond_length $i [expr $i+1]]
    if { $bl > $blmin && $bl < $blmax } then {
	incr blh([expr int(($bl-$blmin)/$blstep)])
    }
}

# print the histogram
puts "Bond length distribution..."
for { set i 0 } { $i < $blbins } { incr i } {
    puts -nonewline [format "%5.3f: " [expr $blmin+$i*$blstep]]

    for { set j 0 } { $j < $blh($i) } { incr j } {
	puts -nonewline "#"
    }
    puts ""
}

# set f [open "gen_harm.vtf" w]
# writevsf $f
# writevcf $f
# close $f

##################################################
# compute and write the testcase data
thermostat off
velocities 0
invalidate_system
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

