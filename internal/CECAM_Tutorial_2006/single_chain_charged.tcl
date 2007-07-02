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
#############################################################
#                                                           #
#  ESPResSo Tutorial                                        #
#                                                           #
#  Created:       25.07.2004 by HL                          #
#  Modified:      26.06.2006 by OL                          #
#  Modified:      28.11.2006 by TS                          #
#                                                           #
#############################################################

#############################################################
#   Preparing the System                                    #
#############################################################

set l_poly  20
set n_ci    $l_poly
set n_part  [expr $l_poly+$n_ci]
set density 0.001
set volume     [expr $n_part/$density] 
set box_length [expr pow($volume,1.0/3.0)]

puts "Simulate Polymer N=$l_poly at density $density"
puts "Simulation box: $box_length"

#############################################################
#   Initializng ESPResSo                                    #
#############################################################

setmd box_l $box_length $box_length $box_length
setmd time_step 0.01
setmd skin 0.4
integrate set nvt
thermostat langevin 1.0 1.0

#############################################################
#   Setting up Interactions                                 #
#############################################################

inter 0 fene 7.0 2.0
inter 0 0 lennard-jones 1.0 1.0 1.12246 0.25 0.0

#############################################################
#   Creating Polymer                                        #
#############################################################

polymer 1 $l_poly 1.0 

#############################################################
#   Visualization                                           #
#############################################################

set vmd "yes"
if { $vmd == "yes" } {
    prepare_vmd_connection tutorial 3000
    exec sleep 4
    imd positions
}

#############################################################
#   Warmup                                                  #
#############################################################

set min 0
set cap 10
while { $min < 0.8 } {
    inter ljforcecap $cap
    integrate 20
    set min [analyze mindist]
    incr cap 10
}
puts "Warmup finished. Minimal distance now $min"
inter ljforcecap 0


#############################################################
#   Integration                                             #
#############################################################

set n_cycle 100000
set n_steps 2
set obs [open "re.dat" "w"]

set i 0 
while { $i<$n_cycle } {
    integrate $n_steps

#############################################################
#   Analysis                                                #
#############################################################

    analyze set chains 0 1 $l_poly
    set re [lindex [analyze re] 0]

#############################################################
#    Store and view Results                                 #
#############################################################
   
    puts $obs "[setmd time] $re"
    set filename [ format "snapshot_%04i.xyz" $i ]
#    writexyz $filename

    if { $vmd == "yes" } { imd positions }

    incr i
}
close $obs

