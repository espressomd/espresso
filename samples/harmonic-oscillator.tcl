#############################################################
#                                                           #
#  Template for A Chain of Harmonic Oscillators             #
#                                                           #
#############################################################
#
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

puts " "
puts "======================================================="
puts "=            harmonic_oscillator.tcl                  ="
puts "======================================================="
puts " "

puts "Program Information: \n[code_info]\n"

#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "harmonic_oscillator"
set ident "_s4"

# System parameters
#############################################################
# l_chain: the length of the chain
set l_chain        2

# Interaction parameters
#############################################################
set harmonic_r        1.0
set harmonic_k        1.0

# Integration parameters
#############################################################
set time_step    0.005
set skin         0.5

set int_steps    10
set int_n_times  1000

# Other parameters
#############################################################
set tcl_precision 6
set mypi          3.141592653589793
set vmd_output    "yes"
set vmd_wait      0

#############################################################
#  Setup System                                             #
#############################################################

set box_l       [expr 4.0*$l_chain*$harmonic_r]
set center      [expr $box_l/2.0]

setmd box_l     $box_l $box_l $box_l
setmd periodic  0 0 0
setmd time_step $time_step
setmd skin      $skin
#thermostat is turned off for pure MD simulation
thermostat off

# Interaction setup
#############################################################

# type 0 harmonic
inter 0 harmonic $harmonic_k $harmonic_r

#############################################################
# Harmonic chain setup
set part_id 0
set posx [expr $center - $l_chain/2.0 ]
set posy $center
set posz $center
for {set j 0} { $j < $l_chain } {incr j} {
    part $part_id pos $posx $posy $posz 
	 if { $j > 0 } {
	     part [expr $part_id-1] bond 0 $part_id
	 }
	 set posx [expr $posx+$harmonic_r]
	 incr part_id
}

# Status report
puts "Simulate $l_chain particles connected by harmonic springs\n"

#  VMD connection                                           #
#############################################################
if { $vmd_output=="yes" } {
    writepsf "$name$ident.psf"
    writepdb "$name$ident.pdb"
    for {set port 10000} { $port < 65000 } { incr port } {
	catch {imd connect $port} res
	if {$res == ""} break
    }
    set HOSTNAME [exec hostname]
    set vmdout_file [open "vmdoutput.script" "w"]
    puts $vmdout_file "mol load psf $name$ident.psf pdb $name$ident.pdb"
    puts $vmdout_file "rotate stop"
    puts $vmdout_file "imd connect 127.0.0.1 $port"
    puts $vmdout_file "mol modstyle 0 0 CPK 1.000000 0.300000 8.000000 6.000000"
    puts $vmdout_file "mol modcolor 0 0 SegName"
    close $vmdout_file
    puts "PSF and PDB written. IMD connection on port $port"
    puts "Start VMD in the same directory on the machine you like :"
    puts "vmd -e vmdoutput.script &"
    imd listen $vmd_wait
    exec vmd -e vmdoutput.script &
}


#############################################################
#      Integration                                          #
#############################################################
puts "Interactions are now: {[inter]}"

#fix the position of part 0
part 0 fix 
#give an initial kick
#part [expr 0] v -1 0 0 
part [expr $l_chain-1] pos 5 4 4 v 0 0 0 
#part [expr $l_chain-1] pos 4 5 4

for {set i 0} { $i < $part_id } {incr i} {
    puts "[part $i]"
}

#analyze set chains 0 $n_chains $l_chain
puts "\nStart integration: $int_n_times times $int_steps steps with time step $time_step"

    set energy_file [open "energy.txt" "w"]
    set part0_file [open "part0.txt" "w"]
    set part1_file [open "part1.txt" "w"]

for {set i 0} { $i < $int_n_times } { incr i} {
    puts -nonewline "Run $i at time=[setmd time]\r"
	 flush stdout
#record the energy and coordinates for particles 0 and 1
    puts  $energy_file "Run $i at time=[setmd time], E = [analyze energy ]"
    puts $part0_file "[part 0]"
    puts $part1_file "[part 1]"
    integrate $int_steps
    if { $vmd_output=="yes" } { catch {imd positions} }
}

puts "\nIntegration done."

#############################################################
puts "\nFinished"

exit
