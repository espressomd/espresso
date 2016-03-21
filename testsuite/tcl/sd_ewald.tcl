# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#  Max-Planck-Institute for Polymer Research, Theory Group
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

source "tests_common.tcl"

require_feature "SD"
require_max_nodes_per_side 1
require_cudadevice
require_feature "EXTERNAL_FORCES"
require_feature "SD_NOT_PERIODIC" "off"

proc set_pos { dx dy dz boxl} {
    setmd box_l [expr $boxl*$dx] [expr $boxl*$dy] [expr $boxl*$dz]
    set p 0;
    for {set x 0} {$x < $dx} {incr x} {
	for {set y 0} {$y < $dy} {incr y} {
	    for {set z 0} {$z < $dz} {incr z} {
		part $p pos [expr $x*$boxl] [expr $y*$boxl] [expr $z*$boxl]  ext_force 1 0 0
		incr p
		part $p pos [expr $x*$boxl] [expr ($y+0.5)*$boxl] [expr $z*$boxl]  ext_force -1 0 0
		incr p
	    }
	}
    }
}


proc check_pos { dx dy dz boxl pos0 pos1} {
    set p 0;
    for {set x 0} {$x < $dx} {incr x} {
	for {set y 0} {$y < $dy} {incr y} {
	    for {set z 0} {$z < $dz} {incr z} {
		set np0 [part $p print pos]
		incr p
		set np1 [part $p print pos]
		incr p
		set diff [concat [expr $x*$boxl] [expr $y*$boxl] [expr $z*$boxl] ]
		for {set d 0} { $d < 3 } {incr d} {
		    set d0 [expr abs( [lindex $np0 $d] - [lindex $diff $d] - [lindex $pos0 $d])]
		    set d1 [expr abs( [lindex $np1 $d] - [lindex $diff $d] - [lindex $pos1 $d])]
		    set prec 1e-6
		    if { $d0 > $prec || $d1 > $prec } {
			#puts "org0: $pos0\nnew0: $np0"
			#puts "org1: $pos1\nnew1: $np1"
			set d [expr max($d0,$d1)]
			puts "Error was $d, required was < $prec. Box was $dx $dy $dz."
			error_exit
		    }
		}
	    }
	}
    }
}


setmd sd_visc 1
setmd sd_radius 0.5
setmd time_step 0.01
set boxl 4
set steps 5
setmd periodic 1 1 1
setmd skin 0.1

thermostat sd 0

#simulate with smallest box
set_pos 1 1 1 $boxl
integrate_sd $steps
set pos0 [part 0 print pos]
set pos1 [part 1 print pos]


puts "check non cubic boxes"
set_pos 1 1 2 $boxl
integrate_sd $steps
check_pos 1 1 2 $boxl $pos0 $pos1

set_pos 1 2 1 $boxl
integrate_sd $steps
check_pos 1 2 1 $boxl $pos0 $pos1

set_pos 2 1 1 $boxl
integrate_sd $steps
check_pos 2 1 1 $boxl $pos0 $pos1

puts "check cubic boxes"
set_pos 2 2 2 $boxl
integrate_sd $steps
check_pos 2 2 2 $boxl $pos0 $pos1

set_pos 4 4 4 $boxl
integrate_sd $steps
check_pos 4 4 4 $boxl $pos0 $pos1
