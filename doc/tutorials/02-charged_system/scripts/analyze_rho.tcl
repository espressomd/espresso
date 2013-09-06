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
set bins 50
set cnt 0
for {set i 0} {$i < $bins} {incr i} {
    lappend avg_rho0 0
    lappend avg_rho1 0
}

foreach filename [lrange $argv 1 end] {
    set f [open $filename "r"]
    while { [blockfile $f read auto] != "eof" } {}
    close $f

    set data0 ""
    set data1 ""
    for {set p 0} {$p <= [setmd max_part]} {incr p} {
	if {[part $p pr type] == 0} {
	    lappend data0 [lindex [part $p pr pos] 2]
	} {
	    lappend data1 [lindex [part $p pr pos] 2]
	}
    }
    set rho0 [bin -linbins 0.5 [expr $box_lz + 0.5] $bins $data0]
    set rho1 [bin -linbins 0.5 [expr $box_lz + 0.5] $bins $data1]

    set avg_rho0 [vecadd $avg_rho0 $rho0]
    set avg_rho1 [vecadd $avg_rho1 $rho1]
    incr cnt
}

set avg_rho0 [vecscale [expr 1.0/$cnt] $avg_rho0]
set avg_rho1 [vecscale [expr 1.0/$cnt] $avg_rho1]
set centers [bin -linbins 0.5 [expr $box_lz + 0.5] $bins -binctrwdth]
set plot [open "rho.data" "w"]
puts $plot "\# z rho0(z) rho1(z)"
foreach z $centers rho0 $avg_rho0 rho1 $avg_rho1 {
    puts $plot "[lindex $z 0] $rho0 $rho1"
}
close $plot
