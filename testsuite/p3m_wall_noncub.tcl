# Copyright (C) 2010,2011,2012 The ESPResSo project
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

# check the charge-charge P3M algorithm
source "tests_common.tcl"

require_feature "ELECTROSTATICS"
require_feature "FFTW"
require_feature "ADRESS" off

puts "---------------------------------------------------------------"
puts "- Testcase p3m_wall_noncub.tcl running on [format %02d [setmd n_nodes]] nodes:"
puts "---------------------------------------------------------------"

set accuracy 1e-3
thermostat off
setmd time_step 0.01
setmd skin 0.05

set file "p3m_wall_system.data"
set f [open $file "r"]

while {![eof $f]} { blockfile $f read auto}
    close $f

for { set i 0 } { $i <= [setmd max_part] } { incr i } {
  set F($i) [part $i pr f]
}

set npart_primary [setmd max_part]

set old_box [setmd box]
set new_box_z [expr 2*[lindex $old_box 2]]

setmd box [lindex $old_box 0] [lindex $old_box 1] $new_box_z

for { set i 0 } { $i <= $npart_primary } { incr i } {
    set old_pos [part $i pr pos]
    set x_pos [lindex $old_pos 0]
    set y_pos [lindex $old_pos 1]
    set z_pos [expr [lindex $old_pos 2] + [lindex $old_box 2]]
    set q [part $i pr q]
    part [expr $npart_primary + $i + 1] pos $x_pos $y_pos $z_pos q $q
}

# write vtf file
# set vtffile [open "p3m_wall.vtf" "w"]
# writevsf $vtffile
# for { set i 0 } { $i < 299 } { incr i } {
#     if { $i % 2 == 0 } {
# 	puts $vtffile "atom $i name S radius 0.1"
#     } else {
# 	puts $vtffile "atom $i name O radius 0.1"
#     }
# }
# writevcf $vtffile
# close $vtffile

# To retune, use
#puts [inter coulomb 1.0 p3m tune accuracy 1.e-3]

# If you want to use specific P3M params, use
# inter coulomb 1.0 p3m 2.000000 32 6 1.322773
# inter coulomb n_interpol 0

puts [inter coulomb 1.0 p3m tune accuracy $accuracy]
puts [inter coulomb]
if { [catch {

    invalidate_system
    integrate 0

    set rmsf 0
    set rmsf_x 0
    set rmsf_y 0
    set rmsf_z 0

    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	set resF [part $i pr f]
	set j $i
	if { $i > $npart_primary } {
	    set j [expr $i - $npart_primary - 1]
	}

	set tgtF $F($j)

	set dx [expr ([lindex $resF 0] - [lindex $tgtF 0])]
	set dy [expr ([lindex $resF 1] - [lindex $tgtF 1])]
	set dz [expr ([lindex $resF 2] - [lindex $tgtF 2])]
    
	set rmsf [expr $rmsf + $dx*$dx + $dy*$dy + $dz*$dz]
	set rmsf_x [expr $rmsf_x + $dx*$dx]
	set rmsf_y [expr $rmsf_y + $dy*$dy]
	set rmsf_z [expr $rmsf_z + $dz*$dz]
    }

    set rmsf [expr sqrt($rmsf / [setmd n_part])]
    set rmsf_x [expr sqrt($rmsf_x / [setmd n_part])]
    set rmsf_y [expr sqrt($rmsf_y / [setmd n_part])]
    set rmsf_z [expr sqrt($rmsf_z / [setmd n_part])]

    puts [format "rms_force_error=%e" $rmsf]
    puts [format "rms_force_error_x=%e" $rmsf_x]
    puts [format "rms_force_error_y=%e" $rmsf_y]
    puts [format "rms_force_error_z=%e" $rmsf_z]
    if { $rmsf > $accuracy } {
	error [format "p3m-charges: rms_force_error=%e larger than accuracy=%e" $rmsf $accuracy]
    }
    
} res ] } {
    error_exit $res
}

exit 0
