# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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


# check the charge-charge P3M  algorithm
source "tests_common.tcl"

require_feature "LENNARD_JONES"
require_feature "ELECTROSTATICS"
require_feature "CUDA"
require_feature "FFTW"
require_feature "ADRESS" off
require_feature "NPT" off

puts "---------------------------------------------------------------"
puts "- Testcase p3m.tcl running on [format %02d [setmd n_nodes]] nodes"
puts "---------------------------------------------------------------"

set epsilon 1e-3
thermostat off
setmd time_step 0.01
setmd skin 0.05

proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

proc write_data {file} {
    set f [open $file "w"]
    blockfile $f write variable box_l
    blockfile $f write tclvariable {energy pressure}
    blockfile $f write interactions
    blockfile $f write particles {id pos q f}
    close $f
}

if { [catch {
    puts "Tests for P3M charge-charge interaction"
    read_data "p3m_system_gpu.data"

    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	set F($i) [part $i pr f]
    }
    ############## P3M-specific part
    # the P3M parameters are stored in p3m_system.data

    # to ensure force recalculation
#    inter coulomb n_interpol 0
    invalidate_system
    integrate 0
#after 99999

    
# here you can create the necessary snapshot
    if { 0 } {
	inter coulomb 1.0 p3m tune accuracy 1e-4
	integrate 0

	write_data "p3m_system.data"
    }

    ############## end

    puts [analyze energy]

    set cureng [lindex [analyze   energy coulomb] 0]
    set curprs 0





    ############## end, here RMS force error for P3M

    set rmsf 0
    set tot 0
    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	set resF [part $i pr f]
	set tgtF $F($i)
	set dx [expr abs(([lindex $resF 0] - [lindex $tgtF 0]))]
	set dy [expr abs(([lindex $resF 1] - [lindex $tgtF 1]))]
	set dz [expr abs(([lindex $resF 2] - [lindex $tgtF 2]))]
        set tot [expr $tot + [lindex $tgtF 0] * [lindex $tgtF 0] + [lindex $tgtF 1] * [lindex $tgtF 1] + [lindex $tgtF 2] * [lindex $tgtF 2] ]

	set rmsf [expr $rmsf + $dx*$dx + $dy*$dy + $dz*$dz]
    }

    set rfe [expr $rmsf]
    set rmsf [expr sqrt($rmsf/[setmd n_part])]
    puts "p3m-charges: rms force deviation $rmsf ($rfe $tot)"
    if { $rmsf > $epsilon } {
	error "p3m-charges: force error too large"
   }
   
   
     #end this part of the p3m-checks by cleaning the system .... 
   part deleteall
   inter coulomb 0.0

} res ] } {
    error_exit $res
}

exit 0
