# Copyright (C) 2011 The ESPResSo project
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

require_feature "ROTATION"
require_feature "ROTATIONAL_INERTIA"
require_feature "GAY_BERNE"

puts "----------------------------------------------"
puts "- Testcase rotation.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------------"

set epsilon 5e-4
thermostat off

setmd time_step 0.01
#setmd time_step 0.011
setmd skin 0.5
thermostat langevin .5  0.0001

proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

proc read_kinetic_energy {site} {
    set kx 0 
    set ky 0 
    set kz 0 
    set kin 0 
    set rx 0 
    set ry 0 
    set rz 0 
    set renergy 0
    for {set i 0} {$i < $site} {incr i} {
    set mass [part $i print mass]	
    foreach {vx vy vz} [part $i print v] break
    set kx  [expr $kx + 0.5*$mass*($vx*$vx)]
    set ky  [expr $ky + 0.5*$mass*($vy*$vy)]
    set kz  [expr $kz + 0.5*$mass*($vz*$vz)]
    set kin [expr $kin+ 0.5*$mass*($vx*$vx+$vy*$vy+$vz*$vz)]
    foreach {vx vy vz} [part $i print omega] break	
    foreach {ix iy iz} [part $i print rinertia] break			
    set rx  [expr $rx + 0.5*$ix*($vx*$vx)]
    set ry  [expr $ry + 0.5*$iy*($vy*$vy)]
    set rz  [expr $rz + 0.5*$iz*($vz*$vz)]
    set renergy  [expr $renergy +0.5*($ix*$vx*$vx+$iy*$vy*$vy+$iz*$vz*$vz)]
    }
    puts "TE: [expr $kin+$renergy] KE: $kin RE: $renergy"
  return "$kx $ky $kz $rx $ry $rz"
}


if { [catch {
    read_data "gb_system.data"
# Use custom masses and rotational inertias
    for {set i 0} {$i <1000} {incr i} {
    part $i mass 1
#   part $i rinertia 1 1 1
#    part $i mass 100
    part $i rinertia 100 200 300
    }

    # to ensure force recalculation
    invalidate_system
    
    inter 0 0 gay-berne 1.0 1.0 4.0 3.0 5.0 2.0 1.0
    
    set GBeng_0 [expr [analyze energy gb 0 0]]
    set toteng_0 [analyze energy total]
    if { [expr abs($toteng_0 - $GBeng_0)] > $epsilon } {
	error "system has unwanted energy contribution, i.e. U_GB != U_total"
    }
    puts "energy before integration: [analyze energy]"
    read_kinetic_energy 1000
    
    set kxt 0; set kyt 0; set kzt 0; set kt 0;
    set rxt 0; set ryt 0; set rzt 0; set rt 0;
    for {set i 0} {$i < 40} {incr i} {
    integrate 400
    set a [read_kinetic_energy 1000]
    foreach {kx ky kz rx ry rz} $a break
    set kxt [expr $kxt+$kx] 
    set kyt [expr $kyt+$ky] 
    set kzt [expr $kzt+$kz] 
    set kt [expr $kt+$kx+$ky+$kz] 
    set rxt [expr $rxt+$rx] 
    set ryt [expr $ryt+$ry] 
    set rzt [expr $rzt+$rz] 
    set rt [expr $rt+$rx+$ry+$rz] 
    set tt [expr $rt+$kt]
     puts "        T: [expr $kxt/$tt] [expr $kyt/$tt] [expr $kzt/$tt] "
     puts "        R: [expr $rxt/$tt] [expr $ryt/$tt] [expr $rzt/$tt] "
     puts "energy after integration: [analyze energy]"
    }

    set kxt 0; set kyt 0; set kzt 0; set kt 0;
    set rxt 0; set ryt 0; set rzt 0; set rt 0;
    for {set i 0} {$i < 40} {incr i} {
    integrate 400
    set a [read_kinetic_energy 1000]
    foreach {kx ky kz rx ry rz} $a break
    set kxt [expr $kxt+$kx] 
    set kyt [expr $kyt+$ky] 
    set kzt [expr $kzt+$kz] 
    set kt [expr $kt+$kx+$ky+$kz] 
    set rxt [expr $rxt+$rx] 
    set ryt [expr $ryt+$ry] 
    set rzt [expr $rzt+$rz] 
    set rt [expr $rt+$rx+$ry+$rz] 
    set tt [expr $rt+$kt]
     puts "        T: [expr $kxt/$tt] [expr $kyt/$tt] [expr $kzt/$tt] "
     puts "        R: [expr $rxt/$tt] [expr $ryt/$tt] [expr $rzt/$tt] "
     puts "energy after integration: [analyze energy]"
    }

    # check the conservation of the total energy
    set toteng [analyze energy total]
    set rel_eng_error [expr abs(($toteng_0 - $toteng)/$toteng)]
    puts "total energy deviation: $rel_eng_error"
    if { $rel_eng_error > $epsilon } {
	error "relative energy error is too large"
    }
    
    # check new GB energy against expected value
    set GB_expected -2971.72
    set GBeng [expr [analyze energy gb 0 0]]
    set rel_eng_error [expr abs(($GBeng - $GB_expected)/$toteng)]
    puts "   GB energy deviation: $rel_eng_error"
    if { $rel_eng_error > $epsilon } {
	error "relative energy error is too large"
    }

} res ] } {
    error_exit $res
}

exit 0
