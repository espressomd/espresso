# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
source "tests_common.tcl"

puts "------------------------------------------------"
puts "- Testcase engine_langevin.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------------"

require_feature "ENGINE"

proc z_f {t z0} {
  global f_swim
  global gamma
  return [expr ($f_swim/$gamma)*(-1./$gamma + $t + (1./$gamma)*exp(-$gamma*$t))+$z0]
}

proc z_v {t z0} {
  global v_swim
  global gamma
  return [expr $v_swim*(-1./$gamma + $t + (1./$gamma)*exp(-$gamma*$t))+$z0]
}

set boxl 12
set sampsteps 2000
set tstep 0.01

setmd box_l $boxl $boxl $boxl 
setmd periodic 1 1 1
setmd skin 0.1
setmd time_step $tstep

set pos_0 "[expr $boxl/2.] [expr $boxl/2.] [expr 1.*$boxl/3.]"
set pos_1 "[expr $boxl/2.] [expr $boxl/2.] [expr 2.*$boxl/3.]"
set v_swim 0.3
set f_swim 0.1

part 0 pos [lindex $pos_0 0] [lindex $pos_0 1] [lindex $pos_0 2] \
   q 0 swimming v_swim $v_swim

part 1 pos [lindex $pos_1 0] [lindex $pos_1 1] [lindex $pos_1 2] \
   q 0 swimming f_swim $f_swim

set temp 0
set gamma 1
thermostat langevin $temp $gamma

integrate $sampsteps

set t [setmd time]
set pos_0_simulate "[part 0 print pos]" 
set pos_0_analytic "[lindex $pos_0 0] [lindex $pos_0 1] [z_v $t [lindex $pos_0 2]] "
set pos_1_simulate "[part 1 print pos]" 
set pos_1_analytic "[lindex $pos_1 0] [lindex $pos_1 1] [z_f $t [lindex $pos_1 2]] "
set delta_pos_0 [veclen [vecsub $pos_0_simulate $pos_0_analytic]]
set delta_pos_1 [veclen [vecsub $pos_1_simulate $pos_1_analytic]]

if { ! (1.4e-3 < $delta_pos_0 && $delta_pos_0 < 1.6e-3) } {
  error_exit "Velocity swimming deviation too large."
}
if { ! (4.9e-4 < $delta_pos_1 && $delta_pos_1 < 5.1e-4) } {
  error_exit "Force swimming deviation too large."
}

ok_exit
