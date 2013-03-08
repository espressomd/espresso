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
setmd box_l 10 10 10
inter 0 fene 1.0 1.0
for {set i 0} {$i < 10000} {incr i} {polymer 1 2 1 pos 5 5 5 start [expr $i*2]}
writepdb "sphere.pdb"
set vmdout_file [open "vmd_start.script" "w"]
puts $vmdout_file "mol load pdb sphere.pdb"
puts $vmdout_file "logfile vmd.log"
puts $vmdout_file "rotate stop"
puts $vmdout_file "logfile off"
puts $vmdout_file "mol modstyle 0 0 Points"
close $vmdout_file
exec vmd -e vmd_start.script &

