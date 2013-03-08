# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
#
#  Offline Analysis
#

set filename "data/sim_info.dat"
set in [open "$filename" "r"]

while {  [set btitle [ blockfile $in read auto] ] != "eof" } {
 if {$btitle == "variable"} {
     set times [setmd time]
    }
 if { $btitle == "particles"}   {
           # 
           # Here We have an access to all simulation information on given time 
	   # 
     set part0 [part 0 print pos]
     set part0v [part 0 print v]
     puts "time=$times position of p0= $part0 velo = $part0v"
    }
}
close $in
exit
