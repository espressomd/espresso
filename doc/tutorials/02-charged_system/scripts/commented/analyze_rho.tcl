#
# Copyright (C) 2010,2012 The ESPResSo project
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
#############################################################
#                                                           #
#  Ion Density Distributions between the Slit-Pore Walls    #
#                                                           # 
#############################################################

# set up Tcl list for binning the ion densities
set bins 50
set cnt 0
for {set i 0} {$i < $bins} {incr i} {
  lappend avg_rho0 0
  lappend avg_rho1 0
}

# now read in all the files passed to the script
# using the CLI
foreach filename [lrange $argv 1 end] {

  # for each file name, open a stream and read all
  # the blocks
  set f [open $filename "r"]
  while { [blockfile $f read auto] != "eof" } {}
  close $f

  # here we loop over all the particles. Note that 
  # we use the ESPResSo internal max_part command
  # which is set to N-1 (hence the <= operator). 
  # for each particle we check its type and then 
  # make a list of all the z-coordinates of the 
  # positive and negative charges, respectively
  set data0 ""
  set data1 ""
  for {set p 0} {$p <= [setmd max_part]} {incr p} {
    if {[part $p pr type] == 0} {
      lappend data0 [lindex [part $p pr pos] 2]
    } {
      lappend data1 [lindex [part $p pr pos] 2]
    }
  }

  # we now use the bin command to bin the two lists
  # we created in the previous step
  set rho0 [bin -linbins 0.5 [expr $box_lz + 0.5] $bins $data0]
  set rho1 [bin -linbins 0.5 [expr $box_lz + 0.5] $bins $data1]

  # next we add this binned list to the total and
  # we increment the counter
  set avg_rho0 [vecadd $avg_rho0 $rho0]
  set avg_rho1 [vecadd $avg_rho1 $rho1]
  incr cnt
}

# after reading all files, we determine the average
set avg_rho0 [vecscale [expr 1.0/$cnt] $avg_rho0]
set avg_rho1 [vecscale [expr 1.0/$cnt] $avg_rho1]

# and we calculate the location of the bins' centres
set centers [bin -linbins 0.5 [expr $box_lz + 0.5] $bins -binctrwdth]

# finally we output this data to the disk
set plot [open "data/rho.dat" "w"]
puts $plot "\# z rho0(z) rho1(z)"
foreach z $centers rho0 $avg_rho0 rho1 $avg_rho1 {
  puts $plot "[lindex $z 0] $rho0 $rho1"
}
close $plot

# and exit ESPResSo
exit
