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

puts "----------------------------------------------"
puts "- Testcase pdb_parser.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------------"

setmd box_l 1. 1. 1.

part 0 pos 0 0 0 type 0

set n_part [readpdb pdb_file pdb_parser.pdb.dat type 10 first_id 1 itp_file pdb_parser.itp.dat rescale_box lj_rel_cutoff 2.5 first_type 3 lj_with 0 1.0 1.0]

if { $n_part != 31 } {
    error_exit "Wrong number of particles"
}

exit 0
