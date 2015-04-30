# Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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
puts "-------------------------------------------"
puts "Testcase h5mdfile.tcl"
puts "-------------------------------------------"


# #Create new dataset and write values
# h5mdfile H5Fcreate "h5mdfile.h5" 
# h5mdfile H5Screate_simple int dims 10 5
# h5mdfile H5Pset_chunk 2 2
# h5mdfile H5Dcreate2 "/dset"
# h5mdfile H5_write_value value 100 index 8 3
# h5mdfile H5Dwrite
# h5mdfile H5Dclose
# h5mdfile H5Sclose
# h5mdfile H5Fclose

# #Write to exisiting dataset
# h5mdfile H5Fopen "h5mdfile.h5"
# h5mdfile H5Dopen2 "/dset"
# h5mdfile H5Dread
# h5mdfile H5_write_value value 111 index 7 2
# set E [expr [h5mdfile H5_read_value value 111 index 7 2]]
# h5mdfile H5Dwrite
# h5mdfile H5Dclose
# h5mdfile H5Fclose

#Extend dataset
h5mdfile H5Fcreate "h5mdfile.h5" 
h5mdfile H5Screate_simple int dims 10 3
h5mdfile H5Pset_chunk 2 5
h5mdfile H5Dcreate2 "/dset"
h5mdfile H5_write_value value 100 index 2 2
h5mdfile H5Dwrite

h5mdfile H5Dextend int dims 20 3
h5mdfile H5Sselect_hyperslab 10 0
h5mdfile H5Screate_simple_ext int dims 10 3
h5mdfile H5_write_value_ext value 200 index 1 2	
h5mdfile H5Dwrite_ext

h5mdfile H5Dextend int dims 30 3
h5mdfile H5Sselect_hyperslab 20 0
h5mdfile H5Screate_simple_ext int dims 10 3
h5mdfile H5_write_value_ext value 200 index 1 2	
h5mdfile H5Dwrite_ext

h5mdfile H5Dclose
h5mdfile H5Sclose
h5mdfile H5Fclose

exit 0
