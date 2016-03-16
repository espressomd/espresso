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

puts "-------------------------------------------"
puts "Testcase h5mdfile.tcl"
puts "-------------------------------------------"

# H5Fopen 		Open file
# H5Gopen2 		Open dataspace"
# H5Dopen2 		Open dataset
# H5Dread		Read dataset from h5-file 
# H5Fcreate 		Create file
# H5Gcreate2		Create group
# H5Dcreate2		Create dataset
# H5Screate_simple	Create dataspace for dataset with dimensions XxYxZ
# H5Pset_chunk		Split dataset in smaller sections XxYxZ 
# H5Dwrite		Write dataset to file
# H5Dextend 		Extend dataset to new dimension XxYxZ
# H5Sselect_hyperslab 	Offset XxYxZ for extendet new dataset
# H5Pclose		Close property
# H5Dclose		Close dataset
# H5Sclose		Close space
# H5Gclose		Close group
# H5Fclose		Close file
# H5_free_memeory	Release allocated memory
# H5_write_value	Write value from Tcl to dataset at position XxYxZ
# H5_read_value		Read value from dataset at position XxYxZ and return to Tcl

#CREATE NEW DATASET AND WRITE VALUES
puts "CREATE NEW DATASET AND WRITE VALUES"
h5mdfile H5Fcreate "h5mdfile.h5" 
h5mdfile H5Gcreate2 "group1"
h5mdfile H5Screate_simple type int dims 10 5 3
h5mdfile H5Pset_chunk dims 2 2 2
h5mdfile H5Dcreate2 "/group1/dset"
h5mdfile H5Dopen2 "/group1/dset"
h5mdfile H5_write_value value 100 index 8 3 0
set E [expr [h5mdfile H5_read_value index 8 3 0]]
puts $E
h5mdfile H5Dwrite

h5mdfile H5Gcreate2 "group2"
h5mdfile H5Screate_simple type int dims 15 5 3
h5mdfile H5Pset_chunk dims 2 2 2
h5mdfile H5Dcreate2 "/group2/dset"
h5mdfile H5Dopen2 "/group2/dset"
h5mdfile H5_write_value value 200 index 8 4 0
set E [expr [h5mdfile H5_read_value index 8 4 0]]
puts $E
h5mdfile H5Dwrite

h5mdfile H5Pclose
h5mdfile H5Dclose
h5mdfile H5Gclose
h5mdfile H5Fclose
h5mdfile H5_free_memory



#WRITE TO EXISITNG DATASET
puts "WRITE TO EXISITNG DATASET"
h5mdfile H5Fopen "h5mdfile.h5"
h5mdfile H5Gopen2 "/group1"
h5mdfile H5Dopen2 "/group1/dset"
h5mdfile H5Dread
h5mdfile H5_write_value value 111 index 9 3 0
set E [expr [h5mdfile H5_read_value index 9 3 0]]
puts $E
h5mdfile H5Dwrite

h5mdfile H5Gopen2 "/group2"
h5mdfile H5Dopen2 "/group2/dset"
h5mdfile H5Dread
h5mdfile H5_write_value value 222 index 9 4 0
set E [expr [h5mdfile H5_read_value index 9 4 0]]
puts $E
h5mdfile H5Dwrite

h5mdfile H5Pclose
h5mdfile H5Dclose
h5mdfile H5Gclose
h5mdfile H5Fclose
h5mdfile H5_free_memory


#EXTEND DATASET
puts "EXTEND DATASET"
h5mdfile H5Fcreate "h5mdfile2.h5" 
h5mdfile H5Gcreate2 "group1"
h5mdfile H5Screate_simple type int dims 10 3 3
h5mdfile H5Pset_chunk dims 2 5 3
h5mdfile H5Dcreate2 "/group1/dset"
h5mdfile H5Dopen2 "/group1/dset"
h5mdfile H5_write_value value 111 index 5 1 0
set E [expr [h5mdfile H5_read_value index 5 1 0]]
puts $E
h5mdfile H5Dwrite

h5mdfile H5Dextend dims 20 3 3
h5mdfile H5Sselect_hyperslab offset 10 0 0
h5mdfile H5Screate_simple type int dims 10 3 3
h5mdfile H5_write_value value 222 index 5 1 0	
set E [expr [h5mdfile H5_read_value index 5 1 0]]
puts $E
h5mdfile H5Dwrite

h5mdfile H5Dextend dims 30 3 3
h5mdfile H5Sselect_hyperslab offset 20 0 0
h5mdfile H5Screate_simple type int dims 10 3 3
h5mdfile H5_write_value value 333 index 5 1 0	
set E [expr [h5mdfile H5_read_value index 5 1 0]]
puts $E
h5mdfile H5Dwrite

h5mdfile H5Pclose
h5mdfile H5Dclose
h5mdfile H5Sclose
h5mdfile H5Gclose
h5mdfile H5Fclose
h5mdfile H5_free_memory

#WRITE STRING
puts "CREATE NEW DATASET AND WRITE STRING"
h5mdfile H5Fcreate "h5mdfile3.h5" 
h5mdfile H5Gcreate2 "group1"
h5mdfile H5Screate_simple type str dims 10
h5mdfile H5Pset_chunk dims 2
h5mdfile H5Dcreate2 "/group1/dset"
h5mdfile H5Dopen2 "/group1/dset"
h5mdfile H5_write_value value "Katze" index 5
set strval [h5mdfile H5_read_value index 5]
puts $strval
h5mdfile H5Dwrite

h5mdfile H5Pclose
h5mdfile H5Dclose
h5mdfile H5Gclose
h5mdfile H5Fclose
h5mdfile H5_free_memory

exit 0