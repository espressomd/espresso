#
# Copyright (C) 2012,2013 The ESPResSo project
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
# Measure the amount of time, by which certain features slow down Espresso
# Use myconfig-basic.h and one of the features from featurelist.

# Path to espresso source tree
esdir=~/espresso-dev
# Path to Espresso binary
es=./Espresso

# Scripts to be run to test the performance
$ Long-running script for timing accurately
scriptLong=../lj_liquid_long.tcl
# Short-running script to get profiling information
scriptShort=../lj_liquid_short.tcl

# Iterate over the features
for f in `cat featurelist`
do
# Generate the myconfig.h with the feature
cp myconfig-basic.h myconfig.h
echo "#define $f" >> myconfig.h

# Prepare the build directory
mkdir build
cd build
cp ../myconfig.h .

# Build for profiling
$esdir/configure CFLAGS='-O5 -pg' 
make -j20

# run
$es $scriptShort > ../$f.out 2> ../$f.err
# Save the profile information
gprof Espresso > ../$f.prof

# Build for optimal performance without profiling
make distclean
$esdir/configure 
make -j20
/usr/bin/time -f %e -o ../$f.time  $es $scriptShort > ../$f.out 2> ../$f.err


cd ..

# Cleanup
rm -rf build
done


