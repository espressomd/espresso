#/usr/bin/env python2
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
### This script extracts particle and interaction data from blockfiles ###
import blockfile as bf
import espressomd
import espressomd._system as es
from espressomd.particle_data import ParticleList
from sys import argv
import pickle

# nagnagnag
if len(argv) < 3:
    print "Usage:  ./pypresso ./blockfile2pickle.py <blockfile> <output-prefix>"
    print "Output: <output-prefix>-<blocktype>.pcl"
    print "Example output: lj_fluid-particle.pcl lj_fluid-interactions.pcl"
    print "needs Espresso/tools/blockfile.py in the same dir"

# initialize blockfile class from blockfile
blocks = bf.open(argv[1])
b_list = []

# create a list of tuples which contains block data
# Each element returned by blocks' __iter__ is a tuple for a tcl block
for i in blocks:
    b_list.append(i)

#loop over every block
for k in b_list:
    # where k[0] is a string of the block type
    # and k[1] is a dict for every block argument.
    
    # if k is a particle block
    if k[0] == 'particles':
        if "pos" not in k[1] or "id" not in k[1]:
            print "Error, id and/or pos missing in particle block"
            exit(1)
        # instance of the internal particle list
        pl = ParticleList()
        # Set appropriate params foreach particle
        pl.add( k[1] )
        # pickle it
        with open(argv[2] + "-" + k[0] + ".pcl", 'w') as outfile:
            pickle.dump(pl, outfile, -1)

    elif k[0] == 'interactions':
        pass
    else:
        print "blockfile type " + k[0] + " not implemented, skipping"


