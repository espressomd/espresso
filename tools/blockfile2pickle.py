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

### This script extracts tcl blockfiles and pickles the appropriate py class.

import blockfile as bf
import espressomd
import espressomd._system as es
from espressomd.particle_data import ParticleList
from espressomd.interactions import NonBondedInteractions, BondedInteractions
from sys import argv
import pickle

def safe_str2intflt(string):
    """ Takes a string argument and exits if the argument is not a string.
        Returns an integer, if conversion is possible, otherwise returns
        a python float, otherwise return the original string """
    if not isinstance(string, str):
        raise TypeError("I want strings, not " + type(string) + "exiting.")
    try: return int(string)
    except ValueError: pass
    try: return float(string)
    except ValueError: return string


# nagnagnag
if len(argv) < 3:
    print "Usage:  ./pypresso ./blockfile2pickle.py <blockfile> <output-prefix>"
    print "Output: create files named <output-prefix>-<classname>.pcl in workdir"
    print "Example output: lj_fluid-ParticleList.pcl lj_fluid-NonBondedInteractions.pcl"
    print "Supports: particle, interaction"
    print "Depends:  Espresso/tools/blockfile.py file in a sys.path dir"

# initialize blockfile class from blockfile
blocks = bf.open(argv[1])
b_list = []

# create a list of tuples which contains block data
# Each element returned by blocks' __iter__ is a tuple for a tcl block
for i in blocks:
    b_list.append(i)

# close the file. for this object, this is handled by __del__
del blocks

#loop over every block
for k in b_list:
    #turn tuple into list
    k = list(k)
    # where k[0] is a string of the block type
    # and k[1] is a dict for every block argument.
    
    # if k is a particle block
    if k[0] == 'particles':
        print "Particle Block found"
        if "pos" not in k[1] or "id" not in k[1]:
            print "Error, id and/or pos missing in particle block, skipping particle"
            continue # the continue statement jumps to the next iteration.

        print "Particle properties found:"
        print k[1].keys()
        # instance of the system internal particle list
        pl = ParticleList()
        # Set appropriate params foreach particle, since pl.add is so nice,
        # we just have to pass the k[1] dict.
        try: pl.add( k[1] )
        except AttributeError: 
            print "Error, the blockfile sets a particle parameter which is not compiled in."
            print "To properly pickle a particle, all necessary features need to be compiled."
            print "For example, if the particle property 'q' is set, we need ELECTROSTATICS."
            print "skipping particle"
            continue
        # pickle it
        with open(argv[2] + "-ParticleList.pcl", 'w') as outfile:
            print "Pickling ParticleList"
            pickle.dump(pl, outfile, -1)
            del pl


    # if k is an interactions block
    elif k[0] == 'interactions' or k[0] == 'inter':
        # here, the blockfile tool returns a string containing ugly escapes.
        # this means we have to do some data processing.
        if not isinstance(k[1], str):
            print "Unexpected error, bf.blockfile did not return string for inter"
            print "skipping interaction block"
            continue
        # split block into strings for every interaction directive and
        # get rid of unnecessary characters
        k[1] = [words.strip().strip('{').strip('}')
                for words in k[1].split('\n')]
        # split interaction strings into list of arguments and convert
        # every numeric string into an integer or float.
        # Then filter every empty string.
        # Do that for every command invocation in block.
        k[1] = [filter(lambda x: x!='',[safe_str2intflt(word) for word in words.split(' ')]) 
                            for words in k[1]] 

        # We need to differentiate between Bonded and Nonbonded
        # The difference is: Nonbonded starts with two integers as part id,
        # and Bonded starts with one integer as bond id.
        # so: Nonbonded looks like this: [<numbers> <numbers> <letters> ...]
        # bonded looks like this: [<numbers> <letters> ...]

        #initialize system internal interaction objects
        bia = BondedInteractions()
        nia = NonBondedInteractions()
        # keep track of whether or not nia or bia where actually used to init
        # an interaction.
        nia_success = False
        bia_success = False
        for arg_list in k[1]:
            n_args = len(arg_list)
            if (isinstance(arg_list[0], int) and isinstance(arg_list[1], int) 
                                            and isinstance(arg_list[2], str)):
                # NonBondedInteractions.
                # figure out what interaction in particular we are looking at.
                if arg_list[2] == 'lennard-jones':
                    # Here, we found lennard-jones interaction. 
                    # In this indentation level we write interaction specific
                    # parsing code.
                    print "Found lennard-jones interaction for particles"
                    # test whether or not there are enough arguments for
                    # the required parameters of lj
                    if n_args < 6:
                        print "Not enough arguments for lennard-jones, skipping"
                        continue # The continue statement jumps to the next
                                 # iteration of the innermost for-loop
                    # create input dict for the lj kwargs
                    in_dict =   {   
                                    "epsilon": arg_list[3],
                                    "sigma"  : arg_list[4],
                                    "cutoff" : arg_list[5]
                                }
                    # optional arguments
                    if 6 < n_args: in_dict["shift"]  = arg_list[6]
                    if 7 < n_args: in_dict["offset"] = arg_list[7]
                    if 8 < n_args: in_dict["min"]    = arg_list[8]
                    if 9 < n_args: 
                        print "Error: too many arguments for lennard-jones, skipping"
                        continue
                    # initialization of the interaction with nia.
                    # Also: since set_params doesn't natively support dicts as *args
                    #       we need to expand the in_dict as **kwargs, this is done
                    #       by passing it to a function as **in_dict
                    nia[arg_list[0], arg_list[1]].lennard_jones.set_params(**in_dict)
                    nia_success = True # set this to true if nia actually did something
                    
                    # be nice to the user
                    print "\nSetting up lennard-jones with the following params:"
                    print in_dict
                    print "for id %d and %d" %(arg_list[0], arg_list[1])
                    print ""
            
                # HERE is the spot where you can implement the initialization
                # for an individual NonBondedInteraction. One elif-statement for every
                # type of interaction. No pickledump, this is already taken care of,
                # just parse arg_list and init with nia. see above (lennard-jones) for
                # an example. Don't forget to set nia_success to True if initialization 
                # actually happens.
                elif arg_list[2] == "MyNonBondedInteraction":
                    pass
                else:
                    print "NonBondedInteraction " + arg_list[2] + " is unknown, skipping"
                    continue
        
            elif isinstance(arg_list[0], int) and isinstance(arg_list[1], str):
            #BondedInteractions

            # HERE is the spot where you can implement the parser for BondedInteractions.
            # But instead of using nia and nia_success, use bia and bia_success
            # accordingly. See NonBondedInteractions above, specifically lennard-jones.
                if arg_list[1] == "MyBondedInteraction":
                    pass
                
                else:
                    print "BondedInteraction" + arg_list[1] + " is unknown, skipping"
                    continue
            else:
                print "Unknown interaction format, skipping"
        
        # end of arg_list for loop
        if nia_success:
            with open(argv[2] + "-NonBondedInteractions.pcl", "w") as outfile:
                pickle.dump(nia, outfile, -1)
                print "pickling NonBondedInteractions"
            del nia, nia_success
        if bia_success:
            with open(argv[2] + "-BondedInteractions.pcl", "w") as outfile:
                pickle.dump(bia, outfile, -1)
                print "pickling NonBondedInteractions"
            del bia, bia_success
    elif k[0] == "MyBlockIdentifier":
        pass

    else:
        print "blocktype " + k[0] + " unknown, skipping"


