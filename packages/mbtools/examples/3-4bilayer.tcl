# Copyright (C) 2010,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
#  
# This file is part of ESPResSo.
#   
# ESPResSo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Parameter file for simulating a two component flat lipid bilayer
# where one of the lipids is a different length from the other at
# fixed box size and temperature.
#
#
# To run this file use $ESPRESSO_SOURCE/$PLATFORM/Espresso ./scripts/main.tcl 3-4bilayer.tcl
#

# For more information about the simulation model on which this
# simulation is based see the following reference
# 
#  Cooke, I. R., Kremer, K. and Deserno, M.(2005): Tuneable, generic
# model for fluid bilayer membranes. Phys. Rev. E. 72 - 011506 
#

::mmsg::send [namespace current] "loading parameter file .. " nonewline
flush stdout

# Specify the name of the job <ident> the location of forcetables and
# the output and scripts directories
set ident "3-4bilayer"
set tabledir "./forcetables/"
set outputdir "./$ident/"
set topofile "$ident.top"

# Specify how we want to use vmd for visualization allowable values
# are "interactive" "offline" "none".  interactive rarely works
set use_vmd "offline"

# --- Specify a global key for molecule types -----#
# Here we specify two different lipid types. Both have head bead type
# 0 but one has tail beads of type 1 and the other has tail beads of
# type 2.  Note that the bonding interactions are also different and
# the lengths of the two lipids are also different.  In both cases
# consecutive atoms are bonded by FENE (type 0 ) but the bending
# interactions are different (type 1 and 2 respectively).  See the
# specification of bonded interactions later in this file.
set moltypes [list { 0 lipid { 0 1 1 1 } { 0 1 } } { 1 lipid { 0 2 2 } { 0 2 } }  ]

# --- Specify the system geometry and composition ----#
# Set the geometry to flat
set geometry { geometry "flat -fixz" }

# In this line we specify that 170 lipids of each type are to be used
set n_molslist { n_molslist {  { 0 170 } { 1 170 } } }

# Now bundle the above info into a list
lappend bilayerspec $geometry
lappend bilayerspec $n_molslist

# Now group the bilayerspec with other specs into a list of such
# systems (we can have multiple systems if we like each with different
# composition of molecule types
lappend system_specs $bilayerspec

set setbox_l  { 15.0 15.0 15.0 }


# Warmup parameters
#----------------------------------------------------------#
set warm_time_step 0.002

set free_warmsteps 0
set free_warmtimes 1 

# ------ Integration parameters -----------------#
set main_time_step 0.01
set verlet_skin 0.4
set langevin_gamma 1.0
set systemtemp 1.1

# The number of steps to integrate with each call to integrate
set int_steps   100
# The number of times to call integrate
set int_n_times 100
# Write a configuration file every <write_frequency> calls to
# integrate
set write_frequency 10
# Write results to analysis files every <analysis_write_frequency>
# calls to integrate
set analysis_write_frequency 1

# Bonded and bending Potentials
#----------------------------------------------------------#
lappend bonded_parms [list 0 FENE 30 1.5 ]
lappend bonded_parms [list 1 harmonic 5.0 4.0 ]
lappend bonded_parms [list 2 harmonic 10.0 4.0 ]


# Non Bonded Potentials
#----------------------------------------------------------#
# Here we specify the full interaction matrix which should be
#
#  00  01  12  
#      11  12  
#          22  
#
# Since 00 , 01 and 12 are all the same repulsive lennard jones we use
lappend nb_interactions [list 0 0 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 1 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 2 tabulated 9_095_11.tab ]

# Then we specify the interactions of bead type 1 with other tail
# beads.  Note that the cross term (23) has a smaller c value
# (potential width) which gives rise to a line tension.
lappend nb_interactions [list 1 1 tabulated n9_c160_22.tab ]
lappend nb_interactions [list 1 2 tabulated n9_c140_22.tab ]

# And finally the interaction of tail beads type 2 with themselves
lappend nb_interactions [list 2 2 tabulated n9_c160_22.tab ]

# We also need to make a list of all the forcetable filenames so that
# the forcetables can be copied to the working directory.

lappend tablenames 9_095_11.tab 
lappend tablenames n9_c160_22.tab 
lappend tablenames n9_c140_22.tab 

# Analysis Parameters
#----------------------------------------------------------# 
# Use these flags to specify which observables to calculate during the
# simulation.  Values are calculated after every call to the espresso
# integrate command and written to files like
# time_vs_parametername. See the module analysis for more
# details
lappend analysis_flags orient_order
lappend analysis_flags pressure 
lappend analysis_flags stress_tensor  
lappend analysis_flags boxl
lappend analysis_flags flipflop
lappend analysis_flags energy
lappend analysis_flags stray
lappend analysis_flags fluctuations
set profile_beadtypes [list  0 1 2 ]
lappend analysis_flags "density_profile -beadtypes \{ $profile_beadtypes \} -nogrid"

::mmsg::send [namespace current] "done"

















