# Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
# Parameter file demonstrating the construction of a complex system
# from multiple geometrical objects.  A spherical vesicle and a flat
# bilayer are shown.
#
#
# To run this file use $ESPRESSO_SOURCE/$PLATFORM/Espresso ./scripts/main.tcl sphere-bilayer.tcl
#



::mmsg::send [namespace current] "loading parameter file .. " nonewline
flush stdout

# Set a bunch of default parameters.
set warmup_temp 0
set warmsteps 100
set warmtimes 20
set free_warmsteps 0 
set free_warmtimes 0 
set area_lipid 1.29
set startj 0    
set startk 0
set startmdtime 0
set npt off


# Specify the name of the job <ident> the location of forcetables and
# the output and scripts directories
set ident "sphere-bilayer"
set tabledir "./forcetables/"
set outputdir "./$ident/"
set topofile "$ident.top"

# Espresso Special Parameters #
setmd max_num_cells 2744

# Specify how we want to use vmd for visualization allowable values
# are "interactive" "offline" "none".  interactive rarely works
set use_vmd "offline"

# Set the system box
set setbox_l  { 40.0 40.0 100.0 }

# --------- Specify multi-object system ------- #

# --- Specify a global key for molecule types -----#

# In this line we specify the details of all molecule types in the
# system.  In this example molecule type 0 is assigned to be a lipid.
# The atom types in this lipid are 0 1 and 1 respectively and the
# bonding interactions are type 0 for short bonds and type 1 for long
# bonds.
set moltypes [list { 0 lipid { 0 1 1 } { 0 1 } } { 1 lipid { 0 2 2 2 } { 0 2 } } ]

# --- Specify the system geometry and composition ----#
# Set the geometry to sphere that is offset to 15 sigma above the box center
set geometry { geometry  "sphere -shuffle -c \{ 0.0 0.0 15.0 \} " }

set n_molslist { n_molslist {  { 0 1000 } } }

# Now bundle the above info into a list
lappend spherespec $geometry
lappend spherespec $n_molslist


# Now group the bilayerspec with other specs into a list of such
# systems (we can have multiple systems if we like each with different
# composition of molecule types
lappend system_specs $spherespec



set geometry { geometry "flat -fixz" }
set n_molslist { n_molslist {  { 1 3000 } } }
lappend bilayerspec $geometry
lappend bilayerspec $n_molslist

lappend system_specs $bilayerspec

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
# #----------------------------------------------------------# 
lappend bonded_parms [list 0 FENE 30 1.5 ]
lappend bonded_parms [list 1 harmonic 10.0 4.0 ]
lappend bonded_parms [list 2 harmonic 5.0 4.0 ]

# Non Bonded Potentials 
#----------------------------------------------------------#
lappend nb_interactions [list 0 0 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 1 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 2 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 3 tabulated 9_095_11.tab ]
lappend nb_interactions [list 1 1 tabulated n9_c160_22.tab ]
lappend nb_interactions [list 1 2 tabulated n9_c160_22.tab ]
lappend nb_interactions [list 1 3 tabulated n9_c160_22.tab ]
lappend nb_interactions [list 2 2 tabulated n9_c160_22.tab ]
lappend nb_interactions [list 2 3 tabulated n9_c160_22.tab ]
lappend nb_interactions [list 3 3 tabulated n9_c160_22.tab ]


lappend tablenames 9_095_11.tab
lappend tablenames n9_c160_22.tab 


# Analysis Parameters
#----------------------------------------------------------# 

#Set the size of the 2d grid for calculating the membrane height
#function.  Used to calculate stray lipids lipid flip-flip rates and
#for fluctuation analyses
set mgrid 8

# Distance from the bilayer beyond which a lipid is considered to be
# stray
set stray_cut_off 3

# Use these flags to specify which observables to calculate during the
# simulation.  Values are calculated after every call to the espresso
# integrate command and written to files like
# time_vs_parametername. See the module ::std_analysis for more
# details
lappend analysis_flag pressure
#lappend analysis_flags stress_tensor
lappend analysis_flags boxl
lappend analysis_flags energy


::mmsg::send [namespace current] "done"

















