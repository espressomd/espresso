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
# Parameter file for simulating a colloidal particle interacting with
# a membrane
#
#
# To run this file use $ESPRESSO_SOURCE/$PLATFORM/Espresso ./scripts/main.tcl readfromfile.tcl
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
set ident "readfromfile"
set tabledir "./forcetables/"
set outputdir "./$ident/"
set topofile "$ident.top"

# Specify how we want to use vmd for visualization allowable values
# are "interactive" "offline" "none".  interactive rarely works
set use_vmd "offline"

# Set the box size to the size in the file
set setbox_l  { 102.852650334 102.852650334 15.0 }

# --- Specify a global key for molecule types -----#

set moltypes [subst { { 1 lipid { 0 1 2 } { 0 1 } } } ]

# --- Specify the system geometry and composition ----#
# Set the geometry to flat
set geometry { geometry "readfile -f ./readfile_cfg.gz -t ./readfile_cfg.top" }

# In this line we specify that 1 molecule of type 0 (ie colloid) and
# 20 of type 1 (ie lipid) are to be used
set n_molslist { n_molslist { { 1 10 } } }


# Now bundle the above info into a list
lappend bilayerspec $geometry
lappend bilayerspec $n_molslist

# Now group the bilayerspec with other specs into a list of such
# systems (we can have multiple systems if we like each with different
# composition of molecule types
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

# ---- NPT parameters --------- #
set npt off
set p_ext 0.000
set piston_mass 0.0005
set gamma_0 1.0
set gamma_v 0.0001

# The number of steps to integrate with each call to integrate
set int_steps   2000
# The number of times to call integrate
set int_n_times 100
# Write a configuration file every <write_frequency> calls to
# integrate
set write_frequency 1
# Write results to analysis files every <analysis_write_frequency>
# calls to integrate
set analysis_write_frequency 1

# Bonded and bending Potentials
#----------------------------------------------------------#
lappend bonded_parms [list 0 FENE 30 1.5 ]
lappend bonded_parms [list 1 harmonic 10.0 4.0 ]


# Non Bonded Potentials 
#----------------------------------------------------------#
lappend nb_interactions [list 0 0 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 1 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 2 tabulated 9_095_11.tab ]
lappend nb_interactions [list 1 1 tabulated n9_c160_22.tab ]
lappend nb_interactions [list 1 2 tabulated n9_c160_22.tab ]
lappend nb_interactions [list 2 2 tabulated n9_c160_22.tab ]

lappend tablenames 9_095_11.tab
lappend tablenames n9_c160_22.tab 


# Analysis Parameters
#----------------------------------------------------------# 

# Use these flags to specify which observables to calculate during the
# simulation.  Values are calculated after every call to the espresso
# integrate command and written to files like
# time_vs_parametername. See the module ::std_analysis for more
# details
lappend analysis_flags pressure
lappend analysis_flags stress_tensor
lappend analysis_flags boxl
lappend analysis_flags energy


















