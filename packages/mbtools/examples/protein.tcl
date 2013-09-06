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
# Parameter file for simulating proteins in
# a membrane
#
#
# To run this file use $ESPRESSO_SOURCE/$PLATFORM/Espresso ./scripts/main.tcl protein.tcl
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
set ident "protein"
set tabledir "./forcetables/"
set outputdir "./$ident/"
set topofile "$ident.top"

# Specify how we want to use vmd for visualization allowable values
# are "interactive" "offline" "none".  interactive rarely works
set use_vmd "offline"

# Set the box size
set setbox_l  { 25.0 25.0 25.0 }

# --- Specify a global key for molecule types -----#

# In this line proteins are assigned type 1 and lipids type
# 0. Proteins are also specified to be constructed from atoms of
# type 0, 1, 2  and lipids 3 atoms of types 3, 4 respectively

set moltypes [subst { { 1 protein  { 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 } { 0 1 2 3 4 5 6 7 8 9 10 11 12 } } { 0 lipid { 4 5 5 5 5 } { 13 14 } } }]

#set moltypes [subst { { 1 protein  { 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 } { 0 1 2 3 4 5 6 7 } } { 0 lipid { 3 4 4 } { 8 9 } } }]

# --- Specify the system geometry and composition ----#
# Set the geometry to single molecule

lappend ccenter [expr [lindex $setbox_l 0]*0.5 ] 
lappend ccenter [expr [lindex $setbox_l 1]*0.5 ] 
lappend ccenter [expr [lindex $setbox_l 2]*0.5 ] 

#set geometry [subst { geometry { singlemol { $ccenter } } }]
set geometry { geometry  "flat -fixz" } 


# In this line we specify that 1 molecule of type 1 (ie protein) and
# 1020 of type 0 (ie lipid) are to be used
set n_molslist { n_molslist { { 1 1 } { 0 1 } } }


# Now bundle the above info into a list
lappend proteinspec $geometry
lappend proteinspec $n_molslist

#set geometry { geometry flat }
#set n_molslist { n_molslist { { 0 1000 } } }

#lappend bilayerspec $geometry
#lappend bilayerspec $n_molslist

# Now group the bilayerspec with other specs into a list of such
# systems (we can have multiple systems if we like each with different
# composition of molecule types
lappend system_specs $proteinspec
#lappend system_specs $bilayerspec


# Warmup parameters
#----------------------------------------------------------#
set warm_time_step 0.002

set free_warmsteps 0
set free_warmtimes 1 


# ---- NPT parameters --------- #
set npt on
set p_ext 0.000
set piston_mass 0.0005
set gamma_0 1.0
set gamma_v 0.0001

# ------ Integration parameters -----------------#
set main_time_step 0.01
set verlet_skin 0.4
set langevin_gamma 1.0
set systemtemp 1.4

# The number of steps to integrate with each call to integrate
set int_steps   1000
# The number of times to call integrate
set int_n_times 100
# Write a configuration file every <write_frequency> calls to
# integrate
set write_frequency 1
# Write results to analysis files every <analysis_write_frequency>
# calls to integrate
set analysis_write_frequency 100

# Bonded and bending Potentials
#----------------------------------------------------------#

lappend bonded_parms [list 0 harmonic 200.0 1.0 ]
lappend bonded_parms [list 1 harmonic 200.0 1.1 ]
lappend bonded_parms [list 2 harmonic 200.0 1.2 ]
lappend bonded_parms [list 3 harmonic 200.0 1.3 ]
lappend bonded_parms [list 4 harmonic 200.0 1.4 ]
lappend bonded_parms [list 5 harmonic 200.0 1.5 ]
lappend bonded_parms [list 6 harmonic 200.0 1.6 ]
lappend bonded_parms [list 7 harmonic 200.0 1.7 ]
lappend bonded_parms [list 8 harmonic 200.0 1.8 ]
lappend bonded_parms [list 9 harmonic 200.0 1.9 ]
lappend bonded_parms [list 10 harmonic 200.0 2.0 ]
lappend bonded_parms [list 11 harmonic 200.0 2.1 ]
lappend bonded_parms [list 12 harmonic 100.0 10 ]
lappend bonded_parms [list 13 FENE 30 1.5 ]
lappend bonded_parms [list 14 harmonic 10 4.0 ]


# Non Bonded Potentials 
#----------------------------------------------------------#
#
lappend nb_interactions [list 0 0 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 1 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 2 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 3 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 4 tabulated 9_095_11.tab ]

lappend nb_interactions [list 1 1 tabulated n9_c160_22.tab ]
lappend nb_interactions [list 1 2 tabulated 9_095_11.tab ]
lappend nb_interactions [list 1 3 tabulated 9_095_11.tab ]
lappend nb_interactions [list 1 4 tabulated n9_c160_22.tab ]

lappend nb_interactions [list 2 2 tabulated 9_095_11.tab ]
lappend nb_interactions [list 2 3 tabulated 9_095_11.tab ]
lappend nb_interactions [list 2 4 tabulated 9_095_11.tab ]

lappend nb_interactions [list 3 3 tabulated 9_095_11.tab ]
lappend nb_interactions [list 3 4 tabulated 9_095_11.tab ]

lappend nb_interactions [list 4 4 tabulated n9_c160_22.tab ]

lappend tablenames 9_095_11.tab 
lappend tablenames n9_c160_22.tab 


# Analysis Parameters
#----------------------------------------------------------# 
# Use these flags to specify which observables to calculate during the
# simulation.  Values are calculated after every call to the espresso
# integrate command and written to files like
# time_vs_parametername. See the module analysis for more
# details
lappend analysis_flags pressure
#lappend analysis_flags stress_tensor
lappend analysis_flags boxl
#lappend analysis_flags flipflop
lappend analysis_flags energy
#lappend analysis_flags orient_order

::mmsg::send [namespace current] "done"

















