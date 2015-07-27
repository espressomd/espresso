# Copyright (C) 2010,2012,2013,2014 The ESPResSo project
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
# To run this file use $ESPRESSO_SOURCE/$PLATFORM/Espresso ./scripts/main.tcl colloid.tcl
#
# For more information about the simulation model on which this
# simulation is based see the following reference
# 
#  Cooke, I. R., Kremer, K. and Deserno, M.(2005): Tuneable, generic
# model for fluid bilayer membranes. Phys. Rev. E. 72 - 011506 
#

::mmsg::send [namespace current] "loading parameter file .. " nonewline
flush stdout

# Specify some extra stuff to be written to our vmd_animation script
lappend vmdcommands "mol color SegName"
lappend vmdcommands "mol representation Lines 1.000000"
lappend vmdcommands "mol selection all"
lappend vmdcommands "mol material Opaque"
lappend vmdcommands "mol addrep 0"
lappend vmdcommands "mol modselect 1 0 segname T003"
lappend vmdcommands "mol modstyle 1 0 CPK 4.000000 0.300000 8.000000 6.000000"

# Specify the name of the job <ident> the location of forcetables and
# the output and scripts directories
set ident "januscolloid"
set tabledir "./forcetables/"
set outputdir "./$ident/"
set topofile "$ident.top"


# Specify how we want to use vmd for visualization allowable values
# are "interactive" "offline" "none".  interactive rarely works
set use_vmd "offline"

# Set the box size
set setbox_l  { 22.0 22.0 30.0 }

# --- Specify a global key for molecule types -----#

# In this line colloids are assigned type 0 and lipids type
# 1. Colloids are also specified to be constructed from 4 atoms of
# type 0 and lipids 3 atoms of types 0 1 and 2 respectively
for { set i 0 } { $i < 212 } { incr i } {
    lappend sphereatoms 4
}
for { set i 212 } { $i < 312 } { incr i } {
    lappend sphereatoms 5
}
# And now include those atoms which fill the sphere
for { set i 0 } { $i < 100 } { incr i } {
    lappend sphereatoms 3
}

set moltypes [subst { { 0 hollowsphere  { $sphereatoms } { 0 } { 100 } } { 1 lipid { 0 1 2 } { 0 1 } } }]

# --- Specify the system geometry and composition ----#
# Set the geometry to flat
set geometry { geometry "flat -fixz" }

# In this line we specify that 1 molecule of type 0 (ie colloid) and
# 20 of type 1 (ie lipid) are to be used
set n_molslist { n_molslist { { 1 1000 } } }


# Now bundle the above info into a list
lappend bilayerspec $geometry
lappend bilayerspec $n_molslist

lappend ccenter [expr [lindex $setbox_l 0]*0.5 ] 
lappend ccenter [expr [lindex $setbox_l 1]*0.5 ] 
lappend ccenter [expr [lindex $setbox_l 2]*0.5 + 13.3] 


# Construct a spec for a colloid
# Since the colloid has two different regions it has an orientation.  
set orient { 0.50 1.0 -1.0 }
set geometry [subst { geometry "singlemol -c \{ $ccenter \} -o \{ $orient \} "  } ]
set n_molslist [list n_molslist { { 0 1 } } ]
lappend colloidspec $geometry
lappend colloidspec $n_molslist


# Now group the bilayerspec with other specs into a list of such
# systems (we can have multiple systems if we like each with different
# composition of molecule types
lappend system_specs $bilayerspec
lappend system_specs $colloidspec


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
set npt on
set p_ext 0.000
set piston_mass 0.0005
set gamma_0 1.0
set gamma_v 0.0001

# The number of steps to integrate with each call to integrate
set int_steps   500
# The number of times to call integrate
set int_n_times 20
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
lappend tablenames 9_095_11.tab 
lappend tablenames n9_c160_22.tab 
lappend tablenames sr_e10_c25.tab

# ------------------------------------------------------------- #
# Lipid Interactions #
lappend nb_interactions [list 0 0 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 1 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 2 tabulated 9_095_11.tab ]
lappend nb_interactions [list 1 1 tabulated n9_c160_22.tab ]
lappend nb_interactions [list 1 2 tabulated n9_c160_22.tab ]
lappend nb_interactions [list 2 2 tabulated n9_c160_22.tab ]

# Filler beads interact with themselves
lappend nb_interactions [list 3 3 tabulated sr_e10_c25.tab ]

# Other non-bonded potentials
# Here we specify non-bonded interactions that are not covered by our
# tabulated potentials.  In this case we specify lennard jones
# repulsions for the beads that make up our spherical colloid.
#
set lj_eps 1.0
# We make lj_sigma a bit larger than the bond length so that there is a strong repulsion
set lj_sigma 1.0
set lj_cutoff 2.5 
set ljshift 0.0
set ljoffset 0.0

# Here are the interactions for the attractive part of the colloid surface
lappend nb_interactions [list 4 4 lennard-jones $lj_eps $lj_sigma [expr 1.125*$lj_sigma] [expr 0.25*$lj_eps] $ljoffset ]
lappend nb_interactions [list 4 3 lennard-jones $lj_eps $lj_sigma [expr 1.125*$lj_sigma] [expr 0.25*$lj_eps] $ljoffset ]
lappend nb_interactions [list 4 0 lennard-jones 1.0 $lj_sigma $lj_cutoff $ljshift $ljoffset ]

# And the non adhesive part
lappend nb_interactions [list 5 5 lennard-jones $lj_eps $lj_sigma [expr 1.125*$lj_sigma] [expr 0.25*$lj_eps] $ljoffset ]
lappend nb_interactions [list 5 3 lennard-jones $lj_eps $lj_sigma [expr 1.125*$lj_sigma] [expr 0.25*$lj_eps] $ljoffset ]

# And the cross interaction between them
lappend nb_interactions [list 4 5 lennard-jones $lj_eps $lj_sigma [expr 1.125*$lj_sigma] [expr 0.25*$lj_eps] $ljoffset ]

# Analysis Parameters
#----------------------------------------------------------# 

# Use these flags to specify which observables to calculate during the
# simulation.  Values are calculated after every call to the espresso
# integrate command and written to files like
# time_vs_parametername. See the module ::std_analysis for more
# details
lappend analysis_flags pressure
#lappend analysis_flags stress_tensor
lappend analysis_flags boxl
lappend analysis_flags energy

::mmsg::send [namespace current] "done"

















