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

# Specify some extra stuff to be written to our vmd_animation script for easier visualisation
lappend vmdcommands "mol color SegName"
lappend vmdcommands "mol representation Lines 1.000000"
lappend vmdcommands "mol selection all"
lappend vmdcommands "mol material Opaque"
lappend vmdcommands "mol addrep 0"
lappend vmdcommands "mol modselect 1 0 segname T003"
lappend vmdcommands "mol modstyle 1 0 CPK 4.000000 0.300000 8.000000 6.000000"

# Specify the name of the job <ident> the location of forcetables and
# the output and scripts directories
set ident "colloid"
set tabledir "./forcetables/"
set outputdir "./$ident/"
set topofile "$ident.top"

# Espresso Special Parameters #
setmd max_num_cells 2744

# Specify how we want to use vmd for visualization allowable values
# are "interactive" "offline" "none".  interactive rarely works
set use_vmd "offline"

# Set the box size
set setbox_l  { 50.0 50.0 50.0 }

# --- Specify a global key for molecule types -----#
# First construct a list of the types of all atoms in the colloid. 
# The atoms forming the shell
for { set i 0 } { $i < 312 } { incr i } {
    lappend sphereatoms 4
}
# And now include those atoms which fill the sphere
for { set i 0 } { $i < 100 } { incr i } {
    lappend sphereatoms 3
}
# Finally put this together to specify the molecule details for a hollowsphere (colloid) and a lipid.
set moltypes [subst { { 0 hollowsphere  { $sphereatoms } { 0 } { 100 } } { 1 lipid { 0 1 2 } { 0 1 } } }]

# --- Specify the system geometry and composition ----#
# Set the geometry to flat
set geometry { geometry "flat -fixz" }
# For the bilayer we use 4000 lipids
set n_molslist { n_molslist { { 1 4000 } } }


# Now bundle the above info into a list
lappend bilayerspec $geometry
lappend bilayerspec $n_molslist

# The center of the colloid
lappend ccenter [expr [lindex $setbox_l 0]*0.5 ] 
lappend ccenter [expr [lindex $setbox_l 1]*0.5 ] 
lappend ccenter [expr [lindex $setbox_l 2]*0.5 + 9.5] 

# Construct a spec for a colloid
set geometry [subst { geometry "singlemol -c \{ $ccenter \}" } ] 
# Use just one colloid
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

# The number of steps to integrate with each call to integrate
set int_steps   1000
# The number of times to call integrate
set int_n_times 10
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
# ------------------------------------------------------------- #
lappend nb_interactions [list 0 0 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 1 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 2 tabulated 9_095_11.tab ]
lappend nb_interactions [list 1 1 tabulated n9_c160_22.tab ]
lappend nb_interactions [list 1 2 tabulated n9_c160_22.tab ]
lappend nb_interactions [list 2 2 tabulated n9_c160_22.tab ]
lappend nb_interactions [list 2 3 tabulated n9_c160_22.tab ]
lappend nb_interactions [list 3 3 tabulated sr_e10_c25.tab ]

# We also need to make a list of all the forcetable filenames so that
# the forcetables can be copied to the working directory.

lappend tablenames 9_095_11.tab 
lappend tablenames n9_c160_22.tab 
lappend tablenames sr_e10_c25.tab

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

# Just for fun we make the head groups attract the membrane while the
# tails repel
set colloidrep [list 4 4 lennard-jones $lj_eps $lj_sigma [expr 1.125*$lj_sigma] [expr 0.25*$lj_eps] $ljoffset ]
set fillshellrep [list 4 3 lennard-jones $lj_eps $lj_sigma [expr 1.125*$lj_sigma] [expr 0.25*$lj_eps] $ljoffset ]
set colloidattr [list 4 0 lennard-jones 4.0 $lj_sigma $lj_cutoff $ljshift $ljoffset ]

lappend nb_interactions $colloidrep 
lappend nb_interactions $fillshellrep 
lappend nb_interactions $colloidattr









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
lappend analysis_flags pressure
#lappend analysis_flags stress_tensor
lappend analysis_flags boxl
#lappend analysis_flags flipflop
lappend analysis_flags energy
#lappend analysis_flags orient_order

# It is not recommended to include fluctuation calculations during the
# simulation since they will crash if a hole appears in the membrane

#lappend analysis_flags { fluctuation_calc 1 1}

::mmsg::send [namespace current] "done"

















