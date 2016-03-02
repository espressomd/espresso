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
# Parameter file for simulating a one component random fluid lipid.
# Useful for studying the self-aggregation procedure using the DPD thermostat.
#
#
# To run this file use $ESPRESSO_SOURCE/$PLATFORM/Espresso ./scripts/main.tcl random_DPD.tcl
#

# For more information about the simulation model on which this
# simulation is based see the following reference
# 
#  Cooke, I. R., Kremer, K. and Deserno, M.(2005): Tuneable, generic
# model for fluid bilayer membranes. Phys. Rev. E. 72 - 011506 
#

puts " "
puts "============================================================"
puts " Parameter Script for the Simulation of a Random Lipid System"
puts "============================================================"

::mmsg::send [namespace current] "loading parameter file .. " nonewline
flush stdout

# Specify the name of the job <ident> the location of forcetables and
# the output and scripts directories
set ident "random_DPD"
set tabledir "./forcetables/"
set outputdir "./$ident/"
set topofile "$ident.top"

# Specify how we want to use vmd for visualization allowable values
# are "interactive" "offline" "none".  interactive rarely works
set use_vmd "offline"

# --- Specify a global key for molecule types -----#
set moltypes [list { 0 lipid { 0 1 1 } { 0 1 } } ]

# --- Specify the system geometry and composition ----#
set geometry { geometry random }

# Set number of molecjules of each type
set n_molslist { n_molslist {  { 0 2000 } } }

# Now bundle the above info into a list
lappend bilayerspec $geometry
lappend bilayerspec $n_molslist

# Now group the bilayerspec with other specs into a list of such
# systems (we can have multiple systems if we like each with different
# composition of molecule types
lappend system_specs $bilayerspec


# Specify the system  size
set setbox_l  { 30.0 30.0 30.0 }


# Warmup parameters
#----------------------------------------------------------#
set warm_time_step 0.002
set free_warmsteps 0
set free_warmtimes 1

# ------ Integration parameters -----------------#
set main_time_step 0.01
set verlet_skin 0.4
set systemtemp 1.1
# NVT parameters
#set thermo "Langevin"		;# "Langevin" or "DPD"
set thermo "DPD"		;# "Langevin" or "DPD"
set langevin_gamma 1.0		;# Langevin thermostat parameter
set dpd_gamma 1.0		;# DPD thermostat parameter "gamma"
set dpd_r_cut [expr 1.12 + 1.8]	;# D_Rc should be > Rc + wc

# The number of steps to integrate with each call to integrate
set int_steps   1000
# The number of times to call integrate
set int_n_times 1000
# Total Integration time (in taus)
set total_int_time [expr $main_time_step*$int_n_times*$int_steps ]
puts "Total Integration Time: $total_int_time"
# Write a configuration file every <write_frequency> calls
set write_frequency 10
# Write results to analysis files every <analysis_write_frequency> calls
set analysis_write_frequency 5

# Bonded and bending Potentials
#----------------------------------------------------------#
lappend bonded_parms [list 0 FENE 30 1.5 ]
lappend bonded_parms [list 1 harmonic 10.0 4.0 ]

# Non Bonded Potentials
#----------------------------------------------------------#
#Define the forcetables to be used for the interaction of each pair of atom types
lappend nb_interactions [list 0 0 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 1 tabulated 9_095_11.tab ]
lappend nb_interactions [list 1 1 tabulated n9_c160_22.tab ]

# We also need to make a list of all the forcetable filenames so that
# the forcetables can be copied to the working directory.
lappend tablenames 9_095_11.tab
lappend tablenames n9_c160_22.tab


# Analysis Parameters
#----------------------------------------------------------#

# Use these flags to specify which observables to calculate during the
# simulation.  Values are calculated after every call to the espresso
# integrate command and written to files like
# time_vs_parametername. See the module ::std_analysis for more
# details
lappend analysis_flags boxl	;# calculate box length
lappend analysis_flags energy ;# calculate energy
lappend analysis_flags pressure ;# pressure calculation
lappend analysis_flags stress_tensor	;# stress tensor calculation
lappend analysis_flags clusters	;# analyzing "clusters"

