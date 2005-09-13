#
# Parameter file for simulating a one component flat lipid bilayer
# at fixed box size and temperature.
#
#
# To run this file use $ESPRESSO_SOURCE/$PLATFORM/Espresso ./scripts/main.tcl simplebilayer.tcl
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
set ident "simplebilayer"
set tabledir "./forcetables/"
set outputdir "./$ident/"
set topofile "$ident.top"


# Specify how we want to use vmd for visualization: allowable values
# are "interactive" "offline" "none".  interactive rarely works
set use_vmd "offline"

# --- Specify a global key for molecule types -----#

# In this line we specify the details of all molecule types in the
# system.  In this example molecule type 0 is assigned to be a lipid.
# The atom types in this lipid are 0 1 and 2 respectively and the
# bonding interactions are type 0 for short bonds and type 1 for long
# bonds.
set moltypes [list { 0 lipid { 0 1 2 } { 0 1 } } ]

# --- Specify the system geometry and composition ----#
# Set the geometry to flat
set geometry { geometry "flat -fixz" }

# In this line we specify that 1 molecule of type 0 (ie colloid) and
# 20 of type 1 (ie lipid) are to be used
set n_molslist { n_molslist {  { 0 360 } } }

# Now bundle the above info into a list
lappend bilayerspec $geometry
lappend bilayerspec $n_molslist


# Now group the bilayerspec with other specs into a list of such
# systems (we can have multiple systems if we like each with different
# composition of molecule types
lappend system_specs $bilayerspec

# Set the box size
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
lappend bonded_parms [list 1 harmonic 10.0 4.0 ]

# Non Bonded Potentials 
#----------------------------------------------------------#
# Define the forcetables to be used for the interaction of each pair of atom types

lappend nb_interactions [list 0 0 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 1 tabulated 9_095_11.tab ]
lappend nb_interactions [list 0 2 tabulated 9_095_11.tab ]
lappend nb_interactions [list 1 1 tabulated n9_c160_22.tab ]
lappend nb_interactions [list 1 2 tabulated n9_c160_22.tab ]
lappend nb_interactions [list 2 2 tabulated n9_c160_22.tab ]

# We also need to make a list of all the forcetable filenames so that
# the forcetables can be copied to the working directory.

lappend tablenames 9_095_11.tab 
lappend tablenames n9_c160_22.tab




# Analysis Parameters
#----------------------------------------------------------# 

# These are are parameters that will be passed to the setup_analysis
# command when it is called in the main.tcl script.

set mgrid 8
set stray_cut_off 3

# Use these flags to specify which observables to calculate during the
# simulation.  Values are calculated after every call to the espresso
# integrate command and written to files like
# time_vs_parametername. See the analysis package for more
# details
lappend analysis_flags orient_order
lappend analysis_flags pressure 
lappend analysis_flags pik1  
lappend analysis_flags boxl
lappend analysis_flags flipflop
lappend analysis_flags energy
lappend analysis_flags stray
lappend analysis_flags fluctuations
set profile_beadtypes [list  0 1 2 ]
lappend analysis_flags "density_profile -beadtypes \{ $profile_beadtypes \} -nogrid"

::mmsg::send [namespace current] "done"

















