#
# Parameter file for simulating a one component flat lipid bilayer
# at fixed box size and temperature.
#
#
# To run this file use $ESPRESSO_SOURCE/$PLATFORM/Espresso ./scripts/main.tcl simplebilayer.tcl
#

# For more information about the simulation model on which this simulation is based see the following reference
# 
#  Cooke, I. R., Kremer, K. and Deserno, M.(2005): Efficient,
# tuneable, generic model for fluid bilayer
# membranes. Phys. Rev. E. (to appear) also see cond-mat/0502418
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
set ident "simplebilayer"
set tabledir "./forcetables/"
set outputdir "./$ident/"
set topofile "$ident.top"

# Espresso Special Parameters #
setmd max_num_cells 2744

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
set geometry { geometry flat }

# In this line we specify that 1 molecule of type 0 (ie colloid) and
# 20 of type 1 (ie lipid) are to be used
set n_molslist { n_molslist {  { 0 340 } } }

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
set int_steps   1000
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
set hb_parms  [list  Harm_bend 10  4.0 ]
set fene_parms { FENE 30 1.5 }

lappend bonded_parms $fene_parms
lappend bonded_parms $hb_parms


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

#Set the size of the 2d grid for calculating the membrane height
#function.  Used to calculate stray lipids lipid flip-flip rates and
#for fluctuation analyses
set mgrid 8

# Distance from the bilayer beyond which a lipid is considered to be
# stray
set stray_cut_off 10000


set profile_beadtypes [list  0 1 2 ]

set profilenogrid 1

# Use these flags to specify which observables to calculate during the
# simulation.  Values are calculated after every call to the espresso
# integrate command and written to files like
# time_vs_parametername. See the module ::std_analysis for more
# details
lappend analysis_flags { pressure_calc  0 1 }
lappend analysis_flags { pik1_calc  0 1 }
lappend analysis_flags { box_len_calc  0 1}
lappend analysis_flags { flipflop_calc  0 1}
lappend analysis_flags { energy_calc  0 1}
lappend analysis_flags { orient_order_calc  0 1}
lappend analysis_flags { density_profile_calc  0 1}

# It is not recommended to include fluctuation calculations during the
# simulation since they will crash if a hole appears in the membrane

#lappend analysis_flags { fluctuation_calc 1 1}

::mmsg::send [namespace current] "done"

















