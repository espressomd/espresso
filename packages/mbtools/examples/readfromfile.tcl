# Parameter file for simulating a colloidal particle interacting with
# a membrane
#
#
# To run this file use $ESPRESSO_SOURCE/$PLATFORM/Espresso ./scripts/main.tcl readfromfile.tcl
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
set warmsteps 1000
set warmtimes 20
set free_warmsteps 1000
set free_warmtimes 10 
set area_lipid 1.29
set startj 0    
set startk 0
set startmdtime 0
set npt off


# Specify the name of the job <ident> the location of forcetables and
# the output and scripts directories
set ident "readfromfile"
set tabledir "./forcetables/"
set outputdir "./$ident/"
set topofile "$ident.top"

# Espresso Special Parameters #
setmd max_num_cells 2744

# Specify how we want to use vmd for visualization allowable values
# are "interactive" "offline" "none".  interactive rarely works
set use_vmd "offline"

# Set the box size to the size in the file
set setbox_l  { 102.852650334 102.852650334 15.0 }

# --- Specify a global key for molecule types -----#

set moltypes [subst { { 1 lipid { 0 1 2 } { 0 1 } } } ]

# --- Specify the system geometry and composition ----#
# Set the geometry to flat
set geometry { geometry { readfile ./readfile_cfg.gz ./readfile_cfg.top } }

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
set hb_parms  [list  Harm_bend 10  4.0 ]
set fene_parms { FENE 30 1.5 }

lappend bonded_parms $fene_parms
lappend bonded_parms $hb_parms


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
lappend analysis_flags { pressure_calc  0 1 }
lappend analysis_flags { pik1_calc  0 1 }
lappend analysis_flags { box_len_calc  0 1}
#lappend analysis_flags { flipflop_calc  0 1}
lappend analysis_flags { energy_calc  0 1}
#lappend analysis_flags { orient_order_calc  0 1}

# It is not recommended to include fluctuation calculations during the
# simulation since they will crash if a hole appears in the membrane

#lappend analysis_flags { fluctuation_calc 1 1}

::mmsg::send [namespace current] "done"

















