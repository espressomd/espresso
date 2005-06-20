# 
# Parameter file for simulating a the condensation of lipids onto an
# attractive spherical constraint.
#
#
# To run this file use $ESPRESSO_SOURCE/$PLATFORM/Espresso ./scripts/main.tcl condense.tcl
#
####################
#DOESN'T WORK YET!!!!
###################

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
set free_warmsteps 1000
set free_warmtimes 10
set area_lipid 1.29
set startj 0    
set startk 0
set startmdtime 0
set npt off
set p_ext 0.0
set piston_mass 0.00005
set gamma_0 1.0
set gamma_v 0.0001

# Specify the name of the job <ident> the location of forcetables and
# the output and scripts directories
set ident "condense"
set tabledir "./forcetables/"
set outputdir "./$ident/"
set topofile "$ident.top"

# Espresso Special Parameters #
setmd max_num_cells 2744

# Specify how we want to use vmd for visualization allowable values
# are "interactive" "offline" "none".  interactive rarely works
set use_vmd "offline"


# --- Specify a global key for molecule types -----#

# Here we specify standard lipids and a spherical constraint
set moltypes [list { 0 lipid { 0 1 2 } { 0 1 } } { 1 sphericalconstraint { 3 } { 5.0 1 } } ]


# Specify the box dimensions
set setbox_l  { 25.0 25.0 25.0 }
lappend boxcenter [expr [lindex $setbox_l 0]*0.5 ] 
lappend boxcenter [expr [lindex $setbox_l 1]*0.5 ] 
lappend boxcenter [expr [lindex $setbox_l 2]*0.5 ] 


# Construct a spec for a spherical constraint
set geometry [list geometry { singlemol { 0 0 0 } } ]
set n_molslist [list n_molslist { { 1 1 } } ]
# Note that because this is a constraint we do not specify any
# molecules list thereby preventing the constraint from being added to
# the topology


lappend colloidspec $geometry
lappend colloidspec $n_molslist

# Construct a parameter spec for a random fluid of lipids all of which
# must lie outside a sphere of radius 8 sigma with center at the
# boxcenter.

set geometry [subst { geometry { random sphere { $boxcenter } 8 } }]
set n_molslist [list n_molslist { { 0 500 } } ] 

lappend lipidsspec $geometry
lappend lipidsspec $n_molslist

# Group the two specs together
lappend system_specs $lipidsspec
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
set int_n_times 100
# Write a configuration file every <write_frequency> calls to
# integrate
set write_frequency 2
# Write results to analysis files every <analysis_write_frequency>
# calls to integrate
set analysis_write_frequency 1

# Bonded and bending Potentials
#----------------------------------------------------------#
set hb_parms  [list  Harm_bend 10  4.0 ]
set fene_parms { FENE 30 1.5 }

lappend bonded_parms $fene_parms
lappend bonded_parms $hb_parms


# Tabulated Non Bonded Potentials 
#----------------------------------------------------------#

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


# Other non-bonded potentials
#
# In this case we have a lj interaction between our constraint and the
# lipids specify the details of all other non-bonded potentials into a
# list
#

set lj_eps 4.0
set lj_sigma 1.0
set lj_cutoff 2.5 
set ljshift 0.0
set ljoffset 0.0

# Just for fun we make the head groups attract the membrane while the
# tails repel
set constr_head [list 3 0 lennard-jones $lj_eps $lj_sigma $lj_cutoff $ljshift $ljoffset ]
set constr_tail1 [list 3 1 lennard-jones $lj_eps $lj_sigma [expr 1.125*$lj_sigma] [ expr 0.25*$lj_eps ] $ljoffset ]
set constr_tail2 [list 3 2 lennard-jones $lj_eps $lj_sigma [expr 1.125*$lj_sigma] [expr 0.25*$lj_eps] $ljoffset ]

lappend nb_interactions $constr_head 
lappend nb_interactions $constr_tail1 
lappend nb_interactions $constr_tail2




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

















