#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    exec $ESPRESSO_SOURCE/Espresso $0 $*
#
#  This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
# 
#############################################################
#                                                           #
#  ESPResSo Tutorial                                        #
#                                                           #
#  Created:       25.07.2004 by HL                          #
#  Modified:      26.06.2006 by OL                          #
#  Modified:      26.11.2006 by TS                          #
#                                                           #
#############################################################

puts " "
puts "======================================================="
puts "=               ESPResSo Tutorial                     ="
puts "=          Learning ESPResSo in 10 steps              ="
puts "======================================================="
puts " "

#############################################################
# Introduction:                                             #
# The aim of this tutorial is to get known to the most      #
# common additional tcl commands which will enable you      #
# to write your own ESPResSo simulation script.             #
# The Tutorial is organized in sections each of which       #
# introduces the usage of one feature of ESPResSo. A        #
# comment block will tell you what has to be done.          #
# It is best if you open the additional commands page of    #
# the ESPResSo documentation to look up the syntax of       #
# the needed commands.                                      #
# On your way through the tutorial you will also learn      #
# some basics about tcl.                                    #
# After succesfully finishing a section remove the 'exit'   #
# command and move on to the next section                   #
#############################################################


#############################################################
# 1 Program Information                                     #
# One of the auxiliary tcl-commands 'code_info' will        #
# return information of the compilation status.             #
# TCL:                                                      #
# Placing a command in brackets '[command]' will            #
# execute the command and results the return value of the   #
# command.                                                  #
# To write a string to the screen use 'puts "<string>"'     #         
#############################################################

puts "[code_info]"

#exit

#############################################################
# 2 Preparing a System                                      #
# In the tutorial we want to simulate a single polymer      #
# chain in an NVT ensemble. The System will contain a       #
# polymer chain of 40 uncharged monomers                    #
# Create variables for the chain length                     #
# and the density. Calculate the volume of the              #
# system and the corresponding simulation box length        #
# assuming a cubic simulation box.                          #
# TIP:                                                      #
# Store calculated values in tcl variables so that they     #
# can be used later on. Control your results with help of   #
# the 'puts' command                                        #
#                                                           #
# TCL:                                                      #
# Variables in tcl are created with the 'set' command:      #
# set <variable_name> <value>                               #
# You can reference to that variable later on by preceeding #
# its name with a $ sign:                                   #
# $<variable_name>                                          #
# Example:                                                  #
# set a 20                                                  #
# puts $a                                                   #
#                                                           #
# More TCL:                                                 #
# Since tcl is a string oriented language calculations have #
# to be done using the 'expr' command:                      #
# [expr <analytical expression>]                            #
# will return the value of <analytical expression>          #
# Example:                                                  #
# set a 20; set b 5                                         #
# set res [expr $a*$b]                                      #
# puts $res                                                 #
#############################################################

set l_poly  40
set n_chains 1
set n_part  [expr $n_chains * $l_poly]
set density 0.001

set volume     [expr $n_part/$density] 
set box_length [expr pow($volume,1.0/3.0)]

puts "Simulate single chain polymer N=$l_poly at density $density"
puts "Simulation box: $box_length"

#exit

#############################################################
# 3 Initialize ESPResSo                                     #
# In the last section you have calculated the box length,   #
# but so far ESPResSo does not know about it. Use the       #
# additional command 'setmd' to set the box length in       #
# ESPResSo.                                                 #
# Initialize the integrator by setting the time step to 0.01#
# the skin to 0.4.                                          #
# To tell ESPResSo to use an NVT ensemble use the           #
# additional command 'integrate'.                           #
# Use the additional command 'thermostat' to enable a       #
# langevin thermostat at temperature 1.0 and friction       #
# constant gamma 1.0                                        #
# TIP:                                                      #
# You find all Global variables in the documentation:       #
# additional commands -> The setmd command -> Global variables #
#                                                           #
# Actually it is now where you have started interacting     #
# with  ESPResSo. Congratulations!                          #
#############################################################


setmd box_l $box_length $box_length $box_length
setmd time_step 0.01
setmd skin 0.4
integrate set nvt
thermostat langevin 1.0 1.0

puts [setmd box_l]
puts [setmd time_step]
puts [setmd skin]
puts [integrate]
puts [thermostat]

#############################################################
# 4 Interactions                                            #
# We need a bond potential for the polyelectrolyte,         #
# excluded volume interactions and a coulomb interaction.   #
# Setting interactions is done with the additional command  #
# 'inter'                                                   #
#                                                           #
# For the bond we use a FENE potential with k_fene=7.0 and  #
# r_fene=2.0. (use <bond_type_number> 0)                    #
#                                                           #
# The excluded volume interaction will be a purely          #
# repulsive lennard-jones interactions. The parameter set is#
# epsilon sigma cutoff  shift offset                        #
# 1.0     1.0   1.12246 0.25  0.0                           #
# Set up lennard-jones interactions for two particle types: #
# 0 and 1                                                   #
#  
# We have to introduce the electrostatic interactions       #
# later for technical reasons.                              #
#############################################################

inter 0 fene 7.0 2.0
inter 0 0 lennard-jones 1.0 1.0 1.12246 0.25 0.0
inter 0 1 lennard-jones 1.0 1.0 1.12246 0.25 0.0
inter 1 1 lennard-jones 1.0 1.0 1.12246 0.25 0.0

#############################################################
# 5 Create Particles                                        #
# Create the Polyelectrolyte and the counterions by using   #
# the additional command 'polymer' and 'counterions'.       #
# Use particle type 0 for the polyelectrolyte and 1 for     #
# the counterions.                                          #
# TIP:                                                      #
# You can control the result with the additional command    #
# 'part', e.g: puts [part 10] prints information about      #
# the particle with id 10                                   #
#############################################################

polymer 1 $l_poly 1.0 start 0 charge 1.0 types 0 0 FENE 0
counterions $n_ci start [setmd n_part] charge -1.0 type 1 
 
puts [part 10] 
puts [part 50] 

#puts "[inter coulomb 1.0 p3m tune accuracy 0.001 mesh 8]"
#############################################################
# 6 Visualization                                           #
# Now lets see what you have done so far.                   #
# Set the variable vmd to "yes" for trying.                 #
# TCL:                                                      #
# Here you can already see how an 'if' condition works in   #
# tcl.                                                      #
#############################################################

set vmd "no"

if { $vmd == "yes" } {
# This calls a small tcl script which starts the program    #
# VMD and opens a socket connection between ESPResSo and    #
# VMD.                                                      #
    prepare_vmd_connection tutorial 3000

# Just wait a moment until VMD has started.                 #
# The 'exec' command is quite useful since with that you can#
# call any other program from within your simulation script.#
    exec sleep 4

# The additional command imd steers the socket connection   #
# to VMD, e.g. sending the actual coordinates               #
    imd positions
}
#############################################################
# 7 Warmup                                                  #
# Since we randomly placed the particles we have to         #
# push the ones that are too close to each other appart.    #
# Otherwise you can run into troubles with the singularity  #
# of the Lennard-Jones interaction.                         #
#                                                           #
# To avoid this singularity ESPResSo offers the             #
# possibility to cap potentials. Look at the 'ljforcecap'   #
# option of the 'inter' command. We start with a value of   #
# 10 for the cap and increase it step wise.                 #
#                                                           #
# The particles will be pushed away from each other in an   #
# integration loop. Below you see the empty structure       #
# for this loop. The tcl 'while' loop has pretty much the   #
# same structure as for example in C. The same holds for    #
# other loop constructions or conditions like 'for' or 'if' #
# Fill in the missing lines after the comments              #
#                                                           #
# To check the status of the sytem we perform the first     #
# analysis of the system. This is done with the             #
# additional command 'analyze <observable>'. In this case   #
# the observable is the minimal distance between any two    #
# particles which is called 'mindist'                       #
#############################################################

set min 0
set cap 10
while { $min < 0.8 } {
    # set ljforcecap
    inter ljforcecap $cap
    # integrate a number of steps, e.g. 20
    integrate 20
    # check the status of the sytem
    set min [analyze mindist]
    # this is a shortcut for 'set cap [expr $cap+10]'
    incr cap 10
}
puts "Warmup finished. Minimal distance now $min"
# turn off the ljforcecap, which is done by setting the 
# force cap value to zero:
inter ljforcecap 0


#############################################################
# 4-2 Electrostatic interaction                             #
# The strength of the electrostatic interaction is given    #
# by the bjerrum length which we will set to 1.0            #
# For the calculation of the electrostatic interaction we   #
# will take the P3M method with the following parameters:   #
# <r_cut> <mesh> <cao> <alpha>                              #
# 13.9629 8      3     0.127621                             #
# (Remark: These values have been obtained with the         #
# internal tuning routine yielding an accuracy of 0.001)    #
# Use the 'inter' command.                                  #
#############################################################
inter coulomb 1 p3m 13.9629 8 3 0.127621 0.000832468

#puts "[inter coulomb 1.0 p3m tune accuracy 0.001 mesh 8]"

#############################################################
# 8 Integration                                             #
# As we know the integration allready from the warmup we    #
# can just construct a very similar structure.              #
# Write a while loop which runs n_cycle times performing    #
# n_steps integration steps per cycle.                      #
#############################################################

set n_cycle 100
set n_steps 100

set obs [open "rg.dat" "w"]

set i 0 
while { $i<$n_cycle } {
    puts $i
    integrate $n_steps

# REMARK:
# Include the last two sections into the integration
# while loop in order to perform analysis, data storage
# repeatedly

#############################################################
# 9 Analysis                                                #
# As a simple analysis we want to observe the chain         #
# extension, namely the radius of gyration 'rg'.            #
# For an analyzation we have first to tell ESPResSo about   #
# molecules. This can be done for simple cases like here    # 
# with 'analyze set chains' (see the documentation) or      #
# more generale with the 'topology' command                 #
# The result of the 'analyze rg' command is a tcl list of   #
# numbers. The first one is the mean value of rg.           #
#                                                           #
# TCL:                                                      #
# To access an element within a tcl list use the command    #
# 'lindex $<list_name> <index>'                             #
#############################################################

    analyze set chains 0 1 $l_poly
    set rg [lindex [analyze rg] 0]

#############################################################
# 10 Store and view Results                                 #
# To store results one has to write data to files           #
# Create a file handle and store it in the tcl variable     #
# 'obs' (Since we want to open up the file only once, place #
# the corresponding line befor the integration loop.        #
# For each integration cycle write the actual simulation    #
# time (access via 'setmd time') and the observable rg into #
# the observable file.                                      #
# Don't forget to close the file handle after the           #
# integration loop.                                         #
# You can look at the result with your favorite graphic     #
# package. At the end of the tutorial you will find         #
# a line which does it for you using gnuplot.               #
# TCL:                                                      #
# You can create a file handle with the command             #
# 'open "<file_name>" "w"'                                  #
# The 'puts' command takes the file handle as an optional   #
# argument: 'puts $<file_handle> "<string>"                 #
# Use 'close $<file_handle>' when finished writing          #
#############################################################
   
    puts $obs "[setmd time] $rg"

    # if you have turned on the vmd option you can now
    # follow what your simulation is doing
    if { $vmd == "yes" } { imd positions }

    incr i
}
close $obs


# Uncommenting the following two lines will show
# you a plot of the rg values
#plotObs "rg.dat" { 1:2 } labels { "time" "rg" } out "rg"
#exec gv rg.ps

#############################################################
#                                                           #
# CONGRATULATIONS !!!                                       #
# You learned how to use ESPResSo and you                   #
# have written your first simulation script.                #
#                                                           #
#                                                           #
# Now everything is working you can play arround with it:   #
# 1) Change the Bjerrum length to observe the exciting chain#
# extension behavior of polyelectrolytes.                   #
# 2) Even more fun can be obtained if you move to poor      #
# solvent conditions by changing the monomer-monomer        #
# interaction (inter 0 0 len...) from purely repulsive to   #
# attractive. Simply Increase the lennard-jones cutoff to   #
# 2.5 and play with the epsilon value.                      #
#                                                           #
# HAVE FUN - VIEL SPASS -                                   #
#                                                           #
#############################################################

