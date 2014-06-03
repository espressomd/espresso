#######################################################################
# Simulation Info #

# Date: 05.05.14
# author: Vera
# description: Drag particle parallel to Cell while holding cell with Ring of fixed nanoparticles

#######################################################################
# Parameters setting on for every simulation 
# -> make changes on parameters here
########################################################################

# setting precision of tcl means decimal?
# set tcl_precision 3

# Box geometry
set boxx  48
set boxy  48
set boxz  48

# Fix physical quantities (SI)

# kinematic viscosity 
# -> if visc is given here time step has to be set later
set nu_SI [expr 1e-6] 
# density of fluid
set rho_SI 1e3
# shear modulus
set ks_SI 5e-6
# bending
set kb_SI 1e-19
# grid size
set gridsize_SI 1e-7
# thermal energy
set kbT_SI [expr 1.38e-23 * 0]

# Particle parameters:

# radius:
# hydrodynamic radius, used to calculate friction constant. Diameter should be around grid constant
set rad_SI [expr 0.5*$gridsize_SI] 
# friction constant: is only set here for consitancy
set zeta_SI [expr 6*3.142*$nu_SI*$rho_SI*$rad_SI]
# # drag force on particle 
set part_force_SI [expr 0.01*$zeta_SI]
# particle volume 
set vol_part_SI [expr 4/3*3.142*$rad_SI*$rad_SI*$rad_SI]
# mass 
# set mass_SI [expr $rho_SI*$vol_part_SI]

########### Four possible unit systems ###################################

# # Unit system I

# # a priori fixed scales
#set rho0 $rho_SI
#set nu0 $nu_SI
#set L0 $deltaX_SI

# # derived scales
#set E0 [expr $rho0*$nu0*$nu0*$L0]
#set ks0 [expr $E0/($L0*$L0)]
#set V0 [expr $nu0/$L0]
#set kb0 [expr $E0]

# # Unit system II

# # a priori fixed scales
#set rho0 $rho_SI
#set V0 $V_SI
#set L0 $deltaX_SI

# # derived scales
#set nu0 [expr $L0*$V0]
#set E0 [expr $rho0*$nu0*$nu0*$L0]
#set ks0 [expr $E0/($L0*$L0)]
#set kb0 [expr $E0]

# Unit system  III

# fixed scales
set rho0 1e3
set L0 1e-7
# be careful energy scale is not set to 0°c!
set E0 [expr 1.38e-23 * 273]

# derived scales
set nu0 [expr sqrt($E0/($rho0*$L0))]
# set V0 $nu0/$L0
set ks0 [expr $E0/($L0*$L0) ]
set kb0 $E0
set m0 [expr $rho0*$L0*$L0*$L0]


# # Unit system  IV
# 
# # fixed scales  
# set rho0 $rho_SI
# set L0 $deltaX_SI
# set ks0 $ks_SI
# 
# # derived scales
# set E0 [expr $ks0*$L0*$L0]
# set nu0 [expr sqrt($E0/($rho0*$L0))]
# set V0 $nu0/$L0 
# set kb0 $E0

############### Calculate quantities in Espresso units ###############################

set rho [expr $rho_SI / $rho0]
set nu [expr $nu_SI / $nu0]
set kbT [expr $kbT_SI / $E0]
set kS [expr $ks_SI / $ks0]
set kB [expr $kb_SI / $kb0]
set kC $kS
# set V [expr $V_SI /  $V0]
set gridsize [expr $gridsize_SI / $L0]

# particle scales
#radius
set rad [expr $rad_SI / $L0]
# volume
set vol_part_SI [expr 4/3*3.142*$rad*$rad*$rad]
# mass 
set mass 1
# friction constant: is calculated new from viscosity and paricle radius
set zeta [expr 6*3.142*$nu*$rho*$rad]
# drag force
set part_force [expr 0.01*$zeta]

# Cell parameters
set cellRad 18
set cellhight [expr $cellRad/6.0]
# distance from Cell surface to particle
set distance 2.0

########################################################################
# Numerical parameters

# maximum strech length for bonds
set maxStretch 3
# timestep calculated from condition that relaxation parameter tau in LB has to be 1 (or similar)
set taurelax 1
set dt_lb [expr ($taurelax - 0.5)*$gridsize*$gridsize/(3*$nu)]
set dt_md [expr $dt_lb]

# Integrations per step
set stepSize 1
# number of steps
set numSteps 10000
# set after which number of steps vtkfile should be written
set distVtk [expr $numSteps/200]
set distCellSeq [expr 1]
set distPartSeq [expr 1]
# set distFlow [expr $numSteps/100]

# Dimensionless numbers following Le
# both in espresso and SI units for consistency check
# particle radius needs to be compatible with node file
# lengthscale of simulation
set L [expr $boxz-1] 
set L_SI [expr $L * $L0]
# dynamic viscosity
set mu [expr $nu * $rho]
set mu_SI [expr $nu_SI * $rho_SI]
# thermal mean squared velocity
set v_mean [expr sqrt(3 * $kbT / $mass)]
set mass_SI [expr $mass*$m0]
set v_mean_SI [expr sqrt(3 * $kbT_SI / $mass_SI)]

# Renolds number calculated with meansquared velocity
set Re [expr $rad*$v_mean / $nu]
set Re_SI [expr $rad_SI*$v_mean_SI / $nu_SI]

# "stiffness ratio" - probably should be about one??
# set Sr [expr $kC/$kS]
# set Sr_SI [expr $kc_SI / $ks_SI]

# estimate final velocity, number of timesteps to get near Cell border
set vel_final [expr $part_force / $zeta]
set N_drag [expr $cellRad/($vel_final * $dt_md) ]

# output
puts "nu = $nu"
puts "rho = $rho"
puts "Re = $Re     Re_SI = $Re_SI"
puts "dt_lb = $dt_lb"
puts "dt_md = $dt_md"
puts "kS = $kS"
puts "kB = $kB"
puts "kC = $kC"
# puts "kbT = $kbT"
puts "zeta = $zeta"
puts "final vel = $vel_final"
puts "N_drag = $N_drag "


######## big info file #########
set simInfo [open "paraInfo_step1e4.dat" w]
puts $simInfo "nu = $nu"
puts $simInfo "rho = $rho"
# puts $simInfo "V = $V"
# puts $simInfo "Ca = $Ca     Ca_SI = $Ca_SI"
puts $simInfo "Re = $Re     Re_SI = $Re_SI"
# puts $simInfo "Sf = $Sf     Sf_SI = $Sf_SI"
# puts $simInfo "Sr = $Sr     Sr_SI = $Sr_SI"
puts $simInfo "dt_md = $dt_md"
puts $simInfo "dt_lb = $dt_lb"
puts $simInfo "part_force = $part_force"
puts $simInfo "kS = $kS"
puts $simInfo "kB = $kB"
puts $simInfo "zeta = $zeta"
puts $simInfo "particle radius = $rad"
puts $simInfo "final vel = $vel_final"
puts $simInfo "N_drag = $N_drag "

close $simInfo

####### parameter file #########
# Write an info File, contains relevant information about the Simulation
set infoFile [open "simInfo_step1e4.dat" "w"] 		
puts $infoFile "# Scalar"
puts $infoFile "# timeStep / stepSize / numSteps / lengthScale"
puts $infoFile "$dt_md $stepSize $numSteps $L0"
close $infoFile

#########################################################################
# SETUP 
#########################################################################

# make directory for files
file mkdir "vtkfiles"
# file mkdir "partSeq"
# file mkdir "cellSeq"

#tools needed for system setup (Juropa or local)
set tclToolsDir "~/Desktop/github/espresso/triel_sim/tclTools"
source "$tclToolsDir/addCell_tools.tcl"
source "$tclToolsDir/writeVtk_folded.tcl"
source "$tclToolsDir/addCell_fixed.tcl"

# seed randumnumber generator:
# t_random seed [clock seconds] [clock seconds] [clock seconds] [clock seconds]

# setting Boxlength
setmd box_l $boxx $boxy $boxz

# setting integration parameters
# skin for verlet list
setmd skin 0.1 
# timestep
setmd time_step $dt_md
# coose how good no slip condition is: 0 non, 1 better
# setmd sequ 1

# setting up the fluid with or without using gpu
lbfluid gpu agrid $gridsize dens $rho visc $nu tau $dt_lb friction $zeta
#setting themostat
thermostat lb $kbT

############## ADD WALL ##################

# walls located at z=1 and z=boxz-1 in x-y-plane
lbboundary wall dist 1 normal 0. 0. 1.
lbboundary wall dist [expr -$boxz+1]  normal 0. 0. -1.

################## ADD CELLS #############

# global variables that will be changed by function to add cells
# counter for all interactions
set interCount 0
# counter for all cells
set cellCount 0
# counter for particles
set nodeCount 0
set numNodesPerCell 0
# counter for printing vtk files
set numType 0

# number of cells
set numCells 1
setmd vescnum $numCells

# write file with cell sequence data
set cellSeqFile [ open "cellSeq_step1e4.dat" "w" ]
puts $cellSeqFile "#CellSequence"
puts $cellSeqFile "#NumCells = $numCells"

# add Cell with addCell_tools script
addCell [expr $boxx*0.5] [expr $boxy*0.5] [expr $boxz*0.5] $numType
incr numType

# add ring around Cell
addFixRing_more [expr $boxx*0.5] [expr $boxy*0.5] [expr $boxz*0.5] $cellRad $cellhight $numType
incr numType

############## ADD NANOPARTICLES ##########
# Anzahl der Partikel für spätere Ausgabe
set numPart 1
#set particle
part [expr $nodeCount] pos [expr ($boxx*0.5)-7 ] [expr $boxy*0.5] [expr ($boxz*0.5)+$cellhight+$distance] type $numType ext_force $part_force 0 0 
# incr nodeCount

# Write file with Particle data
set partSeqFile [open "partSeq_step1e4.dat" "w"]
puts $partSeqFile "#VectorSequence"
puts $partSeqFile "#NumValues = [expr 2*$numPart]"
puts $partSeqFile "#Particle trajetories"

############################################################################
####### SIMULATION ######
############################################################################



# update position of all particles
 integrate 0

puts "Starting simulation"

set startT [clock seconds]


for {set step 0} {$step < $numSteps} {incr step} {
  
  # integration command
  integrate $stepSize
  
  # console output 
  if {fmod($step, 100)==0} { puts "Done $step out of $numSteps steps." }
  
  ############ cell output ###########################
  
  # write cell sequence
  if {fmod($step, $distCellSeq)==0} {
    
    puts $cellSeqFile "#Frame = [expr $stepSize*$step]"
    puts $cellSeqFile "#Time = [expr $stepSize*$step*$dt_md]"
    
    for {set j 0} { $j < $numCells } {incr j} {  
      for {set k 0} {$k < $numNodesPerCell} {incr k} {
	puts $cellSeqFile [part [expr $j*$numNodesPerCell + $k ] print pos]
      }
    }
  
  }
  
  # Ausgabe für paraview nur alle x Schritte
  if {fmod($step, $distVtk)==0} { 
      for { set i 0} { $i < $numCells} {incr i} {writevtkCell "vtkfiles/cell-$i $step.vtk" $i}
  }
  
  ################## particle output ######################
  
  #   write particle sequence  
  if {fmod($step, $distPartSeq)==0} {
    
    puts -nonewline $partSeqFile "[expr $stepSize*($step+1)] [expr $stepSize*($step+1)*$dt_md]"     
    
    for {set j 0} { $j < $numPart} {incr j} {
      set partID [expr $nodeCount + $j]
      
      puts -nonewline $partSeqFile "  [part $partID print pos]"
      puts -nonewline $partSeqFile "  [part $partID print v]"

    }
    puts $partSeqFile ""
  }
  
  # output for paraview
  if {fmod($step, $distVtk)==0} { 
      writevtk "vtkfiles/part_$step.vtk" $numType
  }
  ######################## print ring####################################
  
  if {fmod($step, $distVtk)==0} { 
      writevtk "vtkfiles/ring_$step.vtk" 1
  }
  
  ##################### flow and other output ##########################
  
#   if {fmod($step, $distVtk)==0} { 
#       lbfluid print vtk velocity "vtkfiles/vis_$step.vtk"
#   }
    
#   if {fmod($step, $checkDist) ==0} { set name "check_$step"
# 				     writeCheckpoint $name
# 				   }
  
}

set endT [clock seconds]
set logfile [open "log_step1e4.dat" "w"]
puts $logfile "Needed [expr $endT-$startT] seconds ([expr ($endT-$startT)/60] minutes)."

# puts " seed was [t_random seed]" 
close $logfile
# close $cellSeqFile
close $partSeqFile   
