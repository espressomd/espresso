# description: simulates a spherical cell in an LB fluid with Temperature
# author: Vera
# date: 18.01.14


#######################################################################
# Parameters setting on for every simulation

# Geometry of the box
set boxx  20
set boxy  20
set boxz  20

#######################################################################
# calculations to get quantities in espresso units
# Fix dimensions 

# length scale
set L0 1e-7 
# energy scale
set kbT0 3.77e-21 
# densityscale
set rho0 1e3
# viscosity
set nu0 [expr sqrt( $kbT0/($L0*$rho0) )]
# time scale
set t0 [expr pow($L0,5/2.) * sqrt($rho0) / sqrt( $kbT0 ) ]
# hydrodyn. friction
set zeta0 [expr sqrt( $L0*$rho0*$kbT0) ]
# shear konstant
set ks0 [expr $kbT0 / ($L0*$L0) ]
# bending constant
set kb0 $kbT0

######################################################################
# define Simulation Parameters

# time step MD
set dt 0.001
# time step LB-Fluid 
set tau 0.001 
# Integrations per frame
set stepSize 100
# number of steps ( frames)
set numSteps 100
# gridsize setting
set gridPara 1
# maximum strech length
set maxStretch 3

# number of steps between vtkfiles
set distVtk 1

# calculate other values from this setup

# kinematic viscosity 
set nu [expr 0.5* $gridPara * $gridPara / (3* $tau )]

# density
set rho 1

# hydrodynamic particle radius (used to define zeta)
set r [expr $gridPara * 1]
# friction constant
set zeta [expr 6 * 3.141* $nu* $rho * $r]
# thermal energy (equals 300 K here)
set kbT [expr 300.0 *1.38e-23 / $kbT0] 

# bending 
set kB [expr 1e-19 / $kb0]
# shearing
set kS [expr 5e-6/ $ks0]


#################################################################################
######## SETUP #########################################################

# Write an info File, contains relevant information about the Simulation
set infoFile [open "simInfo.dat" "w"] 		
puts $infoFile "# Scalar"
puts $infoFile "# timeStep / stepSize / numSteps / timeScale / lengthScale"
puts $infoFile "$dt $stepSize $numSteps $t0 $L0"
close $infoFile

# info File that contains relevant Data about the system
set paraFile [open "paraInfo.dat" "w"]
puts $paraFile "# Scalar"
puts $paraFile "# nu / rho / radius / zeta / kT / kS / kB"
puts $paraFile "$nu $rho $r $zeta $kbT $kS $kB"
close $paraFile

exec rm -r /home/btpj/s1mokhal/sphere/vtkfiles
exec mkdir /home/btpj/s1mokhal/sphere/vtkfiles

#tools needed for system setup (Juropa or local)
set tclToolsDir "~/sphere/tclTools"
source "$tclToolsDir/addCell_tools.tcl"
source "$tclToolsDir/writeVtk_folded.tcl"
# source "$tclToolsDir/addCell_fixed.tcl"


# setting Boxlength
setmd box_l $boxx $boxy $boxz

# setting integration parameters
# skin for verlet list
setmd skin 0.1 
# timestep
setmd time_step $dt
# coose how good no slip condition is: 0 non, 1 better
# setmd sequ 1

# setting up the fluid with or without using gpu
lbfluid agrid $gridPara dens $rho visc $nu tau $tau friction $zeta
#setting themostat
thermostat lb $kbT


################## ADD CELLS #############
set interCount 0
set cellCount 0
set nodeCount 0
set numNodesPerCell 0

set numType 0

# number of cells
set numCells 1
# setmd vescnum $numCells

# write file with cell sequence data
set cellSeqFile [ open "cellSeq.dat" "w" ]
puts $cellSeqFile "#CellSequence"
puts $cellSeqFile "#NumCells = $numCells"

addCell [expr $boxx/2.0] [expr $boxy/2.0] [expr $boxz/2.0] $numType
incr numType



####### SIMULATION ######

# update position of all particles
 integrate 0
# Checkpoint only without using gpu
# writeCheckpoint start

set startT [clock seconds]


puts "Starting simulation"

for {set step 0} {$step < $numSteps} {incr step} {
  
  integrate $stepSize
  if {fmod($step, 10)==0} { puts "Done $step out of $numSteps steps." }
   
#   write cell sequence to cellSeq.dat
  puts $cellSeqFile "#Frame = [expr $stepSize*$step]"
  puts $cellSeqFile "#Time = [expr $stepSize*$step*$dt]"
  
  for {set j 0} { $j < $numCells } {incr j} {
   for {set k 0} {$k < $numNodesPerCell} {incr k} {
     puts $cellSeqFile [part [expr $j*$numNodesPerCell + $k ] print pos]
    }
  
  }
  
# output for paraview only every 10th step
  if {fmod($step, $distVtk)==0} { 
      for { set i 0} { $i < $numCells} {incr i}  {writevtkCell "~/sphere/vtkfiles/cell-$i $step.vtk" $i}
     }

}

# store simulation time in log.dat
set endT [clock seconds]
set logfile [open "log.dat" "w"]
puts $logfile "Needed [expr $endT-$startT] seconds ([expr ($endT-$startT)/60] minutes)."

