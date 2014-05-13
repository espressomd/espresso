# description: simulates a lbtracers in an LB fluid with Temperature
# author: Mohamed
# date: 7.05.14


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
set numSteps 50
# gridsize setting
set gridPara 1
# maximum strech length
set maxStretch 3

# No. of particles
set numParts 3

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
#set kbT [expr 300.0 *1.38e-23 / $kbT0] 
set kbT 0

# bending 
set kB [expr 1e-19 / $kb0]
# shearing
set kS [expr 5e-6/ $ks0]


#################################################################################
######## SETUP #########################################################

# make directory for simulation files
file mkdir "simfiles"

# Write an info File, contains relevant information about the Simulation
set infoFile [open "simfiles/simInfo.dat" "w"] 		
puts $infoFile "# Scalar"
puts $infoFile "# timeStep / stepSize / numSteps / timeScale / lengthScale"
puts $infoFile "$dt $stepSize $numSteps $t0 $L0"
close $infoFile

# info File that contains relevant Data about the system
set paraFile [open "simfiles/paraInfo.dat" "w"]
puts $paraFile "# Scalar"
puts $paraFile "# nu / rho / radius / zeta / kT / kS / kB"
puts $paraFile "$nu $rho $r $zeta $kbT $kS $kB"
close $paraFile

# write file with particle sequence data
set partFile [ open "partInfo.dat" "w" ]
puts $partFile "#PartInfo"
puts $partFile "#NumPart = $numParts"

# setting Boxlength
setmd box_l $boxx $boxy $boxz

# setting integration parameters
# skin for verlet list
setmd skin 0.1 
# timestep
setmd time_step $dt

# setting up the fluid with or without using gpu
lbfluid agrid $gridPara dens $rho visc $nu tau $tau friction $zeta ext_force 5 0 0


#setting themostat
thermostat lb $kbT

part 0 pos 10 10 10 virtual 1
part 1 pos 12 10 10 virtual 1
part 2 pos 10 13 10 virtual 1

inter 0 triel 0 1 2 $maxStretch $kS 0
part 0 bond 0 1 2

#part 1 (x) - part 0 (x)
set equDistX 2 
#part 2 (y) - part 0 (y)
set equDistY 3

#equDistY 
set equAngle 90

#prepare_vmd_connection test 1 0
prepare_vmd_connection test 0 3000

# The additional command imd steers the socket connection   #
# to VMD, e.g. sending the actual coordinates               #


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
  
  puts $partFile "----- Integration Step : $step -----" 
  for {set j 0} { $j < $numParts } {incr j} {
     
      puts  $partFile "Force acting on part: [part $j print id f]"
      set partPos [part $j print pos]
      puts  $partFile "Fluid Velocity at Particle: [lbfluid print_interpolated_velocity [lindex $partPos 0] [lindex $partPos 1] [lindex $partPos 2]]"
  }

  set partPos0 [part 0 print pos]
  set partPos1 [part 1 print pos]
  set partPos2 [part 2 print pos]
  
  puts $partFile "DistX: [lindex $partPos1 0] - [lindex $partPos0 0] = [expr [lindex $partPos1 0] - [lindex $partPos0 0]]"
  puts $partFile "DistY: [lindex $partPos2 1] - [lindex $partPos0 1] = [expr [lindex $partPos2 1] - [lindex $partPos0 1]] \n"
  
   imd positions
}

# store simulation time in log.dat
set endT [clock seconds]
set logfile [open "log.dat" "w"]
puts $logfile "Needed [expr $endT-$startT] seconds ([expr ($endT-$startT)/60] minutes)."
