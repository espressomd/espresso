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
set L0 1 
# energy scale
set kbT0 1 
# densityscale
set rho0 1
# viscosity
set nu0 1
# time scale
set t0 1
# hydrodyn. friction
set zeta0 1
# shear konstant
set ks0 1
# bending constant
set kb0 $kbT0

######################################################################
# define Simulation Parameters

# time step MD
set dt 0.001
# time step LB-Fluid 
set tau 0.16666666667
# Integrations per frame
set stepSize 100
# number of steps ( frames)
set numSteps 1000
# gridsize setting
set gridPara 1
# maximum strech length
set maxStretch 5

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
set zeta 1
# thermal energy (equals 300 K here)
#set kbT [expr 300.0 *1.38e-23 / $kbT0] 
set kbT 0

# bending 
set kB 1
# shearing
set kS 1


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

#turning off warnings
setmd warnings 0

part 0 pos 0 0 0 virtual 1
part 1 pos 2 0 0 virtual 1
part 2 pos 0 3 0 virtual 1

inter 0 triel 0 1 2 $maxStretch $kS 0
part 0 bond 0 1 2

#part 1 (x) - part 0 (x)
set equDistX 2 
#part 2 (y) - part 0 (y)
set equDistY 3

#equDistY 
set equAngle 90

#prepare_vmd_connection test 1 0
#prepare_vmd_connection test 0 3000

# The additional command imd steers the socket connection   #
# to VMD, e.g. sending the actual coordinates               #

set partCsv0 [ open "simfiles/partInfo0.csv" "w"]
puts $partCsv0 "posx,posy,posz,vx,vy,vz,fx,fy,fz,lbvx,lbvy,lbvz,t"
set partCsv1 [ open "simfiles/partInfo1.csv" "w"]   
puts $partCsv1 "posx,posy,posz,vx,vy,vz,fx,fy,fz,lbvx,lbvy,lbvz,t"
set partCsv2 [ open "simfiles/partInfo2.csv" "w"]
puts $partCsv2 "posx,posy,posz,vx,vy,vz,fx,fy,fz,lbvx,lbvy,lbvz,t"


####### SIMULATION ######

# update position of all particles
 integrate 0
# Checkpoint only without using gpu
# writeCheckpoint start

set startT [clock seconds]


puts "Starting simulation"

for {set step 0} {$step < $numSteps} {incr step} {
  
  integrate $stepSize
  for {set j 0} { $j < $numParts } {incr j} {
      set partCsv [ open "simfiles/partInfo$j.csv" "a"]
      set dataPart [part $j print pos v f]
      set dataFlow [lbfluid print_interpolated_velocity [lindex $dataPart 1] [lindex $dataPart 2] [lindex $dataPart 3]]
      set data [concat $dataPart $dataFlow]
      #puts $data
      #lappend data $step
      foreach x $data { 
	  puts -nonewline $partCsv "$x,"
      }
      puts -nonewline $partCsv "$step\n"
      close $partCsv
  }
  if {fmod($step, 100)==0} { 
      puts "Done $step out of $numSteps"
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
  }

  if {$step > 3500} { 
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
      
  }
  
  
  
  if 0 { 
  for {set j 0} { $j < $numParts } {incr j} {
     
      puts  $partFile "Force acting on part: [part $j print id f]"
      set partPos [part $j print pos]
      puts  $partFile "Fluid Velocity at Particle: [lbfluid print_interpolated_velocity [lindex $partPos 0] [lindex $partPos 1] [lindex $partPos 2]]"
  }
  }

  
 #  imd positions
}

# store simulation time in log.dat
set endT [clock seconds]
set logfile [open "log.dat" "w"]
puts $logfile "Needed [expr $endT-$startT] seconds ([expr ($endT-$startT)/60] minutes)."
