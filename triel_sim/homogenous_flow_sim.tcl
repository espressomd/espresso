# description: simulates 3 lbtracers bonded by trielastic forces in an LB fluid with homogenous flow
# author: Mohamed
# date: 7.05.14

# Unit system  III

# Box geometry
set boxx  20
set boxy  20
set boxz  20

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
set stepSize 100
# number of steps
set numSteps 10000

# Number of particles
set numParts 3


#################################################################################
######## SETUP #########################################################

exec rm -r simfiles
exec rm -r vtkfiles

exec mkdir simfiles
exec mkdir vtkfiles

file mkdir "vtkfiles"


# setting Boxlength
setmd box_l $boxx $boxy $boxz

# setting integration parameters
# skin for verlet list
setmd skin 0.1 
# timestep
setmd time_step $dt_md

# setting up the fluid with or without using gpu
lbfluid agrid $gridsize dens $rho visc $nu tau $dt_lb friction $zeta ext_force 20 0 0

# walls located at z=1 and z=boxz-1 in x-y-plane
# lbboundary wall dist 1 normal 0. 0. 20.
# lbboundary wall dist [expr -$boxz+1]  normal 0. 0. -20.

#setting themostat
thermostat lb $kbT

#turning off warnings
setmd warnings 0

part 0 pos 2 2 0 virtual 1
part 1 pos 4 2 0 virtual 1
part 2 pos 2 5 0 virtual 1

inter 0 triel 0 1 2 $maxStretch $kS 0
part 0 bond 0 1 2

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
puts "[expr ($dt_lb/$gridsize)]"
set part2pos [part 2 print pos]
# part 2 pos [expr [lindex $part2pos 0] + 0] [expr [lindex $part2pos 1]] [expr [lindex $part2pos 2]]

for {set step 0} {$step < $numSteps} {incr step} {
  
  integrate $stepSize
  puts "Done $step out of $numSteps"
  for {set j 0} { $j < $numParts } {incr j} {
      set partCsv [ open "simfiles/partInfo$j.csv" "a"]
      set dataPart [part $j print pos v f]
      set dataFlow [lbfluid print_interpolated_velocity [lindex $dataPart 0] [lindex $dataPart 1] [lindex $dataPart 2]]
      # for {set i 0} { $i < 2 } {incr i} { 
      # 	  lset $dataFlow $i [expr [lindex $dataFlow $i] * ($dt_lb/$gridsize)]
      # }
      set data [concat $dataPart $dataFlow]
      foreach x $data { 
	  puts -nonewline $partCsv "$x,"
      }
      puts -nonewline $partCsv "$step\n"
      close $partCsv

      if {$step == 200} { 
	  lbfluid ext_force 0 0 0

      }
  }
  
  
# output for paraview only every 10th step
  # if {fmod($step, 10)==0} { 
  #     writevtk "vtkfiles/tri $step.vtk"
  # }
  
 
 #imd positions
}

# store simulation time in log.dat
set endT [clock seconds]
set logfile [open "simfiles/log.dat" "w"]
puts $logfile "Needed [expr $endT-$startT] seconds ([expr ($endT-$startT)/60] minutes)."
