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
set pi 3.1415926535897931
set G 0.2
set r 4.9856099295288183473
set H 2
set rho [expr $rho_SI / $rho0]
set nu [expr $nu_SI / $nu0]
set kbT [expr $kbT_SI / $E0]
set kS [expr $ks_SI / $ks0]
set kB [expr $kb_SI / $kb0]
set kC $kS

set shear_rate [expr ($G * $kS) / ($nu * $r)]
set u [expr $shear_rate * $boxz]
# set u 1

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

exec rm -r /home/btpj/s1mokhal/sphere/simfiles
exec mkdir /home/btpj/s1mokhal/sphere/simfiles

# maximum strech length for bonds
set maxStretch 3
# timestep calculated from condition that relaxation parameter tau in LB has to be 1 (or similar)
set taurelax 1
set dt_lb [expr ($taurelax - 0.5)*$gridsize*$gridsize/(3*$nu)]
set dt_md [expr $dt_lb]

# Integrations per step
set stepSize 1
# number of steps
set numSteps 10
# set after which number of steps vtkfile should be written
set distVtk 1
set distDone 1
set distPartSeq 1
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
puts "Radius = $r"
puts "Shear Rate =  [expr ($G * $kS)] / [expr ($nu * $r)] * $H =  [expr ( ($G * $kS) / ($nu * $r) )]"
puts "u = [expr $shear_rate * ($boxz - ($H/2))]"

######## big info file #########
exec touch /home/btpj/s1mokhal/sphere/simfiles/paraInfo_step1e4.dat
set simInfo [open "~/sphere/simfiles/paraInfo_step1e4.dat" w]
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
puts $simInfo "gamma_dot = $shear_rate/$H"

close $simInfo

####### parameter file #########
# Write an info File, contains relevant information about the Simulation
exec touch /home/btpj/s1mokhal/sphere/simfiles/simInfo_step1e4.dat
set infoFile [open "~/sphere/simfiles/simInfo_step1e4.dat" "w"] 		
puts $infoFile "# Scalar"
puts $infoFile "# timeStep / stepSize / numSteps / lengthScale"
puts $infoFile "$dt_md $stepSize $numSteps $L0"
close $infoFile

#########################################################################
# SETUP 
#########################################################################

# exec rm -r /home/btpj/s1mokhal/sphere/vtkfiles
# exec mkdir /home/btpj/s1mokhal/sphere/vtkfiles



exec touch /home/btpj/s1mokhal/sphere/simfiles/triel_debug.dat
exec rm /home/btpj/s1mokhal/sphere/simfiles/triel_debug.dat

exec touch  /home/btpj/s1mokhal/sphere/simfiles/forces_debug.dat
exec rm /home/btpj/s1mokhal/sphere/simfiles/forces_debug.dat

exec touch  /home/btpj/s1mokhal/sphere/simfiles/lbtracers_debug.dat
exec rm /home/btpj/s1mokhal/sphere/simfiles/lbtracers_debug.dat

exec touch  /home/btpj/s1mokhal/sphere/simfiles/integ_debug.dat
exec rm /home/btpj/s1mokhal/sphere/simfiles/integ_debug.dat

exec touch  /home/btpj/s1mokhal/sphere/simfiles/lb_debug.dat
exec rm /home/btpj/s1mokhal/sphere/simfiles/lb_debug.dat

exec touch  /home/btpj/s1mokhal/sphere/simfiles/lb_cells_debug.dat
exec rm /home/btpj/s1mokhal/sphere/simfiles/lb_cells_debug.dat

exec touch /home/btpj/s1mokhal/sphere/simfiles/verlet_debug.dat
exec rm /home/btpj/s1mokhal/sphere/simfiles/verlet_debug.dat


#tools needed for system setup (Juropa or local)
set tclToolsDir "~/sphere/tclTools"
source "$tclToolsDir/addCell_tools.tcl"
source "$tclToolsDir/writeVtk_folded.tcl"
# source "$tclToolsDir/addCell_fixed.tcl"

set numParts 3
set parts {353 356 358}

set partCsv353 [ open "~/sphere/simfiles/partInfo353.csv" "w"]
puts $partCsv353 "posx,posy,posz,vx,vy,vz,fx,fy,fz,lbvx,lbvy,lbvz,nvx,nvy,nvz,t"
set partCsv356 [ open "~/sphere/simfiles/partInfo356.csv" "w"]   
puts $partCsv356 "posx,posy,posz,vx,vy,vz,fx,fy,fz,lbvx,lbvy,lbvz,nvx,nvy,nvz,t"
set partCsv358 [ open "~/sphere/simfiles/partInfo358.csv" "w"]
puts $partCsv358 "posx,posy,posz,vx,vy,vz,fx,fy,fz,lbvx,lbvy,lbvz,nvx,nvy,nvz,t"

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

setmd warnings 0
# setmd verlet_reuse 0

part 353 pos 11.12083 10 5.127 virtual 1
part 356 pos 11.68887 10.4057 5.3111 virtual 1
part 358 pos 11.68887 9.5942 5.3111 virtual 1

inter 0 triel 353 356 358 $maxStretch $kS 0
part 353 bond 0 356 358

# setting up the fluid with or without using gpu
lbfluid agrid $gridsize dens $rho visc $nu tau $dt_lb friction $zeta ext_force 0 0 0
#setting themostat
thermostat lb $kbT

# cellsystem domain_decomposition -no_verlet_list
# cellsystem nsquare

############## ADD WALL ##################

# walls located at z=1 and z=boxz-1 in x-y-plane
# lbboundary wall dist [expr $H/2] normal 0. 0. 1. velocity [expr $u/2] 0 0 
# lbboundary wall dist [expr -$boxz + ($H/2)] normal 0. 0. -1. velocity [expr -$u/2] 0 0

####### SIMULATION ######

# update position of all particles
 integrate 0
# Checkpoint only without using gpu
# writeCheckpoint start

set startT [clock seconds]

set partPos [part 356 print pos]
part 356 pos [expr [lindex $partPos 0] + 2] [expr [lindex $partPos 1]] [expr [lindex $partPos 2]]

puts "Starting simulation"

for {set step 0} {$step < $numSteps} {incr step} {
    
    integrate $stepSize
    if {fmod($step, $distDone)==0} { puts "Done $step out of $numSteps steps." }
	 
  for {set j 0} { $j < $numParts } {incr j} {
      set partCsv [ open "~/sphere/simfiles/partInfo[lindex $parts $j].csv" "a"]
      set dataPart [part [lindex $parts $j] print pos v f]
      set dataFlow [lbfluid print_interpolated_velocity [lindex $dataPart 1] [lindex $dataPart 2] [lindex $dataPart 3]]
      set data [concat $dataPart $dataFlow]
      foreach x $data { 
	  puts -nonewline $partCsv "$x,"
      }
      puts -nonewline $partCsv "$step\n"
      close $partCsv
  }
    if { $step == 1 || $step == 4} { 
	set nodeCsv [ open "~/sphere/simfiles/nodeInfo$step.csv" "a"]
	puts $nodeCsv "nx ny nz ux uy uz"
	for {set i 0} { $i < $boxx } {incr i} { 
	    for {set j 0} { $j < $boxy} {incr j} { 
		for {set k 0} { $k < $boxz} {incr k} { 
		    puts $nodeCsv "$i $j $k [lbnode $i $j $k print u]"
		}
	    }
	}
    }
		
}

# store simulation time in log.dat
set endT [clock seconds]
set logfile [open "log.dat" "w"]
puts $logfile "Needed [expr $endT-$startT] seconds ([expr ($endT-$startT)/60] minutes)."
