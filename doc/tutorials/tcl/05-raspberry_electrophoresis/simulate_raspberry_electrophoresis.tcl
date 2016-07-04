#set the size of the simulation box
set box_length 48
setmd box_l $box_length $box_length $box_length


#initiate the random number generator (one seed for each node) using the process id.
set n_nodes [setmd n_nodes]
set start_node 1
set pid [pid]
puts "pid = $pid"
#here a simple random number generator is used to create different seed for the different nodes
for {set i 1} {$i <= $n_nodes} {incr i} {
    lappend seeds [expr $i + $n_nodes*$start_node*$pid]
}
#send the created seed to ESPResSo
eval t_random seed $seeds

# a proc that returns pi
proc Pi {} {return 3.1415926535897931}

#set the skin used for creating the Verlet lists
setmd skin 0.3
#set the periodicity of the box
setmd periodic 1 1 1


set radius_col 3.0

#setup the interactions necessary for creating the raspberry
#the subscript c is for colloid and s is for salt (also used for the surface beads)
set c_shift 0.25 ;# LJ shift.
set eps_ss 1.0  ;# LJ epsilon between the colloid's surface particles.
set sig_ss 1.0 ;# LJ sigma between the colloid's surface particles.
set eps_cs 48.0 ;# LJ epsilon between the colloid's centeral particle and surface particles. 
set sig_cs $radius_col ;# LJ sigma between the colloid's centeral particle and surface particles (colloid's radius). 
set a_eff 0.32 ;# effective hydrodynamic radius of a bead due to the descreteness of LB.

set harmonic_radius 3
#the harmonic potential pulls surface beads towards the central colloid bead
inter 0 harmonic 3000. $harmonic_radius
#the LJ potential between surface beads causes them to be roughly equidistant on the colloid surface
inter 1 1 lennard-jones $eps_ss [expr 1.0*$sig_ss] [expr $sig_ss*pow(2.0,1.0/6.0)] $c_shift 0.0
#the LJ potential with the central bead keeps all the beads from simply collapsing into the center
inter 1 0 lennard-jones $eps_cs [expr $sig_cs] [expr $sig_cs*pow(2.0,1.0/6.0)] $c_shift 0.0
set n 0


set nTot 0
set n_cycle 1000
set integ_steps 150

#for the warmup we use a Langevin thermostat with an extremely low temperature and high friction coefficient such that the trajectories roughly follow 
#the gradient of the potential while not accelerating too much
thermostat set langevin 0.00001 40.0

puts "Creating raspberry"
set center [expr $box_length/2.0]

set colPosX $center
set colPosY $center
set colPosZ $center

set q_col 40 ;# note the charge on the colloid is actually negative
set n_col_part [expr int(4*[Pi]*$radius_col**2 + 1)] ;# Number of particles making up the raspberry (surface particles + the centeral particle).
part $nTot pos $colPosX $colPosY $colPosZ fix 1 1 1 q [expr 0-$q_col] type 0
set n0 $nTot
incr nTot
# this loop create n_col_part surface beads randomly placed on a spherical shell Rvec away from the central bead
for {set i 1} {$i < $n_col_part} {incr i} {
    set Rvec [vec_random $radius_col]
    set xpos [expr $colPosX + [lindex $Rvec 0]]
    set ypos [expr $colPosY + [lindex $Rvec 1]]
    set zpos [expr $colPosZ + [lindex $Rvec 2]]
    part $nTot pos $xpos $ypos $zpos type 1
    part $nTot bond 0 $n0
    incr nTot
}

#here we let the bead positions relax. The LJ potential with the central bead combined with the
#harmonic bond keep the monomers roughly radius_col away from the central bead. The LJ
#between the surface beads cause them to distribute more or less evenly on the surface.
inter forcecap 1000
setmd time_step 0.001
for {set i 0} {$i <= $n_cycle} {incr i} {
    puts -nonewline "\r # $i of $n_cycle"
    flush stdout
    integrate $integ_steps
}


inter 0 harmonic 0 0

puts "min dist to the center pre moving of surface beads = [analyze distto 0]"

#this loop moves the surface beads such that they are once again exactly radius_col away from the center
for {set i 1} {$i < $n_col_part} {incr i} {
    set xpos [lindex [part $i print pos] 0]
    set ypos [lindex [part $i print pos] 1]
    set zpos [lindex [part $i print pos] 2]
    set rad_pos [expr sqrt(($xpos-$colPosX)**2+($ypos-$colPosY)**2+($zpos-$colPosZ)**2)]
    set xpos [expr ($xpos-$colPosX)*$radius_col/$rad_pos+$colPosX]
    set ypos [expr ($ypos-$colPosY)*$radius_col/$rad_pos+$colPosY]
    set zpos [expr ($zpos-$colPosZ)*$radius_col/$rad_pos+$colPosZ]
    part $i pos $xpos $ypos $zpos
}
puts "min dist to the center post moving of surface beads = [analyze distto 0]"

#setting min_global_cut is necessary when there is no interaction defined with a range larger than the colloid
#such that the virtual particles are able to communicate their forces to the real particle at the center of the colloid
setmd min_global_cut [expr $radius_col + 1]

#here we calculate the center of mass position (cm) and the moment of inertia (momI)
set cmx 0
set cmy 0
set cmz 0
set momI 0

#Calculate the center of mass position
for {set i 0} {$i<$n_col_part} {incr i} {
    set cmx [expr $cmx + [lindex [part $i print pos] 0]]
    set cmy [expr $cmy + [lindex [part $i print pos] 1]]
    set cmz [expr $cmz + [lindex [part $i print pos] 2]]
}
set cmx [expr $cmx/$n_col_part]
set cmy [expr $cmy/$n_col_part]
set cmz [expr $cmz/$n_col_part]

#calculating the moment of inertia of the colloid
for {set i 1} {$i<$n_col_part} {incr i} {
    set momI [expr $momI + pow(($cmx-[lindex [part $i print pos] 0]),2)]
    set momI [expr $momI + pow(($cmy-[lindex [part $i print pos] 1]),2)]
    set momI [expr $momI + pow(($cmz-[lindex [part $i print pos] 2]),2)]
}
#Note the 2/3 factor is for a spherical shell
set momI [expr  2.0*$momI/3.0]

#note that the real particle must be at the center of mass of the colloid because of the integrator
puts "moving central particle from [part 0 print pos] to $cmx $cmy $cmx"
part 0 fix 0 0 0
part 0 pos $cmx $cmy $cmz fix 0 0 0 mass $n_col_part rinertia $momI $momI $momI
for {set i 1} {$i<$n_col_part} {incr i} {
    part $i virtual 1 vs_auto_relate 0
}
puts "Raspberry made"




puts "# Adding the positive ions"
set salt_rho 0.001
set volume [expr pow($box_length,3)]
set N_counter_ions [expr floor(($volume*$salt_rho) + $q_col + 0.5)]
for {set i 0} {$i < $N_counter_ions} {incr i} {
    set xpos [expr $box_length*[t_random]]
    set ypos [expr $box_length*[t_random]]
    set zpos [expr $box_length*[t_random]]
    #make sure the ion is placed outside of the colloid
    set r2 [expr ($xpos-$center)**2+($ypos-$center)**2+($zpos-$center)**2]
    if {$r2 > [expr pow($radius_col,2)+5]} {
        part $nTot pos $xpos $ypos $zpos type 2 q 1
        incr nTot
    } else {
        set i [expr $i-1]
    }
}
set nEndPosIon $nTot
set N_counter_ions $i
puts "\n# Added $N_counter_ions positive ions"

puts "\n# Adding the negative ions"
set nStartNegIon $nTot
set N_co_ions [expr $N_counter_ions - $q_col]
for {set i 0} {$i < $N_co_ions} {incr i} {
    set xpos [expr $box_length*[t_random]]
    set ypos [expr $box_length*[t_random]]
    set zpos [expr $box_length*[t_random]]
    #make sure the ion is placed outside of the colloid
    set r2 [expr ($xpos-$center)**2 + ($ypos-$center)**2 + ($zpos-$center)**2]
    if {$r2 > [expr pow($radius_col,2)+1]} {
        part $nTot pos $xpos $ypos $zpos type 3 q -1
        incr nTot
    } else {
        set i [expr $i-1]
    }
}

# Exerting the applied electric force on the particles and checking the charge neutrality of the system
set Qtot 0
set Efield 0.1 ;#an electric field of 0.1 is the upper limit of the linear response regime for this model
for {set i 0} {$i < $nTot} {incr i} {
    set chargei [part $i print q]
    part $i ext_force [expr $chargei*$Efield] 0 0
    set Qtot [expr $Qtot + $chargei]
}
#make sure that the system has a net charge of zero, otherwise the system as a whole will slowly accelerate
if {$Qtot > 0.0001 || $Qtot < -0.0001} {
    puts "net charge of $Qtot !!! Exiting"
    puts "# Colloid charge $q_col Positive ions $N_counter_ions Negative ions $N_co_ions"
    exit 0
}

# WCA interactions for the ions, essentially giving them a finite volume
inter 0 2 lennard-jones $eps_ss $sig_ss [expr ($sig_ss)*pow(2.0,1.0/6.0)] $c_shift [expr $sig_cs-1+$a_eff]
inter 0 3 lennard-jones $eps_ss $sig_ss [expr ($sig_ss)*pow(2.0,1.0/6.0)] $c_shift [expr $sig_cs-1+$a_eff]
inter 2 2 lennard-jones $eps_ss $sig_ss [expr $sig_ss*pow(2.0,1.0/6.0)] $c_shift 0.0
inter 2 3 lennard-jones $eps_ss $sig_ss [expr $sig_ss*pow(2.0,1.0/6.0)] $c_shift 0.0
inter 3 3 lennard-jones $eps_ss $sig_ss [expr $sig_ss*pow(2.0,1.0/6.0)] $c_shift 0.0

# Langevin thermostat for warmup before turning on the LB.
set kT 1
thermostat langevin $kT 1.0 ;#the second parameter is the friction coefficient of the beads

#fix the colloid for equilibrating the ions
for {set i $n_col_part} {$i < $n_col_part} {incr i} {
    part $i fix 1 1 1
}

set vmd "yes"
if {$vmd == "yes"} {prepare_vmd_connection "colloid_movie" 100}

set ljcap 100
set CapSteps 1000
for {set i 0} {$i < $CapSteps} {incr i} {
    inter forcecap $ljcap
    puts -nonewline "\r# Equilibrating the ions (without electrostatics): step $i of $CapSteps"
    flush stdout
    integrate 100
    if {$vmd == "yes"} {imd positions}
    set ljcap [expr $ljcap + 5]
}
inter forcecap 0

#let the colloid move now that bad overlaps have been eliminated
for {set i $n_col_part} {$i < $n_col_part} {incr i} {
    part $i fix 0 0 0
}


# Turning on the electrostatics
set errorCoulomb 0.01
puts "\n# p3m starting..."
flush stdout
set bjerrum 2
inter coulomb $bjerrum p3m tunev2 accuracy 0.001
puts "# p3m started!"

setmd time_step 0.01
#note that it is important to set all of the particle velocities to zero just before initializing the LB so that the total system momentum is zero
kill_particle_motion
lbfluid gpu dens 1 visc 3 agrid 1 tau 0.01 friction 20

#replace the Langevin thermostat with the lb thermostat
thermostat off
thermostat lb $kT

set posVsTime [open "posVsTime.dat" "w"] ;# file where the raspberry position will be written
for {set curSample 0} {$curSample < 10000} {incr curSample} {
  integrate 1000
  if {$vmd == "yes"} {imd positions}
  puts $posVsTime "[setmd time] [part 0 print pos]"
  puts "[setmd time] [part 0 print pos]"
}