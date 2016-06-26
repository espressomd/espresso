# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
#  
# This file is part of ESPResSo.
#  
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.   
  
set tcl_precision 15

##########################
set n_particle 100

# number of periodic images
set n_images 2000
#########################


##############################################
# volume fraction
set rho 0.3

set particle_radius 0.5
set dipole_lambda 3.0

#################################################

set external_H_z 0.5

set int_time_step 0.002

set dipole_modulus [expr sqrt($dipole_lambda * pow(2*$particle_radius,3))]

# setting box paramters
set box_l [expr pow(((4 * $n_particle * 3.141592654) / (3*$rho)), 1.0/3.0)*$particle_radius]  
puts "box_lenght $box_l"
set skin 0.5


# give Espresso some parameters
setmd time_step $int_time_step
setmd skin $skin
setmd box_l $box_l $box_l $box_l
setmd periodic 1 1 0
setmd max_num_cells 2500

t_random seed 3

# put the  particles into the box 
puts "\nput the particles in the box..."

for {set i 0} {$i < $n_particle} {incr i} {
    # the for the real particeles
    set posx [expr $box_l *[expr [t_random]]]	
    set posy [expr $box_l *[expr [t_random]]]	
    set posz [expr $box_l *[expr [t_random]]]	
    set costheta [expr 2*[expr [t_random]] - 1]
    set sintheta [expr sin(acos($costheta))]
    set phi [expr 2*3.1415926536*[expr [t_random]]] 
    set dipx [expr $sintheta*cos($phi)*$dipole_modulus]
    set dipy [expr $sintheta*sin($phi)*$dipole_modulus]
    set dipz [expr $costheta*$dipole_modulus]
    
    
    part $i pos $posx $posy $posz type 0  dip $dipx $dipy $dipz fix 0 0 1
    
}      

# Relax particle configuration
inter 0 0 lennard-jones 1 1 1.12 auto 
puts "before relaxation [analyze energy]"
minimize_energy 0 35000 0.1 0.01
puts "after relaxation [analyze energy]"
inter 0 0 lennard-jones 0 0 0 0 0

# Remove particles that are in the mdlc gap region
set gap 2.0
set max_z [expr $box_l -$gap]
for {set i 0} {$i < $n_particle} {incr i} {
  set z [lindex [part $i print folded_pos] 2]
  if { ($z> $max_z) || ($z <0) } {
    part $i delete
  }
}
puts "[setmd n_part] particles remain after clearing gap"





thermostat off 


inter magnetic 1 mdds n_cut $n_images

puts "\n\ncalculating forces and torques using direct summation ..."


integrate 0


set force_torque_data [open "mdlc_reference_data_forces_torques.dat" w]

for {set i 0} {$i < $n_particle} {incr i} {
    if { "[part $i print pos]" != "na" } {
      puts $force_torque_data "$i [part $i print folded_position dip force torque_lab]"
    }
}

close $force_torque_data


set energy_data [open "mdlc_reference_data_energy.dat" w]
puts $energy_data [analyze energy magnetic]
close $energy_data


puts "\n\n\ndone" 

exit 0
