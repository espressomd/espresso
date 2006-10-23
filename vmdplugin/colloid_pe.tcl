puts "[code_info]"

set vtf_filename "colloid_pe.vtf"

# parameters
set id_colloid 0
set r_colloid 2.0
set n_colloid 5

set id_monomer 1
set r_monomer 0.5
set n_polymer 5
set l_polymer 30
set q_monomer +1

set box_length 20.0

set n_cycle 50
set n_steps 1000

set time_step 0.002
set skin 0.4
set temp 1.0
set gamma 1.0

set pressure 0.02
set piston_mass 0.001
set gamma_0 0.5
set gamma_v 0.001

set fene_epsilon 1.0
set lj_epsilon 1.0
set lj_sigma 1.0

# derived parameters
set z_colloid [expr ( - $q_monomer * $l_polymer ) / $n_colloid]
set lj_cutoff [expr $lj_sigma*pow(2.,1./6.)]
set lj_shift [expr 1./4.*$lj_epsilon]

# init simulation
setmd box_l $box_length $box_length $box_length
setmd time_step $time_step
setmd skin $skin
setmd max_num_cells 1000
thermostat langevin $gamma $temp
if {[has_feature NPT]} then {
    integrate set npt_isotropic $pressure $piston_mass 
    thermostat set npt_isotropic $temp $gamma_0 $gamma_v
} else {
    integrate set nvt
}

puts "Simulate PE close to a charged colloid"

puts "Setting global parameters..."
puts "  box_length=[setmd box_l]"
puts "  time_step=[setmd time_step]"
puts "  skin=[setmd skin]"
puts "  integrate=[integrate]"
puts "  thermostat=[thermostat]"

puts "Setting up particles and interactions..."
# steric LJ potentials
inter $id_colloid $id_colloid \
    lennard-jones \
    $lj_epsilon $lj_sigma \
    $lj_cutoff $lj_shift \
    [ expr $r_colloid+$r_colloid - $lj_sigma ]

inter $id_colloid $id_monomer \
    lennard-jones \
    $lj_epsilon $lj_sigma \
    $lj_cutoff $lj_shift \
    [ expr $r_monomer+$r_colloid - $lj_sigma ]

inter $id_monomer $id_monomer \
    lennard-jones \
    $lj_epsilon $lj_sigma \
    $lj_cutoff $lj_shift \
    [ expr 2*$r_monomer - $lj_sigma ]

# bond potential
inter 0 fene $fene_epsilon [expr 4*$r_monomer]

# create colloids
for { set pid 0 } { $pid < $n_colloid } { incr pid } {
    part $pid \
	pos \
	[expr [t_random]*$box_length] \
	[expr [t_random]*$box_length] \
	[expr [t_random]*$box_length] \
	q [expr $z_colloid] \
	type $id_colloid
}

polymer $n_polymer $l_polymer [expr 2*$r_monomer] \
    start $pid \
    charge $q_monomer \
    FENE 0 \
    mode SAW

puts "Writing configurations to $vtf_filename..."
set vtf_file [open $vtf_filename w]
writevsf $vtf_file
writevcf $vtf_file

set min 0
set cap 10
set warmup_finished 0
while { ! $warmup_finished } {
    # set ljforcecap
    inter ljforcecap $cap
    # integrate a number of steps, e.g. 20
    integrate 20
    for { set pid 0 } { $pid < [setmd n_part] } { incr pid } {
	part $pid v 0 0 0
    }
    # check the status of the sytem
    set min_colloid_monomer [analyze mindist $id_colloid $id_monomer]
    set min_colloid_colloid [analyze mindist $id_colloid $id_colloid]
    set min_monomer_monomer [analyze mindist $id_monomer $id_monomer]
    set warmup_finished \
	[expr \
	     $min_colloid_monomer > $r_colloid+$r_monomer && \
	     $min_colloid_colloid > 2*$r_colloid && \
	     $min_monomer_monomer > 2*$r_monomer ]

    # this is a shortcut for 'set cap [expr $cap+10]'
    incr cap 10
}
writevcf $vtf_file

puts "Warmup finished. Minimal distance now $min_colloid_monomer, $min_colloid_colloid, $min_monomer_monomer"
# turn off the ljforcecap, which is done by setting the 
# force cap value to zero:
inter ljforcecap 0
# reset velocities
for { set pid 0 } { $pid < [setmd n_part] } { incr pid } {
    part $pid v 0 0 0
}
# reset simulation time
setmd time 0.0

puts "Tuning P3M..."
puts [inter coulomb 1 p3m tunev2 accuracy 0.01]
#inter coulomb 1 p3m 13.9629 8 3 0.127621 0.000832468

puts "Starting simulation..."
set i 0 
while { $i<$n_cycle } {
    puts $i
    integrate $n_steps

    writevcf $vtf_file

    incr i
}
close $vtf_file

puts "Simulation finished."

