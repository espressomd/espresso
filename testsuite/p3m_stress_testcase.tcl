proc rotate_z {} {
	# rotates the system by 90 degrees around the z axis
	for { set i 0 } { $i < [setmd n_part] } { incr i } {
		set pos_old [part $i print pos]
		set pos_new [list [expr -1.0 * [lindex $pos_old 1]] [lindex $pos_old 0 ] [lindex $pos_old 2]]
		part $i pos [lindex $pos_new 0] [lindex $pos_new 1] [lindex $pos_new 2]
	}
}

require_feature "ELECTROSTATICS"
require_feature "FFTW"

setmd time_step 0.01
setmd skin 0.4


set fin [open "p3m_stress.data" "r"]
blockfile $fin read auto
blockfile $fin read auto
blockfile $fin read auto

inter 0 0 lennard-jones 0.0 0.0 [expr double(1.0)*2**(1.0/6.0)]  0 0
set eps 1e-7

set bjerrum 2.0
#156  6   4.19304e-01 1.03503e+01 9.86970e-09 1036 accuracy 1e-8 for this very system
puts [inter coulomb $bjerrum p3m 4.19304 256 7 1.03503]

set N [setmd n_part]
set vol [expr pow([lindex [setmd box_l] 0], 3)]
puts "results for system with $N particles in vol $vol"
integrate 0
set e_c [analyze energy coulomb]
set ptensor_coul [analyze stress_tensor coulomb]
for { set k 0 } { $k < 9 } { incr k } {
	set pc($k) [lindex $ptensor_coul [expr $k+1]]
}

#puts [format "%1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e" $pc(0) $pc(1) $pc(2) $pc(3) $pc(4) $pc(5) $pc(6) $pc(7) $pc(8)]
#puts "rotate\n"

rotate_z
integrate 0
set ptensor_coul [analyze stress_tensor coulomb]
for { set k 0 } { $k < 9 } { incr k } {
	set pc_2($k) [lindex $ptensor_coul [expr $k+1]]
}
#puts [format "%1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e" $pc_2(0) $pc_2(1) $pc_2(2) $pc_2(3) $pc_2(4) $pc_2(5) $pc_2(6) $pc_2(7) $pc_2(8)]
#6.70541e-01 -8.30360e-03 7.47335e-03 -8.30360e-03 6.94174e-01 1.00401e-02 7.47335e-03 1.00401e-02 -1.06683e+00
#rotate
#
#6.94174e-01 8.30360e-03 -1.00401e-02 8.30360e-03 6.70541e-01 7.47335e-03 -1.00401e-02 7.47335e-03 -1.06683e+00

# consistency check s_xx = s_2_yy
#					s_xy = -s_2_xy
#					s_xz = -s_2_yz
#					s_yx = -s_2_yx
#					s_yy = s_2_xx
#					s_yz = -s_2_xz
#					s_zx = s_2_yz
#					s_zy = -s_2_yz
#					s_zz = s_2_zz

set rotate_map(0) 4
set rotate_map(1) 1
set rotate_map(2) 5
set rotate_map(3) 3
set rotate_map(4) 0
set rotate_map(5) 2
set rotate_map(6) 7
set rotate_map(7) 6
set rotate_map(8) 8

set signs(0) 1
set signs(1) -1
set signs(2) 1
set signs(3) -1
set signs(4) 1
set signs(5) -1
set signs(6) 1
set signs(7) -1
set signs(8) 1

for { set i 0 } { $i < 9 } { incr i } { 
	if { [expr abs($pc($i) - $pc_2($rotate_map($i))*$signs($i)) > $eps] } {
		error "p3m stress tensor broken $i"
		break
	} 
}

# trace check 
set p_trace [expr ($pc(0) + $pc(4) + $pc(8))/3.0]
set p_trace_2 [expr ($pc_2(0) + $pc_2(4) + $pc_2(8))/3.0]

if { [expr abs($p_trace - $e_c/(3.0 * $vol))/abs($e_c/(3.0*$vol)) >10.0 * $eps ] } {
	error "tensor trace and pressure deviate too much [expr abs($p_trace - $e_c/(3.0 * $vol))/abs($e_c/(3.0 *$vol))]"
}

if { [expr abs($p_trace - $p_trace_2) > $eps] } {
	error "trace differs after rotation"
}



exit 

