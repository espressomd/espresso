# ----------------------------------------------------------------------
# Copyright (C) 2016 Dr. Bogdan Tanygin<b.m.tanygin@gmail.com>
# All rights reserved.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ----------------------------------------------------------------------

# Constants - math
set PI 3.14159
#[SI] Constants - physics
set Na 6.022E23
set mu0 [expr 4.0 * $PI * 1E-7]
set muB [expr 9.27400968 * 1E-24]
set kb [expr 1.3806488 * 1E-23]
set Oe_to_Am 79.577

#[SI] Materials

# Medium
#[Pa*s] Kerosene dynamic viscosity
set eta_car 0.00164
#[J] Hamaker constant
set A_H 4.0E-20
set T 300.0

# Nanoparticle
#[m] magnetite unit cell size - a cubic spinel structure with space group Fd3m (above the Verwey temperature) \cite{Fleet1981}
set a0 0.8397E-9
#[A/m] magnetite saturation magnetization
set Ms 478.0E3
#[kg/m^3] magnetite mass density
set rop 5240.0

# Oleic acid
#[m^-2] Surface density of oleic acid at the 50% coating \cite{Fertman}
set N_oa 1.0E18
#[m] the oleic acid molecule length
set delta 1.97E-9;
#[kg] the oleic acid molecule mass
set mass_o [expr 282.47E-3/$Na]
# Level of the coverage of the particle surface by the oleic acid
set k_o 0.5
# Flocculation parameter \cite{Cebers1987}
set k_f 1.1

#[SI] Ferrofluid
# The volume fraction value of a FF dispersed phase
#set phi_v 1.2E-6
#set phi_v 5E-2
#set phi_v 5E-5

# The large particles fraction \cite{Ivanov1996}
#set PHI_l 7.5E-2

set PHI_l 0.9

#set PHI_l 8E-1
#set PHI_l 0.99
#[m] the mean nanoparticle diameter. Hard part: a ferromagnetic core & double nonmagnetic surface
set d_mean_hard 11.5E-9
set d_mean_hydrod [expr $d_mean_hard+2.0*$delta]
set d_mean_mag_core [expr $d_mean_hard-2.0*$a0]
set V_mean_hard [expr $PI*pow($d_mean_hard,3.0)/6.0]
set S_mean_hard [expr $PI*pow($d_mean_hard,2.0)]
set N_oleic_mean [expr $S_mean_hard*$k_o*$N_oa/0.5]
set V_mean_mag_core [expr $PI*pow($d_mean_mag_core,3)/6.0]
set m_mean [expr $Ms*$V_mean_mag_core]
set M_mean [expr $rop*$V_mean_hard+$N_oleic_mean*$mass_o]
#[m] the large nanoparticles diameter \cite{Ivanov1996}
set d_l_hard 21.6E-9
set d_l_hydrod [expr $d_l_hard+2*$delta]
set d_l_mag_core [expr $d_l_hard-2*$a0]
set V_l_hard [expr $PI*pow($d_l_hard,3.0)/6.0]
set S_l_hard [expr $PI*pow($d_l_hard,2.0)]
set N_oleic_l [expr $S_l_hard*$k_o*$N_oa/0.5]
set V_l_mag_core [expr $PI*pow($d_l_mag_core,3)/6.0]
set m_l [expr $Ms*$V_l_mag_core]
set M_l [expr $rop*$V_l_hard+$N_oleic_l*$mass_o]

# Dimensionless model
set SIGMA [expr $d_mean_hydrod]
set M0 [expr $M_mean]
set H0 [expr pow($kb*$T/$mu0,0.5)/pow($SIGMA,1.5)]
set m0 [expr sqrt($kb*$T*pow($SIGMA,3.0)/$mu0)]
set t0 [expr sqrt($M0*pow($SIGMA,2.0)/($kb*$T))]
set gamma0 [expr 3*$PI*$eta_car*pow($SIGMA,2)/sqrt($M0*$kb*$T)]

# Computational experiment setup
# TODO: time_step to increase till the stability loss
#setmd time_step [expr 0.00005]

# max for pure BH or master
# setmd time_step [expr 0.0005]

setmd time_step [expr 0.05]
set n_part 10000

#set H_ext_Oe 500.0
set H_ext_Oe 0.0

set n_part_small [expr $n_part*(1-$PHI_l)]
#set box_l [expr pow($n_part*$PI*($PHI_l*pow($d_l_hard/$SIGMA,3.0)+(1-$PHI_l)*pow($d_mean_hard/$SIGMA,3.0))/(6.0*$phi_v),1./3.)]
set start_lattice_a [expr 2*$d_l_hydrod/$SIGMA]
set buf_l [expr 2*$start_lattice_a]
set box_l [expr $start_lattice_a*pow($n_part,1/3.0)+$buf_l]

set V_carr [expr pow($box_l,3.0)]
set box_l_x [expr $box_l*1.0]
set box_l_y [expr $box_l*1.0]
set box_l_z [expr $box_l/1.0]
setmd box_l $box_l_x $box_l_y $box_l_z
setmd periodic 0 0 0
setmd skin 0.4
set temp 1; set gammat 1; set gammar 1
thermostat langevin $temp $gammat $gammar $gammar $gammar

set coord_shift [expr 0.1*$start_lattice_a]
set posx $coord_shift
set posy $coord_shift
set posz $coord_shift

for {set i 0} { $i < $n_part } {incr i 1} {
set rnd_val_for_p_type [expr [t_random]]
set theta [expr $PI*[t_random]]
set phi [expr 2.0*$PI*[t_random]]

set posx [expr $posx+$start_lattice_a]
if {$posx > $box_l_x-$buf_l} {
    set posx [expr $coord_shift]
    set posy [expr $posy+$start_lattice_a]
    if {$posy > $box_l_y-$buf_l} {
        set posy [expr $coord_shift]
        set posz [expr $posz+$start_lattice_a]
        if {$posz > $box_l_z-$buf_l} {
            puts "Not enough box_l for particles!"
        }
    }
}

if {$rnd_val_for_p_type < $PHI_l} {
    set dx [expr ($m_l/$m0)*sin($theta)*cos($phi)]
	set dy [expr ($m_l/$m0)*sin($theta)*sin($phi)]
	set dz [expr ($m_l/$m0)*cos($theta)]
	set massval [expr $M_l/$M0]
	set I_val [expr 0.1*($M_l/$M0)*pow($d_l_hydrod/$SIGMA,2.0)]
	set gamma_tr_val [expr $gamma0*$d_l_hydrod/$SIGMA]
	set gamma_rot_val [expr $gamma0*pow($d_l_hydrod/$SIGMA,3.0)/3.0]

	part $i pos $posx $posy $posz type 1 mass $massval rinertia $I_val $I_val $I_val dip $dx $dy $dz gamma $gamma_tr_val gamma_rot $gamma_rot_val $gamma_rot_val $gamma_rot_val 
} else {
    set dx [expr ($m_mean/$m0)*sin($theta)*cos($phi)]
	set dy [expr ($m_mean/$m0)*sin($theta)*sin($phi)]
	set dz [expr ($m_mean/$m0)*cos($theta)]
	set massval [expr 1.0]
	set I_val [expr 0.1]
	set gamma_tr_val [expr $gamma0*1.0]
	set gamma_rot_val [expr $gamma0*1.0/3.0]

	part $i pos $posx $posy $posz type 0 mass $massval rinertia $I_val $I_val $I_val dip $dx $dy $dz gamma $gamma_tr_val gamma_rot $gamma_rot_val $gamma_rot_val $gamma_rot_val
}
}

puts "[part 10]"

set sig [expr 1.0]; set cut [expr 7.0*$sig]
set eps 1.5; set shift [expr 0.25]
inter 0 0 lennard-jones $eps $sig $cut $shift 0

set sig [expr ($d_l_hydrod/$SIGMA)]; set cut [expr 2.0*$sig]
set eps 1.5; set shift [expr 0.25]
inter 1 1 lennard-jones $eps $sig $cut $shift 0

set sig [expr ((0.5*$d_l_hydrod+0.5*$SIGMA)/$SIGMA)]; set cut [expr 6.0*$sig]
set eps 1.5; set shift [expr 0.25]
inter 0 1 lennard-jones $eps $sig $cut $shift 0

if { [regexp "ROTATION" [code_info]] } {
set deg_free 6
} else { set deg_free 3 }

for {set cap 1} {$cap < 200} {incr cap 20} {
puts "t=[setmd time] E=[analyze energy total]"
inter forcecap $cap; integrate 20
}
inter forcecap 0

#imd connect 10000
prepare_vmd_connection "vmdfile" 10000

for {set i 0} { $i < 5 } { incr i} {
set temp [expr [analyze energy kinetic]/(($deg_free/2.0)*$n_part)]
puts "t=[setmd time] E=[analyze energy total], T=$temp"
integrate 10
imd positions
}

inter magnetic 1.0 bh-gpu
#inter magnetic 1.0 dds-gpu

set H_demag 0.0
set dipm_obs [ observable new com_dipole_moment all ]

puts "=============================="
puts "t_epoch | t_in_silico | Mz | T"
puts "=============================="

for {set i 0} { $i < 1E20 } { incr i} {
set temp [expr [analyze energy kinetic]/(($deg_free/2.0)*$n_part)]
#puts "t=[setmd time] E=[analyze energy total], T=$temp"

#puts "t=[setmd time] dipm=[observable $dipm_obs print]"
#set H_demag [expr -[lindex [ observable $dipm_obs print ] 2]/$V_carr ]
set Mz [expr [lindex [ observable $dipm_obs print ] 2]/$V_carr]
set H_demag 0.0
set B_mag [expr (($H_ext_Oe*$Oe_to_Am)/$H0+$H_demag)]
constraint ext_magn_field 0 0 $B_mag

integrate 10
imd positions

puts "[clock seconds],[setmd time],$Mz,$temp"

}

imd disconnect

set f [open "config_$i" "w"]
blockfile $f write tclvariable {box_l density}
blockfile $f write variable box_l
blockfile $f write particles {id pos type}
close $f
