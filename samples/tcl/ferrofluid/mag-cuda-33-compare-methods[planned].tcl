# ----------------------------------------------------------------------
# Copyright (C) 2016,2017 Dr. Bogdan Tanygin<b.m.tanygin@gmail.com>
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

proc anykey {{msg "Please set a proper camera position & other VMD configuration. Then, press any key: "}} {
    set stty_settings [exec stty -g]
    exec stty raw -echo
    puts -nonewline $msg
    flush stdout
    read stdin 1
    exec stty $stty_settings
    puts ""
}

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
#[A/m] magnetite saturation magnetization \cite{Kemp2016}
#set Ms 478.0E3
set Ms 412.0E3
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

# The large particles fraction \cite{Ivanov1996}
#set PHI_l 7.5E-2
set PHI_l 0.5

#[m] the mean nanoparticle diameter. Hard part: a ferromagnetic core & double nonmagnetic surface
set d_small_hard 11.5E-9
set d_small_hydrod [expr $d_small_hard+2.0*$delta]
set d_small_mag_core [expr $d_small_hard-2.0*$a0]
set V_small_hard [expr $PI*pow($d_small_hard,3.0)/6.0]
set S_small_hard [expr $PI*pow($d_small_hard,2.0)]
set N_oleic_small [expr $S_small_hard*$k_o*$N_oa/0.5]
set V_small_mag_core [expr $PI*pow($d_small_mag_core,3)/6.0]
set m_small [expr $Ms*$V_small_mag_core]
set M_small [expr $rop*$V_small_hard+$N_oleic_small*$mass_o]
set Vo_small_eff [expr ($PI/6.0)*(pow($d_small_hydrod,3)-pow($d_small_hard,3))]
# Effective mass density of the oleic molecules on the nanoparticle surface
set roo_small_eff [expr $N_oleic_small*$mass_o/($Vo_small_eff)]
set I_small [expr ($PI*pow($d_small_hydrod,5)/60.0)*(($rop-$roo_small_eff)*pow(1-2*$delta/$d_small_hydrod,5)+$roo_small_eff)]
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
set Vo_l_eff [expr ($PI/6.0)*(pow($d_l_hydrod,3)-pow($d_l_hard,3))]
# Effective mass density of the oleic molecules on the nanoparticle surface
set roo_l_eff [expr $N_oleic_l*$mass_o/($Vo_l_eff)]
set I_large [expr ($PI*pow($d_l_hydrod,5)/60.0)*(($rop-$roo_l_eff)*pow(1-2*$delta/$d_l_hydrod,5)+$roo_l_eff)]

# Dimensionless model
set SIGMA [expr $d_small_hydrod]
# Approximate meam mass
set M0 [expr $M_small]
set H0 [expr pow($kb*$T/(4*$PI*$mu0),0.5)/pow($SIGMA,1.5)]
set m0 [expr sqrt(4*$PI*$kb*$T*pow($SIGMA,3.0)/$mu0)]
set t0 [expr sqrt($M0*pow($SIGMA,2.0)/($kb*$T))]
set gamma0 [expr 3*$PI*$eta_car*pow($SIGMA,2)/sqrt($M0*$kb*$T)]

# Computational experiment setup
# time_step should be increased till the stability loss

# => max for pure VV:
# setmd time_step [expr 0.0005]
# => for SEMI_INTEGRATED (AI) method:
# setmd time_step [expr 0.05]

setmd time_step [expr 5E-2]
set n_part 100000

#set H_ext_Oe 500.0
set H_ext_Oe 0.0

set n_part_small [expr $n_part*(1-$PHI_l)]
#set box_l [expr pow($n_part*$PI*($PHI_l*pow($d_l_hard/$SIGMA,3.0)+(1-$PHI_l)*pow($d_small_hard/$SIGMA,3.0))/(6.0*$phi_v),1./3.)]
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
thermostat langevin $temp $gammat $gammar

set coord_shift [expr 0.1*$start_lattice_a]
set posx $coord_shift
set posy $coord_shift
set posz $coord_shift

for {set i 0} { $i < $n_part } {incr i 1} {
set rnd_val_for_p_type [expr [t_random]]
set costheta [expr 2*[expr [t_random]] - 1]
set sintheta [expr sin(acos($costheta))]
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
    set dx 0.0
    set dy 0.0
    set dz [expr 1.0 * ($m_l/$m0)]
	set massval_large [expr $M_l/$M0]
    set I_val_large [expr $I_large/($M0*pow($SIGMA,2))]
	set gamma_tr_val_large [expr $gamma0*$d_l_hydrod/$SIGMA]
	set gamma_rot_val_large [expr $gamma0*pow($d_l_hydrod/$SIGMA,3.0)/3.0]

	part $i pos $posx $posy $posz type 1 mass $massval_large rinertia $I_val_large dip $dx $dy $dz gamma $gamma_tr_val_large gamma_rot $gamma_rot_val_large
} else {
    set dx 0.0
    set dy 0.0
    set dz [expr 1.0 * ($m_small/$m0)]
	set massval_small [expr 1.0]
	set I_val_small [expr $I_small/($M0*pow($SIGMA,2))]
	set gamma_tr_val_small [expr $gamma0*1.0]
	set gamma_rot_val_small [expr $gamma0*1.0/3.0]

	part $i pos $posx $posy $posz type 0 mass $massval_small rinertia $I_val_small dip $dx $dy $dz gamma $gamma_tr_val_small gamma_rot $gamma_rot_val_small
}
}

puts "====================="
puts "Temporal parameters"
puts "====================="
puts "Small nanoparticles"
puts "d_small_hydrod=$d_small_hydrod m"
puts "tau_v_tran=[expr $massval_small/$gamma_tr_val_small] arb. units"
puts "tau_v_rot=[expr $I_val_small/$gamma_rot_val_small] arb. units"
puts ""
puts "Large nanoparticles"
puts "d_l_hydrod=$d_l_hydrod m"
puts "tau_v_tran=[expr $massval_large/$gamma_tr_val_large] arb. units"
puts "tau_v_rot=[expr $I_val_large/$gamma_rot_val_large] arb. units"

puts "[part 10]"

set sig [expr 1.0]; set cut [expr 1.12246*$sig]
set eps 1.5; set shift [expr 0.25]
inter 0 0 lennard-jones $eps $sig $cut $shift 0

set sig [expr ($d_l_hydrod/$SIGMA)]; set cut [expr 1.12246*$sig]
set eps 1.5; set shift [expr 0.25]
inter 1 1 lennard-jones $eps $sig $cut $shift 0

set sig [expr ((0.5*$d_l_hydrod+0.5*$SIGMA)/$SIGMA)]; set cut [expr 1.12246*$sig]
set eps 1.5; set shift [expr 0.25]
inter 0 1 lennard-jones $eps $sig $cut $shift 0

if { [regexp "ROTATION" [code_info]] } {
set deg_free 6
} else { set deg_free 3 }

prepare_vmd_connection "vmdfile" 10000
#anykey

for {set cap 1} {$cap < 200} {incr cap 20} {
puts "t=[setmd time] E=[analyze energy total]"
inter forcecap $cap; integrate 20
}
inter forcecap 0

for {set i 0} { $i < 5 } { incr i} {
set temp [expr [analyze energy kinetic]/(($deg_free/2.0)*$n_part)]
puts "t=[setmd time] E=[analyze energy total], T=$temp"
integrate 10
imd positions
}

#the option from the Tanygin's manunscript:
inter magnetic 1.0 bh-gpu
#inter magnetic 1.0 dds-gpu
#inter magnetic 1.0 dawaanr

set H_demag 0.0
set dipm_obs [ observable new com_dipole_moment all ]

puts "=============================="
puts "t_epoch | t_in_silico | Mz | T"
puts "=============================="

integrate 0 recalc_forces
for {set i 0} { $i < 1E20 } { incr i} {
set temp [expr [analyze energy kinetic]/(($deg_free/2.0)*$n_part)]
#puts "t=[setmd time] E=[analyze energy total], T=$temp"
#puts "t=[setmd time] dipm=[observable $dipm_obs print]"
set H_demag 0.0
#set H_demag [expr -[lindex [ observable $dipm_obs print ] 2]/$V_carr ]
set Mz [expr [lindex [ observable $dipm_obs print ] 2]/$V_carr]
set H_mag [expr (($H_ext_Oe*$Oe_to_Am)/$H0+$H_demag)]
constraint ext_magn_field 0 0 $H_mag

integrate 10 reuse_forces
imd positions

puts "[clock seconds],[setmd time],$Mz,$temp"

}

imd disconnect
