# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
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


# check the charge-charge P3M  algorithm
source "tests_common.tcl"

require_feature "REACTION_ENSEMBLE"

puts "---------------------------------------------------------------"
puts "- Testcase reaction_ensemble.tcl running on [format %02d [setmd n_nodes]] nodes"
puts "---------------------------------------------------------------"

proc ideal_degree_of_association { pK_a pH } {
    return [expr 1-1.0/(1+pow(10,$pK_a - $pH))]
}

if { [catch {
    puts "Tests for reaction ensemble in the ideal case (no interactions)"
    
    set type_A 0
    set type_HA 1
    set type_H 2
    set N0 200
    set c0 0.0001
    set box_l [expr pow($N0/$c0,1.0/3)]
    setmd box_l $box_l $box_l $box_l

    for { set i 0 } { $i <= $N0 } { incr i } {
	    part $i pos [expr $box_l/2] [expr $box_l/2] [expr $box_l*[t_random]] type $type_HA
    }
    
    set temperature 1.0
    set standard_pressure_in_simulation_units 0.00108
    set lj_sig 1.0
    
    set K_HA_diss 0.5; #could be in this test for example anywhere in the range 0.000001 ... 9
    #initialize reaction ensemble with first reaction
    reaction_ensemble add equilibrium_constant $K_HA_diss educt_types $type_HA educt_coefficients 1 product_types $type_A $type_H product_coefficients 1 1
    reaction_ensemble add equilibrium_constant [expr 1.0/$K_HA_diss] educt_types $type_A $type_H educt_coefficients 1 1 product_types $type_HA product_coefficients 1
    reaction_ensemble standard_pressure_in_simulation_units $standard_pressure_in_simulation_units;
    reaction_ensemble set_default_charge_of_type $type_A -1
    reaction_ensemble set_default_charge_of_type $type_HA 0
    reaction_ensemble set_default_charge_of_type $type_H +1
    reaction_ensemble temperature_reaction_ensemble $temperature
    reaction_ensemble exclusion_radius $lj_sig
    puts [reaction_ensemble ]
  
    #chemical warmup in order to get to chemical equilibrium before starting to calculate the observable "degree of association"
    for {set i 0} {$i <[expr 10*$N0] } {incr i} {
        reaction_ensemble do
    }  
    
    
    set volume [expr pow($box_l,3)]
    
    set average_NH 0.0
    set average_degree_of_association 0.0
    set num_samples 1000
    for {set i 0} {$i <$num_samples} {incr i} {
        reaction_ensemble do
        set average_NH [expr $average_NH + [part gc num $type_H]]
        set average_degree_of_association [expr $average_degree_of_association+ [part gc num $type_HA]/double($N0)]
    }
    set average_NH [expr $average_NH/$num_samples]
    set average_degree_of_association [expr $average_degree_of_association/$num_samples]
    set pH [expr -log10($average_NH/$volume)]
    set K_apparent_HA_diss [expr $K_HA_diss*$standard_pressure_in_simulation_units/$temperature]
    set pK_a [expr -log10($K_apparent_HA_diss)]
    set real_error_in_degree_of_association [expr abs($average_degree_of_association-[ideal_degree_of_association $pK_a $pH])/[ideal_degree_of_association $pK_a $pH]]
    
    if { $real_error_in_degree_of_association > 0.07  } {
      error "reaction ensemble: relative error in the degreee of association too large"
    }
    
   #end this part of the p3m-checks by cleaning the system .... 
   part deleteall
   inter coulomb 0.0

} res ] } {
    error_exit $res
}

exit 0
