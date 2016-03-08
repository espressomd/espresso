# Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
# 

set tcl_precision 15

# location of gibbs simulation tcl scripts
source "scripts/gibbs_ghmc.tcl"

# intialization of the simulation box
proc init {id param} {

  global box_l Vbox
  global n_md
  global lj_epsilon lj_sigma lj_cut lj_shift lj_cutoff
  global results
  global temperature
  global systype

  # setup initial box parameters
  set Vbox [lindex $param 0]
  set box_l [expr pow($Vbox,1.0/3.0) * $lj_sigma]
  set n_particle [lindex $param 1]
  puts "initial box volume is $Vbox\ninitial number of particles is $n_particle"

  # give Espresso some parameters
  setmd time 0
  setmd time_step 0.01
  setmd skin 0.4
  setmd box_l $box_l $box_l $box_l
  setmd periodic 1 1 1

  # setup LJ interaction
  inter 0 0 lennard-jones $lj_epsilon $lj_sigma $lj_cutoff $lj_shift 0

  # put the  particles into the box 
  puts "\nput the particles in the box..."  
  part deleteall
  for {set i 0} {$i < $n_particle} {incr i} {
      set posx [expr $box_l*[t_random]]	
      set posy [expr $box_l*[t_random]]	
      set posz [expr $box_l*[t_random]] 	 
      part $i pos $posx $posy $posz type 0
  }      

  # warmup settings
  set warmup_lj_cap 100
  set warmup_min_dist 0.9
  set int_step_warmup 1000
  set min_dist 0
  # set up warmup thermostat
  thermostat langevin $temperature 1.0
  # do the warm up
  puts "\ndoing the warmup..."
  while {$min_dist < $warmup_min_dist} {  
    inter forcecap $warmup_lj_cap
    integrate $int_step_warmup
    set min_dist [analyze mindist] 
  }
  puts "minimanl distance is $min_dist"
  # disable force cap   
  inter forcecap 0
  # disable langevin thermostat
  thermostat off

  # enable the gibbs simulation ghmc thermostat
  thermostat ghmc $temperature $n_md 0

  puts "\nequilibration..."

}

# calculate potential energy + tail correction
proc potential_energy {N} {
  global Vbox
  global lj_epsilon lj_sigma lj_cutoff lj_shift

  if {$N > 0} {
    # potential energy
    set eng [expr [analyze energy total] - [analyze energy kinetic]]

    set dens [expr double($N)/$Vbox]

    # tail correction
    set tail [expr $N*8./3.*[PI]*$dens*$lj_epsilon*pow($lj_sigma, 3)*\
          (1./3.*pow($lj_sigma/$lj_cutoff, 9) - pow($lj_sigma/$lj_cutoff, 3))]
    
    # shift correction
    set shift [expr -$N*8./3.*[PI]*$dens*$lj_epsilon*pow($lj_cutoff, 3)*\
            $lj_shift]
    
    return [expr $eng + $tail + $shift]
  } else {
    return 0
  }
}

# volume exchange between boxes
proc swap_vol {id mvtype V_new} {

  global box_l Vbox
  global eo en mvcancel

  #update energy
  if {$mvtype} {
    if {$mvcancel == 0} {
      set eo $en
    } elseif {$mvcancel != 1} {
      set eo [potential_energy [setmd n_part]]
    }
  } else {
      set mvcancel 1
  }

  #change box volume
  if { $V_new != $Vbox } {

    set Vbox $V_new
    set box_l [expr pow($Vbox,1.0/3.0)]
    change_volume $box_l xyz
  }
  
  # update energy  
  if {$mvtype} {
    set mvcancel 0
    set en [potential_energy [setmd n_part]]
    return "$eo $en"
  } else {
    return "0 0"
  }

}

# particle exchange between boxes
proc swap_part {id mvtype part_id} {

  global cp ncp
  global Vbox
  global temperature
  global eo en mvcancel

  # remove particle of rejected swap
  if {$mvtype == 0} {
    delete_particle $part_id
    if {$mvcancel == 1} {set mvcancel 0} else {set mvcancel 1}
    return "0 0 0"
  }

  # update energy
  if {$mvcancel == 0} {
    set eo $en
  } elseif {$mvcancel != 1} {
    set eo [potential_energy [setmd n_part]]
  }

  if {$mvtype == 1} {
    # add particle
    set part_id [create_particle]
    set en [potential_energy [setmd n_part]]
    # measure chemical potential by widom insertion method
    set cp [expr $cp + $Vbox*exp(-($en-$eo)/$temperature)/double([setmd n_part])]; incr ncp
    set mvcancel 0
    } else  {
      # remove particle
      # check if box is empty
      if {[setmd n_part] < 1} { return "$eo 100 0"}

      # TRICK: we do not really delete the particle, just give it another
      # type, so that it does not interact. This is easier than removing it
      # and potentially restoring it, if the move was rejected.

      # randomly pick a particle, and "remove"
      set part_id [t_random int [expr [setmd max_part]] + 1]
      catch {
        part $part_id type 1
      } errMsg
      if { $errMsg != "" } {
        puts "$errMsg, pid: $part_id, max_part: [setmd max_part]"
        set part_id [t_random int [expr [setmd max_part]] + 1]
        part $part_id type 1
      }

      # CAVEAT! Since we make the particle disappear only through its type,
      # it still counts to n_part, but we must not do the tail correction for it!
      set en [potential_energy [expr [setmd n_part] - 1]]
      # fix the particles type again
      eval part $part_id type 0
      set mvcancel 1
    }
    return "$eo $en $part_id"
}

proc create_particle {} {

    global box_l
    # add particle after last known particle
    set pid [expr [setmd max_part] + 1]
    set posx [expr $box_l*[t_random]]
    set posy [expr $box_l*[t_random]]
    set posz [expr $box_l*[t_random]]

    part $pid pos $posx $posy $posz type 0

    return $pid

}

proc delete_particle {pid} {

  if {$pid == [setmd max_part]} {
    # last particle, just delete
    part $pid delete
  } {
    # otherwise, shift down last particle
    # this avoids that the particle identities get excessive
    eval part $pid pos [part [setmd max_part] pr pos] v [part [setmd max_part] pr v ]
    part [setmd max_part] delete
  }
}

# particle moves within boxes
proc move_parts {id time_step} {

  global n_md
  global eo en mvcancel

  setmd time_step $time_step
  catch {
    integrate $n_md
  } errMsg
  if { $errMsg != "" } {
    puts "Time step too big: dt = $time_step, integration failed."
    return "0 100"
  }

  set stats [ghmc statistics]
  set res [lrange $stats 0 1]
  set mvcancel 2

  return $res
}

# sampling of obsevables
proc sample {id param} {

  global Vbox
  global results
  global cycle n_sample_cycles
  global rho nrho

  puts $results [format "%8.2f %8.5f %5d %8.5f" \
    [setmd time] $Vbox [setmd n_part] [expr double([setmd n_part])/$Vbox]]

  incr cycle -1 
  if {$cycle < $n_sample_cycles} {
    incr nrho
    set rho [expr $rho+double([setmd n_part])/$Vbox]
    # console information
    puts -nonewline [format "\r%5d average density = %8.5f" \
       [expr $nrho] [expr $rho/double($nrho)]]
    flush stdout
  } elseif {$cycle == $n_sample_cycles} {
    puts "\nsampling..."
  }

}

#####################################################
# gibbs ensemble simulation global settings #########
#####################################################

# setup system type (master / slave) and host ip
# the system name is the only script argument,
# name starting with "m" will be master
set sysname [lindex $argv 0]
set systype [string range $sysname 0 0]
set host 127.0.0.1

# set gibbs MC cycles
# total cycles
set n_total_cycles 10000
# sampling cycles
set n_sample_cycles 5000
# cycles counter
set cycle $n_total_cycles
# MC parameters adjusting rate
set adjust_rate 250

# set simulation temperature
set temperature  1.20

# init box variables
variable Vbox
variable box_l

# init chemical potential
set cp  0
set ncp 0

# init density
set rho  0
set nrho 0
      
# lj potential parameters
set lj_epsilon 1.0
set lj_sigma   1.0
set lj_cutoff  3.0
set lj_shift   0.0

# number of particles moves integration steps
set n_md 10

# global energy calulation variables
variable eo
variable en
set mvcancel 2

#set initial box densities
set rho_init 0.3

# set inital box particle number
set n_part_init 256

#set initial box volume
set V_init [expr $n_part_init/$rho_init]
#set maximal logarithm of volume change
set Vmax [expr log(2*$V_init)*0.015]

# inital values struct for gibbs_ghmc script
set init_values [list [list $V_init $n_part_init] [list $V_init $n_part_init]]

# init random number generator
t_random seed [clock clicks]

# setup output file
set filename "lj-gibbs_${temperature}_${sysname}.data"
set results [open $filename w];
puts $results "\# time V N rho"

puts "\nGibbs ensemble simulation on system $sysname at temperature $temperature"
set start_time [clock seconds] 

if {$systype == "m"} {

  gibbs_ghmc::main -values $init_values -rounds $n_total_cycles -sample_rounds $n_sample_cycles -resrate $adjust_rate \
    -temp $temperature -vmax $Vmax -time_step 0.01 -max_time_step 0.02 \
    -nmove 16 -nvol 4 -npart 80 -mc_cycle 100 -adjust_dt 0.5 -adjust_vol 0.5 \
    -init init -swap_vol swap_vol -swap_part swap_part -move_part move_parts -sample sample -info comm -port 15000

} else {

  gibbs_ghmc::main -connect $host \
    -init init -swap_vol swap_vol -swap_part swap_part -move_part move_parts -sample sample -port 15000

}

close $results

puts "\n\ndone!"
puts "on system $sysname at temperature $temperature:\naverage density after $n_sample_cycles cycles was [expr $rho/double($n_sample_cycles)]"
puts "with average excess chemical potential: [expr -log($cp/double($ncp))]"
puts "\nelapsed simulation time was [expr [clock seconds]-$start_time] seconds"
