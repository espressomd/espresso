#############################################################
#
#  tcp/ip socket based gibbs Hybrid MC interface
#
#############################################################
#
# Copyright (C) 2010,2011,2013 The ESPResSo project
# Copyright (C) 2007,2008,2009,2010,2011 Axel Arnold
# Copyright (C) 2011 Sela Samin
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

# TO DO: ##################################


# how to sample an arbitrary probability distribution ? (for various biasing)

# profile and then ################################
# write documentation

namespace eval gibbs_ghmc {
    # common variables
    #################################

    variable rounds 1
    variable values ""
    variable load 1
    variable port 12000
    variable master ""
    variable beta 0
    # variables for # of attempts of volume and particles swapping / particles move 
    variable nvol 1
    variable npart 1
		variable nmove 1
		#variables for total # of attempts of in one MC cycle
		variable mc_cycle 0
    #variable for total # of equllibration MC rounds
    variable sample_rounds 0
		#variables for whthear to adjust volume and time steps
		variable adjust_vol 0
		variable adjust_dt 0
		# max step in log V
    variable Vmax 0
		#MD stepsize array and dt_max
		variable md_dt 0.01
		variable dt_max 0.1

    # function calculating the swap volume probabilities
    variable swap_vol ""
    # function calculating the swap particle probabilities
    variable swap_part ""
    # name of the function moving the particles
    variable move_part ""
    # name of the function running the simulation
    variable sample ""
    # name of the initialization function per run
    variable init ""
    # data common to all processors. Can only bet set from the master node
    variable shareddata ""
    # information level
    # - none: total silence
    # - comm: communication information on stderr (connection of clients)
    # - all:  information on swap tries on stdout
    variable info "all"
    variable reset_acceptance_rate 100

    # name space for the routines used by the slave processes
    #########################################################
    namespace eval slave {
	variable sock

	####
	proc init {id param} {
	    if {$gibbs_ghmc::init != ""} {
		eval ::$gibbs_ghmc::init $id $param
	    }
	}
	####
	proc sample {id param} {
	    eval ::$gibbs_ghmc::sample $id $param
	}
	####
	proc swap_vol {id param1 param2} {
	    return [eval ::$gibbs_ghmc::swap_vol $id $param1 $param2]
	}
	####
	proc swap_part {id param1 param2} {
	    return [eval ::$gibbs_ghmc::swap_part $id $param1 $param2]
	}
	####
	proc move_part {id param1} {
	    return [eval ::$gibbs_ghmc::move_part $id $param1]
	}
	####
	proc connect {} {
	    variable sock
			set ntry 10
			while {$ntry > 0} {
				catch {
					set sock [socket $gibbs_ghmc::master $gibbs_ghmc::port]
					fconfigure $sock -buffering line
					set msg [gets $sock]
					if {$msg != "gibbs_ghmc welcome"} { error "not the right server for me" }
					puts $sock "gibbs_ghmc load $gibbs_ghmc::load"
					} errMsg
				if { $errMsg == "" } {
					break
				} elseif { $errMsg != "couldn't open socket: connection refused"} {
					error $errMsg
				}
				if {$gibbs_ghmc::info != "none"} {puts "\n$errMsg, trying again after 1s"}
				after 1000
				incr ntry -1
			}
}
	####
	proc disconnect {} {
	    variable sock
	    close $sock
	}
	####
	proc main {} {
	    variable sock
	    connect
	    while { [set msg [gets $sock]] != "exit" } {
		switch -- [lindex $msg 0] {
		    "system" {
			init [lindex $msg 1] [list [lrange $msg 2 end]]
			puts $sock "system [lindex $msg 1] ok"
		    }
		    "shareddata" {
			gibbs_ghmc::set_shareddata [lrange $msg 1 end]
		    }
		    "sample" { sample [lindex $msg 1] [list [lrange $msg 2 end]] }
 		    "swapprob_vol" {
			puts $sock [swap_vol [lindex $msg 1] [lindex $msg 2] [lindex $msg 3]]
		    }
 		    "swapprob_part" {
			puts $sock [swap_part [lindex $msg 1] [lindex $msg 2] [lindex $msg 3]]
		    }
				"moveprob_part" {
			puts $sock [move_part [lindex $msg 1] [lindex $msg 2]]
		    }
		    default { error "illegal request from master \"[lindex $msg 0]\"" }
		}
	    }
	    disconnect
	}
    }

    # name space for the routines used by the master process
    ########################################################
    namespace eval master {
	# list of all the channels and ids of the clients, in the same order as the values
	variable clients
	# a unique list of just the slave nodes
	variable slaves
	# whether the shareddata block has to be retransmitted
	variable shareddata_changed
	# tries and accepted tries of swaps in an array
	# acceptance(Ta,Tb)={total_tries successful_tries}
	variable acceptance
	variable acceptance_prev

	####
	proc get_acceptance_rate {id1 id2 rtype} {
	    variable acceptance	    
	if {[array names acceptance $id1,$id2,$rtype ] != ""} {
	    foreach {tries succ} $acceptance($id1,$id2,$rtype) { break }
	    return [expr double($succ)/$tries]
		} else { 
			return 0
		}
	}
	####
	proc reset_acceptance_rates {} {
	    variable acceptance
	    array unset acceptance
	}
	####
	proc adjust_acceptance_rates {id1 id2 rtype mv_param} {
	    variable acceptance
			variable acceptance_prev

			if {[array names acceptance $id1,$id2,$rtype] != ""} {
			
			#set acceptance rate target
			if {$rtype == 0} {
				set succ_target $gibbs_ghmc::adjust_vol
			} elseif {$rtype == 2} {
				set succ_target $gibbs_ghmc::adjust_dt
			} else {
				set succ_target 0.5
			}

			foreach {tries succ} $acceptance($id1,$id2,$rtype) { break }

			if {[array names acceptance_prev $id1,$id2,$rtype ] != ""} {
				foreach {triesp succp} $acceptance_prev($id1,$id2,$rtype) { break }
			} else {
				set triesp 0; set succp 0
			}

      if { $tries == 0 || $triesp >= $tries} {
				set acceptance_prev($id1,$id2,$rtype) "$tries $succ"
      } else {
         set frac [expr double($succ-$succp)/double($tries-$triesp)]
         set mvp_old $mv_param
         set mv_param [expr $mv_param*abs($frac/$succ_target)]
					#limit the change:
         if {[expr $mv_param/$mvp_old]>1.5} { set mv_param [expr $mvp_old*1.5]}
         if {[expr $mv_param/$mvp_old]<0.5} { set mv_param [expr $mvp_old*0.5]}
				puts "gibbs_ghmc: MC move parameter set to : $mv_param (old : $mvp_old), Frac. acc.: $frac, attempts: [expr $tries-$triesp] succes: [expr $succ-$succp]"
				#store current try data for next use
        set acceptance_prev($id1,$id2,$rtype) "$tries $succ"
      }
		}
			return $mv_param
	}
	####
	proc slave_connect {channel clientaddr clientport} {
	    variable clients
	    if {$gibbs_ghmc::info != "none"} {
		puts stderr "connection from $clientaddr $clientport"
	    }
	    fconfigure $channel -buffering line
	    set n_systems [llength $gibbs_ghmc::values]
	    set n_clients [expr [llength $clients]/2]
	    # check for overoccupation
	    if {$n_clients >= $n_systems} {
		puts $channel "go away"
		close $channel
		return
	    }
	    puts $channel "gibbs_ghmc welcome"
	    set response [gets $channel]
	    if {![string match "gibbs_ghmc load *" $response]} {
		error "wrong response from client"
	    }
	    set client_load [lindex $response 2]
	    if {[expr $n_clients + $client_load] > $n_systems} {
		set client_load [expr $n_systems - $n_clients]
	    }
	    # register client
	    for {set i 0} {$i < $client_load} {incr i} { lappend clients $channel $i }
	}
	####
	proc wait_for_clients {} {
	    variable clients

	    set n_systems [llength $gibbs_ghmc::values]
	    while { 1 } {
		set n_clients [expr [llength $clients]/2]
		if {$n_clients >= $n_systems} { break }
		if {$gibbs_ghmc::info != "none"} {
		    puts stderr "waiting for [expr $n_systems - $n_clients] more clients"
		}
		vwait gibbs_ghmc::master::clients
	    }
	    if {$gibbs_ghmc::info != "none"} {
		puts stderr "all there, let's go"
	    }
	}
	####
	proc init_clients {} {
	    variable clients

	    # send initialization command and initialize MD time step array
			set init_dt $gibbs_ghmc::md_dt; unset gibbs_ghmc::md_dt
	    foreach param $gibbs_ghmc::values {channel id} $clients {
				if {$channel != "myself"} {
						puts $channel "system $id $param"
				} {
						lappend myparam $id $param
				}
				set gibbs_ghmc::md_dt($channel,$id) $init_dt
	    }

	    # now, initialize myself
	    foreach {id param} $myparam { gibbs_ghmc::slave::init $id [list $param] }
	    # and wait for all clients to report ready
	    foreach {channel id} $clients {
		if {$channel != "myself"} {
		    set response [gets $channel]
		    if {$response != "system $id ok"} {
			error "client connection error, wrong response"
		    }
		    puts $channel "shareddata $gibbs_ghmc::shareddata"
		}
	    }
	    set shareddata_changed 0
	}
	####
	proc sample_clients {} {
	    variable clients
	    variable shareddata_changed

	    # broadcast current parameters
	    foreach param $gibbs_ghmc::values {channel id} $clients {
		if {$channel == "myself"} {
		    lappend myparam $id $param
		} {
		    if {$shareddata_changed} {
			puts $channel "shareddata $gibbs_ghmc::shareddata"
		    }
		    puts $channel "sample $id $param"
		}
	    }
	    # sample myself
	    foreach {id param} $myparam { gibbs_ghmc::slave::sample $id [list $param] }
	    set shareddata_changed 0
	}
	####
	proc ask_swapprob_vol {channel id param1 param2} {
	    if {$channel == "myself"} {
		set res [gibbs_ghmc::slave::swap_vol $id $param1 $param2]
	    } {
	    	puts $channel "swapprob_vol $id $param1 $param2"
	    	set res [gets $channel]
	    }
	    if {[llength $res] != 2} {
		error "swap volume format wrong, \"$res\" is not a list of two values."
	    }
	    return $res
	}
	####
	proc swap_clients_vol {odd} {
	    variable clients
	    variable acceptance

	    if {$gibbs_ghmc::info == "all"} {
		puts "gibbs_ghmc swapping volume start"
	    }
	    # odd-even rule. We either try to exchange 1-2
	    # or we simply keep 1 and start with 2-3
	    if {$odd} {
		# keep first client
		set v [lrange $gibbs_ghmc::values 1 end]
		set c [lrange $clients 2 end]
	    } {
		# do not keep it
		set v $gibbs_ghmc::values
		set c $clients
	    }

 	    foreach {param1 param2} $v {channel1 id1 channel2 id2} $c {
		# single parameter left over, cannot swap
		if {$param2 == ""} { break }
		
		#random walk in ln V1/V2
		set N1 [lindex $param1 1]
		set N2 [lindex $param2 1]
		set V1_old [lindex $param1 0]
		set V2_old [lindex $param2 0]
		set Vtot [expr $V1_old+$V2_old]

		set lnV_new [expr log($V1_old/$V2_old)+([t_random]-0.5)*$gibbs_ghmc::Vmax]
		set V1_new [expr $Vtot*exp($lnV_new)/(1+exp($lnV_new))]
		set V2_new [expr $Vtot-$V1_new]    

		# ask the clients for swapping volumes
		foreach {E1o E1n} [ask_swapprob_vol $channel1 $id1 1 $V1_new] { break }
		foreach {E2o E2n} [ask_swapprob_vol $channel2 $id2 1 $V2_new] { break }

		# check whether to actually swap by Metropolis
		set exp [expr -$gibbs_ghmc::beta*($E1n - $E1o + $E2n - $E2o)+($N1+1)*log($V1_new/$V1_old)+($N2+1)*log($V2_new/$V2_old)]
		if {$exp > 0} { set boltzmann 1	} elseif {$exp < -30} { set boltzmann 0	} {
		    set boltzmann [expr exp($exp)]
		}
		if {$gibbs_ghmc::info == "all"} {
		    puts "gibbs_ghmc tryswap volume $V1_old to $V1_new and $V2_old to $V2_new : energies $E1o $E1n $E2o $E2n -> $exp -> $boltzmann"
		}

		set pairid "$id1,$id2,0"
		if {[array names acceptance $pairid] == ""} {
		    set tries 1
		    set succ 0
		} {
		    foreach {tries succ} $acceptance($pairid) { break }
		    incr tries
		}

		set myrand [t_random]
		if {[expr $myrand < $boltzmann]} {
		    if {$gibbs_ghmc::info == "all"} {
					puts "gibbs_ghmc accept swap volume to: $V1_new $V2_new"
		    }
		    set gibbs_ghmc::values [lreplace $gibbs_ghmc::values 0 0 [list $V1_new $N1]]
				set gibbs_ghmc::values [lreplace $gibbs_ghmc::values 1 1 [list $V2_new $N2]]
		    incr succ
		} else {
		    if {$gibbs_ghmc::info == "all"} {
					puts "gibbs_ghmc undo swap volume to: $V1_old $V2_old"
		    }
		    # ask the clients for swapping volumes back
		  foreach {E1o E1n} [ask_swapprob_vol $channel1 $id1 0 $V1_old] { break }
		  foreach {E2o E2n} [ask_swapprob_vol $channel2 $id2 0 $V2_old] { break }
		}
		
		set acceptance($pairid) "$tries $succ"
	    }

	    if {$gibbs_ghmc::info == "all"} {
		puts "gibbs_ghmc swapping volume done"
	    }
	}
	####
	proc ask_swapprob_part {channel id param1 param2} {
	    if {$channel == "myself"} {
		set res [gibbs_ghmc::slave::swap_part $id $param1 $param2]
	    } {
	    	puts $channel "swapprob_part $id $param1 $param2"
	    	set res [gets $channel]
	    }
	    if {[llength $res] != 3} {
		error "swap part format wrong, \"$res\" is not a list of three elements."
	    }
	    return $res
	}
	####
	proc ask_moveprob_part {channel id param1} {
	    if {$channel == "myself"} {
		set res [gibbs_ghmc::slave::move_part $id $param1]
	    } {
	    	puts $channel "moveprob_part $id $param1"
	    	set res [gets $channel]
	    }
	    if {[llength $res] != 2} {
		error "move part format wrong, \"$res\" is not a list of two element."
	    }
	    return $res
	}
	####
	proc swap_clients_part {odd} {
	    variable clients
	    variable acceptance

	    if {$gibbs_ghmc::info == "all"} {
		puts "gibbs_ghmc swapping part start"
	    }
	    # odd-even rule. We either try to exchange 1-2
	    # or we simply keep 1 and start with 2-3
	    if {$odd} {
		# keep first client
		set v [lrange $gibbs_ghmc::values 1 end]
		set c [lrange $clients 2 end]
	    } {
		# do not keep it
		set v $gibbs_ghmc::values
		set c $clients
	    }

 	    foreach {param1 param2} $v {channel1 id1 channel2 id2} $c {
		# single parameter left over, cannot swap
		if {$param2 == ""} { break }
		
		
		set V1 [lindex $param1 0]
		set V2 [lindex $param2 0]
		set N1_old [lindex $param1 1]
		set N2_old [lindex $param2 1]

		if {[t_random]<0.5} {
			set mv1 1
			set mv2 -1
			set Nin $N1_old
			set Nout $N2_old
			set Vin $V1
			set Vout $V2
		} else {
			set mv1 -1
			set mv2 1
			set Nin $N2_old
			set Nout $N1_old
			set Vin $V2
			set Vout $V1
		}

		# ask the clients for swapping particles
		foreach {E1o E1n part_id1} [ask_swapprob_part $channel1 $id1 $mv1 0] { break }
		foreach {E2o E2n part_id2} [ask_swapprob_part $channel2 $id2 $mv2 0] { break }

		# check whether to actually swap by Metropolis
		set exp [expr -$gibbs_ghmc::beta*($E1n - $E1o + $E2n - $E2o)+log($Nout*$Vin/(($Nin+1)*$Vout))]
		if {$exp > 0} { set boltzmann 1	} elseif {$exp < -30} { set boltzmann 0	} {
		    set boltzmann [expr exp($exp)]
		}
		if {$gibbs_ghmc::info == "all"} {
		    puts "gibbs_ghmc tryswap part $N1_old $N2_old -> $mv1 $mv2 : energies $E1o $E1n $E2o $E2n -> $exp -> $boltzmann"
		}

		set pairid "$id1,$id2,1"
		if {[array names acceptance $pairid] == ""} {
		    set tries 1
		    set succ 0
		} {
		    foreach {tries succ} $acceptance($pairid) { break }
		    incr tries
		}

		set myrand [t_random]
		if {[expr $myrand < $boltzmann]} {
		    if {$gibbs_ghmc::info == "all"} {
					puts "gibbs_ghmc accept swap part to: [expr $N1_old+$mv1] [expr $N2_old+$mv2]"
		    }
		    set gibbs_ghmc::values [lreplace $gibbs_ghmc::values 0 0 [list $V1 [expr $N1_old+$mv1]]]
				set gibbs_ghmc::values [lreplace $gibbs_ghmc::values 1 1 [list $V2 [expr $N2_old+$mv2]]]
				# ask the clients to really deleting chosen particle
				if {$mv1 == -1} {foreach {E1o E1n part_id1} [ask_swapprob_part $channel1 $id1 0 $part_id1] { break }}
				if {$mv2 == -1} {foreach {E2o E2n part_id2} [ask_swapprob_part $channel2 $id2 0 $part_id2] { break }}
		    incr succ
		} else {
		    if {$gibbs_ghmc::info == "all"} {
					puts "gibbs_ghmc undo swap part to: [expr $N1_old+$mv1] [expr $N2_old+$mv2]"
		    }
		    # ask the clients to deleting newly created particle
				if {$mv1 == 1} {foreach {E1o E1n part_id1} [ask_swapprob_part $channel1 $id1 0 $part_id1] { break }}
				if {$mv2 == 1} {foreach {E2o E2n part_id2} [ask_swapprob_part $channel2 $id2 0 $part_id2] { break }}
		}
		
		set acceptance($pairid) "$tries $succ"
	    }

	    if {$gibbs_ghmc::info == "all"} {
		puts "gibbs_ghmc swapping part done"
	    }
	}
	####
	proc move_clients_part {odd} {
	    variable clients
	    variable acceptance

	    if {$gibbs_ghmc::info == "all"} {
		puts "gibbs_ghmc move part start"
	    }
	    # odd-even rule. We either try to exchange 1-2
	    # or we simply keep 1 and start with 2-3
	    if {$odd} {
		# keep first client
		set v [lrange $gibbs_ghmc::values 1 end]
		set c [lrange $clients 2 end]

	    } {
		# do not keep it
		set v $gibbs_ghmc::values
		set c $clients
	    }

 	    foreach {param1 param2} $v {channel1 id1 channel2 id2} $c {
		# single parameter left over, cannot swap
		if {$param2 == ""} { break }

		set N1 [lindex $param1 1]
		set N2 [lindex $param2 1]

# 		if {[t_random]<0.5} {}
		if {[t_random]<[expr double($N1)/double($N1+$N2)]} {
			set boxchannel $channel1
			set boxid $id1
		} else {
			set boxchannel $channel2
			set boxid $id2
		}


		# ask the client for moving particles
		foreach {matt macc} [ask_moveprob_part $boxchannel $boxid $gibbs_ghmc::md_dt($boxchannel,$boxid)] { break }

			# update move statistics
			set pairid "$boxchannel,$boxid,2"
			if {[array names acceptance $pairid] == ""} {
					set tries 1
					set succ 0
			} {
					foreach {tries succ} $acceptance($pairid) { break }
					set tries [expr $tries+$matt]
					set succ [expr $succ+$macc]
			}
			
			set acceptance($pairid) "$tries $succ"
	    }

	    if {$gibbs_ghmc::info == "all"} {
				puts "gibbs_ghmc moving part done"
	    }
	}
	####
	proc release_clients {} {
	    variable slaves
	    foreach channel $slaves { puts $channel "exit"; close $channel }
	}
	####
	proc send {sparam string} {
	    variable clients
	    if {$sparam == "-all"} {
		foreach channel $slaves {
		    puts $channel "$string"
		}
	    } {
		foreach param $gibbs_ghmc::values {channel id} $clients {
		    if {$sparam == $param} { puts $channel "$id $string" }
		}
	    }
	}
	####
	proc recv {sparam} {
	    variable clients
	    foreach param $gibbs_ghmc::values {channel id} $clients {
		if {$sparam == $param} {
		    set res [gets $channel]
		    set sid [lindex $res 0]
		    if {$sid != $id} { error "send/recv order messed up, $sid != $id" }
		    lappend results $res
		}
	    }
	    if {[llength $results] == 1} { set results [lrange [lindex $results 1] 1 end] }
	    return $results
	}
	####
	proc main {} {
	    variable clients
	    variable slaves
	    variable acceptance
	    variable shareddata_changed

	    set shareddata_changed 1
	    set clients ""
	    set slaves  ""

	    # connect myself
	    for {set i 0} {$i < $gibbs_ghmc::load} {incr i} { lappend clients "myself" $i }

	    set server_socket [socket -server slave_connect $gibbs_ghmc::port]

	    wait_for_clients

	    # setup the slaves
	    set slaves ""
	    foreach {channel id} $clients {
				if {$channel != "myself"} { lappend slaves $channel }
	    }
	    set slaves [lsort -unique $slaves]

	    init_clients

	    reset_acceptance_rates
			
	    # main loop
			while { $gibbs_ghmc::rounds > 0 } {

				#if {$gibbs_ghmc::info == "all"} {}
						if {[expr $gibbs_ghmc::rounds % $gibbs_ghmc::reset_acceptance_rate] == 0 && [expr $gibbs_ghmc::rounds] > $gibbs_ghmc::sample_rounds} {
							if {[array names acceptance] != ""} {
								puts "gibbs_ghmc listing acceptance_rates:"
								foreach {channel1 id1 channel2 id2} $clients {
									if {$id1 != "" && $id2 != ""} {
										puts "gibbs_ghmc: volume acceptance rate $id1 $id2 [get_acceptance_rate $id1 $id2 0]"
										puts "gibbs_ghmc: particle acceptance rate $id1 $id2 [get_acceptance_rate $id1 $id2 1]"
										puts "gibbs_ghmc: move acceptance rate $channel1 $id1 [get_acceptance_rate $channel1 $id1 2]"
										puts "gibbs_ghmc: move acceptance rate $channel2 $id2 [get_acceptance_rate $channel2 $id2 2]"
										if {$gibbs_ghmc::adjust_vol} {set gibbs_ghmc::Vmax [adjust_acceptance_rates $id1 $id2 0 $gibbs_ghmc::Vmax]}
										if {$gibbs_ghmc::adjust_dt} {
											set gibbs_ghmc::md_dt($channel1,$id1) [expr min([adjust_acceptance_rates $channel1 $id1 2 $gibbs_ghmc::md_dt($channel1,$id1)],$gibbs_ghmc::dt_max)]
											set gibbs_ghmc::md_dt($channel2,$id2) [expr min([adjust_acceptance_rates $channel2 $id2 2 $gibbs_ghmc::md_dt($channel2,$id2)],$gibbs_ghmc::dt_max)]
										}
									}
								}
								#reset_acceptance_rates
							}
						}

				sample_clients
				
				for {set k 0} {$k<$gibbs_ghmc::mc_cycle} {incr k} { 
					set myrand [expr [t_random]*($gibbs_ghmc::nmove+$gibbs_ghmc::nvol+$gibbs_ghmc::npart)]
					if {$myrand<$gibbs_ghmc::nmove} {
						move_clients_part [expr ($gibbs_ghmc::rounds) == 0]
					} elseif {$myrand<[expr $gibbs_ghmc::nmove+$gibbs_ghmc::nvol]} {
						swap_clients_vol [expr ($gibbs_ghmc::rounds) == 0]
					} else {
						swap_clients_part [expr ($gibbs_ghmc::rounds) == 0]
					}
				}

				incr gibbs_ghmc::rounds -1
			}

	    release_clients
	    close $server_socket
	}
    
}

# gibbs_ghmc public methods
################################################

proc set_shareddata {data} {
	variable shareddata
	set shareddata $data
	set master::shareddata_changed 1
}
###
proc main {args} {
	variable rounds; variable master; variable port
	variable values; variable load;	variable info
	variable init; variable swap_vol; variable swap_part; variable move_part; variable sample
	variable nvol; variable npart; variable nmove;
	variable beta; variable Vmax; variable md_dt; variable dt_max
	variable mc_cycle; variable adjust_vol; variable adjust_dt; variable sample_rounds
	
	set init ""
	set values ""
	set master ""

	while {$args != ""} {
	    set arg [lindex $args 0]
	    switch -- $arg {
		"-info"           { set info [lindex $args 1]; set skip 2 }
		"-rounds"         { set rounds [lindex $args 1]; set skip 2 }
		"-resrate"        { set gibbs_ghmc::reset_acceptance_rate [lindex $args 1]; set skip 2 }
		"-load"           { set load   [lindex $args 1]; set skip 2 }
		"-connect"        { set master [lindex $args 1]; set skip 2 }
		"-port"           { set port [lindex $args 1]; set skip 2 }
		"-values"         { set values [lindex $args 1]; set skip 2 }
		"-init"           { set init [lindex $args 1]; set skip 2 }
		"-swap_vol"       { set swap_vol [lindex $args 1]; set skip 2 }
		"-swap_part"      { set swap_part [lindex $args 1]; set skip 2 }
		"-move_part"      { set move_part [lindex $args 1]; set skip 2 }
		"-sample"         { set sample [lindex $args 1]; set skip 2 }
		"-nvol"           { set nvol [lindex $args 1]; set skip 2 }
		"-npart"          { set npart [lindex $args 1]; set skip 2 }
		"-nmove"          { set nmove [lindex $args 1]; set skip 2 }
		"-temp"           { set beta [expr 1.0/[lindex $args 1]]; set skip 2 }
		"-vmax"           { set Vmax [lindex $args 1]; set skip 2 }
		"-time_step"      { if {[array names md_dt] != ""} {unset md_dt}; set md_dt [lindex $args 1]; set skip 2 }
		"-max_time_step"  { set dt_max [lindex $args 1]; set skip 2 }
		"-mc_cycle"       { set mc_cycle [lindex $args 1]; set skip 2 }
    "-sample_rounds"  { set sample_rounds [lindex $args 1]; set skip 2 }
		"-adjust_vol"     { set adjust_vol [lindex $args 1]; set skip 2 }
		"-adjust_dt"      { set adjust_dt [lindex $args 1]; set skip 2 }
		default { error "argument \"$arg\" is not a gibbs_ghmc argument" }
	    }
	    set args [lrange $args $skip end]
	}

	if {$swap_vol == "" || $swap_part == "" || $move_part == "" || $sample == ""} {
	    error "you have to specify the swap volume, swap prticle, move particles and sampleing routines"
	}

	if {$values == "" && $master == ""} {
	    error "you either have to specify the master or the value set"
	}

	if {$beta == 0 && $master == ""} {
	    error "you either have to specify master or the simulation temperature"
	}

	if {$Vmax == 0 && $master == ""} { set Vmax [expr log([lindex [lindex $values 0] 0]+[lindex [lindex $values 1] 0])*0.015]}

	if {$mc_cycle == 0 && $master == ""} { set mc_cycle [expr $nvol+$nmove+$npart]}

  if {$sample_rounds == 0 && $master == ""} { set sample_rounds [expr $rounds/2]}

	if {$master == ""} { master::main } { slave::main }
    }
}
