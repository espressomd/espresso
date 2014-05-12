#############################################################
#
#  tcp/ip socket based parallel tempering interface
#
#############################################################
#
# Copyright (C) 2010,2012,2013 The ESPResSo project
# Copyright (C) 2007,2008,2009,2010 Axel Arnold
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

namespace eval parallel_tempering {
    # common variables
    #################################

    variable rounds 1
    variable values ""
    variable load 1
    variable port 12000
    variable master ""
    # function calculating the swap probabilities
    variable swap ""
    # name of the function running the simulation
    variable perform ""
    # name of the initialization function per run
    variable init ""
    # data common to all processors. Can only bet set from the master node
    variable shareddata ""
    # information level
    # - none: total silence
    # - comm: communication information on stderr (connection of clients)
    # - all:  information on swap tries on stdout
    variable info "all"
    variable reset_acceptance_rate 10

    # name space for the routines used by the slave processes
    #########################################################
    namespace eval slave {
	variable sock

	####
	proc init {id param} {
	    if {$parallel_tempering::init != ""} {
		eval ::$parallel_tempering::init $id $param
	    }
	}
	####
	proc perform {id param} {
	    eval ::$parallel_tempering::perform $id $param
	}
	####
	proc swap {id param1 param2} {
	    return [eval ::$parallel_tempering::swap $id $param1 $param2]
	}
	####
	proc connect {} {
	    variable sock

	    set sock [socket $parallel_tempering::master $parallel_tempering::port]
	    fconfigure $sock -buffering line
	    set msg [gets $sock]
	    if {$msg != "parallel_tempering welcome"} { error "not the right server for me" }
	    puts $sock "parallel_tempering load $parallel_tempering::load"
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
			init [lindex $msg 1] [lindex $msg 2]
			puts $sock "system [lindex $msg 1] ok"
		    }
		    "shareddata" {
			parallel_tempering::set_shareddata [lrange $msg 1 end]
		    }
		    "perform" { perform [lindex $msg 1] [lindex $msg 2] }
 		    "swapprob" {
			puts $sock [swap [lindex $msg 1] [lindex $msg 2] [lindex $msg 3]]
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

	####
	proc get_acceptance_rate {id1 id2} {
	    variable acceptance	    
	    foreach {tries succ} $acceptance($id1,$id2) { break }
	    return [expr double($succ)/$tries]
	}
	####
	proc reset_acceptance_rates {} {
	    variable acceptance
	    array unset acceptance
	}
	####
	proc slave_connect {channel clientaddr clientport} {
	    variable clients
	    if {$parallel_tempering::info != "none"} {
		puts stderr "connection from $clientaddr $clientport"
	    }
	    fconfigure $channel -buffering line
	    set n_systems [llength $parallel_tempering::values]
	    set n_clients [expr [llength $clients]/2]
	    # check for overoccupation
	    if {$n_clients >= $n_systems} {
		puts $channel "go away"
		close $channel
		return
	    }
	    puts $channel "parallel_tempering welcome"
	    set response [gets $channel]
	    if {![string match "parallel_tempering load *" $response]} {
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

	    set n_systems [llength $parallel_tempering::values]
	    while { 1 } {
		set n_clients [expr [llength $clients]/2]
		if {$n_clients >= $n_systems} { break }
		if {$parallel_tempering::info != "none"} {
		    puts stderr "waiting for [expr $n_systems - $n_clients] more clients"
		}
		vwait parallel_tempering::master::clients
	    }
	    if {$parallel_tempering::info != "none"} {
		puts stderr "all there, let's go"
	    }
	}
	####
	proc init_clients {} {
	    variable clients
	    # send initialization command
	    foreach param $parallel_tempering::values {channel id} $clients {
		if {$channel != "myself"} {
		    puts $channel "system $id $param"
		} {
		    lappend myparam $id $param
		}
	    }
	    # now, initialize myself
	    foreach {id param} $myparam { parallel_tempering::slave::init $id $param }
	    # and wait for all clients to report ready
	    foreach {channel id} $clients {
		if {$channel != "myself"} {
		    set response [gets $channel]
		    if {$response != "system $id ok"} {
			error "client connection error, wrong response"
		    }
		    puts $channel "shareddata $parallel_tempering::shareddata"
		}
	    }
	    set shareddata_changed 0
	}
	####
	proc perform_clients {} {
	    variable clients
	    variable shareddata_changed

	    # broadcast current parameters
	    foreach param $parallel_tempering::values {channel id} $clients {
		if {$channel == "myself"} {
		    lappend myparam $id $param
		} {
		    if {$shareddata_changed} {
			puts $channel "shareddata $parallel_tempering::shareddata"
		    }
		    puts $channel "perform $id $param"
		}
	    }
	    # perform myself
	    foreach {id param} $myparam { parallel_tempering::slave::perform $id $param }
	    set shareddata_changed 0
	}
	####
	proc ask_swapprob {channel id param1 param2} {
	    if {$channel == "myself"} {
		set res [parallel_tempering::slave::swap $id $param1 $param2]
	    } {
	    	puts $channel "swapprob $id $param1 $param2"
	    	set res [gets $channel]
	    }
	    if {[llength $res] != 2} {
		error "swap energies format wrong, \"$res\" is not a list of two values."
	    }
	    return $res
	}
	####
	proc swap_clients {odd} {
	    variable clients
	    variable acceptance

	    if {$parallel_tempering::info == "all"} {
		puts "parallel_tempering swapping start"
	    }
	    # odd-even rule. We either try to exchange 1-2
	    # or we simply keep 1 and start with 2-3
	    if {$odd} {
		# keep first client
		set v [lrange $parallel_tempering::values 1 end]
		set c [lrange $clients 2 end]
		set newclients [lrange $clients 0 1]
	    } {
		# do not keep it
		set v $parallel_tempering::values
		set c $clients
	    }

 	    foreach {param1 param2} $v {channel1 id1 channel2 id2} $c {
		# single parameter left over, cannot swap
		if {$param2 == ""} { lappend newclients $channel1 $id1; break }

		# ask the clients for the swapping energies
		foreach {E11 E12} [ask_swapprob $channel1 $id1 $param1 $param2] { break }
		foreach {E21 E22} [ask_swapprob $channel2 $id2 $param1 $param2] { break }

		# check whether to actually swap by Metropolis
		set exp [expr -($E12 + $E21 - $E11 - $E22)]
		if {$exp > 0} { set boltzmann 1	} elseif {$exp < -30} { set boltzmann 0	} {
		    set boltzmann [expr exp($exp)]
		}
		if {$parallel_tempering::info == "all"} {
		    puts "parallel_tempering tryswap $param1 $param2 : energies $E11 $E12 $E21 $E22 -> $boltzmann"
		}

		set pairid "$param1,$param2"
		if {[array names acceptance $pairid] == ""} {
		    set tries 1
		    set succ 0
		} {
		    foreach {tries succ} $acceptance($pairid) { break }
		    incr tries
		}

		set myrand [t_random]
		if {[expr $myrand < $boltzmann]} {
		    if {$parallel_tempering::info == "all"} {
			puts "parallel_tempering swap $param1 $param2"
		    }
		    lappend newclients $channel2 $id2 $channel1 $id1
		    incr succ
		} {
		    lappend newclients $channel1 $id1 $channel2 $id2
		}
		set acceptance($pairid) "$tries $succ"
	    }

	    set clients $newclients
	    if {$parallel_tempering::info == "all"} {
		puts "parallel_tempering swapping done"
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
		foreach param $parallel_tempering::values {channel id} $clients {
		    if {$sparam == $param} { puts $channel "$id $string" }
		}
	    }
	}
	####
	proc recv {sparam} {
	    variable clients
	    foreach param $parallel_tempering::values {channel id} $clients {
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
	    for {set i 0} {$i < $parallel_tempering::load} {incr i} { lappend clients "myself" $i }

	    set server_socket [socket -server slave_connect $parallel_tempering::port]

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
	    while { $parallel_tempering::rounds > 0 } {
		if {$parallel_tempering::info == "all"} {
		    if {[expr $parallel_tempering::rounds % $parallel_tempering::reset_acceptance_rate] == 0} {
			if {[array names acceptance] != ""} {
			    puts "parallel_tempering listing acceptance_rates"
			    set pp ""
			    foreach cp $parallel_tempering::values {
				if {$pp != ""} {
				    puts "parallel_tempering acceptance_rate $pp $cp [get_acceptance_rate $pp $cp]"
				}
				set pp $cp
			    }
			    reset_acceptance_rates
			}
		    }
		}
		perform_clients
		swap_clients [expr ($parallel_tempering::rounds % 2) == 0]
		incr parallel_tempering::rounds -1
	    }

	    release_clients
	    close $server_socket
	}
    }

    # parallel tempering public methods
    ################################################

    ####
    proc set_shareddata {data} {
	variable shareddata
	set shareddata $data
	set master::shareddata_changed 1
    }
    ####
    proc main {args} {
	variable rounds; variable master; variable port
	variable values; variable load;	variable info
	variable init; variable swap; variable perform

	set init ""
	set values ""
	set master ""

	while {$args != ""} {
	    set arg [lindex $args 0]
	    switch -- $arg {
		"-info"     { set info [lindex $args 1]; set skip 2 }
		"-rounds"   { set rounds [lindex $args 1]; set skip 2 }
		"-resrate"  { set reset_acceptance_rate [lindex $args 1]; set skip 2 }
		"-load"     { set load   [lindex $args 1]; set skip 2 }
		"-connect"  { set master [lindex $args 1]; set skip 2 }
		"-port"     { set port [lindex $args 1]; set skip 2 }
		"-values"   { set values [lindex $args 1]; set skip 2 }
		"-init"     { set init [lindex $args 1]; set skip 2 }
		"-swap"     { set swap [lindex $args 1]; set skip 2 }
		"-perform"  { set perform [lindex $args 1]; set skip 2 }
		default { error "argument \"$arg\" is not a parallel_tempering argument" }
	    }
	    set args [lrange $args $skip end]
	}
	if {$swap == "" || $perform == ""} {
	    error "you have to specify the swap energy and performing routines"
	}
	if {$values == "" && $master == ""} {
	    error "you either have to specify the master or the value set"
	}

	if {$master == ""} { master::main } { slave::main }
    }
}
