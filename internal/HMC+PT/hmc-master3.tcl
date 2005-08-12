#!/bin/sh
source ./hmc-initialize.tcl
source ./hmc-setup.tcl
set mywaittime [expr [lindex $argv 1] * 1000 ]
puts "MYWAITTIME $mywaittime"

#############################################################
#      Integration                                          #
#############################################################
set time_step    0.005
setmd time_step $time_step
setmd gamma     0.0
setmd temp      0.0
thermostat off

analyze set chains 0 $n_poly $p_length

set max_num [setmd n_part]
analyze remove

#############################################################
#      Read the equilibrated configurations                 #
#############################################################
if { $ensemble_num < 2} { set check_num 1} else {set check_num 2}
for {set ens_i 0} { $ens_i < $check_num } {incr ens_i} {
    set bjerrum [lindex $bjerrum_list $ens_i]
    set f [open "|gzip -cd equilibrized.gz" "r"]
    while { [blockfile $f read auto] != "eof" } {}
    analyze append 
    close $f
}
puts "STORED [analyze stored]"



# particle numbers
set n_mono [expr $n_poly * $p_length]
set n_ci   [expr $n_mono/$c_ion_val]



#############################################################
#      Read the coulomb information
#############################################################

set coulombfile [open "coulomb.txt" "r"]
set test " "
while { $test != "=>"} {
    gets $coulombfile cou_info
    set test [lindex $cou_info 0]
}
close $coulombfile
eval inter coulomb [lindex $cou_info 3] p3m [lrange $cou_info 5 9] 
puts "[inter]"

set f [open "wherewasI" "r"]
gets $f i_loop
close $f

after $mywaittime
set oldness 36000
for {set ens_i [expr $ensemble_num-1]} { $ens_i > -1 } {incr ens_i -1} {
    after $mywaittime
    set bjerrum [lindex $bjerrum_list $ens_i]
    if { [file exists "lb_[format %06.3f $bjerrum]/lock_[format %05d $i_loop]"] } {
        if { ![file exists "lb_[format %06.3f $bjerrum]/done_[format %05d $i_loop]"] } {
            set creattime [file atime "lb_[format %06.3f $bjerrum]/lock_[format %05d $i_loop]" ]
            set currenttime [expr int([lindex [exec date +%s])]];
            if { [expr $currenttime - $creattime] >  $oldness} {
                puts "KILL BILL BJERRUM $bjerrum ILOOP $i_loop"
                exec rm -rf "lb_[format %06.3f $bjerrum]/lock_[format %05d $i_loop]"
            }
        }
    }
}


exec touch "masterisrunning"

set success 0

for {set ens_i [expr $ensemble_num-1]} { $ens_i > -1 } {incr ens_i -1} {
    lappend checklist $ens_i
    lappend orig_checklist $ens_i
}

while { $i_loop < $loops } {

    # FIRST INTEGRATE ALL ENSEMBLES
    set checklist $orig_checklist
    set success [llength $checklist]
    set int_counter 0
    while { $success > 0} {
	puts "INTEGRATION LOOP"
	after $mywaittime
	if { $int_counter > 100 } {puts "NORMAL int_counter exit "; exit 0}
	incr int_counter

	set trac_i -1
	foreach ens_i $checklist {
	    incr trac_i
	    set bjerrum [lindex $bjerrum_list $ens_i]
	    if { [file exists "lb_[format %06.3f $bjerrum]/lock_[format %05d $i_loop]"] } {
		set checklist [lreplace $checklist $trac_i $trac_i]
		incr success -1
		incr trac_i -1
	    } else {
		if { [file exists "lb_[format %06.3f $bjerrum]/swpd_[format %05d [expr $i_loop-1]]"] } { 
		    exec touch "lb_[format %06.3f $bjerrum]/lock_[format %05d $i_loop]"
		    file attributes "lb_[format %06.3f $bjerrum]/lock_[format %05d $i_loop]" -permissions -w
		    set energyfile [open "lb_[format %06.3f $bjerrum]/energy.txt" "a"]
		    
		    #read the last configuration
		    set file [open "|gzip -cd lb_[format %06.3f $bjerrum]/$name.t_[format %05d [expr $i_loop * ($int_n_times +1)]].gz" "r"]
		    while { [blockfile $file read auto] != "eof" } {}
		    if {[catch {close $file} err]} {
			if {![regexp "decompression OK, trailing garbage ignored" $err]} {
			    error $err
			}
		    }
		    analyze replace 0
		    
		    # DO ENSEMBLE INTEGRATION
		    puts " ILOOP $i_loop BJERRUM $bjerrum checklist $checklist $trac_i"
		    #####################################################
		    ### HMC move
		    #####################################################
		    
		    if { $bjerrum == 0.0 } {inter coulomb $bjerrum
		    } else {eval inter coulomb $bjerrum p3m [lrange $cou_info 5 9]}
		    for {set i 0} { $i < $int_n_times } { incr i} {
			#puts "B [part 0 print pos]"
			maxwell_velocities
			set oldener [analyze energy]
			integrate $int_steps
			set ener [analyze energy]
			set boltzman [expr exp(-([lindex [lindex $ener 0] 1]-[lindex [lindex $oldener 0] 1]))]
			set myrand [t_random]
			if { $myrand < $boltzman } {
			    #puts "ACCEPTED BOLTZMAN $boltzman $myrand ENER [lindex [lindex $ener 0] 1] OLDENER [lindex [lindex $oldener 0] 1]"
			    #puts "[part 0 print pos]"
			    puts  $energyfile  "Run [expr $i_loop * ($int_n_times +1) + $i +1 ] Energy = $ener RAND $myrand BOLTZMAN $boltzman ENERGIES OLD: $oldener NEW: $ener"
			    analyze replace 0;
			} else {
			    #puts "REJECTED BOLTZMAN $boltzman $myrand ENER [lindex [lindex $ener 0] 1] OLDENER [lindex [lindex $oldener 0] 1]"
			    analyze activate 0;
			    #puts "[part 0 print pos]"
			    puts  $energyfile  "Run [expr $i_loop * ($int_n_times +1) + $i +1 ] Energy = $oldener RAND $myrand BOLTZMAN $boltzman ENERGIES OLD: $oldener NEW: $ener"
			};
			
			####################################################
			# Write out the output from the sensemble          #
			####################################################
			set out [open "|gzip -c - >lb_[format %06.3f $bjerrum]/$name.t_[format %05d [expr $i_loop * ($int_n_times +1) + $i +1 ] ].gz" "w"]
			blockfile $out write particles "id pos type q " all
			close $out	
		    }
		    flush $energyfile
		    close $energyfile
		    exec touch  "lb_[format %06.3f $bjerrum]/done_[format %05d $i_loop]"

		    set checklist [lreplace $checklist $trac_i $trac_i]
		    incr success -1
		    incr trac_i -1
		}
	    }
	}
    }

    
    ###############################################################
    #### Ensemble exchange 
    ###############################################################
    set checklist $orig_checklist
    set success $ensemble_num
    
    while { $success == $ensemble_num} {
	puts "PARALLEL TEMPERING"
	after $mywaittime
	if { [expr fmod($i_loop,2)] == 0 } {
	    set start_ens 0
	    if {[expr fmod($ensemble_num,2)] == 1 } {
		set bjerrum [lindex $bjerrum_list [expr $ensemble_num -1] ]
		if { [file exists "lb_[format %06.3f $bjerrum]/done_[format %05d $i_loop]"] } {
		    if { ![file exists "lb_[format %06.3f $bjerrum]/swpl_[format %05d $i_loop]"] } {
			exec touch "lb_[format %06.3f $bjerrum]/swpl_[format %05d $i_loop]"
			file attributes  "lb_[format %06.3f $bjerrum]/swpl_[format %05d $i_loop]" -permissions -w
			exec cp "lb_[format %06.3f $bjerrum]/$name.t_[format %05d [expr $i_loop * ($int_n_times +1)+ $int_n_times ]].gz" "lb_[format %06.3f $bjerrum]/$name.t_[format %05d [expr $i_loop * ($int_n_times +1)+ $int_n_times+1 ]].gz"
			exec touch  "lb_[format %06.3f $bjerrum]/swpd_[format %05d $i_loop]"
			puts "SWP $bjerrum"
		    }
		}
	    }
	} else { 
	    set start_ens 1 
	    set bjerrum [lindex $bjerrum_list 0 ]
	    if { [file exists "lb_[format %06.3f $bjerrum]/done_[format %05d $i_loop]"] } {
		if { ![file exists "lb_[format %06.3f $bjerrum]/swpl_[format %05d $i_loop]"] } {
		    exec touch  "lb_[format %06.3f $bjerrum]/swpl_[format %05d $i_loop]"
		    file attributes  "lb_[format %06.3f $bjerrum]/swpl_[format %05d $i_loop]" -permissions -w
		    exec cp "lb_[format %06.3f $bjerrum]/$name.t_[format %05d [expr $i_loop * ($int_n_times +1)+ $int_n_times ]].gz" "lb_[format %06.3f $bjerrum]/$name.t_[format %05d [expr $i_loop * ($int_n_times +1)+ $int_n_times+1 ]].gz"
		    exec touch "lb_[format %06.3f $bjerrum]/swpd_[format %05d $i_loop]"
		    puts "SWP $bjerrum"
		}
	    }
	    if {[expr fmod($ensemble_num,2)] == 0 } {
		set bjerrum [lindex $bjerrum_list [expr $ensemble_num -1] ]
		if { [file exists "lb_[format %06.3f $bjerrum]/done_[format %05d $i_loop]"] } {
		    if { ![file exists "lb_[format %06.3f $bjerrum]/swpl_[format %05d $i_loop]"] } {
			exec touch "lb_[format %06.3f $bjerrum]/swpl_[format %05d $i_loop]"
			file attributes  "lb_[format %06.3f $bjerrum]/swpl_[format %05d $i_loop]" -permissions -w
			exec cp "lb_[format %06.3f $bjerrum]/$name.t_[format %05d [expr $i_loop * ($int_n_times +1)+ $int_n_times ]].gz" "lb_[format %06.3f $bjerrum]/$name.t_[format %05d [expr $i_loop * ($int_n_times +1)+ $int_n_times+1 ]].gz"
			exec touch "lb_[format %06.3f $bjerrum]/swpd_[format %05d $i_loop]"
			puts "SWP $bjerrum"
		    }
		}
	    }
	}

	set finish_ens [expr $ensemble_num - 1]
	
	for {set ens_i $start_ens} { $ens_i <  $finish_ens} {set ens_i [expr $ens_i + 2]} {
	    set ens_ip1 [expr $ens_i + 1]
	    set bjerrum0 [lindex $bjerrum_list $ens_i ]
	    set bjerrum1 [lindex $bjerrum_list $ens_ip1 ]

	    if { [file exists "lb_[format %06.3f $bjerrum0]/done_[format %05d $i_loop]"] && [file exists "lb_[format %06.3f $bjerrum1]/done_[format %05d $i_loop]"]} {
		if { ![file exists "lb_[format %06.3f $bjerrum0]/swpl_[format %05d $i_loop]"] && ![file exists "lb_[format %06.3f $bjerrum1]/swpl_[format %05d $i_loop]"]} {
		    exec touch "lb_[format %06.3f $bjerrum0]/swpl_[format %05d $i_loop]"
		    exec touch "lb_[format %06.3f $bjerrum1]/swpl_[format %05d $i_loop]"
		    file attributes  "lb_[format %06.3f $bjerrum0]/swpl_[format %05d $i_loop]" -permissions -w
		    file attributes  "lb_[format %06.3f $bjerrum1]/swpl_[format %05d $i_loop]" -permissions -w

		    set energyfile0 [open "lb_[format %06.3f $bjerrum0]/energy.txt" "a"]
		    set energyfile1 [open "lb_[format %06.3f $bjerrum1]/energy.txt" "a"]

		    
		    set file [open "|gzip -cd lb_[format %06.3f $bjerrum1]/$name.t_[format %05d [expr $i_loop * ($int_n_times +1)+ $int_n_times ]].gz" "r"]
		    while { [blockfile $file read auto] != "eof" } {}
		    if {[catch {close $file} err]} {
			if {![regexp "decompression OK, trailing garbage ignored" $err]} {
			    error $err
			}
		    }
		    analyze replace 1

		    set file [open "|gzip -cd lb_[format %06.3f $bjerrum0]/$name.t_[format %05d [expr $i_loop * ($int_n_times +1)+ $int_n_times ]].gz" "r"]
		    while { [blockfile $file read auto] != "eof" } {}
		    if {[catch {close $file} err]} {
			if {![regexp "decompression OK, trailing garbage ignored" $err]} {
			    error $err
			}
		    }
		    analyze replace 0
		    		    
		    if { $bjerrum0 == 0.0 } {
			inter coulomb $bjerrum0
		    } else {
			eval inter coulomb $bjerrum0 p3m [lrange $cou_info 5 9]	    
		    }
		    set energyold0 [analyze energy]
		    set oldener0 [lindex [lindex $energyold0 0] 1]
		    if { $bjerrum1 == 0.0 } {
			inter coulomb $bjerrum1
		    } else {
			eval inter coulomb $bjerrum1 p3m [lrange $cou_info 5 9]	    
		    }
		    set energy1 [analyze energy]
		    set newener1 [lindex [lindex  $energy1 0] 1]
		    
		    analyze activate 1
		    set energyold1 [analyze energy]
		    set oldener1 [lindex [lindex $energyold1 0] 1]
		    if { $bjerrum0 == 0.0 } {
			inter coulomb $bjerrum0
		    } else {
			eval inter coulomb $bjerrum0 p3m [lrange $cou_info 5 9]
		    }
		    set energy0 [analyze energy]
		    set newener0 [lindex [lindex $energy0 0] 1]
		    
		    set boltzman [expr exp(-( $newener0 + $newener1 - $oldener0 -$oldener1 ))]
		    set myrand [t_random]
		    if { $myrand < $boltzman } {
			#puts " I AM HERE ACCEPTED"
			puts  $energyfile0  "SWA [expr $i_loop * ($int_n_times +1) + $int_n_times +1 ] Energy = $energy0 RAND $myrand BOLTZMAN $boltzman ENERGIES OLD: $energyold0 NEW: $energy0 "
			puts  $energyfile1  "SWA [expr $i_loop * ($int_n_times +1) + $int_n_times +1 ] Energy = $energy1 RAND $myrand BOLTZMAN $boltzman ENERGIES OLD: $energyold1 NEW: $energy1 "
			set out [open "|gzip -c - >lb_[format %06.3f $bjerrum0]/$name.t_[format %05d [expr $i_loop * ($int_n_times +1) + $int_n_times +1 ]].gz" "w"]
			blockfile $out write particles "id pos type q " all
			close $out	
			
			analyze activate 0
			set out [open "|gzip -c - >lb_[format %06.3f $bjerrum1]/$name.t_[format %05d [expr $i_loop * ($int_n_times +1) + $int_n_times +1 ]].gz" "w"]
			blockfile $out write particles "id pos type q " all
			close $out	
			
		    } else {
			puts  $energyfile0  "SWA [expr $i_loop * ($int_n_times +1) + $int_n_times +1] Energy = $energyold0 RAND $myrand BOLTZMAN $boltzman ENERGIES OLD: $energyold0 NEW: $energy0 "
			puts  $energyfile1  "SWA [expr $i_loop * ($int_n_times +1) + $int_n_times +1 ] Energy = $energyold1 RAND $myrand BOLTZMAN $boltzman ENERGIES OLD: $energyold1 NEW: $energy1 "
			set out [open "|gzip -c - >lb_[format %06.3f $bjerrum1]/$name.t_[format %05d [expr $i_loop * ($int_n_times +1) + $int_n_times +1 ]].gz" "w"]
			blockfile $out write particles "id pos type q " all
			close $out	
			
			analyze activate 0
			set out [open "|gzip -c - >lb_[format %06.3f $bjerrum0]/$name.t_[format %05d [expr $i_loop * ($int_n_times +1) + $int_n_times +1 ]].gz" "w"]
			blockfile $out write particles "id pos type q " all
			close $out	
			
		    }
		    flush $energyfile0
		    flush $energyfile1
		    close $energyfile0
		    close $energyfile1
		    exec touch "lb_[format %06.3f $bjerrum0]/swpd_[format %05d $i_loop]"
		    exec touch "lb_[format %06.3f $bjerrum1]/swpd_[format %05d $i_loop]"
		    puts "SWP $bjerrum0"
		    puts "SWP $bjerrum1"
		}
	    }
	}
	
	set success $ensemble_num
	foreach bjerrum $bjerrum_list {
	    if { [file exists "lb_[format %06.3f $bjerrum]/swpd_[format %05d $i_loop]"] } { 
		incr success -1;
	    }
	}
	puts "success $success"
    }

    
    incr i_loop

    if { $success == 0 } {
	set f [open "wherewasI" "w"]
	puts $f "$i_loop"
	close $f
    }

    if { ![file exists "masterisrunning"]} {puts "NORMAL PROGRAM TERMINATION :)"; exit 0;} 

    set currenttime [expr int([lindex [exec date +%s])]];  
    if { [expr $currenttime - $starttime ] > $runtime } {
	puts "NORMAL PROGRAM TERMINATION :)  " ; 
	exec rm -rf "masterisrunning"
	exit 0
    }     
    
}

puts "NORMAL PROGRAM TERMINATION :)  " ; 
exec rm -rf "masterisrunning"
puts "\nIntegration done."
exit 0
#############################################################
