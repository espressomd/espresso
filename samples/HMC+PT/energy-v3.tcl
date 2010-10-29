#!/bin/sh
source ./hmc-initialize.tcl
source ./hmc-setup.tcl

set start [lindex $argv 1]
set end   [lindex $argv 2]
set step  [lindex $argv 3]

set coulombfile [open "coulomb.txt" "r"]
set test " "
while { $test != "=>"} {
    gets $coulombfile cou_info
    set test [lindex $cou_info 0]
}
close $coulombfile
eval inter coulomb [lindex $cou_info 3] p3m [lrange $cou_info 5 9] 
puts "[inter]"

set f [open "|gzip -cd equilibrized.gz" "r"]
while { [blockfile $f read auto] != "eof" } {}
close $f
    
### if something goes wrong and one has to do it partially
#set bjerrum_list [lrange $bjerrum_list 0 20]
#puts "$bjerrum_list"
#exit 0

#foreach bjerrum $bjerrum_list {
#    exec cp "acceptance_lb_[format %08.3f $bjerrum].txt" lb_[format %06.3f $bjerrum]
#}
#exit 0

foreach bjerrum $bjerrum_list {
    puts "*****************************\n BJERRUM $bjerrum"
    eval inter coulomb $bjerrum p3m [lrange $cou_info 5 9] 

    set f [open "lb_[format %06.3f $bjerrum]/poten.txt" a]
    #########################################################################
    for { set config $start } { $config <= $end } { incr config $step} {
	set file [open "|gzip -cd lb_[format %06.3f $bjerrum]/$name.t_[format %05d $config].gz" "r"]
	while { [blockfile $file read auto] != "eof" } {}
	if {[catch {close $file} err]} {
	    if {![regexp "decompression OK, trailing garbage ignored" $err]} {
		error $err
	    }
	}
	set result [analyze energy]
	puts $f "Run= $config $result"
    }
    close $f
}
#exit 0
