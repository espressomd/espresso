

source ./hmc-initialize.tcl
source ./hmc-setup.tcl

set dirname "/tmp/"
##############################################################
##############################################################
set start [lindex $argv 1]
set end   [lindex $argv 2]
set step  [lindex $argv 3]
# flag wether you want to delete the created psf/pdb files
set delete 1

set name "HMC"
set bjerrum  [lindex $argv 4]
#########################################################################33
for { set config $start } { $config <= $end } { incr config $step} {
    puts  -nonewline "Read checkpoint.[format %05d $config] - create pdb\r" 
    flush stdout
    set file [open "|gzip -cd lb_[format %06.3f $bjerrum]/$name.t_[format %05d $config].gz" "r"]
    while { [blockfile $file read auto] != "eof" } {}
    close $file
    if { $config == $start } { writepsf "$dirname/$name.vmd.psf" }
    writepdb "$dirname/$name.vmd[format %04d [expr $config - $start]].pdb" -folded
}



set vmd_file [open "$dirname/vmd_animation.script" "w"]
puts $vmd_file "loadseries $dirname/$name.vmd $step"
puts $vmd_file "rotate stop"
puts $vmd_file "mol modstyle 0 0 CPK 1.000000 0.300000 8.000000 6.000000"
puts $vmd_file "mol modcolor 0 0 SegName"
puts $vmd_file "display nearclip set 0.01"
puts $vmd_file "logfile $dirname/vmd.log"
puts $vmd_file "scale by 1.7"
puts $vmd_file "animate forward"
puts $vmd_file "logfile off"
#puts $vmd_file "quit"
close $vmd_file

exec vmd -e $dirname/vmd_animation.script &

puts "logfile: [file exists \"vmd.log\" ]"
# wait untill vmd has read in the pdb files
set vmd_finish 0
while { $vmd_finish == 0 } {
    if { [file exists "$dirname/vmd.log" ] } { set vmd_finish 1 }
}

exec sleep 10

if { $delete != 0 } {
    for { set config $start } { $config <= $end } { incr config $step} {
	if { $config == $start } { exec rm "$dirname/$name.vmd.psf" }
	exec rm "$dirname/$name.vmd[format %04d [expr $config - $start]].pdb"
    }
    exec rm $dirname/vmd_animation.script
}

exec rm $dirname/vmd.log

puts "\nFinished\n"


exit









