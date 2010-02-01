setmd box_l 10 10 10
inter 0 fene 1.0 1.0
for {set i 0} {$i < 10000} {incr i} {polymer 1 2 1 pos 5 5 5 start [expr $i*2]}
writepdb "sphere.pdb"
set vmdout_file [open "vmd_start.script" "w"]
puts $vmdout_file "mol load pdb sphere.pdb"
puts $vmdout_file "logfile vmd.log"
puts $vmdout_file "rotate stop"
puts $vmdout_file "logfile off"
puts $vmdout_file "mol modstyle 0 0 Points"
close $vmdout_file
exec vmd -e vmd_start.script &

