set conn 0

foreach filename [lrange $argv 1 end] {
    set f [open $filename "r"]
    while { [blockfile $f read auto] != "eof" } {}
    close $f

    if {!$conn} {
	set filename "replay"
	writepsf "$filename.psf"
	writepdb "$filename.pdb" -folded
	for {set port 10000} { $port < 65000 } { incr port } {
	    catch {imd connect $port} res
	    if {$res == ""} break
	}
	set vmdout_file [open "$filename.vmd.tcl" "w"]
	puts $vmdout_file "mol load psf $filename.psf pdb $filename.pdb"
	puts $vmdout_file "rotate stop"
	puts $vmdout_file "logfile off"
	puts $vmdout_file "mol modstyle 0 0 CPK 2.0 0.3 8 6"
	puts $vmdout_file "mol modcolor 0 0 SegName"
	puts $vmdout_file "imd connect localhost $port"
	puts $vmdout_file "imd transfer 1"
	puts $vmdout_file "imd keep 1"
	close $vmdout_file

	exec vmd -e "$filename.vmd.tcl" &
	imd listen 20000
	set conn 1
    }
    after 1000
    imd positions
}
