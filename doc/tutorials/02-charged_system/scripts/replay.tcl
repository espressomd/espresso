#
# Copyright (C) 2010,2012,2013 The ESPResSo project
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
