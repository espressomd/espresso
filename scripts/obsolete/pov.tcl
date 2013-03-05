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
proc writepov {file {folded "no"} {boxopt "no"} {render "no"} {rotate "-10 22.5 0"} {nofold_list 0} } {
#writes the current configuration in pov-ray format
#file is the name of the output
#use -folded if you want folded coordinates
#use -box if you want to output the box
#use -render if you want to render the pov-ray file
#render contains a vector by which the entire configuration is rotated
#nofold_list contains a list of types that should not be folded if 'folded' is set to -folded
    if {$folded == "-folded"} {
	set de "folded"
    } else {set de "pos"}
    if {[file exists "$file"] } {
        error "ERROR: $file already exists; please use a different file name!\nAborting..." 
    } else {
        set f [open $file "w"]
    }
    set box [setmd box_l]
    
    #write background, camera, lighting
    puts $f "background {color rgb <0.98,0.98,0.98>}"
    puts $f "light_source \{\n<1000,1000,-100>\ncolor rgb 1\n\}"
    puts $f "light_source \{\n<-1000,1000,-100>\ncolor rgb 1\n\}"
    puts $f "camera \{\northographic\nlocation <[expr [lindex $box 0]*0.5],\
             [expr [lindex $box 1]*0.5],-[expr [lindex $box 2]*1.5]>\n\
             look_at <[expr [lindex $box 0]*0.5],[expr [lindex $box 1]*0.5],[expr [lindex $box 2]*0.5]>\n\}"
    
    #write atoms and bonds
    set radius [expr [lindex $box 0]*0.004]
    set b_radius [expr $radius*0.25]
    lappend color  <0.36,0.57,0.31> <0.66,0.0,0.12> <0.20,0.32,0.53> 
    set mp [setmd max_part]
    set texture "texture { finish { ambient 0.000 diffuse 0.650 phong 0.1 phong_size 40.000 specular 0.500 } }"
    puts $f "union \{"
    for {set p 0} { $p <= $mp} {incr p} {
	set type [part $p print type]
	# if {$type==0 || $type==1} {set de "pos"} else {set de "folded"}
	if {$nofold_list != 0} { if {[lsearch $nofold_list $type] != -1 } {set de "pos"} else {set de "folded"} }
	set pos [part $p print $de]
	puts $f "\/\/ $p POS\nsphere \{\n<[lindex $pos 0],[lindex $pos 1],[lindex $pos 2]>,\n$radius\n$texture\n\
                 pigment \{color rgb [lindex $color $type] \}\nno_shadow\n\}"
	if {[llength [lindex [part $p print bond] 0]] > 0} {
	    set n_bond [llength [lindex [part $p print bond] 0]]
	    for {set q 0} {$q < $n_bond} {incr q} {
		set end [part [lindex [part $p print bond] 0 $q 1] print $de]
		if {[veclen [vecsub $pos $end]] < [expr [lindex $box 0]*0.8]} {
		    puts $f "\/\/ $p [lindex [part $p print bond] 0 $q 1] BOND\n\
                             cylinder \{\n<[lindex $pos 0],[lindex $pos 1],[lindex $pos 2]>,\n\
                            <[lindex $end 0],[lindex $end 1],[lindex $end 2]>,\n$b_radius\n\
                            $texture\npigment \{ color rgb [lindex $color $type]\}\nno_shadow\n\}"
		}
	    }
	}	
    }
    #write box
    if {$boxopt == "-box"} {
	set c_radius [expr $radius*1.5]
	puts $f "\/\/ box"
	puts $f "cylinder \{\n<0,0,0>,\n<[lindex $box 0],0,0>,\n$c_radius\n$texture\npigment \{ color rgb [lindex $color 0]\}\
               no_shadow \n \}"
	puts $f "cylinder \{\n<0,0,0>,\n<0,[lindex $box 1],0>,\n$c_radius\n$texture\npigment \{ color rgb [lindex $color 0]\}\
               no_shadow \n \}"
	puts $f "cylinder \{\n<0,0,0>,\n<0,0,[lindex $box 2]>,\n$c_radius\n$texture\npigment \{ color rgb [lindex $color 0]\}\
               no_shadow \n \}"

	puts $f "cylinder \{\n<[lindex $box 0],0,[lindex $box 2]>,\n<[lindex $box 0],0,0>,\n$c_radius\n$texture\npigment\ 
              \{ color rgb [lindex $color 0]\} \nno_shadow \n \}"
	puts $f "cylinder \{\n<[lindex $box 0],0,[lindex $box 2]>,\n<[lindex $box 0],[lindex $box 1],[lindex $box 2]>,\n$c_radius\n$texture\npigment\
              \{ color rgb [lindex $color 0]\} \nno_shadow \n \}"
	puts $f "cylinder \{\n<[lindex $box 0],0,[lindex $box 2]>,\n<0,0,[lindex $box 2]>,\n$c_radius\n$texture\npigment\
              \{ color rgb [lindex $color 0]\} \nno_shadow \n \}"

	puts $f "cylinder \{\n<[lindex $box 0],[lindex $box 1],0>,\n<[lindex $box 0],0,0>,\n$c_radius\n$texture\npigment\
              \{ color rgb [lindex $color 0]\} \nno_shadow \n \}"
	puts $f "cylinder \{\n<[lindex $box 0],[lindex $box 1],0>,\n<0,[lindex $box 1],0>,\n$c_radius\n$texture\npigment\
              \{ color rgb [lindex $color 0]\} \nno_shadow \n \}"
	puts $f "cylinder \{\n<[lindex $box 0],[lindex $box 1],0>,\n<[lindex $box 0],[lindex $box 1],[lindex $box 2]>,\n$c_radius\n$texture\npigment\
              \{ color rgb [lindex $color 0]\} \nno_shadow \n \}"

	puts $f "cylinder \{\n<0,[lindex $box 1],[lindex $box 2]>,\n<0,[lindex $box 1],0>,\n$c_radius\n$texture\npigment\
              \{ color rgb [lindex $color 0]\} \nno_shadow \n \}"
	puts $f "cylinder \{\n<0,[lindex $box 1],[lindex $box 2]>,\n<0,0,[lindex $box 2]>,\n$c_radius\n$texture\npigment\
              \{ color rgb [lindex $color 0]\} \nno_shadow \n \}"
	puts $f "cylinder \{\n<0,[lindex $box 1],[lindex $box 2]>,\n<[lindex $box 0],[lindex $box 1],[lindex $box 2]>,\n$c_radius\n$texture\npigment\
              \{ color rgb [lindex $color 0]\} \nno_shadow \n \}"
    }
    #write rotation
    puts $f "rotate <[lindex $rotate 0],[lindex $rotate 1],[lindex $rotate 2]>"
    puts $f "\}"
    close $f
    
    #Render pov-file
    if {$render=="-render"} {
	catch { exec povray +FT +W800 +H600 +A0.2 -I$file -O$file.tga }
    }
}

proc morph {f1 f2 f3 {arg 100} {render ""}} {
# morphs 2 povray-files in arg steps
#f1 is the start of the morph
#f2 is the end of the morph
#f3 is the name of the output

    set inp2 [open $f2 r]
    set list ""
    while {[eof $inp2]==0} {
    	set l2 [gets $inp2]
	set wl2 [string map { "\{" " "} $l2]
	if {[lindex $wl2 2]=="POS"} {
	    gets $inp2
	    set l2 [gets $inp2]
	    lappend list $l2
	}
    }
    close $inp2

    for {set j 0} {$j < [expr $arg+1]} {incr j} {
	set outp [open $f3-[format %05d $j].pov w]
	set inp1 [open $f1 r]
	
	while {[eof $inp1]==0} {

	    set l1 [gets $inp1]
	    set wl1 [string map { "\{" " "} $l1]
	    set wl1 [string map { "\}" " "} $wl1]
	
	    if {[lindex $wl1 2]=="POS"} {
		puts $outp $l1
		set part [lindex $l1 1]
		set l1 [gets $inp1]
		puts $outp $l1	    
		set l1 [gets $inp1]		
		set l2 [lindex $list $part]
		set pos1 [string map { "<" " "} $l1]
		set pos2 [string map { "<" " "} $l2]
		set pos1 [string map { ">" " "} $pos1]
		set pos2 [string map { ">" " "} $pos2]
		set pos1 [string map { "," " "} $pos1]
		set pos2 [string map { "," " "} $pos2]
		set morph_pos [vecadd $pos1 [vecscale [expr $j.0/$arg.0] [vecsub $pos2 $pos1]]]
		puts $outp "<[lindex $morph_pos 0],[lindex $morph_pos 1],[lindex $morph_pos 2]>"
		set l1 [gets $inp1]
		puts $outp $l1
		set l1 [gets $inp1]
		puts $outp $l1
		set l1 [gets $inp1]
		set rgb [string map { "<" " "} $l1]
		set rgb [string map { ">" " "} $rgb]
		set rgb [string map { "," " "} $rgb]
		set color "[lindex $rgb 1 2] [lindex $rgb 1 3] [lindex $rgb 1 4]"
		puts $outp "pigment \{ color rgbt <[lindex $color 0],[lindex $color 1],[lindex $color 2],0.6> \}"
		set l1 [gets $inp1]
		puts $outp $l1
		set l1 [gets $inp1]
		puts $outp $l1
	    } elseif {[lindex $wl1 3]=="BOND"} {
 
		for {set cnt 0} {$cnt < 8} {incr cnt} {
		    gets $inp1
		}
	    } else {
		puts $outp $l1
	    }
	}
	close $outp
	close $inp1
       
	if {$render=="-render"} {
	    puts "rendering"
	    catch {exec povray +FT +W800 +H600 +A0.2 -I$f3-[format %05d $j].pov -O$f3-[format %05d $j].tga}
	    puts "rendering done"
	}
    } 
}
