#!/bin/sh
# tricking...\
PLATFORM=$(uname -s); exec $PLATFORM/tcl_md $0 $*

if {$argc < 2} {
    mdvar N       500
    mdvar Nchains 2
    
    for {set i 0} {$i < 500} {incr i} {
	mdvar x $i [expr $i / 10000.0]
	mdvar y $i 0
	mdvar z $i 0
    }
    
    mdchain 0   0-200 1
    mdchain 1 200-400 2
} {
    puts "reading $argv[1]"
    source $argv[1]
}

# write data in rereadable form
puts "mdvar N [mdvar N]"
puts "mdvar Nchains [mdvar Nchains]"

puts "mdchain 0 [mdchain 0]"
puts "mdchain 1 [mdchain 1]"

for {set i 0} {$i < 500} {incr i} {
    puts "mdvar x $i [mdvar x $i]"
    puts "mdvar y $i [mdvar y $i]"
    puts "mdvar z $i [mdvar z $i]"
}

mdvar N 0
exit
