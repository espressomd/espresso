# general parameters
############################################################################

# whether to view the simulation on VMD
set vmd "y"

# integrator settings for the simulation
setmd time_step 0.1    
setmd skin 0.2
thermostat off

# rectangular channel box geometry
setmd box_l 100 20 20

# parameters for the object-in-fluid code
############################################################################

# stretching of the cell
set ks 0.05
# bending
set kb 0.01
# local area preservation
set kal 0.01
# global area preservation
set kag 0.01
# volume preservation
set kv 10.
# mass of mess nodes
set massPart 1.0

# the initial coordinates of the center
# of the immersed object
set originX 10
set originY 10
set originZ 10; 

# the immersed object will be scaled 
# by these factors in correspomding directions
set stretchX 1.0
set stretchY 1.0
set stretchZ 1.0;

# Rotation by specified angles around 
# X axis, Y axis and Z axis. Angles given
# in radians. rotateX=Pi/2 rotates the object 
# by 90 degrees with the axes of rotation x 
# such that vector 0,1,0 changes to 0,0,1 and
# 0,0,1 changes to 0,-1,0
set rotateX 0.
set rotateY 0.
set rotateZ 0.
                      
# each immersed object must have 	
# different type and mol ID
set typePart 0
set molPart 0

# make sure that when adding several objects, that these increase accordingly
set firstBondId 0
set firstPartId 0

# what files to read/generate where
############################################################################

# input files, describe the object shape
set inputdir "input"
# change these to one of the example files
# currently there are CELL and TETRA
set type [lindex $argv 0]
if {$type == ""} {
    set type "TETRA"
}
set fileNodes  "$inputdir/${type}mesh-nodes.dat"
set fileTriangles "$inputdir/${type}mesh-triangles.dat"

# temporary mesh descriptions
set tmpdir   "tmp"
set createPart  "$tmpdir/createPart"
set bondS       "$tmpdir/bondsStretching"
set bondB       "$tmpdir/bondsBending"
set bondAlocal  "$tmpdir/bondsAreaLocal"
set bondAglobal "$tmpdir/bondsAreaGlobal"
set bondV       "$tmpdir/bondsVolume"
set bondVA      "$tmpdir/bondsVolumeAreaGlobal"
set partS       "$tmpdir/partStretching"
set partB       "$tmpdir/partBending"
set partAlocal  "$tmpdir/partAreaLocal"
set partAglobal "$tmpdir/partAreaGlobal"
set partV       "$tmpdir/partVolume"
set partVA      "$tmpdir/partVolumeAreaGlobal"
    
# get the number of nodes and triangles
############################################################################

proc getnlines {file} {
    set f [open $file "r"]
    set lines 0
    while {![eof $f]} { gets $f; incr lines }
    close $f
    incr lines -1
    return $lines
}

# number of IB points on the surface
# of the immersed object
set nnode [getnlines $fileNodes]
# number of triangles at the surface
# of the immersed object
set ntriangle [getnlines $fileTriangles]

puts "nnode $nnode ntriangle $ntriangle"
file mkdir $tmpdir

# the following script generates the source files
# for particles, interactions and bonds
############################################################################

exec ./bondGenerator $fileNodes $fileTriangles $nnode \
    $ntriangle $rotateX $rotateY $rotateZ $originX $originY $originZ \
    $stretchX $stretchY $stretchZ $createPart$molPart \
    $typePart $molPart $firstPartId  $bondS$molPart \
    $bondB$molPart $bondAlocal$molPart $bondAglobal$molPart \
    $bondV$molPart $bondVA$molPart $partS$molPart $partB$molPart \
    $partAlocal$molPart $partAglobal$molPart $partV$molPart $partVA$molPart\
    $ks $kb $kal $kag $kv $firstBondId $massPart;

# create particles and bonds from source files

source $createPart$molPart
source $bondS$molPart
source $partS$molPart 
source $bondAlocal$molPart
source $partAlocal$molPart 
source $bondAglobal$molPart 
source $partAglobal$molPart 
source $bondB$molPart 
source $partB$molPart 
source $bondV$molPart 
source $partV$molPart

# run it!
############################################################################

lbfluid grid 1 dens 1.0 visc 1.5 tau 0.1 friction 0.5

if { $vmd == "y" } {
    prepare_vmd_connection simEspresso 3000 1
    exec sleep 2   
    imd positions
}

# main iteration loop

set cycle 0 
while { $cycle<1000 } {
    puts "$cycle"
    if { $vmd == "y"} { imd positions }

    # setting the constant velocity
    # of the fluid on the left side of the md_box
    for { set i 0 } { $i < 1} { incr i } {
        for { set j 0 } { $j < 20 } { incr j } {
            for { set k 0 } { $k < 20 } { incr k } {
                lbnode $i $j $k set u 0.5 0.0 0.0
            }
        }
    }
    integrate 1
    incr cycle
}
