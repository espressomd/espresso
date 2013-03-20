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

# initialization of the object-in-fluid mechanisms
############################################################################
init_objects_in_fluid	


# adding one object in fluid - oif. Some parameters are self-explanatory, some not:
############################################################################
add_oif_object origin 10 10 10 nodesfile $fileNodes trianglesfile $fileTriangles stretch 1.0 1.0 1.0 ks 0.05 kb 0.01 kal 0.01 kag 0.01 kv 10.0 type 0 mol 0 rotate 0.0 0.0 0.0


# origin 			sets the coordinates of the center of the object
# stretch 			the immersed object will be scaled 
# 					by these factors in correspomding directions
# nodesfile 		meshfile for vertices
# trianglesfile 	meshfile for triangles
#
# elastic parameters of the object:
# ks	stretching of the cell
# kb	bending
# kal	local area preservation
# kag	global area preservation
# kv	volume preservation
#
# rotate	Rotation by specified angles around 
# 			X axis, Y axis and Z axis. Angles given
# 			in radians. rotateX=Pi/2 rotates the object 
# 			by 90 degrees with the axes of rotation x 
# 			such that vector 0,1,0 changes to 0,0,1 and
# 			0,0,1 changes to 0,-1,0

# type		each immersed object must have 	
# mol			different type and mol ID


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
while { $cycle<200 } {
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
