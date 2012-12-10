# the filenames must be set
set inputdir INPUT
set tmpdir TMP


set createPart "$tmpdir/createPart"
set fileNodes  "$inputdir/TETRAmesh-nodes.dat"
set fileTriangles "$inputdir/TETRAmesh-triangles.dat"
                     # for example with an red blood cell
                     # use CELLmesh-nodes.dat 
                     # and CELLmesh-triangles.dat 

set bondS "$tmpdir/bondsStretching"
set bondB "$tmpdir/bondsBending"
set bondAlocal "$tmpdir/bondsAreaLocal"
set bondAglobal "$tmpdir/bondsAreaGlobal"
set bondV "$tmpdir/bondsVolume"
set bondVA "$tmpdir/bondsVolumeAreaGlobal"
set partS "$tmpdir/partStretching"
set partB "$tmpdir/partBending"
set partAlocal "$tmpdir/partAreaLocal"
set partAglobal "$tmpdir/partAreaGlobal"
set partV "$tmpdir/partVolume"
set partVA "$tmpdir/partVolumeAreaGlobal"
set vmd "y"





setmd time_step 0.1    
	# timestep for the integrator
setmd skin 0.2
thermostat off
setmd box_l 100 20 20
    # rectangular channel is defined


set ks 0.05
    # stretching of the cell
set kb 0.01
    # bending
set kal 0.01
    # local area preservation
set kag 0.01
    # global area preservation
set kv 10.
    # volume preservation

set massPart 1.0
    # setting the mass of praticles

set nnode 4
    # number of IB points on the surface
                      # of the immersed object
                      # for example with an red blood cell
                      # set 400
set ntriangle 4
    # number of triangles at the surface
    # of the immersed object
    # for example with an red blood cell
    # set 796
                      
set originX 10
set originY 10
set originZ 10; 
    # the initial coordinates of the center
    # of the immersed object

set stretchX 1.0
set stretchY 1.0
set stretchZ 1.0;
     # the immersed object will be scaled 
     # by these factors in correspomding directions

set rotateX 0.
set rotateY 0.
set rotateZ 0.
       # Rotation by specified angles around 
       # X axis, Y axis and Z axis. Angles given
       # in radians. rotateX=Pi/2 rotates the object 
       # by 90 degrees with the axes of rotation x 
       # such that vector 0,1,0 changes to 0,0,1 and
       # 0,0,1 changes to 0,-1,0
                      
set typePart 0
set molPart 0
       # each immersed object must have 	
       # different type and mol ID

set firstBondId 0
set firstPartId 0
    # when adding more immersed objects, 
    # these parameters need to be set properly
    
# the following script generates the source files
# for particles, interactions and bonds
file mkdir TMP

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

cellsystem domain_decomposition -no_verlet_list
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
