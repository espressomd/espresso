if { $argc != 1 } {
        puts "1 argument is expected."
        puts "desired velocity of the sphere (in x direction)"
        puts " "
        exit;
    } else {
	set velocity [expr [lindex $argv 0]]; # velocity of the sphere   
        }

# the filenames must be set
set createPart TMP/createPart;
set fileNodes  INPUT/SPHEREmesh-nodes.dat;
set fileTriangles INPUT/SPHEREmesh-triangles.dat;

set nnode 624;
set ntriangle 1244;

set createPart TMP/createPart;

set bondS TMP/bondsStretching
set bondB TMP/bondsBending
set bondAlocal TMP/bondsAreaLocal
set bondAglobal TMP/bondsAreaGlobal
set bondV TMP/bondsVolume
set bondVA TMP/bondsVolumeAreaGlobal
set partS TMP/partStretching
set partB TMP/partBending
set partAlocal TMP/partAreaLocal
set partAglobal TMP/partAreaGlobal
set partV TMP/partVolume
set partVA TMP/partVolumeAreaGlobal
set vmd "n"

set export_fluid "n"
set export_cell "n"
set export_cyl "n"

setmd time_step 0.1;    # timestep for the integrator
setmd skin 0.2;
thermostat off;
set pipeX 100;
set pipeY 20;
set pipeZ 20;
setmd box_l $pipeX $pipeY $pipeZ;  # rectangular channel is defined

set massPart 1.0;     # setting the mass of particles
set mytau 0.1;
set force 1;
set beginning 1;

set ks 1.0;          # stretching of the cell
set kb 0.0;          # bending
set kal 0.0;         # local area preservation
set kag 0.0;         # global area preservation
set kv 0.0;          # volume preservation
                      
set originX 10;       # the initial coordinates of the center
set originY 10;       # of the immersed object
set originZ 10; 

set stretchX 1.0;     # the immersed object will be scaled 
set stretchY 1.0;     # by these factors in corresponding directions
set stretchZ 1.0;

set rotateX 0.;       # Rotation by specified angles around 
set rotateY 0.;       # X axis, Y axis and Z axis. Angles given
set rotateZ 0.;       # in radians. rotateX=Pi/2 rotates the object 
                      # by 90 degrees with the axes of rotation x 
                      # such that vector 0,1,0 changes to 0,0,1 and
                      # 0,0,1 changes to 0,-1,0
                      
set typePart 0;       # each immersed object must have 
set molPart 0;        # different type and mol ID
	
set firstBondId 0;    # when adding more immersed objects, 
set firstPartId 0;    # these parameters need to be set properly

# the following script generates the source files
# for particles, interactions and bonds

exec ./bondGenerator $fileNodes $fileTriangles $nnode \
$ntriangle $rotateX $rotateY $rotateZ $originX $originY $originZ \
$stretchX $stretchY $stretchZ $createPart$molPart \
$typePart $molPart $firstPartId  $bondS$molPart \
$bondB$molPart $bondAlocal$molPart $bondAglobal$molPart \
$bondV$molPart $bondVA$molPart $partS$molPart $partB$molPart \
$partAlocal$molPart $partAglobal$molPart $partV$molPart $partVA$molPart\
$ks $kb $kal $kag $kv $firstBondId $massPart;

# create particles and bonds from source files

source $createPart$molPart; 
source $bondS$molPart; 
source $partS$molPart; 
source $bondAlocal$molPart; 
source $partAlocal$molPart; 
source $bondAglobal$molPart; 
source $partAglobal$molPart; 
source $bondB$molPart; 
source $partB$molPart; 
source $bondV$molPart; 
source $partV$molPart;

cellsystem domain_decomposition -no_verlet_list;      
lbfluid grid 1 dens 1.0 visc 1.5 tau $mytau friction 0.5;
                           
set upscale 12.;      # for scaling in the visualization

if { $vmd == "y" } {
    prepare_vmd_connection simEspresso 3000 1 $upscale;
                     #visualization
    exec sleep 0.5   
    imd positions $upscale
}

set out_file [open "OUT/sphere1.txt" "w"]
close $out_file;

set cycle 0 
set cur_time 0;
while { $cycle<100 } {
	set cur_time [expr $cur_time + $mytau];
	#puts "$cur_time";

        #if { $vmd == "y"} { imd positions $upscale; };
	set centerX 0.;
	set centerY 0.;
	set centerZ 0.;
        set v_sphereX 0.;
	set v_sphereY 0.;
	set v_sphereZ 0.;
	
	for { set pts 0} { $pts < $nnode } {incr pts} {
	    set aa [part $pts print pos];
	    set xcoord [lindex $aa 0];
	    set ycoord [lindex $aa 1];
	    set zcoord [lindex $aa 2];
	    set aa [part $pts print v];
	    set xvel [lindex $aa 0];
	    set yvel [lindex $aa 1];
	    set zvel [lindex $aa 2];
	    set centerX [expr $centerX + $xcoord];
	    set centerY [expr $centerY + $ycoord];
	    set centerZ [expr $centerZ + $zcoord];
	    set v_sphereX [expr $v_sphereX + $xvel];
	    set v_sphereY [expr $v_sphereY + $yvel];
	    set v_sphereZ [expr $v_sphereZ + $zvel];
		}
	set centerX [expr $centerX / (1.0*$nnode)];
	set centerY [expr $centerY / (1.0*$nnode)];
	set centerZ [expr $centerZ / (1.0*$nnode)];
	set v_sphereX [expr $v_sphereX / (1.0*$nnode)];
	set v_sphereY [expr $v_sphereY / (1.0*$nnode)];
	set v_sphereZ [expr $v_sphereZ / (1.0*$nnode)];
    
    if {$v_sphereX<$velocity} {
	if {$beginning == 1} {
	    for { set pts 0} { $pts < $nnode } {incr pts} {
		part $pts ext_force [expr $force*$cur_time] 0.0 0.0;
	    }
	}
    } else {
	for { set pts 0} { $pts < $nnode } {incr pts} {
	    part $pts ext_force 0.0 0.0 0.0;
	    set beginning 0;
	}
    }
 
    # setting the velocity
    # of the fluid inside the sphere
    for { set i [expr int(floor($centerX-4))+1] } { $i < [expr int(floor($centerX+4))] } { incr i } {
	for { set j [expr int(floor($centerY-4))+1] } { $j < [expr int(floor($centerY+4))] } { incr j } {
	    for { set k [expr int(floor($centerZ-4))+1] } { $k < [expr int(floor($centerZ+4))] } { incr k } {
		lbnode $i $j $k set u $v_sphereX $v_sphereY $v_sphereZ;
	    }
	}
   }

  if { $cycle % 5 == 0} {
    set out_file [open "OUT/sphere1.txt" "a"]
    puts $out_file "$cur_time $centerX $v_sphereX";
    close $out_file;
  }

  integrate 1;
  incr cycle;
}
