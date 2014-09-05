# description: all tcl-Tools to add a cell in LB simulations
# date: 17.02.14
# author: Vera
# functions:
# 	- addCell: add a normal cell 
# 	- addCellRotated: cell that can be rotated by 2 angles
# 	- add2ndCell: if you want more than 1 cell and not define all interactions new

# the following global variables are needed for all functions:

# global boxx boxy boxz
# global interCount cellCount nodeCount
# global kB kS maxStretch
# global numNodesPerCell
# global cellSeqFile

# ---------------------------------------------------------------
# add a Cell orientated in z direction 
# ---------------------------------------------------------------

proc addCell {originX originY originZ {num "0"} {conserve_vol "0"}} {
    cd ..
    # edited: 13.9.13
    # contains variable num to set type for virtual partikles of cell (needed to show more cells with paraview)  

    #declare global variables
    global boxx boxy boxz
    global interCount cellCount nodeCount
    global kB kS maxStretch
    global numNodesPerCell
    global cellSeqFile

    incr cellCount

    # read cell nodes
    set nodefile [open "immersed_boundary_system-nodes.data" "r"]
    gets $nodefile numNodes
    set numNodesPerCell $numNodes

    for {set i 0} {$i < $numNodes} {incr i} {
	gets $nodefile point
	set x [expr [lindex $point 0] + $originX ]
	set y [expr [lindex $point 1] + $originY ]
	set z [expr [lindex $point 2] + $originZ ]
	if {$conserve_vol == 1} { 
	    part [expr $i + $nodeCount] pos $x $y $z type $num virtual 1 molecule_id $cellCount} else { 
		part [expr $i + $nodeCount] pos $x $y $z type $num virtual 1
	    }
    }

}
close $nodefile

# open triangle file
set trianfile [open "immersed_boundary_system-triangles.data" "r"]
set intercount 0
gets $trianfile numTri

# read cell triangles & store in output file          
for {set i 0} {$i < $numTri } {incr i} {
    gets $trianfile tri
    lappend indices $tri
    # set second ks parameter to rubbish. Since we use Neo-Hookean, this does not matter  
    set ind1 [expr $nodeCount + [lindex $tri 0]]
    set ind2 [expr $nodeCount + [lindex $tri 1]]  
    set ind3 [expr $nodeCount + [lindex $tri 2]]  
    
    inter $intercount triel $ind1 $ind2 $ind3 $maxStretch $kS -1e6
    part $ind1 bond $intercount $ind2 $ind3
    incr intercount
    incr interCount
}

close $trianfile

# read bending information
set bendfile [open "immersed_boundary_system-bend.data" "r"]
gets $bendfile numBend

for {set i 0} {$i< $numBend } {incr i} {
    gets $bendfile bend

    set ind1 [expr $nodeCount + [lindex $bend 0]]
    set ind2 [expr $nodeCount + [lindex $bend 1]]
    set ind3 [expr $nodeCount + [lindex $bend 2]]
    set ind4 [expr $nodeCount + [lindex $bend 3]]
    
    inter $intercount tribend $ind1 $ind2 $ind3 $ind4 0 $kB $maxStretch
    part $ind1 bond $intercount $ind2 $ind3 $ind4
    incr intercount
    incr interCount
}

close $bendfile

# increase counter
incr nodeCount $numNodes

}


# --------------------------------------------------------------------------------
# add cells that can be rotated with 2 angles
# --------------------------------------------------------------------------------

proc addCellRotated {originX originY originZ {alpha "0"} {beta "0"}  {num "0"} } {

  # edited: 08.01.14
  # contains variable num to set type for virtual partikles of cell (needed to show more cells with paraview)  
  # contains variables alpha (rotation x-Axis) and beta (rotation new y-axis. Angles are entered in radial measurement. 

  #declare global variables
  global boxx boxy boxz
  global interCount cellCount nodeCount
  global kB kS maxStretch
  global numNodesPerCell
  global cellSeqFile

  incr cellCount

  # read cell nodes
  set nodefile [open "~/sphere/nodes" "r"]
  gets $nodefile numNodes
  set numNodesPerCell $numNodes
  
  # calculate Sin and cos values
  set alphaSin [expr sin($alpha)]
  set alphaCos [expr cos($alpha)]
  set betaCos [expr cos($beta)]
  set betaSin [expr sin($beta)]

  for {set i 0} {$i < $numNodes} {incr i} {
    gets $nodefile point
    # before the point is set to absolute coordinates in the box the cell will be rotated. 
    set xNode [lindex $point 0] 
    set yNode [lindex $point 1]
    set zNode [lindex $point 2]
  
    set x [expr $xNode*$betaCos + $yNode*$betaSin*$alphaSin + $zNode*$betaSin*$alphaCos + $originX ]
    set y [expr $yNode*$alphaCos - $zNode*$alphaSin + $originY ]
    set z [expr - $xNode*$betaSin + $yNode*$betaCos*$alphaSin + $zNode*$betaCos*$alphaCos + $originZ ]
    part [expr $i + $nodeCount] pos $x $y $z type $num virtual 1 molecule_id $cellCount
  }
  close $nodefile

  # open triangle file
  set trianfile [open "~/sphere/triangles" "r"]
  set intercount 0
  gets $trianfile numTri

  #open cell sequence file for output
  puts $cellSeqFile "#NumNodes = $numNodes" 
  puts $cellSeqFile "#NumTriangles = $numTri"
          
  # read cell triangles & store in output file          
  for {set i 0} {$i < $numTri } {incr i} {
    gets $trianfile tri
    lappend indices $tri
    # set second ks parameter to rubbish. Since we use Neo-Hookean, this does not matter  
    set ind1 [expr $nodeCount + [lindex $tri 0]]
    set ind2 [expr $nodeCount + [lindex $tri 1]]  
    set ind3 [expr $nodeCount + [lindex $tri 2]]  
  
    inter $intercount triel $ind1 $ind2 $ind3 $maxStretch $kS -1e6
    part $ind1 bond $intercount $ind2 $ind3
    incr intercount
    incr interCount
    puts $cellSeqFile "[lindex $tri 0] [lindex $tri 1] [lindex $tri 2]"
  }
  close $trianfile

  # read bending information
  set bendfile [open "~/sphere/bend" "r"]
  gets $bendfile numBend

  for {set i 0} {$i< $numBend } {incr i} {
    gets $bendfile bend

    set ind1 [expr $nodeCount + [lindex $bend 0]]
    set ind2 [expr $nodeCount + [lindex $bend 1]]
    set ind3 [expr $nodeCount + [lindex $bend 2]]
    set ind4 [expr $nodeCount + [lindex $bend 3]]
  
    inter $intercount tribend $ind1 $ind2 $ind3 $ind4 0 $kB $maxStretch
    part $ind1 bond $intercount $ind2 $ind3 $ind4
    incr intercount
    incr interCount
  }
        
  close $bendfile

  # increase counter
  incr nodeCount $numNodes

} 
  
# ---------------------------------------------------------------------------------
# add a second Cell if you already have a first. No rotated cells at the moment.
# ----------------------------------------------------------------------------------

proc add2ndCell {originX originY originZ {num "0"}} {

  # edited: 21.10.13
  # sets equal cells with different mass center  

  #declare global variables
  global boxx boxy boxz
  global interCount cellCount nodeCount
  global kB kS maxStretch
  global numNodesPerCell
  global cellSeqFile

  incr cellCount

  # read cell nodes
  set nodefile [open "~/sphere/nodes" "r"]
  gets $nodefile numNodes
  set numNodesPerCell $numNodes

  for {set i 0} {$i < $numNodes} {incr i} {
    gets $nodefile point
    set x [expr [lindex $point 0] + $originX ]
    set y [expr [lindex $point 1] + $originY ]
    set z [expr [lindex $point 2] + $originZ ]
        
    part [expr $i + $nodeCount] pos $x $y $z type $num virtual 1 molecule_id $cellCount
  }
  close $nodefile

  # open triangle file
  set trianfile [open "~/sphere/triangles" "r"]
  set intercount 0
  gets $trianfile numTri

  #open cell sequence file for output
  puts $cellSeqFile "#NumNodes = $numNodes" 
  puts $cellSeqFile "#NumTriangles = $numTri"
          
  # read cell triangles & store in output file          
  for {set i 0} {$i < $numTri } {incr i} {
    gets $trianfile tri
    lappend indices $tri
    # set second ks parameter to rubbish. Since we use Neo-Hookean, this does not matter  
    set ind1 [expr $nodeCount + [lindex $tri 0]]
    set ind2 [expr $nodeCount + [lindex $tri 1]]  
    set ind3 [expr $nodeCount + [lindex $tri 2]]  
    # here no interaction has to be set.
    part $ind1 bond $intercount $ind2 $ind3
    incr intercount
  
    puts $cellSeqFile "[lindex $tri 0] [lindex $tri 1] [lindex $tri 2]"
  }
  close $trianfile

  # read bending information
  set bendfile [open "~/sphere/bend" "r"]
  gets $bendfile numBend

  for {set i 0} {$i< $numBend } {incr i} {
    gets $bendfile bend

    set ind1 [expr $nodeCount + [lindex $bend 0]]
    set ind2 [expr $nodeCount + [lindex $bend 1]]
    set ind3 [expr $nodeCount + [lindex $bend 2]]
    set ind4 [expr $nodeCount + [lindex $bend 3]]
  
    #  inter $intercount tribend $ind1 $ind2 $ind3 $ind4 0 $kB $maxStretch
    part $ind1 bond $intercount $ind2 $ind3 $ind4
    incr intercount
    #  incr interCount
  }
        
  close $bendfile

  # increase counter
  incr nodeCount $numNodes

}                    
