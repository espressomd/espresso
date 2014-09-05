### output flow field

proc printFlow {filename} {

  set flow [open $filename  "w"]
  global boxx
  global boxy
  global boxz
  # output at each grid point
  set numX $boxx
  set numY $boxy
  set numZ $boxz
  
  # header
  puts $flow "#VectorGridSingle"
  
  # number of points and coordinates
  set stepX [expr $boxx/double($numX)]
  puts $flow "#NumX = $numX"
  for {set i 0} {$i < $numX} {incr i} {
    puts $flow [expr $i*$stepX]
  }

  set stepY [expr $boxy/double($numY)]
  puts $flow "#NumY = $numY"
  for {set i 0} {$i < $numY} {incr i} {
    puts $flow [expr $i*$stepY]
  }
  
  set stepZ [expr $boxz/double($numZ)]
  puts $flow "#NumZ = $numZ"
  for {set i 0} {$i < $numZ} {incr i} {
    puts $flow [expr $i*$stepZ]
  }

  for {set i 0} {$i < $numX} {incr i} {
    for {set j 0} {$j < $numY} {incr j} {
      for {set k 0} {$k < $numZ} {incr k} {
        set xPos [expr $i*$stepX]
        set yPos [expr $j*$stepY]
        set zPos [expr $k*$stepZ]
        puts $flow [lbfluid print_interpolated_velocity $xPos $yPos $zPos]
      }
    }
  }
  
  close $flow
}