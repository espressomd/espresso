# Copyright (C) 2014,2015,2016 The ESPResSo project
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

#-----------------------------------------------------------------------------------------------------
# define your own walls and boundaries here
#
# remember that 
# output_vtk_* writes boundary for visualisation later
# constraint sets up boundary for objects 
# and lbboundary sets up boundary for fluid

# wall - bottom
set corX 0; set corY 0; set corZ 0;
set aX $boxX; set aY 0; set aZ 0;
set bX 0; set bY $boxY; set bZ 0;
set cX 0; set cY 0; set cZ 1;
set rhomFile "output/wallbottom.vtk"
output_vtk_rhomboid $corX $corY $corZ $aX $aY $aZ $bX $bY $bZ $cX $cY $cZ $rhomFile 
constraint rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1 type 10 
lbboundary rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1

# wall - top
set corX 0; set corY 0; set corZ [expr $boxZ - 1];
set aX $boxX; set aY 0; set aZ 0;
set bX 0; set bY $boxY; set bZ 0;
set cX 0; set cY 0; set cZ 1;
set rhomFile "output/walltop.vtk"
output_vtk_rhomboid $corX $corY $corZ $aX $aY $aZ $bX $bY $bZ $cX $cY $cZ $rhomFile 
constraint rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1 type 10 
lbboundary rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1

# wall - front side
set corX 0; set corY 0; set corZ 1;
set aX $boxX; set aY 0; set aZ 0;
set bX 0; set bY 1; set bZ 0;
set cX 0; set cY 0; set cZ [expr $boxZ - 2];
set rhomFile "output/wallfront.vtk"
output_vtk_rhomboid $corX $corY $corZ $aX $aY $aZ $bX $bY $bZ $cX $cY $cZ $rhomFile 
constraint rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1 type 10 
lbboundary rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1

# wall - back side
set corX 0; set corY [expr $boxY - 1]; set corZ 1;
set aX $boxX; set aY 0; set aZ 0;
set bX 0; set bY 1; set bZ 0;
set cX 0; set cY 0; set cZ [expr $boxZ - 2];
set rhomFile "output/wallback.vtk"
output_vtk_rhomboid $corX $corY $corZ $aX $aY $aZ $bX $bY $bZ $cX $cY $cZ $rhomFile 
constraint rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1 type 10 
lbboundary rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1

# obstacle rhomboid1
set corX 10; set corY 1; set corZ 1;
set aX 8; set aY 0; set aZ 0;
set bX 0; set bY 4; set bZ 0;
set cX 0; set cY 0; set cZ 18;
set rhomFile "output/rhomboid1.vtk"
output_vtk_rhomboid $corX $corY $corZ $aX $aY $aZ $bX $bY $bZ $cX $cY $cZ $rhomFile 
constraint rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1 type 10 
lbboundary rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1

# obstacle cylinder1 
set cX 16; set cY 17; set cZ 10;
set nX 0; set nY 0; set nZ 1;
set L 9
set r 3
set cylFile "output/cylinder1.vtk"
set n 20
output_vtk_cylinder $cX $cY $cZ $nX $nY $nZ $r $L $n $cylFile
constraint cylinder center $cX $cY $cZ axis $nX $nY $nZ radius $r length $L direction 1 type 10 
lbboundary cylinder center $cX $cY $cZ axis $nX $nY $nZ radius $r length $L direction 1

# obstacle rhomboid2
set corX 25; set corY 1; set corZ 1;
set aX 5; set aY 0; set aZ 0;
set bX 0; set bY 20; set bZ 0;
set cX 0; set cY 0; set cZ 10;
set rhomFile "output/rhomboid2.vtk"
output_vtk_rhomboid $corX $corY $corZ $aX $aY $aZ $bX $bY $bZ $cX $cY $cZ $rhomFile 
constraint rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1 type 10 
lbboundary rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1

# obstacle cylinder2 
set cX 37; set cY 10; set cZ 10;
set nX 0; set nY 0; set nZ 1;
set L 9
set r 3
set cylFile "output/cylinder2.vtk"
set n 20
output_vtk_cylinder $cX $cY $cZ $nX $nY $nZ $r $L $n $cylFile
constraint cylinder center $cX $cY $cZ axis $nX $nY $nZ radius $r length $L direction 1 type 10 
lbboundary cylinder center $cX $cY $cZ axis $nX $nY $nZ radius $r length $L direction 1
