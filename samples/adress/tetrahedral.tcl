#  This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#  
#############################################################
#                                                           #
# tetrahedral.tcl                                           #
# ===========                                               #
#                                                           #
# Script to create the tetrahedral molecules                #
#                                                           #
#                                                           #
#############################################################



proc create_pyramide { i posx posy posz bond_length n_part_per_mol mass_ex } { 
# Creates the ith  molecule with the center of mass  at coordinates (<posx>, <posy>, <posz>) 
# Parameters:  bond_lentgth       =  bond length between beads
#              n_part_per_mol     =  number of particles per molecule
    
    set j [expr ($n_part_per_mol+1)*$i]
    
    
    set posx1 [expr $posx-$bond_length*0.5]
    set posy1 [expr $posy-$bond_length*1.73205080757/6.0]
    set posz1 [expr $posz-$bond_length/4.89897948557] 
    
    part $j pos $posx1 $posy1 $posz1 mass $mass_ex type 0 virtual 0

    set posx2 [expr $posx+$bond_length*0.5]
    set posy2 [expr $posy-$bond_length*1.73205080757/6.0]
    set posz2 [expr $posz-$bond_length/4.89897948557]

    part [expr $j+1] pos $posx2 $posy2 $posz2 mass $mass_ex type 0 virtual 0

    set posx3 [expr $posx+$bond_length*0.0]
    set posy3 [expr $posy+$bond_length*1.73205080757/3.0]
    set posz3 [expr $posz-$bond_length/4.89897948557] 

    part [expr $j+2] pos $posx3 $posy3 $posz3 mass $mass_ex type 0 virtual 0

    set posx4 [expr $posx+$bond_length*0.0]
    set posy4 [expr $posy+$bond_length*0.0]
    set posz4 [expr $posz+$bond_length*0.612372435688] 

    part [expr $j+3] pos $posx4 $posy4 $posz4 mass $mass_ex type 0 virtual 0
 
    

   for {set ii 0} { $ii < [expr $n_part_per_mol-1]} {incr ii} { 
       for {set k [expr $ii+1]} { $k < $n_part_per_mol} {incr k} { 
           part [expr $j+$ii] bond 0 [expr $j+$k]

       }
   }
}
