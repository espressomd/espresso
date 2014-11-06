#
# Copyright (C) 2013,2014 The ESPResSo project
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

def changeVolume(dNew, dir="xyz"):
    if dNew<0:
        raise ValueError("No negative lengths")
    if dir=="xyz":
        dNew=dNew**(1./3.)
        rescale_boxl(3, dNew)
    elif dir=="x":
        rescale_boxl(0, dNew)
    elif dir=="y":
        rescale_boxl(1, dNew)
    elif dir=="z":
        rescale_boxl(2, dNew)
    else:
        raise ValueError('Usage: changeVolume { <V_new> | <L_new> { "x" | "y" | "z" | "xyz" } }')


    
