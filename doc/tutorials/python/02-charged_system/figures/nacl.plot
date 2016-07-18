# Copyright (C) 2012,2013,2014,2015,2016 The ESPResSo project
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
set term post eps enh "Times-Roman" 25
set out "nacl-rdf.eps"
set xlabel "r"
set ylabel "g(r)"
plot [.9:0.5*6.75] \
	"data/rdf_from_melt_00.data" notitle w linesp pt 4, \
	"data/rdf_from_melt_10.data" notitle w linesp pt 6, \
	"data/rdf_lj_00.data"        notitle w linesp pt 1
unset out

set out "neutral-rho.eps"
set xlabel "z"
set ylabel "{/Symbol r}(z)"
plot [0.8:3.5] \
	"data/neutral-rho.data" u 1:2 notitle w linesp pt 4, \
	"data/neutral-rho.data" u 1:3 notitle w linesp pt 5
unset out

set out "nonneutral-rho.eps"
set xlabel "z"
set ylabel "{/Symbol r}(z)"
plot [0.8:3.5] \
	"data/nonneutral-rho.data" u 1:2 notitle w linesp pt 4, \
	"data/nonneutral-rho.data" u 1:3 notitle w linesp pt 5
unset out
