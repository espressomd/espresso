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
include "myconfig.pxi"
IF ELECTROSTATICS == 1:
    def setParams(kappa, rCut):
        if rCut < 0:
            raise ValueError("rCut must be > 0")
        dh_set_params(kappa, rCut)

    def setRcut(rCut):
        if rCut < 0:
            raise ValueError("rCut must be > 0")
        dh_set_params(dh_params.kappa, rCut)

    def setKappa(kappa):
        dh_set_params(kappa, dh_params.r_cut)

    def getRcut():
        return dh_params.r_cut

    def getKappa():
        return dh_params.kappa
