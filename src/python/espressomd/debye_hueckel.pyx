#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
    def set_params(kappa, r_cut):
        if r_cut < 0:
            raise ValueError("r_cut must be > 0")
        dh_set_params(kappa, r_cut)

    def set_rcut(r_cut):
        if r_cut < 0:
            raise ValueError("r_cut must be > 0")
        dh_set_params(dh_params.kappa, r_cut)

    def set_kappa(kappa):
        dh_set_params(kappa, dh_params.r_cut)

    def get_rcut():
        return dh_params.r_cut

    def get_kappa():
        return dh_params.kappa
