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
from globals cimport *
cimport thermostat

cdef class Thermostat:

    def __init__(self):
        pass

    def get_status(self):
        """Returns the thermostat status."""
        thermo_list = []
        if temperature == -1:
            return "{ not initialized }"
        if thermo_switch == THERMO_OFF:
            return "{ off }"
        if thermo_switch and THERMO_LANGEVIN:
            lang_dict = {}
            lang_dict["type"] = "langevin"
            lang_dict["T"] = temperature
            lang_dict["gamma"] = langevin_gamma
            thermo_list.append(lang_dict)
        if thermo_switch and THERMO_NPT_ISO:
            npt_dict = {}
            npt_dict["type"] = "NPT_ISO"
            npt_dict["T"] = temperature
            npt_dict.update(nptiso)
            # thermo_dict["gamma0"] = nptiso_gamma0
            # thermo_dict["gammav"] = nptiso_gammav
            # thermo_dict["p_ext"] = nptiso.p_ext
            # thermo_dict["p_inst"] = nptiso.p_inst
            # thermo_dict["p_inst_av"] = nptiso.p_inst_av
            # thermo_dict["piston"] = nptiso.piston
            # thermo_dict["p_diff"] = nptiso.p_diff
            thermo_list.append(npt_dict)
        if (thermo_switch and THERMO_DPD) or (thermo_switch and THERMO_INTER_DPD):
            dpd_dict = {}
            if (thermo_switch and THERMO_DPD):
                dpd_dict["type"] = "DPD"
            else:
                dpd_dict["type"] = "INTER_DPD"
            dpd_dict["T"] = temperature
            dpd_dict["gamma"] = dpd_gamma
            dpd_dict["r_cut"] = dpd_r_cut
            dpd_dict["wf"] = dpd_wf
            dpd_dict["tgamma"] = dpd_tgamma
            dpd_dict["tr_cut"] = dpd_tr_cut
            dpd_dict["twf"] = dpd_twf
            thermo_list.append(dpd_dict)
        return thermo_list


    def turn_off(self):
        """Turns off all the thermostat and sets all the thermostat variables to zero"""

        global temperature
        temperature = 0.
        mpi_bcast_parameter(FIELD_TEMPERATURE)
        global langevin_gamma
        langevin_gamma = 0.
        mpi_bcast_parameter(FIELD_LANGEVIN_GAMMA)
        global thermo_switch
        thermo_switch = THERMO_OFF
        mpi_bcast_parameter(FIELD_THERMO_SWITCH)
        # here other thermostats stuff
        return True

    def set_langevin(self, kT="", gamma=""):
        """Sets the Langevin thermostat with required parameters 'temperature' 'gamma'"""

        if kT == "" or gamma == "":
            raise ValueError(
                "Both, kT and gamma have to be given as keyword args")
        if not isinstance(kT, float) or not isinstance(gamma, float) or float(kT) < 0. or float(gamma) < 0.:
            raise ValueError("temperature and gamma must be positive numbers")
        global temperature
        temperature = float(kT)
        global langevin_gamma
        langevin_gamma = float(gamma)
        global thermo_switch
        thermo_switch = (thermo_switch or THERMO_LANGEVIN)
        mpi_bcast_parameter(FIELD_THERMO_SWITCH)
        mpi_bcast_parameter(FIELD_TEMPERATURE)
        mpi_bcast_parameter(FIELD_LANGEVIN_GAMMA)
        return True

    def set_npt(self, kT="", p_diff="", piston=""):
        """Sets the NPT thermostat with required parameters 'temperature' 'p_diff' 'piston'"""
        if kT == "" or p_diff == "" or piston == "":
            raise ValueError(
                "kT, p_diff and piston have to be given as keyword args")
        if not isinstance(kT, float) or not isinstance(p_diff, piston) or float(kT) < 0. or float(piston) < 0.:
            raise ValueError("temperature and piston must be positive numbers")
        global temperature
        temperature = float(kT)
        global thermo_switch
        thermo_switch = (thermo_switch or THERMO_NPT_ISO)
        global npt_piston
        nptiso.piston = piston
        global npt_p_diff
        nptiso.p_diff = p_diff
        mpi_bcast_parameter(FIELD_THERMO_SWITCH)
        mpi_bcast_parameter(FIELD_TEMPERATURE)
        mpi_bcast_parameter(FIELD_NPTISO_PISTON)
        mpi_bcast_parameter(FIELD_NPTISO_PDIFF)
        
    def set_dpd(self, **kwargs):
        req = ["kT","gamma","r_cut"]
        valid = ["kT","gamma","r_cut","tgamma","tr_cut","wf","twf"]
        raise Exception("Not implemented yet.")


