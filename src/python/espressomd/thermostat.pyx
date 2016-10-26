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
from __future__ import print_function, absolute_import
from . cimport thermostat
include "myconfig.pxi"
from globals cimport *
import numpy as np
from . cimport utils

cdef class Thermostat:

    def __init__(self):
        pass


    # __getstate__ and __setstate__ define the pickle interaction
    def __getstate__(self):
        # Attributes to pickle.
        thermolist = self.get_state()
        return thermolist

    def __setstate__(self, thermolist):
        if thermolist == []:
            return

        for thmst in thermolist:
            if thmst["type"] == "OFF":
                self.turn_off()
            if thmst["type"] == "LANGEVIN":
                self.set_langevin(kT=thmst["kT"], gamma=thmst["gamma"], gamma_rotation=thmst["gamma_rotation"])
            if thmst["type"] == "LB":
                self.set_lb(kT=thmst["kT"])
            if thmst["type"] == "NPT_ISO":
                self.set_npt(kT=thmst["kT"], p_diff=thmst["p_diff"], piston=thmst["piston"])
            if thmst["type"] == "DPD" or thmst["type"] == "INTER_DPD":
                pass

    def get_ts(self):
        return thermo_switch


    def get_state(self):
        """Returns the thermostat status."""
        thermo_list = []
        if temperature == -1:
            raise Exception("Thermostat is not initialized")
        if thermo_switch == THERMO_OFF:
            return [{"type": "OFF"}]
        if thermo_switch & THERMO_LANGEVIN:
            lang_dict = {}
            lang_dict["type"] = "LANGEVIN"
            lang_dict["kT"] = temperature
            lang_dict["gamma"] = langevin_gamma
            IF ROTATIONAL_INERTIA:
                lang_dict["gamma_rotation"] = [langevin_gamma_rotation[0],
                                               langevin_gamma_rotation[1],
                                               langevin_gamma_rotation[2]]
            ELSE:
                lang_dict["gamma_rotation"] = langevin_gamma_rotation

            thermo_list.append(lang_dict)
        if thermo_switch & THERMO_LB:
            lb_dict = {}
            lb_dict["type"] = "LB"
            lb_dict["kT"] = temperature
            thermo_list.append(lb_dict)
        if thermo_switch & THERMO_NPT_ISO:
            npt_dict = {}
            npt_dict["type"] = "NPT_ISO"
            npt_dict["kT"] = temperature
            npt_dict.update(nptiso)
            # thermo_dict["gamma0"] = nptiso_gamma0
            # thermo_dict["gammav"] = nptiso_gammav
            # thermo_dict["p_ext"] = nptiso.p_ext
            # thermo_dict["p_inst"] = nptiso.p_inst
            # thermo_dict["p_inst_av"] = nptiso.p_inst_av
            # thermo_dict["piston"] = nptiso.piston
            # thermo_dict["p_diff"] = nptiso.p_diff
            thermo_list.append(npt_dict)
        if (thermo_switch & THERMO_DPD) or (thermo_switch & THERMO_INTER_DPD):
            dpd_dict = {}
            if (thermo_switch & THERMO_DPD):
                dpd_dict["type"] = "DPD"
            else:
                dpd_dict["type"] = "INTER_DPD"
            dpd_dict["kT"] = temperature
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
        global langevin_gamma_rotation
        IF ROTATION:
            IF ROTATIONAL_INERTIA:
                for i in range(3):
                    langevin_gamma_rotation[i] = 0.
            ELSE:
                langevin_gamma_rotation = 0.
            mpi_bcast_parameter(FIELD_LANGEVIN_GAMMA_ROTATION)

        global thermo_switch
        thermo_switch = THERMO_OFF
        mpi_bcast_parameter(FIELD_THERMO_SWITCH)
        # here other thermostats stuff
        return True

    def set_langevin(self, kT="", gamma="", gamma_rotation=""):
        """Sets the Langevin thermostat with required parameters 'temperature' 'gamma'
        and optional parameter 'gamma_rotation'"""

        if kT == "" or gamma == "":
            raise ValueError(
                "Both, kT and gamma have to be given as keyword args")
        utils.check_type_or_throw_except(kT,1,float,"kT must be a number")
        utils.check_type_or_throw_except(gamma,1,float,"gamma must be a number")
        if float(kT) < 0. or float(gamma) < 0.:
            raise ValueError("temperature and gamma must be positive numbers")

        global langevin_gamma_rotation
        if gamma_rotation != "":
            IF ROTATIONAL_INERTIA:
                langevin_gamma_rotation[0] = gamma_rotation[0]
                langevin_gamma_rotation[1] = gamma_rotation[1]
                langevin_gamma_rotation[2] = gamma_rotation[2]
            ELSE:
                langevin_gamma_rotation = gamma_rotation

        global temperature
        temperature = float(kT)
        global langevin_gamma
        langevin_gamma = float(gamma)
        global thermo_switch
        thermo_switch = (thermo_switch | THERMO_LANGEVIN)
        mpi_bcast_parameter(FIELD_THERMO_SWITCH)
        mpi_bcast_parameter(FIELD_TEMPERATURE)
        mpi_bcast_parameter(FIELD_LANGEVIN_GAMMA)
        return True

    IF LB_GPU or LB:
        def set_lb(self, kT=""):
            """Sets the LB thermostat with required parameter 'temperature'"""

            if kT == "":
                raise ValueError(
                    "kT has to be given as keyword arg")
            utils.check_type_or_throw_except(kT,1,float,"kT must be a number")
            if float(kT) < 0.:
                raise ValueError("temperature must be non-negative")
            global temperature
            temperature = float(kT)
            global thermo_switch
            thermo_switch = (thermo_switch or THERMO_LB)
            mpi_bcast_parameter(FIELD_THERMO_SWITCH)
            mpi_bcast_parameter(FIELD_TEMPERATURE)
            return True

    IF NPT:
        def set_npt(self, kT="", gamma0="", gammav=""):
            """Sets the NPT thermostat with required parameters 'temperature' 'gamma0' 'gammav'"""
            if kT == "" or gamma0 == "" or gammav == "":
                raise ValueError(
                    "kT, gamma0 and gammav have to be given as keyword args")
            if not isinstance(kT, float):
                raise ValueError("temperature must be a positive number")
            global temperature
            temperature = float(kT)
            global thermo_switch
            thermo_switch = (thermo_switch | THERMO_NPT_ISO)
            global nptiso_gamma0
            nptiso_gamma0 = gamma0
            global nptiso_gammav
            nptiso_gammav = gammav
            mpi_bcast_parameter(FIELD_THERMO_SWITCH)
            mpi_bcast_parameter(FIELD_TEMPERATURE)
            mpi_bcast_parameter(FIELD_NPTISO_G0)
            mpi_bcast_parameter(FIELD_NPTISO_GV)
        
            
    IF DPD or INTER_DPD:
        def set_dpd(self, **kwargs):
            """Sets the DPD thermostat with required parameters 'kT' 'gamma' 'r_cut'"""
            req = ["kT","gamma","r_cut"]
            valid = ["kT","gamma","r_cut","tgamma","tr_cut","wf","twf"]
            raise Exception("Not implemented yet.")


