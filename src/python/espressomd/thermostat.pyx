#
# Copyright (C) 2013-2018 The ESPResSo project
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
from functools import wraps
from . cimport thermostat
include "myconfig.pxi"
from globals cimport *
import numpy as np
from . cimport utils
from .lb cimport *
from .lb import HydrodynamicInteraction
from .lb cimport lb_lbcoupling_set_gamma
from .lb cimport lb_lbcoupling_get_gamma


def AssertThermostatType(*allowedthermostats):
    """Assert that only a certain thermostat is active

    Decorator class to assure that only a given thermostat is active
    at a time.  Usage:

        @AssertThermostatType(THERMO_LANGEVIN)
        def set_langevin(self, kT=None, gamma=None, gamma_rotation=None):

    This will prefix an assertion for THERMO_LANGEVIN to the call.

    """
    def decoratorfunction(function):
        @wraps(function, assigned=('__name__', '__doc__'))
        def wrapper(*args, **kwargs):
            if (not (thermo_switch in allowedthermostats) and
                    (thermo_switch != THERMO_OFF)):
                raise Exception(
                    "This combination of thermostats is not allowed!")
            function(*args, **kwargs)
        return wrapper
    return decoratorfunction


cdef class Thermostat(object):

    # We have to cdef the state variable because it is a cdef class
    cdef _state

    def __init__(self):
        self._state = None
        pass

    def suspend(self):
        """Suspend the thermostat

        The thermostat can be suspended, e.g. to perform an energy
        minimization.

        """
        self._state = self.__getstate__()
        self.turn_off()

    def recover(self):
        """Recover a suspended thermostat

        If the thermostat had been suspended using .suspend(), it can
        be recovered with this method.

        """
        if self._state is not None:
            self.__setstate__(self._state)

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
                self.set_langevin(kT=thmst["kT"], gamma=thmst[
                                  "gamma"], gamma_rotation=thmst["gamma_rotation"], act_on_virtual=thmst["act_on_virtual"], seed=thmst["seed"])
            if thmst["type"] == "LB":
                self.set_lb(
                    act_on_virtual=thmst["act_on_virtual"],
                    seed=thmst["rng_counter_fluid"])
            if thmst["type"] == "NPT_ISO":
                self.set_npt(kT=thmst["kT"], p_diff=thmst[
                             "p_diff"], piston=thmst["piston"])
            if thmst["type"] == "DPD":
                self.set_dpd(kT=thmst["kT"])
            IF BROWNIAN_DYNAMICS:
                if thmst["type"] == "BROWNIAN":
                    self.set_brownian(kT=thmst["kT"], gamma=thmst[
                                      "gamma"], gamma_rotation=thmst["gamma_rotation"],act_on_virtual=thmst["act_on_virtual"], seed=thmst["seed"])

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
            lang_dict["act_on_virtual"] = thermo_virtual
            lang_dict["seed"] = int(langevin_get_rng_state())
            IF PARTICLE_ANISOTROPY:
                lang_dict["gamma"] = [langevin_gamma[0],
                                      langevin_gamma[1],
                                      langevin_gamma[2]]
            ELSE:
                lang_dict["gamma"] = langevin_gamma
            IF ROTATION:
                IF PARTICLE_ANISOTROPY:
                    lang_dict["gamma_rotation"] = [langevin_gamma_rotation[0],
                                                   langevin_gamma_rotation[1],
                                                   langevin_gamma_rotation[2]]
                ELSE:
                    lang_dict["gamma_rotation"] = langevin_gamma_rotation
            ELSE:
                lang_dict["gamma_rotation"] = None

            thermo_list.append(lang_dict)
        IF BROWNIAN_DYNAMICS:
            if thermo_switch & THERMO_BROWNIAN:
                lang_dict = {}
                lang_dict["type"] = "BROWNIAN"
                lang_dict["kT"] = temperature
                lang_dict["act_on_virtual"] = thermo_virtual
                lang_dict["seed"] = int(langevin_get_rng_state())
                IF PARTICLE_ANISOTROPY:
                    lang_dict["gamma"] = [langevin_gamma[0],
                                          langevin_gamma[1],
                                          langevin_gamma[2]]
                ELSE:
                    lang_dict["gamma"] = langevin_gamma
                IF ROTATION:
                    IF PARTICLE_ANISOTROPY:
                        lang_dict["gamma_rotation"] = [langevin_gamma_rotation[0],
                                                       langevin_gamma_rotation[1],
                                                       langevin_gamma_rotation[2]]
                    ELSE:
                        lang_dict["gamma_rotation"] = langevin_gamma_rotation
                ELSE:
                    lang_dict["gamma_rotation"] = None

                thermo_list.append(lang_dict)
        if thermo_switch & THERMO_LB:
            lb_dict = {}
            lb_dict["gamma"] = lb_lbcoupling_get_gamma()
            lb_dict["type"] = "LB"
            lb_dict["act_on_virtual"] = thermo_virtual
            lb_dict["rng_counter_fluid"] = lb_lbcoupling_get_rng_state()
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
        if (thermo_switch & THERMO_DPD):
            dpd_dict = {}
            dpd_dict["type"] = "DPD"
            dpd_dict["kT"] = temperature
            thermo_list.append(dpd_dict)
        return thermo_list

    def turn_off(self):
        """
        Turns off all the thermostat and sets all the thermostat variables to zero.

        """

        global temperature
        global thermo_virtual
        thermo_virtual = True
        temperature = 0.
        mpi_bcast_parameter(FIELD_THERMO_VIRTUAL)
        mpi_bcast_parameter(FIELD_TEMPERATURE)
        global langevin_gamma
        IF PARTICLE_ANISOTROPY:
            for i in range(3):
                langevin_gamma[i] = 0.
        ELSE:
            langevin_gamma = 0.
        mpi_bcast_parameter(FIELD_LANGEVIN_GAMMA)
        global langevin_gamma_rotation
        IF ROTATION:
            IF PARTICLE_ANISOTROPY:
                for i in range(3):
                    langevin_gamma_rotation[i] = 0.
            ELSE:
                langevin_gamma_rotation = 0.
            mpi_bcast_parameter(FIELD_LANGEVIN_GAMMA_ROTATION)

        global thermo_switch
        thermo_switch = THERMO_OFF
        mpi_bcast_parameter(FIELD_THERMO_SWITCH)
        lb_lbcoupling_set_gamma(0.0)
        return True

    @AssertThermostatType(THERMO_LANGEVIN)
    def set_langevin(self, kT=None, gamma=None, gamma_rotation=None,
                     act_on_virtual=False, seed=None):
        """
        Sets the Langevin thermostat with required parameters 'kT' 'gamma'
        and optional parameter 'gamma_rotation'.

        Parameters
        -----------
        kT : :obj:`float`
             Thermal energy of the simulated heat bath.
        gamma : :obj:`float`
                Contains the friction coefficient of the bath. If the feature 'PARTICLE_ANISOTROPY'
                is compiled in then 'gamma' can be a list of three positive floats, for the friction
                coefficient in each cardinal direction.
        gamma_rotation : :obj:`float`, optional
                         The same applies to 'gamma_rotation', which requires the feature
                         'ROTATION' to work properly. But also accepts three floating point numbers
                         if 'PARTICLE_ANISOTROPY' is also compiled in.
        act_on_virtual : :obj:`bool`, optional
                If true the thermostat will act on virtual sites, default is off.
        seed : :obj:`int`, required
                Initial counter value (or seed) of the philox RNG.
                Required on first activation of the langevin thermostat.

        """

        scalar_gamma_def = True
        scalar_gamma_rot_def = True
        IF PARTICLE_ANISOTROPY:
            if hasattr(gamma, "__iter__"):
                scalar_gamma_def = False
            else:
                scalar_gamma_def = True

        IF PARTICLE_ANISOTROPY:
            if hasattr(gamma_rotation, "__iter__"):
                scalar_gamma_rot_def = False
            else:
                scalar_gamma_rot_def = True

        if kT is None or gamma is None:
            raise ValueError(
                "Both, kT and gamma have to be given as keyword args")
        utils.check_type_or_throw_except(kT, 1, float, "kT must be a number")
        if scalar_gamma_def:
            utils.check_type_or_throw_except(
                gamma, 1, float, "gamma must be a number")
        else:
            utils.check_type_or_throw_except(
                gamma, 3, float, "diagonal elements of the gamma tensor must be numbers")
        if gamma_rotation is not None:
            if scalar_gamma_rot_def:
                utils.check_type_or_throw_except(
                    gamma_rotation, 1, float, "gamma_rotation must be a number")
            else:
                utils.check_type_or_throw_except(
                    gamma_rotation, 3, float, "diagonal elements of the gamma_rotation tensor must be numbers")

        if scalar_gamma_def:
            if float(kT) < 0. or float(gamma) < 0.:
                raise ValueError(
                    "temperature and gamma must be positive numbers")
        else:
            if float(kT) < 0. or float(gamma[0]) < 0. or float(gamma[1]) < 0. or float(gamma[2]) < 0.:
                raise ValueError(
                    "temperature and diagonal elements of the gamma tensor must be positive numbers")
        if gamma_rotation is not None:
            if scalar_gamma_rot_def:
                if float(gamma_rotation) < 0.:
                    raise ValueError(
                        "gamma_rotation must be positive number")
            else:
                if float(gamma_rotation[0]) < 0. or float(gamma_rotation[1]) < 0. or float(gamma_rotation[2]) < 0.:
                    raise ValueError(
                        "diagonal elements of the gamma_rotation tensor must be positive numbers")

        #Seed is required if the rng is not initialized
        if not seed and langevin_is_seed_required():
            raise ValueError(
                "A seed has to be given as keyword argument on first activation of the thermostat")

        if seed:
            utils.check_type_or_throw_except(
                seed, 1, int, "seed must be a positive integer")
            langevin_set_rng_state(seed)

        global temperature
        temperature = float(kT)
        global langevin_gamma
        IF PARTICLE_ANISOTROPY:
            if scalar_gamma_def:
                langevin_gamma[0] = gamma
                langevin_gamma[1] = gamma
                langevin_gamma[2] = gamma
            else:
                langevin_gamma[0] = gamma[0]
                langevin_gamma[1] = gamma[1]
                langevin_gamma[2] = gamma[2]
        ELSE:
            langevin_gamma = float(gamma)

        global langevin_gamma_rotation
        IF ROTATION:
            if gamma_rotation is not None:
                IF PARTICLE_ANISOTROPY:
                    if scalar_gamma_rot_def:
                        langevin_gamma_rotation[0] = gamma_rotation
                        langevin_gamma_rotation[1] = gamma_rotation
                        langevin_gamma_rotation[2] = gamma_rotation
                    else:
                        langevin_gamma_rotation[0] = gamma_rotation[0]
                        langevin_gamma_rotation[1] = gamma_rotation[1]
                        langevin_gamma_rotation[2] = gamma_rotation[2]
                ELSE:
                    if scalar_gamma_rot_def:
                        langevin_gamma_rotation = gamma_rotation
                    else:
                        raise ValueError(
                            "gamma_rotation must be a scalar since feature PARTICLE_ANISOTROPY is disabled")
            else:
                IF PARTICLE_ANISOTROPY:
                    if scalar_gamma_def:
                        langevin_gamma_rotation[0] = gamma
                        langevin_gamma_rotation[1] = gamma
                        langevin_gamma_rotation[2] = gamma
                    else:
                        langevin_gamma_rotation[0] = gamma[0]
                        langevin_gamma_rotation[1] = gamma[1]
                        langevin_gamma_rotation[2] = gamma[2]
                ELSE:
                    langevin_gamma_rotation = langevin_gamma

        global thermo_switch
        thermo_switch = (thermo_switch | THERMO_LANGEVIN)
        mpi_bcast_parameter(FIELD_THERMO_SWITCH)
        mpi_bcast_parameter(FIELD_TEMPERATURE)
        mpi_bcast_parameter(FIELD_LANGEVIN_GAMMA)

        global thermo_virtual
        thermo_virtual = act_on_virtual
        mpi_bcast_parameter(FIELD_THERMO_VIRTUAL)

        IF ROTATION:
            mpi_bcast_parameter(FIELD_LANGEVIN_GAMMA_ROTATION)
        return True

    @AssertThermostatType(THERMO_LB)
    def set_lb(
        self,
        seed=None,
        act_on_virtual=True,
        LB_fluid=None,
            gamma=0.0):
        """
        Sets the LB thermostat.

        This thermostat requires the feature LBFluid or LBFluidGPU.

        Parameters
        ----------
        LB_fluid : instance of :class:`espressomd.lb.LBFluid` or :class:`espressomd.lb.LBFluidGPU`
        seed : :obj:`int`
             Seed for the random number generator, required
             if kT > 0.
        act_on_virtual : :obj:`bool`, optional
            If true the thermostat will act on virtual sites, default is on.
        gamma : :obj:`float`
            Frictional coupling constant for the MD particle coupling.

        """
        if not isinstance(LB_fluid, HydrodynamicInteraction):
            raise ValueError(
                "The LB thermostat requires a LB / LBGPU instance as a keyword arg.")

        if lb_lbfluid_get_kT() > 0.:
            if not seed and lb_lbcoupling_is_seed_required():
                raise ValueError(
                    "seed has to be given as keyword arg")
            elif seed:
                lb_lbcoupling_set_rng_state(seed)

        global thermo_switch
        thermo_switch = (thermo_switch or THERMO_LB)
        mpi_bcast_parameter(FIELD_THERMO_SWITCH)

        global thermo_virtual
        thermo_virtual = act_on_virtual
        mpi_bcast_parameter(FIELD_THERMO_VIRTUAL)
        lb_lbcoupling_set_gamma(gamma)
        return True

    IF NPT:
        @AssertThermostatType(THERMO_NPT_ISO)
        def set_npt(self, kT=None, gamma0=None, gammav=None):
            """
            Sets the NPT thermostat with required parameters 'temperature', 'gamma0', 'gammav'.

            Parameters
            ----------
            kT : :obj:`float`
                 Thermal energy of the heat bath
            gamma0 : :obj:`float`
                     Friction coefficient of the bath
            gammav : :obj:`float`
                     Artificial friction coefficient for the volume
                     fluctuations. Mass of the artificial piston.

            """

            if kT is None or gamma0 is None or gammav is None:
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

    IF DPD:
        def set_dpd(self, kT=None):
            """
            Sets the DPD thermostat with required parameters 'kT'.
            This also activates the DPD interactions.

            Parameters
            ----------
            'kT' : float
                Thermal energy of the heat bath, floating point number

            """
            if kT is None:
                raise ValueError(
                    "kT has to be given as keyword args")
            if not isinstance(kT, float):
                raise ValueError("temperature must be a positive number")
            global temperature
            temperature = float(kT)
            global thermo_switch
            thermo_switch = (thermo_switch | THERMO_DPD)

            mpi_bcast_parameter(FIELD_THERMO_SWITCH)
            mpi_bcast_parameter(FIELD_TEMPERATURE)

    IF BROWNIAN_DYNAMICS:
        @AssertThermostatType(THERMO_BROWNIAN)
        def set_brownian(self, kT=None, gamma=None, gamma_rotation=None,
                         act_on_virtual=False, seed=None):
            """Sets the Brownian Dynamics thermostat with required parameters 'kT' 'gamma'
            and optional parameter 'gamma_rotation'.
    
            Parameters
            -----------
            kT : :obj:`float`
             Thermal energy of the simulated heat bath.
            gamma : :obj:`float`
                    Contains the friction coefficient of the bath. If the feature 'PARTICLE_ANISOTROPY'
                    is compiled in then 'gamma' can be a list of three positive floats, for the friction
                    coefficient in each cardinal direction.
            gamma_rotation : :obj:`float`, optional
                             The same applies to 'gamma_rotation', which requires the feature
                             'ROTATION' to work properly. But also accepts three floating point numbers
                             if 'PARTICLE_ANISOTROPY' is also compiled in.
            act_on_virtual : :obj:`bool`, optional
                    If true the thermostat will act on virtual sites, default is off.
            seed : :obj:`int`, required
                    Initial counter value (or seed) of the philox RNG.
                    Required on first activation of the langevin thermostat.
    
            """

            self.set_langevin(kT, gamma, gamma_rotation, act_on_virtual, seed)
            global thermo_switch
            # this is safe because this combination of thermostats is not allowed
            thermo_switch = (thermo_switch & (~THERMO_LANGEVIN))
            thermo_switch = (thermo_switch | THERMO_BROWNIAN)
            mpi_bcast_parameter(FIELD_THERMO_SWITCH)
            return True
