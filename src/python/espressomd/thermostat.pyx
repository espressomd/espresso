#
# Copyright (C) 2013-2022 The ESPResSo project
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
import functools
include "myconfig.pxi"
from . cimport utils
from .lb import HydrodynamicInteraction


def AssertThermostatType(*allowed_thermostats):
    """Assert that only a certain group of thermostats is active at a time.

    Decorator class to ensure that only specific combinations of thermostats
    can be activated together by the user. Usage:

    .. code-block:: cython

        cdef class Thermostat:
            @AssertThermostatType(THERMO_LANGEVIN, THERMO_DPD)
            def set_langevin(self, kT=None, gamma=None, gamma_rotation=None,
                     act_on_virtual=False, seed=None):
                ...

    This will prefix an assertion that prevents setting up the Langevin
    thermostat if the list of active thermostats contains anything other
    than the DPD and Langevin thermostats.

    Parameters
    ----------
    allowed_thermostats : :obj:`str`
        Allowed list of thermostats which are known to be compatible
        with one another.

    """
    def decoratorfunction(function):
        @functools.wraps(function, assigned=('__name__', '__doc__'))
        def wrapper(*args, **kwargs):
            if (not (thermo_switch in allowed_thermostats) and
                    (thermo_switch != THERMO_OFF)):
                raise Exception(
                    "This combination of thermostats is not allowed!")
            function(*args, **kwargs)
        return wrapper
    return decoratorfunction


cdef class Thermostat:

    # We have to cdef the state variable because it is a cdef class
    cdef _state
    cdef _LB_fluid

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

        If the thermostat had been suspended using :meth:`suspend`, it can
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
                self.set_langevin(kT=thmst["kT"], gamma=thmst["gamma"],
                                  gamma_rotation=thmst["gamma_rotation"],
                                  act_on_virtual=thmst["act_on_virtual"],
                                  seed=thmst["seed"])
                langevin_set_rng_counter(thmst["counter"])
            if thmst["type"] == "LB":
                self.set_lb(
                    LB_fluid=thmst["LB_fluid"],
                    act_on_virtual=thmst["act_on_virtual"],
                    gamma=thmst["gamma"],
                    seed=thmst["rng_counter_fluid"])
            if thmst["type"] == "NPT_ISO":
                if NPT:
                    self.set_npt(kT=thmst["kT"], gamma0=thmst["gamma0"],
                                 gammav=thmst["gammav"], seed=thmst["seed"])
                    npt_iso_set_rng_counter(thmst["counter"])
            if thmst["type"] == "DPD":
                if DPD:
                    self.set_dpd(kT=thmst["kT"], seed=thmst["seed"])
                    dpd_set_rng_counter(thmst["counter"])
            if thmst["type"] == "BROWNIAN":
                self.set_brownian(kT=thmst["kT"], gamma=thmst["gamma"],
                                  gamma_rotation=thmst["gamma_rotation"],
                                  act_on_virtual=thmst["act_on_virtual"],
                                  seed=thmst["seed"])
                brownian_set_rng_counter(thmst["counter"])
            if thmst["type"] == "SD":
                IF STOKESIAN_DYNAMICS:
                    self.set_stokesian(kT=thmst["kT"], seed=thmst["seed"])
                    stokesian_set_rng_counter(thmst["counter"])

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
            lang_dict["seed"] = langevin.rng_seed()
            lang_dict["counter"] = langevin.rng_counter()
            IF PARTICLE_ANISOTROPY:
                lang_dict["gamma"] = [langevin.gamma[0],
                                      langevin.gamma[1],
                                      langevin.gamma[2]]
            ELSE:
                lang_dict["gamma"] = langevin.gamma
            IF ROTATION:
                IF PARTICLE_ANISOTROPY:
                    lang_dict["gamma_rotation"] = [langevin.gamma_rotation[0],
                                                   langevin.gamma_rotation[1],
                                                   langevin.gamma_rotation[2]]
                ELSE:
                    lang_dict["gamma_rotation"] = langevin.gamma_rotation
            ELSE:
                lang_dict["gamma_rotation"] = None

            thermo_list.append(lang_dict)
        if thermo_switch & THERMO_BROWNIAN:
            lang_dict = {}
            lang_dict["type"] = "BROWNIAN"
            lang_dict["kT"] = temperature
            lang_dict["act_on_virtual"] = thermo_virtual
            lang_dict["seed"] = brownian.rng_seed()
            lang_dict["counter"] = brownian.rng_counter()
            IF PARTICLE_ANISOTROPY:
                lang_dict["gamma"] = [brownian.gamma[0],
                                      brownian.gamma[1],
                                      brownian.gamma[2]]
            ELSE:
                lang_dict["gamma"] = brownian.gamma
            IF ROTATION:
                IF PARTICLE_ANISOTROPY:
                    lang_dict["gamma_rotation"] = [brownian.gamma_rotation[0],
                                                   brownian.gamma_rotation[1],
                                                   brownian.gamma_rotation[2]]
                ELSE:
                    lang_dict["gamma_rotation"] = brownian.gamma_rotation
            ELSE:
                lang_dict["gamma_rotation"] = None

            thermo_list.append(lang_dict)
        if thermo_switch & THERMO_LB:
            lb_dict = {}
            lb_dict["LB_fluid"] = self._LB_fluid
            lb_dict["gamma"] = lb_lbcoupling_get_gamma()
            lb_dict["type"] = "LB"
            lb_dict["act_on_virtual"] = thermo_virtual
            lb_dict["rng_counter_fluid"] = lb_lbcoupling_get_rng_state()
            thermo_list.append(lb_dict)
        if thermo_switch & THERMO_NPT_ISO:
            if NPT:
                npt_dict = {}
                npt_dict["type"] = "NPT_ISO"
                npt_dict["kT"] = temperature
                npt_dict["seed"] = npt_iso.rng_seed()
                npt_dict["counter"] = npt_iso.rng_counter()
                npt_dict["gamma0"] = npt_iso.gamma0
                npt_dict["gammav"] = npt_iso.gammav
                thermo_list.append(npt_dict)
        if thermo_switch & THERMO_DPD:
            IF DPD:
                dpd_dict = {}
                dpd_dict["type"] = "DPD"
                dpd_dict["kT"] = temperature
                dpd_dict["seed"] = dpd.rng_seed()
                dpd_dict["counter"] = dpd.rng_counter()
                thermo_list.append(dpd_dict)
        if thermo_switch & THERMO_SD:
            IF STOKESIAN_DYNAMICS:
                sd_dict = {}
                sd_dict["type"] = "SD"
                sd_dict["kT"] = get_sd_kT()
                sd_dict["seed"] = stokesian.rng_seed()
                sd_dict["counter"] = stokesian.rng_counter()
                thermo_list.append(sd_dict)
        return thermo_list

    def _set_temperature(self, kT):
        mpi_set_temperature(kT)
        utils.handle_errors("Temperature change failed")

    def turn_off(self):
        """
        Turns off all the thermostat and sets all the thermostat variables to zero.

        """

        self._set_temperature(0.)
        mpi_set_thermo_virtual(True)
        IF PARTICLE_ANISOTROPY:
            mpi_set_langevin_gamma(utils.make_Vector3d((0., 0., 0.)))
            mpi_set_brownian_gamma(utils.make_Vector3d((0., 0., 0.)))
            IF ROTATION:
                mpi_set_langevin_gamma_rot(utils.make_Vector3d((0., 0., 0.)))
                mpi_set_brownian_gamma_rot(utils.make_Vector3d((0., 0., 0.)))
        ELSE:
            mpi_set_langevin_gamma(0.)
            mpi_set_brownian_gamma(0.)
            IF ROTATION:
                mpi_set_langevin_gamma_rot(0.)
                mpi_set_brownian_gamma_rot(0.)

        mpi_set_thermo_switch(THERMO_OFF)
        lb_lbcoupling_set_gamma(0.0)
        mpi_bcast_lb_particle_coupling()

    @AssertThermostatType(THERMO_LANGEVIN, THERMO_DPD)
    def set_langevin(self, kT, gamma, gamma_rotation=None,
                     act_on_virtual=False, seed=None):
        """
        Sets the Langevin thermostat.

        Parameters
        -----------
        kT : :obj:`float`
            Thermal energy of the simulated heat bath.
        gamma : :obj:`float`
            Contains the friction coefficient of the bath. If the feature
            ``PARTICLE_ANISOTROPY`` is compiled in, then ``gamma`` can be a list
            of three positive floats, for the friction coefficient in each
            cardinal direction.
        gamma_rotation : :obj:`float`, optional
            The same applies to ``gamma_rotation``, which requires the feature
            ``ROTATION`` to work properly. But also accepts three floats
            if ``PARTICLE_ANISOTROPY`` is also compiled in.
        act_on_virtual : :obj:`bool`, optional
            If ``True`` the thermostat will act on virtual sites, default is
            ``False``.
        seed : :obj:`int`
            Initial counter value (or seed) of the philox RNG.
            Required on first activation of the Langevin thermostat.
            Must be positive.

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
            if float(kT) < 0. or float(gamma[0]) < 0. or float(
                    gamma[1]) < 0. or float(gamma[2]) < 0.:
                raise ValueError(
                    "temperature and diagonal elements of the gamma tensor must be positive numbers")
        if gamma_rotation is not None:
            if scalar_gamma_rot_def:
                if float(gamma_rotation) < 0.:
                    raise ValueError(
                        "gamma_rotation must be positive number")
            else:
                if float(gamma_rotation[0]) < 0. or float(
                        gamma_rotation[1]) < 0. or float(gamma_rotation[2]) < 0.:
                    raise ValueError(
                        "diagonal elements of the gamma_rotation tensor must be positive numbers")

        # Seed is required if the RNG is not initialized
        if seed is None and langevin.is_seed_required():
            raise ValueError(
                "A seed has to be given as keyword argument on first activation of the thermostat")

        if seed is not None:
            utils.check_type_or_throw_except(
                seed, 1, int, "seed must be a positive integer")
            if seed < 0:
                raise ValueError("seed must be a positive integer")
            langevin_set_rng_seed(seed)

        self._set_temperature(kT)
        IF PARTICLE_ANISOTROPY:
            cdef utils.Vector3d gamma_vec
            if scalar_gamma_def:
                for i in range(3):
                    gamma_vec[i] = gamma
            else:
                gamma_vec = utils.make_Vector3d(gamma)
        IF ROTATION:
            IF PARTICLE_ANISOTROPY:
                cdef utils.Vector3d gamma_rot_vec
                if gamma_rotation is None:
                    # rotational gamma is translational gamma
                    gamma_rot_vec = gamma_vec
                else:
                    if scalar_gamma_rot_def:
                        for i in range(3):
                            gamma_rot_vec[i] = gamma_rotation
                    else:
                        gamma_rot_vec = utils.make_Vector3d(gamma_rotation)
            ELSE:
                if gamma_rotation is None:
                    # rotational gamma is translational gamma
                    gamma_rotation = gamma

        global thermo_switch
        mpi_set_thermo_switch(thermo_switch | THERMO_LANGEVIN)
        IF PARTICLE_ANISOTROPY:
            mpi_set_langevin_gamma(gamma_vec)
            IF ROTATION:
                mpi_set_langevin_gamma_rot(gamma_rot_vec)
        ELSE:
            mpi_set_langevin_gamma(gamma)
            IF ROTATION:
                mpi_set_langevin_gamma_rot(gamma_rotation)

        mpi_set_thermo_virtual(act_on_virtual)

    @AssertThermostatType(THERMO_BROWNIAN)
    def set_brownian(self, kT, gamma, gamma_rotation=None,
                     act_on_virtual=False, seed=None):
        """Sets the Brownian thermostat.

        Parameters
        -----------
        kT : :obj:`float`
            Thermal energy of the simulated heat bath.
        gamma : :obj:`float`
            Contains the friction coefficient of the bath. If the feature
            ``PARTICLE_ANISOTROPY`` is compiled in, then ``gamma`` can be a list
            of three positive floats, for the friction coefficient in each
            cardinal direction.
        gamma_rotation : :obj:`float`, optional
            The same applies to ``gamma_rotation``, which requires the feature
            ``ROTATION`` to work properly. But also accepts three floats
            if ``PARTICLE_ANISOTROPY`` is also compiled in.
        act_on_virtual : :obj:`bool`, optional
            If ``True`` the thermostat will act on virtual sites, default is
            ``False``.
        seed : :obj:`int`
            Initial counter value (or seed) of the philox RNG.
            Required on first activation of the Brownian thermostat.
            Must be positive.

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
            if float(kT) < 0. or float(gamma[0]) < 0. or float(
                    gamma[1]) < 0. or float(gamma[2]) < 0.:
                raise ValueError(
                    "temperature and diagonal elements of the gamma tensor must be positive numbers")
        if gamma_rotation is not None:
            if scalar_gamma_rot_def:
                if float(gamma_rotation) < 0.:
                    raise ValueError(
                        "gamma_rotation must be positive number")
            else:
                if float(gamma_rotation[0]) < 0. or float(
                        gamma_rotation[1]) < 0. or float(gamma_rotation[2]) < 0.:
                    raise ValueError(
                        "diagonal elements of the gamma_rotation tensor must be positive numbers")

        # Seed is required if the RNG is not initialized
        if seed is None and brownian.is_seed_required():
            raise ValueError(
                "A seed has to be given as keyword argument on first activation of the thermostat")

        if seed is not None:
            utils.check_type_or_throw_except(
                seed, 1, int, "seed must be a positive integer")
            if seed < 0:
                raise ValueError("seed must be a positive integer")
            brownian_set_rng_seed(seed)

        self._set_temperature(kT)
        IF PARTICLE_ANISOTROPY:
            cdef utils.Vector3d gamma_vec
            if scalar_gamma_def:
                for i in range(3):
                    gamma_vec[i] = gamma
            else:
                gamma_vec = utils.make_Vector3d(gamma)
        IF ROTATION:
            IF PARTICLE_ANISOTROPY:
                cdef utils.Vector3d gamma_rot_vec
                if gamma_rotation is None:
                    # rotational gamma is translational gamma
                    gamma_rot_vec = gamma_vec
                else:
                    if scalar_gamma_rot_def:
                        for i in range(3):
                            gamma_rot_vec[i] = gamma_rotation
                    else:
                        gamma_rot_vec = utils.make_Vector3d(gamma_rotation)
            ELSE:
                if gamma_rotation is None:
                    # rotational gamma is translational gamma
                    gamma_rotation = gamma

        global thermo_switch
        mpi_set_thermo_switch(thermo_switch | THERMO_BROWNIAN)

        IF PARTICLE_ANISOTROPY:
            mpi_set_brownian_gamma(gamma_vec)
            IF ROTATION:
                mpi_set_brownian_gamma_rot(gamma_rot_vec)
        ELSE:
            mpi_set_brownian_gamma(gamma)
            IF ROTATION:
                mpi_set_brownian_gamma_rot(gamma_rotation)

        mpi_set_thermo_virtual(act_on_virtual)

    @AssertThermostatType(THERMO_LB, THERMO_DPD)
    def set_lb(
        self,
        seed=None,
        act_on_virtual=True,
        LB_fluid=None,
            gamma=0.0):
        """
        Sets the LB thermostat.

        This thermostat requires the feature ``WALBERLA``.

        Parameters
        ----------
        LB_fluid : :class:`~espressomd.lb.LBFluidWalberla`
        seed : :obj:`int`
            Seed for the random number generator, required if kT > 0.
            Must be positive.
        act_on_virtual : :obj:`bool`, optional
            If ``True`` the thermostat will act on virtual sites (default).
        gamma : :obj:`float`
            Frictional coupling constant for the MD particle coupling.

        """
        if not isinstance(LB_fluid, HydrodynamicInteraction):
            raise ValueError(
                "The LB thermostat requires a LB / LBGPU instance as a keyword arg.")

        self._LB_fluid = LB_fluid
        if lb_lbfluid_get_kT() > 0.:
            if seed is None and lb_lbcoupling_is_seed_required():
                raise ValueError(
                    "seed has to be given as keyword arg")
            if seed is not None:
                utils.check_type_or_throw_except(
                    seed, 1, int, "seed must be a positive integer")
                if seed < 0:
                    raise ValueError("seed must be a positive integer")
                lb_lbcoupling_set_rng_state(seed)
                mpi_bcast_lb_particle_coupling()
        else:
            lb_lbcoupling_set_rng_state(0)
            mpi_bcast_lb_particle_coupling()

        global thermo_switch
        mpi_set_thermo_switch(thermo_switch | THERMO_LB)

        mpi_set_thermo_virtual(act_on_virtual)
        lb_lbcoupling_set_gamma(gamma)
        mpi_bcast_lb_particle_coupling()

    IF NPT:
        @AssertThermostatType(THERMO_NPT_ISO)
        def set_npt(self, kT, gamma0, gammav, seed=None):
            """
            Sets the NPT thermostat.

            Parameters
            ----------
            kT : :obj:`float`
                Thermal energy of the heat bath
            gamma0 : :obj:`float`
                Friction coefficient of the bath
            gammav : :obj:`float`
                Artificial friction coefficient for the volume fluctuations.
            seed : :obj:`int`
                Initial counter value (or seed) of the philox RNG.
                Required on first activation of the Langevin thermostat.
                Must be positive.

            """

            utils.check_type_or_throw_except(
                kT, 1, float, "kT must be a number")

            # Seed is required if the RNG is not initialized
            if seed is None and npt_iso.is_seed_required():
                raise ValueError(
                    "A seed has to be given as keyword argument on first activation of the thermostat")

            if seed is not None:
                utils.check_type_or_throw_except(
                    seed, 1, int, "seed must be a positive integer")
                if seed < 0:
                    raise ValueError("seed must be a positive integer")
                npt_iso_set_rng_seed(seed)

            self._set_temperature(kT)
            global thermo_switch
            mpi_set_thermo_switch(thermo_switch | THERMO_NPT_ISO)
            mpi_set_nptiso_gammas(gamma0, gammav)

    IF DPD:
        @AssertThermostatType(THERMO_DPD, THERMO_LANGEVIN, THERMO_LB)
        def set_dpd(self, kT, seed=None):
            """
            Sets the DPD thermostat with required parameters 'kT'.
            This also activates the DPD interactions.

            Parameters
            ----------
            kT : :obj:`float`
                Thermal energy of the heat bath.
            seed : :obj:`int`
                Initial counter value (or seed) of the philox RNG.
                Required on first activation of the DPD thermostat.
                Must be positive.

            """
            utils.check_type_or_throw_except(
                kT, 1, float, "kT must be a number")

            # Seed is required if the RNG is not initialized
            if seed is None and dpd.is_seed_required():
                raise ValueError(
                    "A seed has to be given as keyword argument on first activation of the thermostat")

            if seed is not None:
                utils.check_type_or_throw_except(
                    seed, 1, int, "seed must be a positive integer")
                if seed < 0:
                    raise ValueError("seed must be a positive integer")
                dpd_set_rng_seed(seed)

            self._set_temperature(kT)
            global thermo_switch
            mpi_set_thermo_switch(thermo_switch | THERMO_DPD)

    IF STOKESIAN_DYNAMICS:
        @AssertThermostatType(THERMO_SD)
        def set_stokesian(self, kT=None, seed=None):
            """
            Sets the SD thermostat with required parameters.

            This thermostat requires the feature ``STOKESIAN_DYNAMICS``.

            Parameters
            ----------
            kT : :obj:`float`, optional
                Temperature.
            seed : :obj:`int`, optional
                Seed for the random number generator

            """

            if (kT is None) or (kT == 0):
                set_sd_kT(0.0)
                if stokesian.is_seed_required():
                    stokesian_set_rng_seed(0)
            else:
                utils.check_type_or_throw_except(
                    kT, 1, float, "kT must be a float")
                set_sd_kT(kT)

                # Seed is required if the RNG is not initialized
                if seed is None and stokesian.is_seed_required():
                    raise ValueError(
                        "A seed has to be given as keyword argument on first activation of the thermostat")
                if seed is not None:
                    utils.check_type_or_throw_except(
                        seed, 1, int, "seed must be an integer")
                    if seed < 0:
                        raise ValueError("seed must be a positive integer")
                    stokesian_set_rng_seed(seed)

            global thermo_switch
            mpi_set_thermo_switch(thermo_switch | THERMO_SD)
