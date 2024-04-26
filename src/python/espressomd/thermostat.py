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

from .script_interface import script_interface_register, ScriptInterfaceHelper


@script_interface_register
class Thermostat(ScriptInterfaceHelper):
    """
    Container for the system thermostat. Multiple thermostats can be active
    simultaneously. When setting up a thermostat raises an exception, e.g.
    due to invalid parameters, the thermostat gets disabled. Therefore,
    updating an already active thermostat with invalid parameters will
    deactivate it.

    Methods
    -------
    turn_off()
        Turn off all thermostats.

    set_langevin()
        Set the Langevin thermostat.

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
        seed : :obj:`int`
            Initial counter value (or seed) of the philox RNG.
            Required on first activation of the Langevin thermostat.
            Must be positive.

    set_brownian()
        Set the Brownian thermostat.

        Parameters
        ----------
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
        seed : :obj:`int`
            Initial counter value (or seed) of the philox RNG.
            Required on first activation of the Brownian thermostat.
            Must be positive.

    set_lb()
        Set the LB thermostat. The kT parameter is automatically extracted
        from the currently active LB solver.

        This thermostat requires the feature ``WALBERLA``.

        Parameters
        ----------
        seed : :obj:`int`
            Seed for the random number generator, required if kT > 0.
            Must be positive.
        gamma : :obj:`float`
            Frictional coupling constant for the MD particle coupling.

    set_npt()
        Set the NPT thermostat.

        Parameters
        ----------
        kT : :obj:`float`
            Thermal energy of the heat bath.
        gamma0 : :obj:`float`
            Friction coefficient of the bath.
        gammav : :obj:`float`
            Artificial friction coefficient for the volume fluctuations.
        seed : :obj:`int`
            Initial counter value (or seed) of the philox RNG.
            Required on first activation of the Langevin thermostat.
            Must be positive.

    set_dpd()
        Set the DPD thermostat with required parameters 'kT'.
        This also activates the DPD interactions.

        Parameters
        ----------
        kT : :obj:`float`
            Thermal energy of the heat bath.
        seed : :obj:`int`
            Initial counter value (or seed) of the philox RNG.
            Required on first activation of the DPD thermostat.
            Must be positive.

    set_stokesian()
        Set the SD thermostat with required parameters.

        This thermostat requires the feature ``STOKESIAN_DYNAMICS``.

        Parameters
        ----------
        kT : :obj:`float`, optional
            Temperature.
        seed : :obj:`int`, optional
            Seed for the random number generator.

    """
    _so_name = "Thermostat::Thermostat"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = (
        "set_langevin",
        "set_brownian",
        "set_npt",
        "set_dpd",
        "set_lb",
        "set_stokesian",
        "set_thermalized_bond",
        "turn_off")


@script_interface_register
class Langevin(ScriptInterfaceHelper):
    _so_name = "Thermostat::Langevin"
    _so_creation_policy = "GLOBAL"


@script_interface_register
class Brownian(ScriptInterfaceHelper):
    _so_name = "Thermostat::Brownian"
    _so_creation_policy = "GLOBAL"


@script_interface_register
class IsotropicNpt(ScriptInterfaceHelper):
    _so_name = "Thermostat::IsotropicNpt"
    _so_creation_policy = "GLOBAL"


@script_interface_register
class Stokesian(ScriptInterfaceHelper):
    _so_name = "Thermostat::Stokesian"
    _so_creation_policy = "GLOBAL"


@script_interface_register
class LBThermostat(ScriptInterfaceHelper):
    _so_name = "Thermostat::LB"
    _so_creation_policy = "GLOBAL"


@script_interface_register
class DPDThermostat(ScriptInterfaceHelper):
    _so_name = "Thermostat::DPD"
    _so_creation_policy = "GLOBAL"


@script_interface_register
class ThermalizedBond(ScriptInterfaceHelper):
    _so_name = "Thermostat::ThermalizedBond"
    _so_creation_policy = "GLOBAL"
