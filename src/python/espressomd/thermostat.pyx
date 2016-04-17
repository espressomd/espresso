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
cimport thermostat

cdef class Thermostat:

    def __init__(self):
        pass

    def get_status(self):
        """Returns the thermostat status. Equivalent to the tcl command 'thermostat' with no arguments"""

        if temperature == -1:
            return "{ not initialized }"
        if thermo_switch == THERMO_OFF:
            return "{ off }"
        if thermo_switch and THERMO_LANGEVIN:
            return "{ langevin " + str(temperature) + " " + str(langevin_gamma) + " }"

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
