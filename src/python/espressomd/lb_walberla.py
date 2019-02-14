# Copyright (C) 2010-2018 The ESPResSo project
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
from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper, script_interface_register
import numpy as np
import array
from . import utils
from espressomd.utils import array_locked, is_valid_type


@script_interface_register
class LbWalberla(ScriptInterfaceHelper):

    """Interface to a Lattice-Boltzmann fluid backed by Walberla

    """
    _so_name = "LbWalberla"

    def __getitem__(self, key):
        if isinstance(key, tuple) or(key, list) or isinstance(key, np.ndarray):
            if len(key) == 3:
                return LbWalberlaFluidRoutines(np.array(key), self)
        else:
            raise Exception(
                    "%s is not a valid key. Should be a point on the nodegrid e.g. lbf[0,0,0]," % key)

class LbWalberlaFluidRoutines(object):
    def __init__(self, key, lb):
        utils.check_type_or_throw_except(
            key, 3, int, "The index of an lb fluid node consists of three integers.")
        self.node = np.int_([0,0,0])
        self.node[0] = int(key[0])
        self.node[1] = int(key[1])
        self.node[2] = int(key[2])
        self.lbWalberla = lb
# TODO: check for valid input range
#        if not lb_lbnode_is_index_valid(self.node):
#            raise ValueError("LB node index out of bounds")

    @property 
    def velocity(self):
        double_return = self.lbWalberla.call_method("node_get_v", node=self.node.tolist())
        return array_locked(double_return)

    @velocity.setter
    def velocity(self, value):
        if all(is_valid_type(v, float) for v in value) and len(value) == 3:
            host_velocity = np.zeros(3)
            host_velocity[0] = value[0]
            host_velocity[1] = value[1]
            host_velocity[2] = value[2]
            self.lbWalberla.call_method("node_set_v", node=self.node.tolist(), v=host_velocity.tolist())
        else:
            raise ValueError(
                "Velocity has to be of shape 3 and type float.")





