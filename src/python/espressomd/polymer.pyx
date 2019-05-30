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
include "myconfig.pxi"
from . cimport polymer
import numpy as np
from espressomd.utils import is_valid_type, array_locked
from espressomd.utils cimport make_Vector3d


def validate_params(_params, default):
    if _params["n_polymers"] <= 0:
        raise ValueError(
            "n_polymers has to be a positive integer")
    if _params["beads_per_chain"] <= 1:
        raise ValueError(
            "beads_per_chain has to be a positive integer larger than 1")
    if _params["bond_length"] < 0:
        raise ValueError(
            "bond_length has to be a positive float")
    if not ((np.size(_params["start_positions"]) == 0)
            or np.shape(_params["start_positions"]) == (_params["n_polymers"], 3)):
        raise ValueError(
            "start_positions has to be a numpy array with shape (n_polymers, 3)")
    if _params["min_distance"] < 0:
        raise ValueError(
            "min_distance has to be a positive float")
    if _params["max_tries"] < 0:
        raise ValueError(
            "max_tries has to be a positive integer")
    if _params["use_bond_angle"] and np.isnan(_params["bond_angle"]):
        raise ValueError(
            "bond_angle has to be a positive float")
    if type(_params["respect_constraints"]) != bool:
        raise ValueError(
            "respect_constraints has to be either True or False")
    if type(_params["seed"]) != int:
        raise ValueError(
            "seed has to be an integer")

# wrapper function to expose to the user interface


def positions(**kwargs):
    """
    Generates particle positions for polymer creation.

    Parameters
    ----------
    n_polymers : :obj:`int`, required
        Number of polymer chains
    beads_per_chain : :obj:`int`, required
        Number of monomers per chain
    bond_length : :obj:`float`, required
        distance between adjacent monomers in a chain
    seed : :obj:`int`, required
        Seed for the RNG used to generate the particle positions.
    bond_angle : :obj:`float`, optional
        If set, this parameter defines the angle between adjacent bonds
        withing a polymer.
    start_positions : array_like :obj:`float`.
        If set, this vector defines the start positions for the polymers, i.e.,
        the position of each polymer's first monomer bead.
        Here, a numpy array of shape (n_polymers, 3) is expected.
    min_distance : :obj:`float`, optional
        Minimum distance between all generated positions. Defaults to 0
    respect_constraints : :obj:`bool`, optional
        If True, the particle setup tries to obey previously defined constraints.
        Default value is False.
    max_tries : :obj:`int`, optional
        Maximal number of attempts to generate every monomer position,
        as well as maximal number of retries per polymer, if choosing
        suitable monomer positions fails. Default value is 1000.
        Depending on the total number of beads and constraints,
        this value needs to be adapted.

    Returns
    -------
    array_like :obj:`float`
        Three-dimensional numpy array, namely a list of polymers containing the
        coordinates of the respective monomers.

    """
    params = dict()
    default_params = dict()
    default_params["n_polymers"] = 0
    default_params["beads_per_chain"] = 0
    default_params["bond_length"] = 0
    default_params["start_positions"] = np.array([])
    default_params["min_distance"] = 0
    default_params["max_tries"] = 1000
    default_params["bond_angle"] = -1
    default_params["respect_constraints"] = False
    default_params["seed"] = None

    params = default_params

    # use bond_angle if set via kwarg
    params["use_bond_angle"] = "bond_angle" in kwargs

    valid_keys = [
        "n_polymers",
        "beads_per_chain",
     "bond_length",
     "start_positions",
     "min_distance",
     "max_tries",
     "bond_angle",
     "respect_constraints",
     "seed"]

    required_keys = ["n_polymers", "beads_per_chain", "bond_length", "seed"]

    for k in kwargs:
        if not k in valid_keys:
            raise ValueError("Unknown parameter '%s'" % k)
        params[k] = kwargs[k]

    for k in required_keys:
        if k not in kwargs:
            print(k)
            raise ValueError(
                "At least the following keys have to be given as keyword arguments: " + required_keys.__str__())

    validate_params(params, default_params)

    cdef vector[Vector3d] start_positions
    if (params["start_positions"] != []):
        for i in range(len(params["start_positions"])):
            start_positions.push_back(
                make_Vector3d(params["start_positions"][i]))

    data = draw_polymer_positions(
        partCfg(),
        params["n_polymers"],
     params["beads_per_chain"],
     params["bond_length"],
     start_positions,
     params["min_distance"],
     params["max_tries"],
     int(params["use_bond_angle"]),
     params["bond_angle"],
     int(params["respect_constraints"]),
     params["seed"])
    positions = []
    for polymer in data:
        p = []
        for monomer in polymer:
            m = array_locked([monomer[0], monomer[1], monomer[2]])
            p.append(m)
        positions.append(p)
    return np.array(positions)
