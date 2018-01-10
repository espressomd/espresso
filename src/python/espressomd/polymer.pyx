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
include "myconfig.pxi"
from . cimport polymer
import numpy as np
from espressomd.utils import is_valid_type

def validate_params(_params, default):
    if _params["N_P"] <= 0:
        raise ValueError(
                "N_P has to be a positive Integer" )
    if _params["MPC"] <= 1:
        raise ValueError(
                "MPC has to be a positive Integer larger than 1" )
    if _params["bond_length"] < 0 :
        raise ValueError(
                "bond_length has to be a positive float" )
    if _params["start_id"] < 0:
        raise ValueError(
                "start_id has to be a positive Integer")
    if not isinstance(_params["start_pos"], np.ndarray) or len(_params["start_pos"]) != 3:
        raise ValueError(
                "start_pos has to be an numpy array with 3 Elements" )
    if not is_valid_type(_params["mode"], int):
        raise ValueError(
                "mode has to be a positive Integer" )
    if _params["shield"] < 0 and default["shield"] != _params["shield"]:
        raise ValueError(
                "shield has to be a positive float")
    if _params["max_tries"] < 0 and default["max_tries"] != _params["max_tries"]:
        raise ValueError(
                "max_tries has to be a positive Integer")
    if not is_valid_type(_params["val_poly"], float) and default["val_poly"] != _params["val_poly"]:
        raise ValueError(
                "val_poly has to be a float")
    if _params["charge_distance"] < 0:
        raise ValueError(
                "charge_distance has to be a positive Integer")
    if _params["type_poly_neutral"] < 0:
        raise ValueError(
                "type_poly_neutral has to be a nonnegative Integer")
    if _params["type_poly_charged"] < 0:
        raise ValueError(
                "type_poly_charged has to be a nonnegative Integer")
    if _params["angle"] < 0 and default["angle"] != _params["angle"]:
        raise ValueError(
                "angle has to be a positive float")
    if _params["angle2"] < 0 and default["angle2"] != _params["angle2"]:
        raise ValueError(
                "angle2 has to be a positive float")
    if _params["constraints"] < 0 :
        raise ValueError(
                "constraint has to be either 0 or 1" )

# wrapper function to expose to the user interface
def create_polymer(**kwargs):
    """ Generators have a ``Yields`` section instead of a ``Returns`` section.

    Parameters
    ----------
    n : :obj:`intd`
        The upper limit of the range to generate, from 0 to `n` - 1.
    N_P : :obj:`int`
        Number of polymer chains
    MPC : :obj:`int`
        Number of monomers per chain
    bond_length : :obj:`float`
        distance between adjacent monomers in a chain
    bond : :obj:`espressomd.interactions.BondedInteraction`
        The bonded interaction to be set up between the monomers. 
    start_id : :obj:`int`, optional
        Particle ID of the first monomer, all other particles will have larger IDs. Defaults to 0
    start_pos : array_like :obj:`float`. Defaults to numpy.array([0, 0, 0])
        Position of the first monomer
    mode : :obj:`int`, optional
        Selects a specific random walk procedure for the
        polymer setup mode = 1 uses a common random walk,
        mode = 2 produces a pruned self-avoiding random walk,
        and mode = 0 a self-avoiding random walk. Note that
        mode = 2 does not produce a true self-avoiding
        random walk distribution but is much faster than mode = 0. Defaults to 1
    shield : :obj:`float`, optional
        Shielding radius for the pruned self-avoiding walk mode. Defaults to 0
    max_tries : :obj:`int`, optional
        Maximal number of attempts to set up a polymer,
        default value is 1,000. Depending on the random walk
        mode and the polymer length this value needs to be
        adapted. 
    val_poly : :obj:`float`, optional
        Valency of the monomers, default is 0.0
    charge_distance : :obj:`int`, optional
        Distance between charged monomers along the chain. Default is 1
    type_poly_neutral : :obj:`int`, optional
        Particle type of neutal monomers, default is 0.
    type_poly_charged : :obj:`int`, optional
        Particle type for charged monomers, default is 1
    angle : :obj:`float`, optional
    angle2 : :obj:`float`, optional
        The both angles angle and angle2 allow to set up
        planar or helical polymers, they fix the angles
        between adjacent bonds.
    pos2 : array_like, optional
        Sets the position of the second monomer. Defaults to numpy.array([0, 0, 0]).
    constraints : :obj:`int`, optional
        Either 0 or 1, default is 0. If 1, the particle setup-up tries to obey previously defined constraints.
        
    Examples
    --------
    This example sets 2 polyelectrolyte chains of the length 10. Beads are connected by FENE potential.

    >>> fene = interactions.FeneBond(k=10, d_r_max=2)
    >>> polymer.create_polymer(
            N_P = 2, 
            MPC = 10, 
            bond_length = 1, 
            bond = fene, 
            val_poly = -1.0)
    """

    params=dict()
    default_params=dict()
    default_params["N_P"] = 0 
    default_params["MPC"] = 0
    default_params["bond_length"] = 0 
    default_params["start_id"] = 0
    default_params["start_pos"] = np.array([0, 0, 0])
    default_params["mode"] = 1 
    default_params["shield"] = 0
    default_params["max_tries"] = 1000
    default_params["val_poly" ] = 0.0
    default_params["charge_distance"] = 1
    default_params["type_poly_neutral"] = 0 
    default_params["type_poly_charged"] = 1
    default_params["angle"] = -1.0
    default_params["angle2"] = -1.0
    default_params["constraints"]=0 
    default_params["pos2"] = np.array([0, 0, 0])

    params = default_params 

    valid_keys=["N_P", "MPC", "bond_length", "bond", "start_id", "start_pos", "mode", "shield", "max_tries", "val_poly", "charge_distance", "type_poly_neutral", "type_poly_charged", "angle", "angle2", "constraints"]

    required_keys=["N_P", "MPC", "bond_length", "bond"]

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

    bond_id = params["bond"]._bond_id

    cdef double start_pos[3];
    cdef double start_pos2[3];
    for i in range(3):
        start_pos[i] = params["start_pos"][i]
        start_pos2[i] =params["pos2"][i]

    polymerC(partCfg(), params["N_P"], params["MPC"], params["bond_length"], params["start_id"], \
             start_pos, params["mode"], params["shield"], params["max_tries"], \
             params["val_poly"], params["charge_distance"], params["type_poly_neutral"], \
             params["type_poly_charged"], bond_id, params["angle"], \
             params["angle2"], start_pos2, params["constraints"])
