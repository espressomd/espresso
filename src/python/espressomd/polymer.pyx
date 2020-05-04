#
# Copyright (C) 2013-2019 The ESPResSo project
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
from . cimport polymer
import numpy as np
from .system import System
from .interactions import BondedInteraction
from .utils cimport make_Vector3d, check_type_or_throw_except
from .utils import array_locked


def validate_params(_params, default):
    if _params["n_polymers"] <= 0:
        raise ValueError(
            "n_polymers has to be a positive integer")
    if _params["beads_per_chain"] <= 0:
        raise ValueError(
            "beads_per_chain has to be a positive integer")
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


def linear_polymer_positions(**kwargs):
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
        within a polymer.
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
    :obj:`ndarray`
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

    # use bond_angle if set via kwargs
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
        if k not in valid_keys:
            raise ValueError("Unknown parameter '%s'" % k)
        params[k] = kwargs[k]

    for k in required_keys:
        if k not in kwargs:
            print(k)
            raise ValueError(
                "At least the following keys have to be given as keyword arguments: " + required_keys.__str__())

    validate_params(params, default_params)

    cdef vector[Vector3d] start_positions
    if (params["start_positions"].size > 0):
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


def setup_diamond_polymer(system=None, bond=None, MPC=0, 
                          dist_cM=1, val_cM=0.0, val_nodes=0.0, 
                          start_id='auto', no_bonds=False, 
                          type_nodes=0, type_nM=1, type_cM=2):
    """
    Places particles to form a diamond lattice shaped polymer.
    Can also assign charges and bonds at the appropriate places.

    Parameters
    ----------
    system : :class:`espressomd.system.System`, required
        System to which the particles will be added.
    bond : :class:`espressomd.interactions.BondedInteraction`, required if ``no_bonds == False``
        The bond to be created between monomers. Should be compatible with the 
        spacing ``system.box_l[0]*(0.25 * sqrt(3))/(MPC + 1)`` between monomers.
    no_bonds : :obj:`bool`, optional
        If True, the particles will only be placed in the system but not connected by bonds. 
        In that case, the ``bond`` argument can be omitted. Defaults to ``False``.
    MPC : :obj:`int`, optional
        Monomers per chain, where chain refers to the connection 
        between the 8 lattice nodes of the diamond lattice.
        Defaults to 0.
    dist_cM : :obj:`int`, optional
        Distance between charged monomers in the chains. Defaults to 1.
    val_cM : :obj:`float`, optional
        Valence of the charged monomers in the chains. Defaults to 0.
    val_nodes : :obj:`float`, optional
        Valence of the node particles. Defaults to 0.
    start_id : :obj:`int` or ``'auto'``, optional
        Start id for particle creation. Subsequent ids will be contiguous integers.
        If ``'auto'``, particle ids will start after the highest id of particles already in the system. 
    type_nodes : :obj:`int`, optional
        Type assigned to the node particles. Defaults to 0.
    type_nM : :obj:`int`, optional
        Type assigned to the neutral monomers in the chains. Defaults to 1.
    type_cM : :obj:`int`, optional
        Type assigned to the charged monomers in the chains. Defaults to 2.
    """

    if start_id == 'auto':
        start_id = system.part.highest_particle_id + 1

    check_type_or_throw_except(
        no_bonds, 1, bool, "no_bonds must be one bool")
    if not no_bonds and not isinstance(bond, BondedInteraction):
        raise TypeError(
            "bond argument must be an instance of espressomd.interaction.BondedInteraction")
    if not isinstance(system, System):
        raise TypeError(
            "System argument must be an instance of an espressomd System")

    check_type_or_throw_except(
        MPC, 1, int, "MPC must be one int")
    check_type_or_throw_except(
        val_cM, 1, float, "val_cM must be one float")
    check_type_or_throw_except(
        val_nodes, 1, float, "val_nodes must be one float")
    check_type_or_throw_except(
        type_nodes, 1, int, "type_nodes must be one int")
    check_type_or_throw_except(
        type_nM, 1, int, "type_nM must be one int")
    check_type_or_throw_except(
        type_cM, 1, int, "type_cM must be one int")
    check_type_or_throw_except(
        dist_cM, 1, int, "dist_cM must be one int")
    check_type_or_throw_except(
        start_id, 1, int, "start_id must be one int or 'auto'")

    box = system.box_l
    if not box[0] == box[1] == box[2]:
        raise Exception("Simulation box must be cubic but is {}".format(box))
    box_length = box[0]

    node_positions = box_length / 4. * np.array([[0, 0, 0], [1, 1, 1],
                                                 [2, 2, 0], [0, 2, 2],
                                                 [2, 0, 2], [3, 3, 1],
                                                 [1, 3, 3], [3, 1, 3]])
    connected_nodes = [(0, 1), (1, 2), (1, 3), (1, 4),
                       (2, 5), (3, 6), (4, 7), (5, 0),
                       (5, 3), (5, 4), (6, 0), (6, 2),
                       (6, 4), (7, 0), (7, 2), (7, 3)]

    # place nodes
    node_ids = []
    current_id = start_id
    for node_pos in node_positions:
        system.part.add(pos=node_pos, id=current_id, 
                        type=type_nodes, q=val_nodes)
        node_ids.append(current_id)
        current_id += 1

    # place monomers inbetween
    if MPC > 0:
        for start_node, end_node in connected_nodes:
            node_connection_vec = (
                node_positions[end_node, :] - node_positions[start_node, :])
            # find minimum image of neighbour node
            node_connection_vec -= np.rint(node_connection_vec /
                                           box_length) * box_length

            for j in range(1, MPC + 1):
                if np.mod(j, dist_cM) == 0:
                    mono_q = val_cM
                    mono_type = type_cM
                else:
                    mono_q = 0
                    mono_type = type_nM
                pos = node_positions[start_node, :] + \
                    j / (MPC + 1) * node_connection_vec
                system.part.add(
                    pos=pos,
                    id=current_id,
                    type=mono_type,
                    q=mono_q)

                # add bonds along the chain
                if not no_bonds:
                    if j == 1:
                        system.part[node_ids[start_node]].add_bond(
                            (bond, current_id))
                    else:
                        system.part[current_id].add_bond(
                            (bond, current_id - 1))
                    if j == MPC:
                        system.part[node_ids[end_node]].add_bond(
                            (bond, current_id))
                current_id += 1
