#
# Copyright (C) 2010-2024 The ESPResSo project
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

from .script_interface import ScriptInterfaceHelper, script_interface_register


@script_interface_register
class CollisionDetection(ScriptInterfaceHelper):

    """
    Interface to the collision detection / dynamic binding.

    See :ref:`Creating bonds when particles collide` for detailed instructions.

    This class should not be instantiated by the user. Instead, use
    the :attr:`~espressomd.system.System.collision_detection` attribute
    of the system class to access the collision detection.

    Attributes
    ----------
    protocol :
        Protocol instance, or ``None`` if not set.

    """

    _so_name = "CollisionDetection::CollisionDetection"
    _so_features = ("COLLISION_DETECTION",)


@script_interface_register
class Off(ScriptInterfaceHelper):

    """
    Disable collision detection.
    """

    _so_name = "CollisionDetection::Off"
    _so_features = ("COLLISION_DETECTION",)


@script_interface_register
class BindCenters(ScriptInterfaceHelper):

    """
    Add a pair bond between two particles upon collision.
    Particles can still slide around their contact point.

    Parameters
    ----------
    distance : :obj:`float`
        Distance below which two particles are considered to have collided.

    bond_centers : :obj:`espressomd.interactions.BondedInteraction`
        Bond to add between the colliding particles.

    """
    _so_name = "CollisionDetection::BindCenters"
    _so_features = ("COLLISION_DETECTION",)


@script_interface_register
class BindAtPointOfCollision(ScriptInterfaceHelper):

    """
    Add a pair bond between two particles upon collision, and two extra
    pair bonds or angle bonds with two automatically-generated virtual sites.
    This protocol prevents sliding of the particles at the contact point.

    Parameters
    ----------
    distance : :obj:`float`
        Distance below which two particles are considered to have collided.

    bond_centers : :obj:`espressomd.interactions.BondedInteraction`
        Bond to add between the colliding particles.

    bond_vs : :obj:`espressomd.interactions.BondedInteraction`
        Bond to add between virtual sites.

    part_type_vs : :obj:`int`
        Particle type of the virtual sites created on collision.

    vs_placement : :obj:`float`
        Barycenter of the virtual sites. A value of 0 means that the virtual sites
        are placed at the same position as the colliding particles on which they are based.
        A value of 0.5 will result in the virtual sites being placed at the mid-point between
        the two colliding particles. A value of 1 will result the virtual site associated to
        the first colliding particle to be placed at the position of the second colliding particle.
        In most cases, 0.5, is a good choice. Then, the bond connecting the virtual sites
        should have an equilibrium length of zero.

    """
    _so_name = "CollisionDetection::BindAtPointOfCollision"
    _so_features = ("COLLISION_DETECTION", "VIRTUAL_SITES_RELATIVE")


@script_interface_register
class GlueToSurface(ScriptInterfaceHelper):

    """
    Attach small particles to the surface of a larger particle.

    It is asymmetric: several small particles can be bound to a large particle
    but not vice versa. It can be made irreversible: the small particles can
    change type after collision to become *inert*.

    On collision, a single virtual site is placed and related to the large
    particle. Then a bond (``bond_centers``) connects the large and the small
    particle. A second bond (``bond_vs``) connects the virtual site and
    the small particle.

    Parameters
    ----------
    distance : :obj:`float`
        Distance below which two particles are considered to have collided.

    bond_centers : :obj:`espressomd.interactions.BondedInteraction`
        Bond to add between the colliding particles.

    bond_vs : :obj:`espressomd.interactions.BondedInteraction`
        Bond to add between virtual sites.

    part_type_vs : :obj:`int`
        Particle type of the virtual sites created on collision.

    part_type_to_attach_vs_to : :obj:`int`
        Type of the large particle.

    part_type_to_be_glued : :obj:`int`
        Type of the small particle.

    part_type_after_glueing : :obj:`int`
        Type of the small particle after collision. If different from
        ``part_type_to_be_glued``, the bond is irreversible.

    distance_glued_particle_to_vs : :obj:`float`
        Distance of the virtual site to the small particle,
        as a fraction of the pair distance.

    """
    _so_name = "CollisionDetection::GlueToSurface"
    _so_features = ("COLLISION_DETECTION", "VIRTUAL_SITES_RELATIVE")
