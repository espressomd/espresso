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

import numpy as np
import collections
import functools
from .interactions import BondedInteraction
from .interactions import BondedInteractions
from .utils import nesting_level, array_locked, is_valid_type, handle_errors
from .utils import check_type_or_throw_except
from .code_features import assert_features, has_features
from .script_interface import script_interface_register, ScriptInterfaceHelper
import itertools


@script_interface_register
class ParticleHandle(ScriptInterfaceHelper):
    """
    Attributes
    ----------
    id: :obj:`int`
        Particle identifier.

    type: :obj:`int`
        The particle type for non-bonded interactions.

        .. note::
           The value of ``type`` has to be an integer >= 0.

    mol_id: :obj:`int`
        The molecule id of the particle.

        The particle ``mol_id`` is used to differentiate between
        particles belonging to different molecules, e.g. when virtual
        sites are used, or object-in-fluid cells. The default
        ``mol_id`` for all particles is 0.

        .. note::
           The value of ``mol_id`` has to be an integer >= 0.

    pos: (3,) array_like of :obj:`float`
        The unwrapped (not folded into central box) particle position.

    pos_folded: (3,) array_like of :obj:`float`
        The wrapped (folded into central box) position vector of a particle.

        .. note::
           Setting the folded position is ambiguous and is thus not possible,
           please use ``pos`` instead.

        Examples
        --------
        >>> import espressomd
        >>> system = espressomd.System(box_l=[10, 10, 10])
        >>> system.part.add(pos=(5, 0, 0))
        >>> system.part.add(pos=(10, 0, 0))
        >>> system.part.add(pos=(25, 0, 0))
        >>> for p in system.part:
        ...     print(p.pos)
        [ 5.  0.  0.]
        [ 10.   0.   0.]
        [ 25.   0.   0.]
        >>> for p in system.part:
        ...     print(p.pos_folded)
        [5.0, 0.0, 0.0]
        [0.0, 0.0, 0.0]
        [5.0, 0.0, 0.0]

    image_box: (3,) array_like of :obj:`int`
        The image box the particles is in.

        This is the number of times the particle position has been folded by
        the box length in each direction.

    lees_edwards_offset: :obj:`float`
        The accumulated Lees-Edwards offset.
        Can be used to reconstruct continuous trajectories.

    lees_edwards_flag: :obj:`int`
        The Lees-Edwards flag that indicate if the particle crossed
        the upper or lower boundary.

    v: (3,) array_like of :obj:`float`
        The particle velocity in the lab frame.

        .. note::
           The velocity will be updated during integration.

    f: (3,) array_like of :obj:`float`
        The instantaneous force acting on this particle.

        .. note::
           The force is recomputed during the integration step and any force
           set in this way is immediately lost at the next integration step.

    node: (3,) array_like of :obj:`int`
        The node the particle is on, identified by its MPI rank.

    mass: :obj:`float`
        Particle mass.

    omega_lab: (3,) array_like of :obj:`float`
        The particle angular velocity the lab frame.

        .. note::
           This needs the feature ``ROTATION``.

           If you set the angular velocity of the particle in the lab
           frame, the orientation of the particle
           (:attr:`~espressomd.particle_data.ParticleHandle.quat`) must be
           set before setting ``omega_lab``, otherwise the conversion from
           lab to body frame will not be handled properly.

        See Also
        ---------
        :attr:`~espressomd.particle_data.ParticleHandle.omega_body`

    quat: (4,) array_like of :obj:`float`
        Quaternion representation of the particle rotational position.

        .. note::
           This needs the feature ``ROTATION``.

    director: (3,) array_like of :obj:`float`
        The particle director.

        The ``director`` defines the the z-axis in the body-fixed frame.
        If particle rotations happen, the director, i.e., the body-fixed
        coordinate system co-rotates. Properties such as the angular
        velocity :attr:`espressomd.particle_data.ParticleHandle.omega_body`
        are evaluated in this body-fixed coordinate system.
        When using particle dipoles, the dipole moment is co-aligned with
        the particle director. Setting the director thus modifies the
        dipole moment orientation (:attr:`espressomd.particle_data.ParticleHandle.dip`)
        and vice versa.
        See also :ref:`Rotational degrees of freedom and particle anisotropy`.

        .. note::
           This needs the feature ``ROTATION``.

    omega_body: (3,) array_like of :obj:`float`
        The particle angular velocity in body frame.

        This property sets the angular momentum of this particle in the
        particles co-rotating frame (or body frame).

        .. note::
           This needs the feature ``ROTATION``.

    torque_lab: (3,) array_like of :obj:`float`
        The particle torque in the lab frame.

        This property defines the torque of this particle
        in the fixed frame (or laboratory frame).

        .. note::
           The orientation of the particle
           (:attr:`~espressomd.particle_data.ParticleHandle.quat`) must be
           set before setting this property, otherwise the conversion from
           lab to body frame will not be handled properly.

    rinertia: (3,) array_like of :obj:`float`
        The particle rotational inertia.

        Sets the diagonal elements of this particles rotational inertia
        tensor. These correspond with the inertial moments along the
        coordinate axes in the particle's co-rotating coordinate system.
        When the particle's quaternions are set to ``[1, 0, 0, 0,]``, the
        co-rotating and the fixed (lab) frames are co-aligned.

        .. note::
           This needs the feature ``ROTATIONAL_INERTIA``.

    q: :obj:`float`
        Particle charge.

        .. note::
           This needs the feature ``ELECTROSTATICS``.

    mu_E: :obj:`float`
        Particle electrophoretic velocity.

        This effectively acts as a velocity offset between
        a lattice-Boltzmann fluid and the particle. Has only
        an effect if LB is turned on.

        .. note::
           This needs the feature ``LB_ELECTROHYDRODYNAMICS``.

    virtual: :obj:`bool`
        Virtual flag.

        Declares the particles as virtual (``True``) or non-virtual
        (``False``, default).

        .. note::
           This needs the feature ``VIRTUAL_SITES``

    vs_quat: (4,) array_like of :obj:`float`
        Virtual site quaternion.

        This quaternion describes the virtual particles orientation in the
        body fixed frame of the related real particle.

        .. note::
           This needs the feature ``VIRTUAL_SITES_RELATIVE``.

    vs_relative: :obj:`tuple`
        Virtual sites relative parameters.

        Allows for manual access to the attributes of virtual sites in the
        "relative" implementation. Format: ``(PID, distance, quaternion)``.
        PID denotes the id of the particle to which this virtual site is
        related and distance the distance between non-virtual and virtual particle.
        The relative orientation is specified as a quaternion.

        .. note::
           This needs the feature ``VIRTUAL_SITES_RELATIVE``

    dip: (3,) array_like of :obj:`float`
        The orientation of the dipole axis.

        .. note::
           This needs the feature ``DIPOLES``.

    dipm: :obj:`float`
        The magnitude of the dipole moment.

        .. note::
           This needs the feature ``DIPOLES``.

    dip_fld: (3,) array_like of :obj:`float`
        Total dipole field value at the position of the particle.

        .. note::
           This needs the feature ``DIPOLE_FIELD_TRACKING``.

    ext_force: (3,) array_like of :obj:`float`
        An additional external force applied to the particle.

        .. note::
           This needs the feature ``EXTERNAL_FORCES``.

    fix: (3,) array_like of :obj:`bool`
        Fixes the particle motion in the specified cartesian directions.

        Fixes the particle in space. It is possible to fix motion in the
        x-, y-, or z-direction independently. For example::

            part.by_id(1).fix = [False, False, True]

        will fix motion for particle with index 1 only in the z-direction.

        .. note::
           This needs the feature ``EXTERNAL_FORCES``.

    ext_torque: (3,) array_like of :obj:`float`
        An additional external torque is applied to the particle.

        ..  note::
            * This torque is specified in the laboratory frame!
            * This needs features ``EXTERNAL_FORCES`` and ``ROTATION``.

    gamma: :obj:`float` or (3,) array_like of :obj:`float`
        The translational frictional coefficient used in the Langevin
        and Brownian thermostats.

        .. note::
            This needs feature ``THERMOSTAT_PER_PARTICLE`` and
            optionally ``PARTICLE_ANISOTROPY``.

        See Also
        ----------
        :meth:`espressomd.thermostat.Thermostat.set_langevin` : Setting the parameters of the Langevin thermostat

    gamma_rot: :obj:`float` or (3,) array_like of :obj:`float`
        The particle rotational frictional coefficient used in
        the Langevin and Brownian thermostats.

        gamma_rot : :obj:`float` or (3,) array_like of :obj:`float`

        .. note::
            This needs features ``THERMOSTAT_PER_PARTICLE``, ``ROTATION`` and
            optionally ``PARTICLE_ANISOTROPY``.

    rotation: (3,) array_like of :obj:`bool`
        Switches the particle's rotational degrees of freedom in the
        Cartesian axes in the body-fixed frame. The content of the torque
        and omega variables are meaningless for the co-ordinates for which
        rotation is disabled.

        The default is not to integrate any rotational degrees of freedom.

        rotation : (3,) array_like of :obj:`bool`

        .. note::
            This needs the feature ``ROTATION``.

    swimming:
        Set swimming parameters.

        This property takes a dictionary with a different number of entries
        depending whether there is an implicit fluid (i.e. with the Langevin
        thermostat) of an explicit fluid (with lattice-Boltzmann).

        Swimming enables particle self-propulsion in the direction
        determined by its quaternion. For setting the quaternion of the
        particle see :attr:`~espressomd.particle_data.ParticleHandle.quat`.
        Self-propulsion is achieved by imposing a constant force term
        ``f_swim`` along the particle direction. The steady-state propulsion
        speed (``v_swim``) can be calculated from the friction (``gamma``)
        of a thermostat: ``v_swim = f_swim / gamma``.
        When resolving hydrodynamics via lattice-Boltzmann, the swimming
        attribute can be used to create the typical dipolar flowfield of
        self-propelled particles: setting ``is_engine_force_on_fluid``
        to ``True`` will make the particle not experience any friction
        or noise, but instead apply the swim force ``f_swim`` to the fluid.
        Use :func:`espressomd.swimmer_helpers.add_dipole_particle`
        to automate such a setup.

        Parameters
        ----------
        f_swim : :obj:`float`
            Magnitude of the self-propulsion force.
        is_engine_force_on_fluid : :obj:`bool`
            Default: ``False``.
            If ``True``, the particle will apply the swimming force to the fluid
            instead of experiencing drag.

        Notes
        -----
        This needs feature ``ENGINE``, and optionally ``VIRTUAL_SITES_RELATIVE``
        to add the propulsion force on a lattice-Boltzmann fluid.

        Examples
        --------
        >>> import espressomd
        >>> # swimming withut hydrodynamics
        >>> system = espressomd.System(box_l=[10, 10, 10])
        >>> partcl = system.part.add(pos=[1, 0, 0], swimming={'f_swim': 0.03})
        >>> # swimming with hydrodynamics
        >>> import espressomd.swimmer_helpers.add_dipole_particle as add_dip
        >>> dipole_partcl = add_dip(system, partcl, 2., 0)

    Methods
    -------
    delete_all_bonds()
        Delete all bonds from the particle.

        See Also
        ----------
        delete_bond : Delete an unverified bond held by the particle.
        bonds : ``Particle`` property containing a list of all current bonds held by ``Particle``.

    """

    _so_name = "Particles::ParticleHandle"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = (
        "delete_all_bonds",
    )

    # here we must redefine the script interface setters

    def set_params(self, **kwargs):
        for name, value in kwargs.items():
            self.set_parameter(name, value)

    def set_parameter(self, name, value):
        return self.call_method("set_param_parallel", name=name, value=value)

    def remove(self):
        """
        Delete the particle.

        See Also
        --------
        espressomd.particle_data.ParticleList.add
        espressomd.particle_data.ParticleList.clear

        """
        self.call_method("remove_particle")
        del self

    def to_dict(self):
        """
        Returns the particle's attributes as a dictionary.

        It includes the content of ``particle_attributes``, minus a few exceptions:

        - :attr:`~ParticleHandle.dip`, :attr:`~ParticleHandle.director`:
          Setting only the director will overwrite the orientation of the
          particle around the axis parallel to dipole moment/director.
          Quaternions contain the full info.
        - :attr:`~ParticleHandle.image_box`, :attr:`~ParticleHandle.node`

        """

        pdict = self.get_params()
        for k in ["director", "dip", "pos_folded",
                  "image_box", "node", "lees_edwards_flag"]:
            if k in pdict:
                del pdict[k]
        if has_features("EXCLUSIONS"):
            pdict["exclusions"] = self.exclusions
        pdict["bonds"] = self.bonds
        return pdict

    def __str__(self):
        res = collections.OrderedDict()
        # Id and pos first, then the rest
        res["id"] = self.id
        res["pos"] = self.pos
        for attr in particle_attributes:
            tmp = getattr(self, attr)
            # Remove array type names from output
            if isinstance(tmp, array_locked):
                res[attr] = tuple(tmp)
            else:
                res[attr] = tmp

        # Get rid of OrderedDict in output
        return str(res).replace("OrderedDict(", "ParticleHandle(")

    def add_exclusion(self, partner):
        """
        Exclude non-bonded interactions with the given partner.

        .. note::
            This needs the feature ``EXCLUSIONS``.

        Parameters
        -----------
        partner : :class:`~espressomd.particle_data.ParticleHandle` or :obj:`int`
            Particle to exclude.

        """
        if isinstance(partner, ParticleHandle):
            p_id = partner.id
        else:
            p_id = partner
        check_type_or_throw_except(
            p_id, 1, int, "Argument 'partner' has to be a ParticleHandle or int.")
        if self.call_method("has_exclusion", pid=p_id):
            raise RuntimeError(
                f"Particle with id {p_id} is already in exclusion list of particle with id {self.id}")
        self.call_method("add_exclusion", pid=p_id)

    def delete_exclusion(self, partner):
        """
        Remove exclusion of non-bonded interactions with the given partner.

        .. note::
            This needs the feature ``EXCLUSIONS``.

        Parameters
        -----------
        partner : :class:`~espressomd.particle_data.ParticleHandle` or :obj:`int`
            Particle to remove from exclusions.

        """
        if isinstance(partner, ParticleHandle):
            p_id = partner.id
        else:
            p_id = partner
        check_type_or_throw_except(
            p_id, 1, int, "Argument 'partner' has to be a ParticleHandle or int.")
        if not self.call_method("has_exclusion", pid=p_id):
            raise RuntimeError(
                f"Particle with id {p_id} is not in exclusion list of particle with id {self.id}")
        self.call_method("del_exclusion", pid=p_id)

    @property
    def bonds(self):
        """
        The bonds stored by this particle. Note that bonds are only stored by
        one partner. You need to define a bonded interaction.

        A bond tuple is specified as a bond identifier associated with
        a particle ``(bond_ID, (*part_ID,))``. A single particle may contain
        multiple bonds.

        Type: Ragged array.

        .. note::
           Bond ids have to be an integer >= 0.

        See Also
        --------
        espressomd.particle_data.ParticleHandle.add_bond : Method to add bonds to a ``Particle``
        espressomd.particle_data.ParticleHandle.delete_bond : Method to remove bonds from a ``Particle``

        """
        bonds = []
        for bond_view in self.call_method("get_bonds_view"):
            bond_id = bond_view[0]
            partner_ids = bond_view[1:]
            bonds.append((BondedInteractions()[bond_id], *partner_ids))

        return tuple(bonds)

    @bonds.setter
    def bonds(self, bonds):
        # Assigning to the bond property means replacing the existing value
        # i.e., we delete all existing bonds
        self.delete_all_bonds()

        if bonds:
            nlvl = nesting_level(bonds)
            if nlvl == 1:
                self.add_bond(bonds)
            elif nlvl == 2:
                for bond in bonds:
                    self.add_bond(bond)
            else:
                raise ValueError(
                    "Bonds have to specified as lists of tuples/lists or a single list.")

    @property
    def exclusions(self):
        """
        The exclusion list of particles where non-bonded interactions are ignored.

        .. note::
            This needs the feature ``EXCLUSIONS``.

        Type: (N,) array_like of :obj:`int`

        """
        assert_features("EXCLUSIONS")
        return array_locked(
            np.array(self.call_method("get_exclusions"), dtype=int))

    @exclusions.setter
    def exclusions(self, p_ids):
        assert_features("EXCLUSIONS")
        self.call_method("set_exclusions", p_ids=p_ids)

    def vs_auto_relate_to(self, rel_to):
        """
        Setup this particle as virtual site relative to the particle
        in argument ``rel_to``. A particle cannot relate to itself.

        Parameters
        -----------
        rel_to : :obj:`int` or :obj:`ParticleHandle`
            Particle to relate to (either particle id or particle object).

        """
        if isinstance(rel_to, ParticleHandle):
            rel_to = rel_to.id
        else:
            check_type_or_throw_except(
                rel_to, 1, int, "Argument of 'vs_auto_relate_to' has to be of type ParticleHandle or int")
        self.call_method("vs_relate_to", pid=rel_to)
        handle_errors("vs_auto_relate_to")

    def add_verified_bond(self, bond):
        """
        Add a bond, the validity of which has already been verified.

        See Also
        --------
        add_bond : Add an unverified bond to the ``Particle``.
        bonds : ``Particle`` property containing a list of all current bonds held by ``Particle``.

        """
        if self.id in bond[1:]:
            raise Exception(
                f"Bond partners {bond[1:]} include the particle {self.id} itself")
        self.call_method("add_bond",
                         bond_id=bond[0]._bond_id,
                         part_id=bond[1:])

    def delete_verified_bond(self, bond):
        """
        Delete a single bond from the particle. The validity of which has already been verified.

        Parameters
        ----------
        bond : :obj:`tuple`
            tuple where the first element is either a bond ID of a bond type,
            and the last element is the ID of the partner particle to be bonded
            to.

        See Also
        --------
        delete_bond : Delete an unverified bond held by the ``Particle``.
        bonds : ``Particle`` property containing a list of all current bonds held by ``Particle``.

        """
        self.call_method("del_bond",
                         bond_id=bond[0]._bond_id,
                         part_id=bond[1:])

    def normalize_and_check_bond_or_throw_exception(self, bond):
        """
        Checks the validity of the given bond:

        - If the bondtype is given as an object or a numerical id
        - If all partners are of type :obj:`int`
        - If the number of partners satisfies the bond
        - If the bond type used exists (is lower than ``n_bonded_ia``)
        - If the number of bond partners fits the bond type

        Throws an exception if any of these are not met.

        Normalize the bond, i.e. replace bond ids by bond objects and particle
        objects by particle ids.

        """
        # Has it []-access
        if not hasattr(bond, "__getitem__"):
            raise ValueError(
                "Bond needs to be a tuple or list containing bond type and partners.")

        bond = list(bond)
        # Bond type or numerical bond id
        if is_valid_type(bond[0], int):
            bond[0] = BondedInteractions()[bond[0]]
        elif not isinstance(bond[0], BondedInteraction):
            raise Exception(
                f"1st element of Bond has to be of type BondedInteraction or int, got {type(bond[0])}")
        # Check the bond is in the list of active bonded interactions
        if bond[0]._bond_id == -1:
            raise Exception(
                "The bonded interaction has not yet been added to the list of active bonds in ESPResSo")
        # Validity of the numeric id
        if not self.call_method("is_valid_bond_id", bond_id=bond[0]._bond_id):
            raise ValueError(
                f"The bond type {bond[0]._bond_id} does not exist.")
        # Number of partners
        expected_num_partners = bond[0].call_method("get_num_partners")
        if len(bond) - 1 != expected_num_partners:
            raise ValueError(
                f"Bond {bond[0]} needs {expected_num_partners} partners")
        # Type check on partners
        for i in range(1, len(bond)):
            if isinstance(bond[i], ParticleHandle):
                # Put the particle id instead of the particle handle
                bond[i] = bond[i].id
            elif not is_valid_type(bond[i], int):
                raise ValueError(
                    "Bond partners have to be of type integer or ParticleHandle.")
        return tuple(bond)

    def add_bond(self, bond):
        """
        Add a single bond to the particle.

        Parameters
        ----------
        bond : :obj:`tuple`
            tuple where the first element is either a bond ID or a bond object,
            and the next elements are particle ids or particle objects to be
            bonded to.

        See Also
        --------
        bonds : ``Particle`` property containing a list of all current bonds held by ``Particle``.

        Examples
        --------
        >>> import espressomd.interactions
        >>>
        >>> system = espressomd.System(box_l=3 * [10])
        >>>
        >>> # define a harmonic potential and add it to the system
        >>> harm_bond = espressomd.interactions.HarmonicBond(r_0=1, k=5)
        >>> system.bonded_inter.add(harm_bond)
        >>>
        >>> # add two particles
        >>> p1 = system.part.add(pos=(1, 0, 0))
        >>> p2 = system.part.add(pos=(2, 0, 0))
        >>>
        >>> # bond them via the bond type
        >>> p1.add_bond((harm_bond, p2))
        >>> # or via the bond index (zero in this case since it is the first one added)
        >>> p1.add_bond((0, p2))

        """
        _bond = self.normalize_and_check_bond_or_throw_exception(bond)
        if _bond in self.bonds:
            raise RuntimeError(
                f"Bond {_bond} already exists on particle {self.id}")
        self.add_verified_bond(_bond)

    def delete_bond(self, bond):
        """
        Delete a single bond from the particle.

        Parameters
        ----------
        bond : :obj:`tuple`
            tuple where the first element is either a bond ID or a bond object,
            and the next elements are particle ids or particle objects that are
            bonded to.

        See Also
        --------
        bonds : ``Particle`` property containing a list of all bonds currently held by ``Particle``.


        Examples
        --------

        >>> import espressomd.interactions
        >>>
        >>> system = espressomd.System(box_l=3 * [10])
        >>>
        >>> # define a harmonic potential and add it to the system
        >>> harm_bond = espressomd.interactions.HarmonicBond(r_0=1, k=5)
        >>> system.bonded_inter.add(harm_bond)
        >>>
        >>> # bond two particles to the first one
        >>> p0 = system.part.add(pos=(1, 0, 0))
        >>> p1 = system.part.add(pos=(2, 0, 0))
        >>> p2 = system.part.add(pos=(1, 1, 0))
        >>> p0.add_bond((harm_bond, p1))
        >>> p0.add_bond((harm_bond, p2))
        >>>
        >>> print(p0.bonds)
        ((HarmonicBond(0): {'r_0': 1.0, 'k': 5.0, 'r_cut': 0.0}, 1),
         (HarmonicBond(0): {'r_0': 1.0, 'k': 5.0, 'r_cut': 0.0}, 2))
        >>> # delete the first bond
        >>> p0.delete_bond(p0.bonds[0])
        >>> print(p0.bonds)
        ((HarmonicBond(0): {'r_0': 1.0, 'k': 5.0, 'r_cut': 0.0}, 2),)

        """

        _bond = self.normalize_and_check_bond_or_throw_exception(bond)
        if _bond not in self.bonds:
            raise RuntimeError(
                f"Bond {_bond} doesn't exist on particle {self.id}")
        self.delete_verified_bond(_bond)

    def update(self, new_properties):
        """
        Update properties of a particle.

        Parameters
        ----------
        new_properties : :obj:`dict`
            Map particle property names to values. All properties except
            for the particle id can be changed.

        Examples
        --------

        >>> import espressomd
        >>> system = espressomd.System(box_l=[10, 10, 10])
        >>> p = system.part.add(pos=[1, 2, 3], q=1, virtual=True)
        >>> print(p.pos, p.q, p.virtual)
        [1. 2. 3.] 1.0 True
        >>> p.update({'pos': [4, 5, 6], 'virtual': False, 'q': 0})
        >>> print(p.pos, p.q, p.virtual)
        [4. 5. 6.] 0.0 False

        """
        if "id" in new_properties:
            raise RuntimeError("Cannot change particle id.")

        for k, v in new_properties.items():
            setattr(self, k, v)

    def convert_vector_body_to_space(self, vec):
        """
        Convert the given vector from the particle's body frame to the space frame.
        """
        assert_features("ROTATION")
        return self.call_method("convert_vector_body_to_space", vec=vec)

    def convert_vector_space_to_body(self, vec):
        """
        Convert the given vector from the space frame to the particle's body frame.
        """
        assert_features("ROTATION")
        return self.call_method("convert_vector_space_to_body", vec=vec)

    def rotate(self, axis, angle):
        """
        Rotate the particle around the given axis.

        Parameters
        ----------
        axis : (3,) array_like of :obj:`float`

        angle : :obj:`float`

        """
        assert_features("ROTATION")
        self.call_method("rotate_particle", axis=axis, angle=angle)


particle_attributes = set(ParticleHandle(id=0)._valid_parameters())
if has_features("EXCLUSIONS"):
    particle_attributes.add("exclusions")
particle_attributes.add("bonds")


@script_interface_register
class ParticleSlice(ScriptInterfaceHelper):
    """
    Handle slice inputs. Set values for selected slices or
    return values as a single list.

    """
    _so_name = "Particles::ParticleSlice"
    _so_creation_policy = "LOCAL"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if "sip" not in kwargs:
            for p_id in self.id_selection:
                if not self.call_method("particle_exists", p_id=p_id):
                    raise IndexError(f"Particle does not exist: {p_id}")

    def __iter__(self):
        return self._id_gen()

    def _id_gen(self):
        """
        Generator for chunked and prefetched iteration of particles.
        """
        for chunk in self.chunks(self.id_selection, self.chunk_size):
            self.call_method("prefetch_particle_data", chunk=chunk)
            for i in chunk:
                yield ParticleHandle(id=i)

    def chunks(self, l, n):
        """
        Generator returning chunks of length n from l.
        """
        for i in range(0, len(l), n):
            yield l[i:i + n]

    def __len__(self):
        return len(self.id_selection)

    @property
    def pos_folded(self):
        """
        Particle position (folded into central image).

        """
        pos_array = np.zeros((len(self.id_selection), 3))
        for i in range(len(self.id_selection)):
            pos_array[i, :] = ParticleHandle(
                id=self.id_selection[i]).pos_folded
        return pos_array

    @pos_folded.setter
    def pos_folded(self, value):
        raise RuntimeError("Parameter 'pos_folded' is read-only.")

    def add_exclusion(self, _partner):
        assert_features(["EXCLUSIONS"])
        for p_id in self.id_selection:
            ParticleHandle(id=p_id).add_exclusion(_partner)

    def delete_exclusion(self, _partner):
        assert_features(["EXCLUSIONS"])
        for p_id in self.id_selection:
            ParticleHandle(id=p_id).delete_exclusion(_partner)

    def __str__(self):
        return "ParticleSlice([" + \
            ", ".join(str(ParticleHandle(id=i))
                      for i in self.id_selection) + "])"

    def update(self, new_properties):
        if "id" in new_properties:
            raise RuntimeError("Cannot change particle id.")

        for k, v in new_properties.items():
            setattr(self, k, v)

    # Bond related methods
    def add_bond(self, _bond):
        """
        Add a single bond to the particles.

        """
        for p_id in self.id_selection:
            ParticleHandle(id=p_id).add_bond(_bond)

    def delete_bond(self, _bond):
        """
        Delete a single bond from the particles.

        """
        for p_id in self.id_selection:
            ParticleHandle(id=p_id).delete_bond(_bond)

    def delete_all_bonds(self):
        for p_id in self.id_selection:
            ParticleHandle(id=p_id).delete_all_bonds()

    def remove(self):
        """
        Delete the particles.

        See Also
        --------
        :meth:`espressomd.particle_data.ParticleList.add`

        """
        for p_id in self.id_selection:
            ParticleHandle(id=p_id).remove()

    def __setattr__(self, name, value):
        if name != "chunk_size" and name != "id_selection" and name not in particle_attributes:
            raise AttributeError(
                f"ParticleHandle does not have the attribute {name}.")
        super().__setattr__(name, value)

    def to_dict(self):
        """
        Returns the particles attributes as a dictionary.

        It can be used to save the particle data and recover it by using

        >>> p = system.part.add(...)
        >>> particle_dict = p.to_dict()
        >>> system.part.add(particle_dict)

        It includes the content of ``particle_attributes``, minus a few exceptions:

        - :attr:`~ParticleHandle.dip`, :attr:`~ParticleHandle.director`:
          Setting only the director will overwrite the orientation of the
          particle around the axis parallel to dipole moment/director.
          Quaternions contain the full info.
        - :attr:`~ParticleHandle.image_box`, :attr:`~ParticleHandle.node`

        """

        odict = {}
        for p in self:
            pdict = ParticleHandle(id=p.id).to_dict()
            for p_key, p_value in pdict.items():
                if p_key in odict:
                    odict[p_key].append(p_value)
                else:
                    odict[p_key] = [p_value]
        return odict


@script_interface_register
class ParticleList(ScriptInterfaceHelper):
    """
    Provides access to the particles.

    Methods
    -------
    clear()
        Remove all particles.

        See Also
        --------
        :meth:`espressomd.particle_data.ParticleHandle.remove`

    auto_exclusions()
        Add exclusions between particles that are connected by pair bonds,
        including virtual bonds. Angle and dihedral bonds are ignored. The most
        common use case for this method is to auto-exclude virtual sites.

        Another use case is to exclude 1-2, 1-3 and optionally 1-4 non-nonded
        interactions on polymer chains. This technique is commonly used in
        atomistic molecular dynamics engines such as NAMD, AMBER or GROMACS,
        where the short-range part of the potential energy surface is better
        approximated with Fourier sums (using dihedral bonds) than with pair
        potentials. Linear, branched and circular topologies are supported.

        Requires feature ``EXCLUSIONS``.

        Parameters
        ----------
        distance : :obj:`int`
            Maximal length of a chain in unit of bonds. The topology
            will be traversed recursively until the bond chain either
            terminates or reaches that distance.

    """
    _so_name = "Particles::ParticleList"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = (
        "clear", "auto_exclusions"
    )

    def by_id(self, p_id):
        """
        Access a particle by its integer id.
        """
        return ParticleHandle(id=p_id)

    def by_ids(self, ids):
        """
        Get a slice of particles by their integer ids.
        """
        return ParticleSlice(id_selection=ids)

    def all(self):
        """
        Get a slice containing all particles.
        """
        all_ids = self.call_method("get_particle_ids")
        return self.by_ids(all_ids)

    def __len__(self):
        return self.call_method("get_n_part")

    @property
    def highest_particle_id(self):
        """
        Largest particle id.
        """
        return self.call_method("get_highest_particle_id")

    def add(self, *args, **kwargs):
        """
        Adds one or several particles to the system

        Parameters
        ----------
        Either a dictionary or a bunch of keyword args.

        Returns
        -------
        Returns an instance of :class:`espressomd.particle_data.ParticleHandle` for each added particle.

        See Also
        --------
        :meth:`espressomd.particle_data.ParticleHandle.remove`

        Examples
        --------

        >>> import espressomd
        >>> system = espressomd.System(box_l=[10, 10, 10])
        >>> # add two particles
        >>> system.part.add(id=0, pos=(1, 0, 0))
        >>> system.part.add(id=1, pos=(2, 0, 0))

        ``pos`` is mandatory, ``id`` can be omitted, in which case it is assigned automatically.
        Several particles can be added by passing one value per particle to each property::

            system.part.add(pos=((1, 2, 3), (4, 5, 6)), q=(1, -1))

        """

        # Did we get a dictionary
        if len(args) == 1 and isinstance(
                args[0], (dict, collections.OrderedDict)):
            particles_dict = args[0]
        else:
            if len(args) == 0 and len(kwargs) != 0:
                particles_dict = kwargs
            else:
                raise ValueError(
                    "add() takes either a dictionary or a bunch of keyword args.")

        # Check for presence of pos attribute
        if "pos" not in particles_dict:
            raise ValueError(
                "pos attribute must be specified for new particle")

        if len(np.array(particles_dict["pos"]).shape) == 2:
            return self._place_new_particles(particles_dict)
        else:
            return self._place_new_particle(particles_dict)

    def _place_new_particle(self, p_dict):
        bonds = []
        if "bonds" in p_dict:
            bonds = p_dict.pop("bonds")
            if nesting_level(bonds) == 1:
                bonds = [bonds]
        p_id = self.call_method("add_particle", **p_dict)
        p = self.by_id(p_id)
        for bond in bonds:
            if len(bond):
                bond = p.normalize_and_check_bond_or_throw_exception(bond)
                p.add_verified_bond(bond)
        return p

    def _place_new_particles(self, p_list_dict):
        # Check if all entries have the same length
        n_parts = len(p_list_dict["pos"])
        if not all(np.array(v, dtype=object).shape and len(v) ==
                   n_parts for v in p_list_dict.values()):
            raise ValueError(
                "When adding several particles at once, all lists of attributes have to have the same size")

        # If particle ids haven't been provided, use free ones
        # beyond the highest existing one
        if "id" not in p_list_dict:
            first_id = self.highest_particle_id + 1
            p_list_dict["id"] = np.arange(first_id, first_id + n_parts)

        # Place the particles
        for i in range(n_parts):
            p_dict = {k: v[i] for k, v in p_list_dict.items()}
            self._place_new_particle(p_dict)

        # Return slice of added particles
        return self.by_ids(p_list_dict["id"])

    # Iteration over all existing particles
    def __iter__(self):
        for p_id in self.call_method("get_particle_ids"):
            yield self.by_id(p_id)

    def exists(self, idx):
        if is_valid_type(idx, int):
            return self.call_method("particle_exists", p_id=idx)
        if isinstance(idx, (slice, tuple, list, np.ndarray)):
            tf_array = np.zeros(len(idx), dtype=type(True))
            for i in range(len(idx)):
                tf_array[i] = self.call_method("particle_exists", p_id=idx[i])
            return tf_array

    def __str__(self):
        return "ParticleList([" + \
            ",".join(map(str, self.call_method("get_particle_ids"))) + "])"

    def writevtk(self, fname, types='all'):
        """
        Write the positions and velocities of particles with specified
        types to a VTK file.

        Parameters
        ----------
        fname: :obj:`str`
            Filename of the target output file
        types: list of :obj:`int` or the string 'all', optional (default: 'all')
            A list of particle types which should be output to 'fname'

        Examples
        --------

        >>> import espressomd
        >>> system = espressomd.System(box_l=[10, 10, 10])
        >>> # add several particles
        >>> system.part.add(pos=0.5 * system.box_l, v=[1, 0, 0], type=0)
        >>> system.part.add(pos=0.4 * system.box_l, v=[0, 2, 0], type=1)
        >>> system.part.add(pos=0.7 * system.box_l, v=[2, 0, 1], type=1)
        >>> system.part.add(pos=0.1 * system.box_l, v=[0, 0, 1], type=2)
        >>> # write to VTK
        >>> system.part.writevtk("part_type_0_1.vtk", types=[0, 1])
        >>> system.part.writevtk("part_type_2.vtk", types=[2])
        >>> system.part.writevtk("part_all.vtk")

        .. todo:: move to ``./io/writer/``

        """

        if not hasattr(types, '__iter__'):
            types = [types]

        n = 0
        for p in self:
            if types == 'all' or p.type in types:
                n += 1

        with open(fname, "w") as vtk:
            vtk.write("# vtk DataFile Version 2.0\n")
            vtk.write("particles\n")
            vtk.write("ASCII\n")
            vtk.write("DATASET UNSTRUCTURED_GRID\n")
            vtk.write("POINTS {} floats\n".format(n))
            for p in self:
                if types == 'all' or p.type in types:
                    vtk.write("{} {} {}\n".format(*(p.pos_folded)))

            vtk.write("POINT_DATA {}\n".format(n))
            vtk.write("SCALARS velocity float 3\n")
            vtk.write("LOOKUP_TABLE default\n")
            for p in self:
                if types == 'all' or p.type in types:
                    vtk.write("{} {} {}\n".format(*p.v))

    def pairs(self):
        """
        Generate all pairs of particles.

        """

        ids = self.call_method("get_particle_ids")
        id_pairs = itertools.combinations(ids, 2)
        for id_pair in id_pairs:
            yield (self.by_id(id_pair[0]), self.by_id(id_pair[1]))

    def select(self, *args, **kwargs):
        """
        Generate a particle slice by filtering particles via a user-defined
        criterion.

        Parameters:

        Either: a keyword arguments in which the keys are names of particle
        properties and the values are the values to filter for. E.g.,::

            system.part.select(type=0, q=1)

        Or: a function taking a ParticleHandle as argument and returning True if
        the particle is to be filtered for. E.g.,::

            system.part.select(lambda p: p.pos[0] < 0.5)

        Returns
        -------
        :class:`ParticleSlice` :
            An instance of :class:`ParticleSlice` containing the selected particles

        """

        # Ids of the selected particles
        ids = []
        # Did we get a function as argument?
        if len(args) == 1 and len(kwargs) == 0 and callable(args[0]):
            # Go over all particles and pass them to the user-provided function
            for p in self:
                if args[0](p):
                    ids.append(p.id)
            return ParticleSlice(id_selection=ids)

        # Did we get a set of keyword args?
        elif len(args) == 0:
            for p in self:
                select = True
                # Check, if the particle fails any required criteria
                for k in kwargs:
                    # Fetch user-provided value and value in particle
                    val1 = kwargs[k]
                    val2 = getattr(p, k)
                    # Get tolerance from numerical accuracy limits
                    tol = max(
                        np.amax(np.spacing(val1)), np.amax(np.spacing(val2)))

                    # Compare
                    if not np.allclose(val1, val2, atol=tol):
                        select = False
                        break
                if select:
                    ids.append(p.id)
            return ParticleSlice(id_selection=ids)
        else:
            raise Exception(
                "select() takes either selection function as positional argument or a set of keyword arguments.")


def set_slice_one_for_all(particle_slice, attribute, values):
    for i in particle_slice.id_selection:
        setattr(ParticleHandle(id=i), attribute, values)


def set_slice_one_for_each(particle_slice, attribute, values):
    for i, v in zip(particle_slice.id_selection, values):
        setattr(ParticleHandle(id=i), attribute, v)


def _add_particle_slice_properties():
    """
    Automatically add all of ParticleHandle's properties to ParticleSlice.

    """

    def set_attribute(particle_slice, values, attribute):
        """
        Setter function that sets attribute on every member of particle_slice.
        If values contains only one element, all members are set to it. If it
        contains as many elements as there are members, each of them gets set
        to the corresponding one. For attributes that are lists of various length,
        (bonds, exclusions) the nesting level decides if it is one-for-all or one-for-each.

        """

        N = len(particle_slice.id_selection)

        if N == 0:
            raise AttributeError(
                "Cannot set properties of an empty ParticleSlice")

        # Special attributes
        if attribute == "bonds":
            nlvl = nesting_level(values)
            if nlvl == 1 or nlvl == 2:
                set_slice_one_for_all(particle_slice, attribute, values)
            elif nlvl == 3 and len(values) == N:
                set_slice_one_for_each(particle_slice, attribute, values)
            else:
                raise Exception("Failed to set bonds for particle slice.")

            return

        elif attribute == "exclusions":
            nlvl = nesting_level(values)
            if nlvl == 0 or nlvl == 1:
                set_slice_one_for_all(particle_slice, attribute, values)
            elif nlvl == 2 and len(values) == N:
                set_slice_one_for_each(particle_slice, attribute, values)
            else:
                raise Exception("Failed to set exclusions for particle slice.")

            return

        elif attribute == "vs_relative":
            nlvl = nesting_level(values)
            if nlvl in [1, 2]:
                set_slice_one_for_all(particle_slice, attribute, values)
            elif nlvl == 3 and len(values) == N:
                set_slice_one_for_each(particle_slice, attribute, values)
            else:
                raise Exception(
                    "Failed to set vs_relative for particle slice.")

            return

        else:
            target = getattr(
                ParticleHandle(id=particle_slice.id_selection[0]), attribute)
            target_shape = np.shape(target)

            if not target_shape:  # scalar quantity
                if not np.shape(values):
                    set_slice_one_for_all(particle_slice, attribute, values)
                elif np.shape(values)[0] == N:
                    set_slice_one_for_each(particle_slice, attribute, values)
                else:
                    raise Exception(
                        f"Value shape {np.shape(values)} does not broadcast to attribute shape {target_shape}.")

                return

            else:  # fixed length vector quantity
                if target_shape == np.shape(values):
                    set_slice_one_for_all(particle_slice, attribute, values)
                elif target_shape == tuple(np.shape(values)[1:]) and np.shape(values)[0] == N:
                    set_slice_one_for_each(particle_slice, attribute, values)
                else:
                    raise Exception(
                        f"Value shape {np.shape(values)} does not broadcast to attribute shape {target_shape}.")

                return

    def get_attribute(particle_slice, attribute):
        """
        Getter function that copies attribute from every member of
        particle_slice into an array (if possible).
        For special properties, a tuple of tuples is used.

        """

        N = len(particle_slice.id_selection)
        if N == 0:
            return np.empty(0, dtype=type(None))

        # get first slice member to determine its type
        target = getattr(ParticleHandle(
            id=particle_slice.id_selection[0]), attribute)
        if isinstance(target, array_locked):  # vectorial quantity
            target_type = target.dtype
        else:  # scalar quantity
            target_type = type(target)

        if attribute in ["exclusions", "bonds", "vs_relative", "swimming"]:
            values = []
            for part in particle_slice._id_gen():
                values.append(getattr(part, attribute))
        else:
            values = np.empty((N,) + np.shape(target), dtype=target_type)
            i = 0
            for part in particle_slice._id_gen():
                values[i] = getattr(part, attribute)
                i += 1

        return values

    for attribute_name in sorted(particle_attributes):
        if attribute_name in dir(ParticleSlice):
            continue

        # synthesize a new property
        new_property = property(
            functools.partial(get_attribute, attribute=attribute_name),
            functools.partial(set_attribute, attribute=attribute_name),
            doc="")
        # attach the property to ParticleSlice
        setattr(ParticleSlice, attribute_name, new_property)


_add_particle_slice_properties()
