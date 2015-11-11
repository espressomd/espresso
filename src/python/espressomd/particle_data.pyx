#
# Copyright (C) 2013,2014 The ESPResSo project
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
cimport numpy as np
import numpy as np
cimport utils
from utils cimport *
cimport particle_data
from interactions import BondedInteraction

PARTICLE_EXT_FORCE = 1


def COORD_FIXED(coord):
    return 2L << coord
COORDS_FIX_MASK = COORD_FIXED(0) | COORD_FIXED(1) | COORD_FIXED(2)
COORDS_ALL_FIXED = COORD_FIXED(0) & COORD_FIXED(1) & COORD_FIXED(2)
PARTICLE_EXT_TORQUE = 16

cdef class ParticleHandle:
    def __cinit__(self, _id):
        #    utils.init_intlist(self.particleData.el)
        utils.init_intlist( & (self.particleData.bl))
        self.id = _id

    cdef int updateParticleData(self) except -1:
        #    utils.realloc_intlist(self.particleData.el, 0)
        utils.realloc_intlist( & (self.particleData.bl), 0)

        if get_particle_data(self.id, & self.particleData):
            raise Exception("Error updating particle data")
        else:
            return 0

    # The individual attributes of a particle are implemented as properties.

    # Particle Type
    property type:
        """Particle type"""

        def __set__(self, _type):
            if isinstance(_type, int) and _type >= 0:
                if set_particle_type(self.id, _type) == 1:
                    raise Exception("set particle position first")
            else:
                raise ValueError("type must be an integer >= 0")

        def __get__(self):
            self.updateParticleData()
            return self.particleData.p.type

    # Position
    property pos:
        """Particle position (not folded into central image)."""

        def __set__(self, _pos):
            cdef double mypos[3]
            checkTypeOrExcept(_pos, 3, float, "Postion must be 3 floats")
            for i in range(3):
                mypos[i] = _pos[i]
            if place_particle(self.id, mypos) == -1:
                raise Exception("particle could not be set")

        def __get__(self):
            self.updateParticleData()
            return np.array([self.particleData.r.p[0],
                             self.particleData.r.p[1],
                             self.particleData.r.p[2]])

    # Velocity
    property v:
        """Particle velocity"""

        def __set__(self, _v):
            cdef double myv[3]
            checkTypeOrExcept(_v, 3, float, "Velocity has to be floats")
            for i in range(3):
                myv[i] = _v[i]
            if set_particle_v(self.id, myv) == 1:
                raise Exception("set particle position first")

        def __get__(self):
            self.updateParticleData()
            return np.array([self.particleData.m.v[0],
                             self.particleData.m.v[1],
                             self.particleData.m.v[2]])

    # Force
    property f:
        """Particle force"""

        def __set__(self, _f):
            cdef double myf[3]
            checkTypeOrExcept(_f, 3, float, "Force has to be floats")
            for i in range(3):
                myf[i] = _f[i]
            if set_particle_f(self.id, myf) == 1:
                raise Exception("set particle position first")

        def __get__(self):
            self.updateParticleData()
            return np.array([self.particleData.f.f[0],
                             self.particleData.f.f[1],
                             self.particleData.f.f[2]])

    # Bonds
    property bonds:
        """Bond partners with respect to bonded interactions."""

        def __set__(self, _bonds):
            # First, we check that we got a list/tuple.
            if not hasattr(_bonds, "__getitem__"):
                raise ValueError(
                    "bonds have to specified as a tuple of tuples. (Lists can also be used)")
            # Check individual bonds
            for bond in _bonds:
                self.checkBondOrThrowException(bond)

            # Assigning to the bond property means replacing the existing value
            # i.e., we delete all existing bonds
            if change_particle_bond(self.id, NULL, 1):
                raise Exception("Deleting existing bonds failed.")

            # And add the new ones
            for bond in _bonds:
                self.addVerifiedBond(bond)

        def __get__(self):
            self.updateParticleData()
            bonds = []
            # Go through the bond list of the particle
            i = 0
            while i < self.particleData.bl.n:
                bond = []
                # Bond type:
                bondId = self.particleData.bl.e[i]
                bond.append(bondId)
                # Number of partners
                nPartners = bonded_ia_params[bondId].num

                i += 1

                # Copy bond partners
                for j in range(nPartners):
                    bond.append(self.particleData.bl.e[i])
                    i += 1
                bonds.append(tuple(bond))

            return tuple(bonds)

    # Properties that exist only when certain features are activated
    # MASS
    IF MASS == 1:
        property mass:
            """Particle mass"""

            def __set__(self, _mass):
                checkTypeOrExcept(_mass, 1, float, "Mass has to be 1 floats")
                if set_particle_mass(self.id, _mass) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.updateParticleData()
                cdef double * x = NULL
                pointer_to_mass( & (self.particleData), x)
                return x[0]

    IF ROTATION == 1:
        # Omega (angular velocity) lab frame
        property omega_lab:
            """Angular velocity in lab frame"""

            def __set__(self, _o):
                cdef double myo[3]
                checkTypeOrExcept(_o, 3, float, "Omega_lab has to be 3 floats")
                for i in range(3):
                    myo[i] = _o[i]
                if set_particle_omega_lab(self.id, myo) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.updateParticleData()
                cdef double o[3]
                convert_omega_body_to_space(& (self.particleData), o)
                return np.array([o[0], o[1], o[2]])

    # ROTATIONAL_INERTIA
    IF ROTATIONAL_INERTIA == 1:
        property rinertia:
            """Rotational inertia"""

            def __set__(self, _rinertia):
                cdef double rinertia[3]
                checkTypeOrExcept(
                    _rinertia, 3, float, "Rotation_inertia has to be 3 floats")
                for i in range(3):
                    rinertia[i] = _rinertia[i]
                if set_particle_rotational_inertia(self.id, rinertia) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.updateParticleData()
                cdef double * rinertia = NULL
                pointer_to_rotational_inertia( & (self.particleData), rinertia)
                return np.array([rinertia[0], rinertia[1], rinertia[2]])

# Omega (angular velocity) body frame
        property omega_body:
            """Angular velocity in body frame"""

            def __set__(self, _o):
                cdef double myo[3]
                checkTypeOrExcept(
                    _o, 3, float, "Omega_body has to be 3 floats")
                for i in range(3):
                    myo[i] = _o[i]
                if set_particle_omega_body(self.id, myo) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.updateParticleData()
                cdef double * o = NULL
                pointer_to_omega_body( & (self.particleData), o)
                return np.array([o[0], o[1], o[2]])


# Torque in lab frame
        property torque_lab:
            """Torque in lab frame"""

            def __set__(self, _t):
                cdef double myt[3]
                checkTypeOrExcept(_t, 3, float, "Torque has to be 3 floats")
                for i in range(3):
                    myt[i] = _t[i]
                if set_particle_torque_lab(self.id, myt) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.updateParticleData()
                cdef double x[3]
                convert_torques_body_to_space(& (self.particleData), x)
                return np.array([x[0], x[1], x[2]])

# Quaternion
        property quat:
            """Quaternions"""

            def __set__(self, _q):
                cdef double myq[4]
                checkTypeOrExcept(
                    _q, 4, float, "Quaternions has to be 4 floats")
                for i in range(4):
                    myq[i] = _q[i]
                if set_particle_quat(self.id, myq) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.updateParticleData()
                cdef double * x = NULL
                pointer_to_quat(& (self.particleData), x)
                return np.array([x[0], x[1], x[2], x[3]])
# Director ( z-axis in body fixed frame)
        property director:
            """Director"""

            def __set__(self, _q):
                raise Exception(
                    "Setting the director is not implemented in the c++-core of Espresso")
#        cdef double myq[3]
#        checkTypeOrExcept(_q,3,float,"Director has to be 3 floats")
#        for i in range(3):
#            myq[i]=_q[i]
#        if set_particle_quatu(self.id, myq) == 1:
#          raise Exception("set particle position first")

            def __get__(self):
                self.updateParticleData()
                cdef double * x = NULL
                pointer_to_quatu(& (self.particleData), x)
                return np.array([x[0], x[1], x[2]])

# Charge
    IF ELECTROSTATICS == 1:
        property q:
            """particle charge"""

            def __set__(self, _q):
                cdef double myq
                checkTypeOrExcept(_q, 1, float, "Charge has to be floats")
                myq = _q
                if set_particle_q(self.id, myq) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.updateParticleData()
                cdef double * x = NULL
                pointer_to_q(& (self.particleData), x)
                return x[0]

    def delete(self):
        """Delete the particle"""
        if remove_particle(self.id):
            raise Exception("Could not delete particle")
        del self

    IF VIRTUAL_SITES == 1:
        # virtual flag

        property virtual:
            """virtual flag"""

            def __set__(self, _v):
                if isinstance(_v, int):
                    if set_particle_virtual(self.id, _v) == 1:
                        raise Exception("set particle position first")
                else:
                    raise ValueError("virtual must be an integer >= 0")

            def __get__(self):
                self.updateParticleData()
                cdef int * x = NULL
                pointer_to_virtual(& (self.particleData), x)
                return x[0]

    IF VIRTUAL_SITES_RELATIVE == 1:
        property vs_relative:
            """virtual sites relative parameters"""

            def __set__(self, x):
                if len(x) != 3:
                    raise ValueError("vs_relative needs six args")
                _relto = x[0]
                _dist = x[1]
                q = x[3]
                checkTypeOrExcept(
                    q, 4, float, "The relative orientation has to be specified as quaternion with 4 floats.")
                cdef double _q[4]
                for i in range(4):
                    _q[i] = q[i]

                if isinstance(_relto, int) and isinstance(_dist, float):
                    if set_particle_vs_relative(self.id, _relto, _dist, _q) == 1:
                        raise Exception("set particle position first")
                else:
                    raise ValueError(
                        "vs_relative takes one int and one float as parameters.")

            def __get__(self):
                self.updateParticleData()
                cdef int * rel_to = NULL
                cdef double * dist = NULL
                cdef double * q = NULL
                pointer_to_vs_relative( & (self.particleData), rel_to, dist, q)
                return (rel_to[0], dist[0], np.array((q[0], q[1], q[2], q[3])))

        # vs_auto_relate_to
        def vs_auto_relate_to(self, _relto):
            """Setup this particle as virtual site relative to the particle with the given id"""
            if isinstance(_relto, int):
                if vs_relate_to(self.id, _relto):
                    raise Exception("Vs_relative setup failed.")
            else:
                raise ValueError(
                    "Argument of vs_auto_relate_to has to be of type int")

                # Virtual sites relative parameters

        # vs_auto_relate_to
        def vs_auto_relate_to(self, _relto):
            """Setup this particle as virtual site relative to the particle with the given id"""
            checkTypeOrExcept(
                _relto, 1, int, "Argument of vs_auto_relate_to has to be of type int")
            if vs_relate_to(self.id, _relto):
                raise Exception("vs_relative setup failed.")

    IF DIPOLES:
        # Vector dipole moment
        property dip:
            """Dipole moment as vector"""

            def __set__(self, _q):
                cdef double myq[3]
                checkTypeOrExcept(
                    _q, 3, float, "Dipole moment vector has to be 3 floats")
                for i in range(3):
                    myq[i] = _q[i]
                if set_particle_dip(self.id, myq) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.updateParticleData()
                cdef double * x = NULL
                pointer_to_dip(& (self.particleData), x)
                return np.array([x[0], x[1], x[2]])

        # Scalar magnitude of dipole moment
        property dipm:
            """Dipole moment (magnitude)"""

            def __set__(self, _q):
                checkTypeOrExcept(
                    _q, 1, float, "Magnitude of dipole moment has to be 1 floats")
                if set_particle_dipm(self.id, _q) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.updateParticleData()
                cdef double * x = NULL
                pointer_to_dipm(& (self.particleData), x)
                return x[0]

    IF EXTERNAL_FORCES:
        property ext_force:
            """External force on a particle defined by a vector"""

            def __set__(self, _ext_f):
                cdef double ext_f[3]
                cdef int ext_flag
                checkTypeOrExcept(
                    _ext_f, 3, float, "External force vector has to be 3 floats")
                for i in range(3):
                    ext_f[i] = _ext_f[i]
                if (ext_f[0] == 0 and ext_f[1] == 0 and ext_f[2] == 0):
                    ext_flag = 0
                else:
                    ext_flag = PARTICLE_EXT_FORCE
                if set_particle_ext_force(self.id, ext_flag, ext_f) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.updateParticleData()
                cdef double * ext_f = NULL
                cdef int * ext_flag = NULL
                pointer_to_ext_force( & (self.particleData), ext_flag, ext_f)
                if (ext_flag[0] & PARTICLE_EXT_FORCE):
                    return np.array([ext_f[0], ext_f[1], ext_f[2]])
                else:
                    return np.array([0.0, 0.0, 0.0])

        property fix:
            """Fix the particle at current position"""

            def __set__(self, _fixed_coord_flag):
                cdef int ext_flag
                checkTypeOrExcept(
                    _fixed_coord_flag, 3, int, "Fix has to be 3 ints")
                for i in map(long, range(3)):
                    if (_fixed_coord_flag[i]):
                        ext_flag |= COORD_FIXED(i)
                if set_particle_fix(self.id, ext_flag) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.updateParticleData()
                fixed_coord_flag = np.array([0, 0, 0], dtype=int)
                cdef int * ext_flag = NULL
                pointer_to_fix(& (self.particleData), ext_flag)
                for i in map(long, range(3)):
                    if (ext_flag[0] & COORD_FIXED(i)):
                        fixed_coord_flag[i] = 1
                return fixed_coord_flag

        IF ROTATION:
            property ext_torque:
                """External torque on a particle defined by a vector"""

                def __set__(self, _ext_t):
                    cdef double ext_t[3]
                    cdef int ext_flag
                    checkTypeOrExcept(
                        _ext_t, 3, float, "External force vector has to be 3 floats")
                    for i in range(3):
                        ext_t[i] = _ext_t[i]
                    if (ext_t[0] == 0 and ext_t[1] == 0 and ext_t[2] == 0):
                        ext_flag = 0
                    else:
                        ext_flag = PARTICLE_EXT_TORQUE
                    if set_particle_ext_torque(self.id, ext_flag, ext_t) == 1:
                        raise Exception("set particle position first")

                def __get__(self):
                    self.updateParticleData()
                    cdef double * ext_t = NULL
                    cdef int * ext_flag = NULL
                    pointer_to_ext_torque( & (self.particleData), ext_flag, ext_t)
                    if (ext_flag[0] & PARTICLE_EXT_TORQUE):
                        return np.array([ext_t[0], ext_t[1], ext_t[2]])
                    else:
                        return np.array([0.0, 0.0, 0.0])

    IF LANGEVIN_PER_PARTICLE:
        property gamma:
            """Friction coefficient per particle in Langevin"""

            def __set__(self, _gamma):
                checkTypeOrExcept(_gamma, 1, float, "gamma has to be a float")
                if set_particle_gamma(self.id, _gamma) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.updateParticleData()
                cdef double * gamma = NULL
                pointer_to_gamma(& (self.particleData), gamma)
                return gamma[0]

        property temp:
            """Temperature per particle in Langevin"""

            def __set__(self, _temp):
                checkTypeOrExcept(_temp, 1, float, "temp has to be a float")
                if set_particle_temperature(self.id, _temp) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.updateParticleData()
                cdef double * temp = NULL
                pointer_to_temperature(& (self.particleData), temp)
                return temp[0]

    IF ROTATION_PER_PARTICLE:
        property rotation:
            """Friction coefficient per particle in Langevin"""

            def __set__(self, _rot):
                cdef int rot
                if _rot:
                    rot = 1
                else:
                    rot = 0
                if set_particle_rotation(self.id, rot) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.updateParticleData()
                cdef short int * _rot = NULL
                pointer_to_rotation(& (self.particleData), _rot)
                if _rot[0] == 1:
                    rot = True
                else:
                    rot = False
                return rot

    IF EXCLUSIONS:
        property exclude:
            """Exclude particle from interaction"""

            def __set__(self, _partners):
                delete = 0
                if type(_partners[0]) == str:
                    if _partners.pop(0) == "delete":
                        delete = 1
                for partner in _partners:
                    checkTypeOrExcept(
                        partner, 1, int, "PID of partner has to be an int")
                    if change_exclusion(self.id, partner, delete) == 1:
                        raise Exception("set particle position first")

            def __get__(self):
                self.updateParticleData()
                cdef int * num_partners = NULL
                cdef int * partners = NULL
                py_partners = []
                pointer_to_exclusions( & (self.particleData), num_partners, partners)
                for i in range(num_partners[0]):
                    py_partners.append(partners[i])
                return np.array(py_partners)

    IF ENGINE:
        property swimming:
            """Set swimming parameters"""

            def __set__(self, _params):
                cdef ParticleParametersSwimming swim
                swim.swimming = True
                swim.v_swim = 0.0
                swim.f_swim = 0.0
                IF LB or LB_GPU:
                    swim.push_pull = 0
                    swim.dipole_length = 0.0
                    swim.rotational_friction = 0.0

                if type(_params) == type(True):
                    if _params == True:
                        raise Exception(
                            "To enable swimming supply a dictionary of parameters")
                else:
                    if 'f_swim' in _params and 'v_swim' in _params:
                        raise Exception(
                            "You can't set v_swim and f_swim at the same time")
                    if 'f_swim' in _params:
                        checkTypeOrExcept(
                            _params['f_swim'], 1, float, "f_swim has to be a float")
                        swim.f_swim = _params['f_swim']
                    if 'v_swim' in _params:
                        checkTypeOrExcept(
                            _params['v_swim'], 1, float, "v_swim has to be a float")
                        swim.v_swim = _params['v_swim']

                    IF LB or LB_GPU:
                        if 'mode' in _params:
                            if _params['mode'] == "pusher":
                                swim.push_pull = -1
                            elif _params['mode'] == "puller":
                                swim.push_pull = 1
                            elif _params['mode'] == "N/A":
                                swim.push_pull = 0
                            else:
                                raise Exception(
                                    "'mode' has to be either 'pusher' or 'puller'")

                        if 'dipole_length' in _params:
                            checkTypeOrExcept(
                                _params['dipole_length'], 1, float, "dipole_length has to be a float")
                            swim.dipole_length = _params['dipole_length']

                        if 'rotational_friction' in _params:
                            checkTypeOrExcept(
                                _params['rotational_friction'], 1, float, "rotational_friction has to be a float")
                            swim.rotational_friction = _params[
                                'rotational_friction']

                if set_particle_swimming(self.id, swim) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.updateParticleData()
                swim = {}
                mode = "N/A"
                cdef ParticleParametersSwimming * _swim = NULL
                pointer_to_swimming( & (self.particleData), _swim)
                IF LB or LB_GPU:
                    if _swim.push_pull == -1:
                        mode = 'pusher'
                    elif _swim.push_pull == 1:
                        mode = 'puller'
                    swim = {
                        'v_swim': _swim.v_swim,
                        'f_swim': _swim.f_swim,
                        'mode': mode,
                        'dipole_length': _swim.dipole_length,
                        'rotational_friction': _swim.rotational_friction
                    }
                ELSE:
                    swim = {
                        'v_swim': _swim.v_swim,
                        'f_swim': _swim.f_swim,
                    }
                return swim

    def delete(self):
        """Delete the particle"""
        if remove_particle(self.id):
            raise Exception("Could not delete particle")
        del self

    # Bond related methods
    def addVerifiedBond(self, bond, partner):
        """Add a bond, the validity of which has already been verified"""
        # If someone adds bond types with more than four partners, this has to
        # be changed
        cdef int bondInfo[5]
        bondInfo[0] = bond._bondId
#    for i in range(len(bond)):
#       bondInfo[i]=bond[i]
        if change_particle_bond(self.id, bondInfo, 0):
            raise Exception("Adding the bond failed.")

    def deleteVerifiedBond(self, bond, partner):
        cdef int bondInfo[5]
        bondInfo[0] = bond._bondId
#    for i in range(len(bond)):
#      bondInfo[i]=bond[i]
        if change_particle_bond(self.id, bondInfo, 1):
            raise Exception("Deleting the bond failed.")

    def checkBondOrThrowException(self, bond, partner):
        """Checks the validity of the given bond:
        * if the bond is given as an object
        * if all partners are of type int
        * if the number of partners satisfies the bond
        * If the bond type used exists (is lower than n_bonded_ia)
        * If the number of bond partners fits the bond type
        Throw an exception if any of these are not met"""

        if not isinstance(bond, BondedInteraction):
            raise Exception(
                "Bond argument has to be of type BondedInteraction.")
        if bond._bondId >= n_bonded_ia:
            raise ValueError("The bond type", bond._bondId, "does not exist.")
        if not hasattr(partner, "__getitem"):
            partner = (partner,)
        if bonded_ia_params[bond._bondId].num != len(partner):
            raise ValueError("Bond of type", bond._bondId, "needs", bonded_ia_params[
                             bond._bondId], "partners.")

        for y in partner:
            if not isinstance(y, int):
                raise ValueError("Partners have to be integer.")

    def addBond(self, bond, partner):
        """Add a single bond to the particle"""
        self.checkBondOrThrowException(bond, partner)
        self.addVerifiedBond(bond, partner)

    def deleteBond(self, bond, partner):
        """Delete a single bond from the particle"""
        self.checkBondOrThrowException(bond, partner)
        self.deleteVerifiedBond(bond, partner)

    def deleteAllBonds(self):
        if change_particle_bond(self.id, NULL, 1):
            raise Exception("Deleting all bonds failed.")

    def deleteAllBonds(self):
        if change_particle_bond(self.id, NULL, 1):
            raise Exception("Deleting all bonds failed.")

cdef class particleList:
    """Provides access to the particles via [i], where i is the particle id. Returns a ParticleHandle object """

    def __getitem__(self, key):
        return ParticleHandle(key)
