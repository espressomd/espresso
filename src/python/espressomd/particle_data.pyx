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
include "myconfig.pxi"

cimport numpy as np
import numpy as np
cimport utils
from utils cimport *
cimport particle_data
from interactions import BondedInteraction
from interactions import BondedInteractions
from copy import copy
from globals cimport max_seen_particle, time_step, smaller_time_step

PARTICLE_EXT_FORCE = 1

def COORD_FIXED(coord):
    return 2L << coord
COORDS_FIX_MASK = COORD_FIXED(0) | COORD_FIXED(1) | COORD_FIXED(2)
COORDS_ALL_FIXED = COORD_FIXED(0) & COORD_FIXED(1) & COORD_FIXED(2)
PARTICLE_EXT_TORQUE = 16


particle_attributes = ["type","pos", "v", "f", "bonds"]

IF MULTI_TIMESTEP:
    particle_attributes.append("smaller_timestep")
                    
IF MASS:
    particle_attributes.append("mass")

IF ROTATION:
    particle_attributes.append("omega_lab")
    particle_attributes.append("ext_torque")

IF ROTATIONAL_INERTIA:
    particle_attributes.append("rinertia")
    particle_attributes.append("omega_body")
    particle_attributes.append("torque_lab")
    particle_attributes.append("quat")

IF ELECTROSTATICS:
    particle_attributes.append("q")

IF VIRTUAL_SITES:
    particle_attributes.append("virtual")
    
IF VIRTUAL_SITES_RELATIVE:
    particle_attributes.append("vs_relative")

IF DIPOLES:
    particle_attributes.append("dip")
    particle_attributes.append("dipm")

IF EXTERNAL_FORCES:
    particle_attributes.append("ext_force")
    particle_attributes.append("fix")

IF LANGEVIN_PER_PARTICLE:
    particle_attributes.append("gamma")
    particle_attributes.append("temp")

IF ROTATION_PER_PARTICLE:
    particle_attributes.append("rotation")

IF EXCLUSIONS:
    particle_attributes.append("exclude")

IF ENGINE:
    particle_attributes.append("swimming")


cdef class ParticleHandle:
    def __cinit__(self, _id):
        #    utils.init_intlist(self.particle_data.el)
        utils.init_intlist( & (self.particle_data.bl))
        self.id = _id

    cdef int update_particle_data(self) except -1:
        #    utils.realloc_intlist(self.particle_data.el, 0)
        utils.realloc_intlist( & (self.particle_data.bl), 0)

        if get_particle_data(self.id, & self.particle_data):
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
            self.update_particle_data()
            return self.particle_data.p.type

    # Position
    property pos:
        """Particle position (not folded into central image)."""

        def __set__(self, _pos):
            cdef double mypos[3]
            check_type_or_throw_except(
                _pos, 3, float, "Postion must be 3 floats")
            for i in range(3):
                mypos[i] = _pos[i]
            if place_particle(self.id, mypos) == -1:
                raise Exception("particle could not be set")

        def __get__(self):
            self.update_particle_data()
            return np.array([self.particle_data.r.p[0],
                             self.particle_data.r.p[1],
                             self.particle_data.r.p[2]])

    # Velocity
    property v:
        """Particle velocity"""

        def __set__(self, _v):
            global time_step
            cdef double myv[3]
            check_type_or_throw_except(
                _v, 3, float, "Velocity has to be floats")
            for i in range(3):
                myv[i] = _v[i]
                myv[i] *= time_step
            if set_particle_v(self.id, myv) == 1:
                raise Exception("set particle position first")

        def __get__(self):
            global time_step, smaller_time_step
            self.update_particle_data()
            IF MULTI_TIMESTEP:
                if smaller_time_step > 0. and self.smaller_timestep:
                    return np.array([self.particle_data.m.v[0]/smaller_time_step,
                                     self.particle_data.m.v[1]/smaller_time_step,
                                     self.particle_data.m.v[2]/smaller_time_step])
                else:
                    return np.array([self.particle_data.m.v[0]/time_step,
                                     self.particle_data.m.v[1]/time_step,
                                     self.particle_data.m.v[2]/time_step])
            ELSE:
                return np.array([self.particle_data.m.v[0]/time_step,
                                 self.particle_data.m.v[1]/time_step,
                                 self.particle_data.m.v[2]/time_step])

    # Force
    property f:
        """Particle force"""

        def __set__(self, _f):
            global time_step
            cdef double myf[3]
            check_type_or_throw_except(_f, 3, float, "Force has to be floats")
            for i in range(3):
                myf[i] = _f[i]
                myf[i] *= (0.5*time_step**2)
            if set_particle_f(self.id, myf) == 1:
                raise Exception("set particle position first")

        def __get__(self):
            global time_step
            self.update_particle_data()
            return np.array([self.particle_data.f.f[0]/(0.5*time_step**2),
                             self.particle_data.f.f[1]/(0.5*time_step**2),
                             self.particle_data.f.f[2]/(0.5*time_step**2)])

    # Bonds
    property bonds:
        """Bond partners with respect to bonded interactions."""

        def __set__(self, _bonds):
            # First, we check that we got a list/tuple.
            if not hasattr(_bonds, "__getitem__"):
                raise ValueError(
                    "bonds have to specified as a tuple of tuples. (Lists can also be used)")
            bonds = list(_bonds)  # as we modify it

            # Assigning to the bond property means replacing the existing value
            # i.e., we delete all existing bonds
            if change_particle_bond(self.id, NULL, 1):
                raise Exception("Deleting existing bonds failed.")

            # And add the new ones
            for bond in bonds:
                self.add_bond(bond)

        def __get__(self):
            self.update_particle_data()
            bonds = []
            # Go through the bond list of the particle
            i = 0
            while i < self.particle_data.bl.n:
                bond = []
                # Bond type:
                bond_id = self.particle_data.bl.e[i]
                bond.append(BondedInteractions()[bond_id])
                # Number of partners
                nPartners = bonded_ia_params[bond_id].num

                i += 1

                # Copy bond partners
                for j in range(nPartners):
                    bond.append(self.particle_data.bl.e[i])
                    i += 1
                bonds.append(tuple(bond))

            return tuple(bonds)

    # Properties that exist only when certain features are activated
    # MULTI_TIMESTEP
    IF MULTI_TIMESTEP == 1:
        property smaller_timestep:
            """Particle flag specifying whether particle trajectory should be integrated with time_step of small_time_step"""

            def __set__(self, _smaller_timestep):
                check_type_or_throw_except(
                    _smaller_timestep, 1, int, "Smaller time step flag has to be 1 ints")
                if set_particle_smaller_timestep(self.id, _smaller_timestep) == 1:
                    raise Exception("error setting particle smaller_timestep")

            def __get__(self):
                self.update_particle_data()
                cdef int * x = NULL
                pointer_to_smaller_timestep( & (self.particle_data), x)
                return x[0]

    # MASS
    IF MASS == 1:
        property mass:
            """Particle mass"""

            def __set__(self, _mass):
                check_type_or_throw_except(
                    _mass, 1, float, "Mass has to be 1 floats")
                if set_particle_mass(self.id, _mass) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.update_particle_data()
                cdef double * x = NULL
                pointer_to_mass( & (self.particle_data), x)
                return x[0]

    IF ROTATION == 1:
        # Omega (angular velocity) lab frame
        property omega_lab:
            """Angular velocity in lab frame"""

            def __set__(self, _o):
                cdef double myo[3]
                check_type_or_throw_except(
                    _o, 3, float, "Omega_lab has to be 3 floats")
                for i in range(3):
                    myo[i] = _o[i]
                if set_particle_omega_lab(self.id, myo) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.update_particle_data()
                cdef double o[3]
                convert_omega_body_to_space(& (self.particle_data), o)
                return np.array([o[0], o[1], o[2]])

    # ROTATIONAL_INERTIA
    IF ROTATIONAL_INERTIA == 1:
        property rinertia:
            """Rotational inertia"""

            def __set__(self, _rinertia):
                cdef double rinertia[3]
                check_type_or_throw_except(
                    _rinertia, 3, float, "Rotation_inertia has to be 3 floats")
                for i in range(3):
                    rinertia[i] = _rinertia[i]
                if set_particle_rotational_inertia(self.id, rinertia) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.update_particle_data()
                cdef double * rinertia = NULL
                pointer_to_rotational_inertia( & (self.particle_data), rinertia)
                return np.array([rinertia[0], rinertia[1], rinertia[2]])

# Omega (angular velocity) body frame
        property omega_body:
            """Angular velocity in body frame"""

            def __set__(self, _o):
                cdef double myo[3]
                check_type_or_throw_except(
                    _o, 3, float, "Omega_body has to be 3 floats")
                for i in range(3):
                    myo[i] = _o[i]
                if set_particle_omega_body(self.id, myo) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.update_particle_data()
                cdef double * o = NULL
                pointer_to_omega_body( & (self.particle_data), o)
                return np.array([o[0], o[1], o[2]])


# Torque in lab frame
        property torque_lab:
            """Torque in lab frame"""

            def __set__(self, _t):
                cdef double myt[3]
                check_type_or_throw_except(
                    _t, 3, float, "Torque has to be 3 floats")
                for i in range(3):
                    myt[i] = _t[i]
                if set_particle_torque_lab(self.id, myt) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.update_particle_data()
                cdef double x[3]
                convert_torques_body_to_space(& (self.particle_data), x)
                return np.array([x[0], x[1], x[2]])

# Quaternion
        property quat:
            """Quaternions"""

            def __set__(self, _q):
                cdef double myq[4]
                check_type_or_throw_except(
                    _q, 4, float, "Quaternions has to be 4 floats")
                for i in range(4):
                    myq[i] = _q[i]
                if set_particle_quat(self.id, myq) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.update_particle_data()
                cdef double * x = NULL
                pointer_to_quat(& (self.particle_data), x)
                return np.array([x[0], x[1], x[2], x[3]])
# Director ( z-axis in body fixed frame)
        property director:
            """Director"""

            def __set__(self, _q):
                raise Exception(
                    "Setting the director is not implemented in the c++-core of Espresso")
#        cdef double myq[3]
#        check_type_or_throw_except(_q,3,float,"Director has to be 3 floats")
#        for i in range(3):
#            myq[i]=_q[i]
#        if set_particle_quatu(self.id, myq) == 1:
#          raise Exception("set particle position first")

            def __get__(self):
                self.update_particle_data()
                cdef double * x = NULL
                pointer_to_quatu(& (self.particle_data), x)
                return np.array([x[0], x[1], x[2]])

# Charge
    IF ELECTROSTATICS == 1:
        property q:
            """particle charge"""

            def __set__(self, _q):
                cdef double myq
                check_type_or_throw_except(
                    _q, 1, float, "Charge has to be floats")
                myq = _q
                if set_particle_q(self.id, myq) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.update_particle_data()
                cdef double * x = NULL
                pointer_to_q(& (self.particle_data), x)
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
                self.update_particle_data()
                cdef int * x = NULL
                pointer_to_virtual(& (self.particle_data), x)
                return x[0]

    IF VIRTUAL_SITES_RELATIVE == 1:
        property vs_relative:
            """virtual sites relative parameters"""

            def __set__(self, x):
                if len(x) != 3:
                    raise ValueError("vs_relative needs six args")
                _relto = x[0]
                _dist = x[1]
                q = x[2]
                check_type_or_throw_except(
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
                self.update_particle_data()
                cdef int * rel_to = NULL
                cdef double * dist = NULL
                cdef double * q = NULL
                pointer_to_vs_relative( & (self.particle_data), rel_to, dist, q)
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
            check_type_or_throw_except(
                _relto, 1, int, "Argument of vs_auto_relate_to has to be of type int")
            if vs_relate_to(self.id, _relto):
                raise Exception("vs_relative setup failed.")

    IF DIPOLES:
        # Vector dipole moment
        property dip:
            """Dipole moment as vector"""

            def __set__(self, _q):
                cdef double myq[3]
                check_type_or_throw_except(
                    _q, 3, float, "Dipole moment vector has to be 3 floats")
                for i in range(3):
                    myq[i] = _q[i]
                if set_particle_dip(self.id, myq) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.update_particle_data()
                cdef double * x = NULL
                pointer_to_dip(& (self.particle_data), x)
                return np.array([x[0], x[1], x[2]])

        # Scalar magnitude of dipole moment
        property dipm:
            """Dipole moment (magnitude)"""

            def __set__(self, _q):
                check_type_or_throw_except(
                    _q, 1, float, "Magnitude of dipole moment has to be 1 floats")
                if set_particle_dipm(self.id, _q) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.update_particle_data()
                cdef double * x = NULL
                pointer_to_dipm(& (self.particle_data), x)
                return x[0]

    IF EXTERNAL_FORCES:
        property ext_force:
            """External force on a particle defined by a vector"""

            def __set__(self, _ext_f):
                cdef double ext_f[3]
                cdef int ext_flag
                check_type_or_throw_except(
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
                self.update_particle_data()
                cdef double * ext_f = NULL
                cdef int * ext_flag = NULL
                pointer_to_ext_force( & (self.particle_data), ext_flag, ext_f)
                if (ext_flag[0] & PARTICLE_EXT_FORCE):
                    return np.array([ext_f[0], ext_f[1], ext_f[2]])
                else:
                    return np.array([0.0, 0.0, 0.0])

        property fix:
            """Fix the particle at current position"""

            def __set__(self, _fixed_coord_flag):
                cdef int ext_flag
                check_type_or_throw_except(
                    _fixed_coord_flag, 3, int, "Fix has to be 3 ints")
                for i in map(long, range(3)):
                    if (_fixed_coord_flag[i]):
                        ext_flag |= COORD_FIXED(i)
                if set_particle_fix(self.id, ext_flag) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.update_particle_data()
                fixed_coord_flag = np.array([0, 0, 0], dtype=int)
                cdef int * ext_flag = NULL
                pointer_to_fix(& (self.particle_data), ext_flag)
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
                    check_type_or_throw_except(
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
                    self.update_particle_data()
                    cdef double * ext_t = NULL
                    cdef int * ext_flag = NULL
                    pointer_to_ext_torque( & (self.particle_data), ext_flag, ext_t)
                    if (ext_flag[0] & PARTICLE_EXT_TORQUE):
                        return np.array([ext_t[0], ext_t[1], ext_t[2]])
                    else:
                        return np.array([0.0, 0.0, 0.0])

    IF LANGEVIN_PER_PARTICLE:
        property gamma:
            """Friction coefficient per particle in Langevin"""

            def __set__(self, _gamma):
                check_type_or_throw_except(
                    _gamma, 1, float, "gamma has to be a float")
                if set_particle_gamma(self.id, _gamma) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.update_particle_data()
                cdef double * gamma = NULL
                pointer_to_gamma(& (self.particle_data), gamma)
                return gamma[0]

        property temp:
            """Temperature per particle in Langevin"""

            def __set__(self, _temp):
                check_type_or_throw_except(
                    _temp, 1, float, "temp has to be a float")
                if set_particle_temperature(self.id, _temp) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.update_particle_data()
                cdef double * temp = NULL
                pointer_to_temperature(& (self.particle_data), temp)
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
                self.update_particle_data()
                cdef short int * _rot = NULL
                pointer_to_rotation(& (self.particle_data), _rot)
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
                if len(_partners) == 0:
                    return
                if type(_partners[0]) == str:
                    if _partners.pop(0) == "delete":
                        delete = 1
                for partner in _partners:
                    check_type_or_throw_except(
                        partner, 1, int, "PID of partner has to be an int")
                    if change_exclusion(self.id, partner, delete) == 1:
                        raise Exception("set particle position first")

            def __get__(self):
                self.update_particle_data()
                cdef int * num_partners = NULL
                cdef int * partners = NULL
                py_partners = []
                pointer_to_exclusions( & (self.particle_data), num_partners, partners)
                for i in range(num_partners[0]):
                    py_partners.append(partners[i])
                return np.array(py_partners)

    IF ENGINE:
        property swimming:
            """Set swimming parameters"""

            def __set__(self, _params):
                cdef particle_parameters_swimming swim
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
                        if _params["f_swim"] == 0 or _params["v_swim"] == 0:
                            pass
                        else:
                            raise Exception(
                                "You can't set v_swim and f_swim at the same time")
                    if 'f_swim' in _params:
                        check_type_or_throw_except(
                            _params['f_swim'], 1, float, "f_swim has to be a float")
                        swim.f_swim = _params['f_swim']
                    if 'v_swim' in _params:
                        check_type_or_throw_except(
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
                            check_type_or_throw_except(
                                _params['dipole_length'], 1, float, "dipole_length has to be a float")
                            swim.dipole_length = _params['dipole_length']

                        if 'rotational_friction' in _params:
                            check_type_or_throw_except(
                                _params['rotational_friction'], 1, float, "rotational_friction has to be a float")
                            swim.rotational_friction = _params[
                                'rotational_friction']

                if set_particle_swimming(self.id, swim) == 1:
                    raise Exception("set particle position first")

            def __get__(self):
                self.update_particle_data()
                swim = {}
                mode = "N/A"
                cdef particle_parameters_swimming * _swim = NULL
                pointer_to_swimming( & (self.particle_data), _swim)
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

    def remove(self):
        """Delete the particle"""
        if remove_particle(self.id):
            raise Exception("Could not delete particle")
        del ParticleList.key_dict["%i"%self.id]
        del self

    # Bond related methods
    # Does not work properly with 3 or more partner bonds!
    def add_verified_bond(self, bond):
        """Add a bond, the validity of which has already been verified"""
        # If someone adds bond types with more than four partners, this has to
        # be changed
        cdef int bond_info[5]
        bond_info[0] = bond[0]._bond_id
        for i in range(1, len(bond)):
            bond_info[i] = bond[i]
        if change_particle_bond(self.id, bond_info, 0):
            raise Exception("Adding the bond failed.")

    def delete_verified_bond(self, bond):
        cdef int bond_info[5]
        bond_info[0] = bond[0]._bond_id
        for i in range(1, len(bond)):
            bond_info[i] = bond[i]
        if change_particle_bond(self.id, bond_info, 1):
            raise Exception("Deleting the bond failed.")

    def check_bond_or_throw_exception(self, bond):
        """Checks the validity of the given bond:
        * if the bondtype is given as an object or a numerical id
        * if all partners are of type int
        * if the number of partners satisfies the bond
        * If the bond type used exists (is lower than n_bonded_ia)
        * If the number of bond partners fits the bond type
        Throw an exception if any of these are not met"""

        # Has it []-access
        if not hasattr(bond, "__getitem__"):
            raise ValueError(
                "Bond needs to be a tuple or list containing bond type and partners")

        # Bond type or numerical bond id
        if not isinstance(bond[0], BondedInteraction):
            if isinstance(bond[0], int):
                bond[0] = BondedInteractions()[bond[0]]
            else:
                raise Exception(
                    "1st element of Bond has to be of type BondedInteraction or int.")

        # Validity of the numeric id
        if bond[0]._bond_id >= n_bonded_ia:
            raise ValueError("The bond type", bond._bond_id, "does not exist.")

        # Number of partners
        if bonded_ia_params[bond[0]._bond_id].num != len(bond) - 1:
            raise ValueError("Bond of type", bond._bond_id, "needs", bonded_ia_params[
                             bond[0]._bond_id], "partners.")

        # Type check on partners
        for y in bond[1:]:
            if not isinstance(y, int):
                raise ValueError("Partners have to be integer.")

    def add_bond(self, _bond):
        """Add a single bond to the particle"""
        bond = list(_bond)  # As we will modify it
        self.check_bond_or_throw_exception(bond)
        self.add_verified_bond(bond)

    def delete_bond(self, _bond):
        """Delete a single bond from the particle"""
        bond = list(_bond)  # as we modify it
        self.check_bond_or_throw_exception(bond)
        self.delete_verified_bond(bond)

    def delete_all_bonds(self):
        if change_particle_bond(self.id, NULL, 1):
            raise Exception("Deleting all bonds failed.")

    def update(self, P):

        if "id" in P:
            raise Exception("Cannot change particle id.")

        for k in P.keys():
            setattr(self, k, P[k])



cdef class ParticleSlice:
    """Handles slice inputs e.g. part[0:2]. Sets values for selected slices or returns values as a single list."""


    def __cinit__(self,slice_):
        id_list=np.arange(max_seen_particle+1)
        self.id_selection=id_list[slice_]

    cdef int update_particle_data(self, id) except -1:
        utils.realloc_intlist( & (self.particle_data.bl), 0)

        if get_particle_data(id, & self.particle_data):
            raise Exception("Error updating particle data")
        else:
            return 0


    # Particle Type
    property type:

        def __get__(self):
            type_list = []
            for id in self.id_selection:
                type_list.append(ParticleHandle(id).type)
            return type_list

        def __set__(self, _type_list):
            if isinstance(_type_list,int):
                for id in self.id_selection:
                    ParticleHandle(id).type = _type_list
                return
            if len(self.id_selection) != len(_type_list):
                raise Exception("Input list size (%i) does not match slice size (%i)"%(len(_type_list),len(self.id_selection)))
            for i in range(len(self.id_selection)):
                ParticleHandle(self.id_selection[i]).type = _type_list[i]
                    

    # Position
    property pos:
        """Particle position (not folded into central image)."""

        def __set__(self, _pos_array):
            if len(self.id_selection) != len(_pos_array):
                raise Exception("Input list size (%i) does not match slice size (%i)"%(len(_pos_array),len(self.id_selection)))

            cdef double mypos[3]
            for i in range(len(_pos_array)):
                ParticleHandle(self.id_selection[i]).pos = _pos_array[i]

        def __get__(self):
            pos_array = np.zeros((len(self.id_selection),3))
            for i in range(len(self.id_selection)):
                pos_array[i,:] = ParticleHandle(self.id_selection[i]).pos
            return pos_array


    # Velocity
    property v:
        """Particle velocity"""

        def __set__(self, _v_array):
            if len(np.array(_v_array).shape) == 1:
                for id in self.id_selection:
                    ParticleHandle(id).v = _v_array
                return

            if len(self.id_selection) != len(_v_array):
                raise Exception("Input list size (%i) does not match slice size (%i)"%(len(_v_array),len(self.id_selection)))
                
            for i in range(len(self.id_selection)):
                ParticleHandle(self.id_selection[i]).v = _v_array[i]

        def __get__(self):
            v_array = np.zeros((len(self.id_selection),3))
            for i in range(len(self.id_selection)):
                v_array[i,:] = ParticleHandle(self.id_selection[i]).v
            return v_array


    # Force
    property f:
        """Particle force"""

        def __set__(self, _f_array):
            if len(np.array(_f_array).shape) == 1:
                for id in self.id_selection:
                    ParticleHandle(id).f = _f_array
                return

            if len(self.id_selection) != len(_f_array):
                raise Exception("Input list size (%i) does not match slice size (%i)"%(len(_f_array),len(self.id_selection)))
            for i in range(len(_f_array)):
                ParticleHandle(self.id_selection[i]).f = _f_array[i]

        def __get__(self):
            f_array = np.zeros((len(self.id_selection),3))
            for i in range(len(self.id_selection)):
                f_array[i,:] = ParticleHandle(self.id_selection[i]).f
            return f_array


    IF MASS:
        property mass:
            """Particle mass"""

            def __set__(self, _mass_array):
                if isinstance(_mass_array, int) or isinstance(_mass_array, float):
                    for i in range(len(self.id_selection)):
                        ParticleHandle(self.id_selection[i]).mass = _mass_array
                    return
                if len(self.id_selection) != len(_mass_array):
                    raise Exception("Input list size (%i) does not match slice size (%i)"%(len(_mass_array),len(self.id_selection)))
                for i in range(len(_mass_array)):
                    ParticleHandle(self.id_selection[i]).mass = _mass_array[i]

            def __get__(self):
                mass_array = np.zeros_like(self.id_selection)
                for i in range(len(self.id_selection)):
                    mass_array[i] = ParticleHandle(self.id_selection[i]).mass
                return mass_array


    IF ELECTROSTATICS == 1:
        property q:
            """particle charge"""

            def __set__(self, _q_array):
                if isinstance(_q_array, int) or isinstance(_q_array, float):
                    for i in range(len(self.id_selection)):
                        ParticleHandle(self.id_selection[i]).q = _q_array
                    return

                if len(self.id_selection) != len(_q_array):
                    raise Exception("Input list size (%i) does not match slice size (%i)"%(len(_q_array),len(self.id_selection)))
                for i in range(len(self.id_selection)):
                    ParticleHandle(self.id_selection[i]).q = _q_array[i]

            def __get__(self):
                q_array = np.zeros_like(self.id_selection)
                for i in range(len(self.id_selection)):
                    q_array[i] = ParticleHandle(self.id_selection[i]).q
                return q_array


    IF EXTERNAL_FORCES:
        property ext_force:
            """External force on a particle defined by a vector"""

            def __set__(self, _ext_f_array):
                if len(np.array(_ext_f_array).shape) == 1:
                    for i in range(len(self.id_selection)):
                        ParticleHandle(self.id_selection[i]).ext_force = _ext_f_array
                    return

                if len(self.id_selection) != len(_ext_f_array):
                    raise Exception("Input list size (%i) does not match slice size (%i)"%(len(_ext_f_array),len(self.id_selection)))

                for i in range(len(self.id_selection)):
                    ParticleHandle(self.id_selection[i]).ext_force = _ext_f_array[i]

            def __get__(self):
                ext_f_array = np.zeros((len(self.id_selection),3))
                for i in range(len(self.id_selection)):
                    ext_f_array[i,:] = ParticleHandle(self.id_selection[i]).ext_force
                    
                return ext_f_array


    def update(self, P):
        if "id" in P:
            raise Exception("Cannot change particle id.")
            
        for k in P.keys():
            setattr(self, k, P[k])

    # Bond related methods
    def add_bond(self, _bond):
        """Add a single bond to the particles"""
        bond = list(_bond)  # As we will modify it
        for i in range(len(self.id_selection)):
            partners = []
            for j in range(1,len(bond)):
                partners.append(bond[j][i])
            ParticleHandle(self.id_selection[i]).add_bond((bond[0],*partners))

    def delete_bond(self, _bond):
        """Delete a single bond from the particles"""
        bond = list(_bond)  # as we modify it
        for i in range(len(self.id_selection)):
            partners = []
            for j in range(1,len(bond)):
                partners.append(bond[j][i])
            ParticleHandle(self.id_selection[i]).delete_bond((bond[0],*partners))

    def delete_all_bonds(self):
        for i in range(len(self.id_selection)):
            ParticleHandle(self.id_selection[i]).delete_all_bonds()

    def remove(self):
        """Delete the particles"""
        for id in self.id_selection:
            ParticleHandle(id).remove()


cdef class ParticleList:
    """Provides access to the particles via [i], where i is the particle id. Returns a ParticleHandle object """
    key_dict={}
    # Retrieve a particle
    def __getitem__(self, key):
        if isinstance(key, slice):
            return ParticleSlice(key)
            
        if not np.all(self.exists(key)):
            if isinstance(key, int):
                non_existing = key
            else:
                non_existing =np.trim_zeros((np.array(key)*np.invert(self.exists(key)))) 
            raise Exception("Particle(s) %s does not exist." % non_existing)

        if isinstance(key, tuple) or isinstance(key, list) or isinstance(key, np.ndarray):
            return ParticleSlice(np.array(key))

        return ParticleHandle(key)


    # __getstate__ and __setstate__ define the pickle interaction
    def __getstate__(self):
        odict={}
        key_list = sorted(ParticleList.key_dict.values())
        for particle_number in key_list:
            pdict={}
            for property_ in particle_attributes:
                pdict[property_] = ParticleHandle(particle_number).__getattribute__(property_)
            odict[particle_number] = pdict
        return odict

    def __setstate__(self,params):
        for particle_number in params.keys():
            params[particle_number]["id"] = particle_number
            self.add(params[particle_number])


    def add(self, *args, **kwargs):

        # Did we get a dictionary
        if len(args) == 1:
            if hasattr(args[0], "__getitem__"):
                # Check for presence of pos attribute
                if not "pos" in args[0]:
                    raise ValueError(
                        "pos attribute must be specified for new particle")

                if len(np.array(args[0]["pos"]).shape) == 2:
                    self._place_new_particles(args[0])
                else:
                    self._place_new_particle(args[0])
        else:
            if len(args) == 0 and len(kwargs.keys()) != 0:
                # Check for presence of pos attribute
                if not "pos" in kwargs:
                    raise ValueError(
                        "pos attribute must be specified for new particle")

                if len(np.array(kwargs["pos"]).shape) == 2:
                    self._place_new_particles(kwargs)
                else:
                    self._place_new_particle(kwargs)
            else:
                raise ValueError(
                    "add() takes either a dictionary or a bunch of keyword args")


    def _place_new_particle(self, P):

        # Handling of particle id
        if not "id" in P:
            # Generate particle id
            P["id"] = max_seen_particle + 1
        else:
            if particle_exists(P["id"]):
                raise Exception("Particle %d already exists." % P["id"])

        # The ParticleList[]-getter ist not valid yet, as the particle
        # doesn't yet exist. Hence, the setting of position has to be
        # done here. the code is from the pos:property of ParticleHandle
        cdef double mypos[3]
        check_type_or_throw_except(
            P["pos"], 3, float, "Postion must be 3 floats")
        for i in range(3):
            mypos[i] = P["pos"][i]
        if place_particle(P["id"], mypos) == -1:
            raise Exception("particle could not be set")
        # Pos is taken care of
        del P["pos"]
        id = P["id"]
        del P["id"]

        if P != {}:
            self[id].update(P)
        ParticleList.key_dict["%i"%id] = id


    def _place_new_particles(self, P):

        if not "id" in P:
            # Generate particle ids
            ids = np.arange(np.array(P["pos"]).shape[0]) + max_seen_particle + 1
        else:
            ids = P["id"]
            del P["id"]

        # Place particles
        cdef double mypos[3]
        for j in range(len(P["pos"])):
            for i in range(3):
                mypos[i] = P["pos"][j][i]
            if place_particle(ids[j], mypos) == -1:
                raise Exception("particle could not be set")
            ParticleList.key_dict["%i"%ids[j]] = ids[j]


        del P["pos"]
        
        if P!= {}:
            self[ids].update(P)


    # Iteration over all existing particles
    def __iter__(self):
        for i in range(max_seen_particle + 1):
            if particle_exists(i):
                yield self[i]

    def exists(self, idx):
        if isinstance(idx,int):
            return particle_exists(idx)
        if isinstance(idx,slice) or isinstance(idx,tuple) or isinstance(idx,list) or isinstance(idx,np.ndarray):
            tf_array=np.zeros(len(idx), dtype=np.bool)
            for i in range(len(idx)):
                tf_array[i]=particle_exists(idx[i])
            return tf_array
    


