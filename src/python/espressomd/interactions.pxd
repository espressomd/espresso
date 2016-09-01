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
# Handling of interactions

from __future__ import print_function, absolute_import
include "myconfig.pxi"
from espressomd._system cimport *
cimport numpy as np
from espressomd.utils cimport *

cdef extern from "interaction_data.hpp":
    ctypedef struct ia_parameters "IA_parameters":
        double LJ_eps
        double LJ_sig
        double LJ_cut
        double LJ_shift
        double LJ_offset
        double LJ_capradius
        double LJ_min

        double LJGEN_eps
        double LJGEN_sig
        double LJGEN_cut
        double LJGEN_shift
        double LJGEN_offset
        double LJGEN_capradius
        int LJGEN_a1
        int LJGEN_a2
        double LJGEN_b1
        double LJGEN_b2
        double LJGEN_lambda
        double LJGEN_softrad
        int TAB_npoints;
        int TAB_startindex;
        double TAB_minval;
        double TAB_minval2;
        double TAB_maxval;
        double TAB_stepsize;
        char TAB_filename[256]; 

    cdef ia_parameters * get_ia_param(int i, int j)
    cdef ia_parameters * get_ia_param_safe(int i, int j)
    cdef void make_bond_type_exist(int type)

cdef extern from "lj.hpp":
    cdef int lennard_jones_set_params(int part_type_a, int part_type_b,
                                      double eps, double sig, double cut,
                                      double shift, double offset,
                                      double cap_radius, double min)

cdef extern from "forcecap.hpp":
    double force_cap
    int forcecap_set_params(double forcecap)

cdef extern from "ljgen.hpp":
    IF LJGEN_SOFTCORE:
        cdef int ljgen_set_params(int part_type_a, int part_type_b,
                                  double eps, double sig, double cut,
                                  double shift, double offset,
                                  int a1, int a2, double b1, double b2,
                                  double cap_radius,
                                  double genlj_lambda, double softrad)
    ELSE:
        cdef int ljgen_set_params(int part_type_a, int part_type_b,
                                  double eps, double sig, double cut,
                                  double shift, double offset,
                                  int a1, int a2, double b1, double b2,
                                  double cap_radius)


IF TABULATED==1:
    cdef extern from "tab.hpp":
        int tabulated_set_params(int part_type_a, int part_type_b, char* filename);


cdef extern from "interaction_data.hpp":
    ctypedef struct Fene_bond_parameters:
        double k
        double drmax
        double r0
        double drmax2
        double drmax2i


#* Parameters for oif_global_forces */
    ctypedef struct Oif_global_forces_bond_parameters:
        double A0_g
        double ka_g
        double V0
        double kv

#* Parameters for oif_local_forces */
    ctypedef struct Oif_local_forces_bond_parameters:
        double r0
        double ks
        double kslin
        double phi0
        double kb
        double A01
        double A02
        double kal

#* Parameters for harmonic bond Potential */
    ctypedef struct Harmonic_bond_parameters:
        double k
        double r
        double r_cut

#* Parameters for the harmonic dumbbell bond potential */
    ctypedef struct Harmonic_dumbbell_bond_parameters:
        double k1
        double k2
        double r
        double r_cut

#* Parameters for three body angular potential (bond-angle potentials).
    ctypedef struct Angle_bond_parameters:
        double bend
        double phi0
        double cos_phi0
        double sin_phi0

#* Parameters for three body angular potential (bond_angle_harmonic).
    ctypedef struct Angle_harmonic_bond_parameters:
        double bend
        double phi0


#* Parameters for three body angular potential (bond_angle_cosine).
    ctypedef struct Angle_cosine_bond_parameters:
        double bend
        double phi0
        double cos_phi0
        double sin_phi0


#* Parameters for three body angular potential (bond_angle_cossquare).
    ctypedef struct Angle_cossquare_bond_parameters:
        double bend
        double phi0
        double cos_phi0

#* Parameters for four body angular potential (dihedral-angle potentials). */
    ctypedef struct Dihedral_bond_parameters:
        double mult
        double bend
        double phase


#* Parameters for n-body tabulated potential (n=2,3,4). */
    ctypedef struct Tabulated_bond_parameters:
        char * filename
        int    type
        int    npoints
        double minval
        double maxval
        double invstepsize
        double * f
        double * e

#* Parameters for n-body overlapped potential (n=2,3,4). */
    ctypedef struct Overlap_bond_parameters:
        char * filename
        int    type
        double maxval
        int    noverlaps
        double * para_a
        double * para_b
        double * para_c

#* Parameters for one-directional harmonic potential */
    ctypedef struct Umbrella_bond_parameters:
        double k
        int    dir
        double r

#* Dummy parameters for -LJ Potential */
    ctypedef struct Subt_lj_bond_parameters:
        double k
        double r
        double r2

#*Parameters for the rigid_bond/SHAKE/RATTLE ALGORITHM*/
    ctypedef struct Rigid_bond_parameters:
        #*Length of rigid bond/Constrained Bond*/
        # double d
        #*Square of the length of Constrained Bond*/
        double d2
        #*Positional Tolerance/Accuracy value for termination of RATTLE/SHAKE iterations during position corrections*/
        double p_tol
        #*Velocity Tolerance/Accuracy for termination of RATTLE/SHAKE iterations during velocity corrections */
        double v_tol

#* Parameters for three body angular potential (bond-angle potentials) that
    ctypedef struct Angledist_bond_parameters:
        double bend
        double phimin
        double distmin
        double phimax
        double distmax
        double cos_phi0
        double sin_phi0


#* Parameters for chainend angular potential with wall  */
    ctypedef struct Endangledist_bond_parameters:
        double bend
        double phi0
        double distmin
        double distmax

#* Union in which to store the parameters of an individual bonded interaction */
    ctypedef union bond_parameters "Bond_parameters":
        Fene_bond_parameters fene
        Oif_global_forces_bond_parameters oif_global_forces
        Oif_local_forces_bond_parameters oif_local_forces
        Harmonic_bond_parameters harmonic
        Harmonic_dumbbell_bond_parameters harmonic_dumbbell
        Angle_bond_parameters angle
        Angle_harmonic_bond_parameters angle_harmonic
        Angle_cosine_bond_parameters angle_cosine
        Angle_cossquare_bond_parameters angle_cossquare
        Dihedral_bond_parameters dihedral
        Tabulated_bond_parameters tab
        Overlap_bond_parameters overlap
        Subt_lj_bond_parameters subt_lj
        Rigid_bond_parameters rigid_bond
        Angledist_bond_parameters angledist
        Endangledist_bond_parameters endangledist

    ctypedef struct bonded_ia_parameters:
        int type
        int num
        #* union to store the different bonded interaction parameters. */
        bond_parameters p

    bonded_ia_parameters * bonded_ia_params
    cdef int n_bonded_ia

cdef extern from "fene.hpp":
    int fene_set_params(int bond_type, double k, double drmax, double r0)
cdef extern from "harmonic.hpp":
    int harmonic_set_params(int bond_type, double k, double r, double r_cut)
cdef extern from "dihedral.hpp":
    int dihedral_set_params(int bond_type, int mult, double bend, double phase)
cdef extern from "angle_harmonic.hpp":
    int angle_harmonic_set_params(int bond_type, double bend, double phi0)
cdef extern from "rattle.hpp":
    int rigid_bond_set_params(int bond_type, double d, double p_tol, double v_tol)
cdef extern from "angle_cosine.hpp":
    int angle_cosine_set_params(int bond_type, double bend, double phi0)
cdef extern from "angle_cossquare.hpp":
    int angle_cossquare_set_params(int bond_type, double bend, double phi0)
cdef extern from "subt_lj.hpp":
    int subt_lj_set_params(int bond_type, double k, double r)
cdef extern from "object-in-fluid/oif_global_forces.hpp":
    int oif_global_forces_set_params(int bond_type, double A0_g, double ka_g, double V0, double kv)
cdef extern from "object-in-fluid/oif_local_forces.hpp":
    int oif_local_forces_set_params(int bond_type, double r0, double ks, double kslin, double phi0, double kb, double A01, double A02, double kal)

IF ROTATION:
    cdef extern from "harmonic_dumbbell.hpp":
        int harmonic_dumbbell_set_params(int bond_type, double k1, double k2, double r, double r_cut)

IF TABULATED == 1:
    cdef extern from "interaction_data.hpp":
        cdef enum TabulatedBondedInteraction:
            TAB_UNKNOWN = 0, TAB_BOND_LENGTH, TAB_BOND_ANGLE, TAB_BOND_DIHEDRAL
    cdef extern from "tab.hpp":
        int tabulated_bonded_set_params(int bond_type, TabulatedBondedInteraction tab_type, char * filename)

IF BOND_ENDANGLEDIST == 1:
    cdef extern from "endangledist.hpp":
        int endangledist_set_params(int bond_type, double bend, double phi0, double distmin, double distmax)

IF OVERLAPPED == 1:
    cdef extern from "interaction_data.hpp":
        cdef enum OverlappedBondedInteraction:
            OVERLAP_UNKNOWN = 0, OVERLAP_BOND_LENGTH, OVERLAP_BOND_ANGLE,\
                OVERLAP_BOND_DIHEDRAL
    cdef extern from "overlap.hpp":
        int overlapped_bonded_set_params(int bond_type, OverlappedBondedInteraction overlap_type,
                                         char * filename)

IF BOND_VIRTUAL == 1:
    cdef extern from "interaction_data.hpp":
        int virtual_set_params(int bond_type)


cdef extern from "interaction_data.hpp":
    cdef enum enum_bonded_interaction "BondedInteraction":
        BONDED_IA_NONE = -1,
        BONDED_IA_FENE,
        BONDED_IA_HARMONIC,
        BONDED_IA_HARMONIC_DUMBBELL,
        BONDED_IA_QUARTIC,
        BONDED_IA_BONDED_COULOMB,
        BONDED_IA_ANGLE_OLD,
        BONDED_IA_DIHEDRAL,
        BONDED_IA_TABULATED,
        BONDED_IA_SUBT_LJ,
        BONDED_IA_RIGID_BOND,
        BONDED_IA_VIRTUAL_BOND,
        BONDED_IA_ANGLEDIST,
        BONDED_IA_ENDANGLEDIST,
        BONDED_IA_OVERLAPPED,
        BONDED_IA_ANGLE_HARMONIC,
        BONDED_IA_ANGLE_COSINE,
        BONDED_IA_ANGLE_COSSQUARE,
        BONDED_IA_OIF_LOCAL_FORCES,
        BONDED_IA_OIF_GLOBAL_FORCES,
        BONDED_IA_CG_DNA_BASEPAIR,
        BONDED_IA_CG_DNA_STACKING,
        BONDED_IA_CG_DNA_BACKBONE,
        BONDED_IA_IBM_TRIEL,
        BONDED_IA_IBM_VOLUME_CONSERVATION,
        BONDED_IA_IBM_TRIBEND,
        BONDED_IA_UMBRELLA
