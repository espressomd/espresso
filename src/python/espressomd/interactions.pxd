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
from espressomd.system cimport *
cimport numpy as np
from espressomd.utils cimport *

cdef extern from "TabulatedPotential.hpp":
    struct TabulatedPotential:
        double maxval
        double minval
        vector[double] energy_tab
        vector[double] force_tab

cdef extern from "interaction_data.hpp":
    ctypedef struct ia_parameters "IA_parameters":
        double LJ_eps
        double LJ_sig
        double LJ_cut
        double LJ_shift
        double LJ_offset
        double LJ_min

        double LJCOS_eps
        double LJCOS_sig
        double LJCOS_cut
        double LJCOS_offset

        double LJCOS2_eps
        double LJCOS2_sig
        double LJCOS2_offset
        double LJCOS2_w

        double LJGEN_eps
        double LJGEN_sig
        double LJGEN_cut
        double LJGEN_shift
        double LJGEN_offset
        double LJGEN_a1
        double LJGEN_a2
        double LJGEN_b1
        double LJGEN_b2
        double LJGEN_lambda
        double LJGEN_softrad

        TabulatedPotential TAB

        double GB_eps
        double GB_sig
        double GB_cut
        double GB_k1
        double GB_k2
        double GB_mu
        double GB_nu

        double SmSt_eps
        double SmSt_sig
        double SmSt_cut
        double SmSt_d
        int SmSt_n
        double SmSt_k0

        double BMHTF_A;
        double BMHTF_B;
        double BMHTF_C;
        double BMHTF_D;
        double BMHTF_sig;
        double BMHTF_cut;

        double MORSE_eps
        double MORSE_alpha
        double MORSE_rmin
        double MORSE_cut
        double MORSE_rest

        double BUCK_A
        double BUCK_B
        double BUCK_C
        double BUCK_D
        double BUCK_cut
        double BUCK_discont
        double BUCK_shift

        double soft_a
        double soft_n
        double soft_cut
        double soft_offset

        double Hertzian_eps
        double Hertzian_sig

        double Gaussian_eps
        double Gaussian_sig
        double Gaussian_cut

        int dpd_wf
        int dpd_twf
        double dpd_gamma
        double dpd_r_cut
        double dpd_pref1
        double dpd_pref2
        double dpd_tgamma
        double dpd_tr_cut
        double dpd_pref3
        double dpd_pref4

        double HAT_Fmax
        double HAT_r

    cdef ia_parameters * get_ia_param(int i, int j)
    cdef ia_parameters * get_ia_param_safe(int i, int j)
    cdef void make_bond_type_exist(int type)

cdef extern from "lj.hpp":
    cdef int lennard_jones_set_params(int part_type_a, int part_type_b,
                                      double eps, double sig, double cut,
                                      double shift, double offset,
                                      double min)
IF LJCOS:
    cdef extern from "ljcos.hpp":
        cdef int ljcos_set_params(int part_type_a, int part_type_b,
                                  double eps, double sig,
                                  double cut, double offset);

IF LJCOS2:
    cdef extern from "ljcos2.hpp":
        cdef int ljcos2_set_params(int part_type_a, int part_type_b,
                                   double eps, double sig, double offset,
                                   double w)

IF GAY_BERNE:
    cdef extern from "gb.hpp":
        int gay_berne_set_params(int part_type_a, int part_type_b,
                                 double eps, double sig, double cut,
                                 double k1, double k2,
                                 double mu, double nu);

cdef extern from "ljgen.hpp":
    IF LJGEN_SOFTCORE:
        cdef int ljgen_set_params(int part_type_a, int part_type_b,
                                  double eps, double sig, double cut,
                                  double shift, double offset,
                                  double a1, double a2, double b1, double b2,
                                  double genlj_lambda, double softrad)
    ELSE:
        cdef int ljgen_set_params(int part_type_a, int part_type_b,
                                  double eps, double sig, double cut,
                                  double shift, double offset,
                                  double a1, double a2, double b1, double b2)

IF SMOOTH_STEP:
    cdef extern from "steppot.hpp":
        int smooth_step_set_params(int part_type_a, int part_type_b,
                                   double d, int n, double eps,
                                   double k0, double sig,
                                   double cut);
IF BMHTF_NACL:
    cdef extern from "bmhtf-nacl.hpp":
        int BMHTF_set_params(int part_type_a, int part_type_b,
                             double A, double B, double C,
                             double D, double sig, double cut);

IF MORSE:
    cdef extern from "morse.hpp":
        int morse_set_params(int part_type_a, int part_type_b,
                             double eps, double alpha,
                             double rmin, double cut);

IF BUCKINGHAM:
    cdef extern from "buckingham.hpp":
        int buckingham_set_params(int part_type_a, int part_type_b,
                                  double A, double B, double C, double D, double cut,
                                  double discont, double shift);

IF SOFT_SPHERE:
    cdef extern from "soft_sphere.hpp":
        int soft_sphere_set_params(int part_type_a, int part_type_b,
                                   double a, double n, double cut, double offset);

IF HERTZIAN:
    cdef extern from "hertzian.hpp":
        int hertzian_set_params(int part_type_a, int part_type_b,
                                double eps, double sig);

IF GAUSSIAN:
    cdef extern from "gaussian.hpp":
        int gaussian_set_params(int part_type_a, int part_type_b,
                                double eps, double sig, double cut);

IF DPD:
    cdef extern from "dpd.hpp":
        int dpd_set_params(int part_type_a, int part_type_b,
                           double gamma, double r_c, int wf,
                           double tgamma, double tr_c, int twf)

IF HAT:
    cdef extern from "hat.hpp":
        int hat_set_params(int part_type_a, int part_type_b,
                           double Fmax, double r)

IF TABULATED==1:
    cdef extern from "tab.hpp":
        int tabulated_set_params(int part_type_a, int part_type_b,
                                 double min, double max,
                                 vector[double] energy,
                                 vector[double] force);

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
        int    type
        TabulatedPotential * pot

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
    int subt_lj_set_params(int bond_type)
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
        int tabulated_bonded_set_params(int bond_type, TabulatedBondedInteraction tab_type, double min, double max, vector[double] energy, vector[double] force)

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
