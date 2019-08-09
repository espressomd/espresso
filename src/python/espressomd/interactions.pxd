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

# Handling of interactions

from libcpp.string cimport string
from libc cimport stdint

include "myconfig.pxi"
from espressomd.system cimport *
cimport numpy as np
#force include of config.hpp
cdef extern from "config.hpp":
    pass

from espressomd.utils cimport *

cdef extern from "TabulatedPotential.hpp":
    struct TabulatedPotential:
        double maxval
        double minval
        vector[double] energy_tab
        vector[double] force_tab

cdef extern from "dpd.hpp":
    cdef struct DPDParameters:
        double gamma
        double cutoff
        int wf
        double pref

cdef extern from "nonbonded_interactions/nonbonded_interaction_data.hpp":
    cdef struct LJ_Parameters:
        double eps
        double sig
        double cut
        double shift
        double offset
        double min

    cdef struct WCA_Parameters:
        double eps
        double sig
        double cut

    cdef struct LJGen_Parameters:
        double eps
        double sig
        double cut
        double shift
        double offset
        double a1
        double a2
        double b1
        double b2
        double lambda1
        double softrad

    cdef struct SmoothStep_Parameters:
        double eps
        double sig
        double cut
        double d
        int n
        double k0

    cdef struct Hertzian_Parameters:
        double eps
        double sig

    cdef struct Gaussian_Parameters:
        double eps
        double sig
        double cut

    cdef struct BMHTF_Parameters:
        double A
        double B
        double C
        double D
        double sig
        double cut
        double computed_shift

    cdef struct Morse_Parameters:
        double eps
        double alpha
        double rmin
        double cut
        double rest

    cdef struct Buckingham_Parameters:
        double A
        double B
        double C
        double D
        double cut
        double discont
        double shift
        double F1
        double F2

    cdef struct SoftSphere_Parameters:
        double a
        double n
        double cut
        double offset

    cdef struct Affinity_Parameters:
        int type
        double kappa
        double r0
        double Kon
        double Koff
        double maxBond
        double cut

    cdef struct Membrane_Parameters:
        double a
        double n
        double cut
        double offset

    cdef struct Hat_Parameters:
        double Fmax
        double r

    cdef struct LJcos_Parameters:
        double eps
        double sig
        double cut
        double offset

    cdef struct LJcos2_Parameters:
        double eps
        double sig
        double offset
        double w

    cdef struct GayBerne_Parameters:
        double eps
        double sig
        double cut
        double k1
        double k2
        double mu
        double nu

    cdef struct Thole_Parameters:
        double scaling_coeff
        double q1q2

    cdef struct IA_parameters:
        LJ_Parameters lj

        WCA_Parameters wca

        LJcos_Parameters ljcos

        LJcos2_Parameters ljcos2

        LJGen_Parameters ljgen

        Affinity_Parameters affinity

        Membrane_Parameters membrane

        SoftSphere_Parameters soft_sphere

        TabulatedPotential tab

        GayBerne_Parameters gay_berne

        SmoothStep_Parameters smooth_step

        BMHTF_Parameters bmhtf

        Morse_Parameters morse

        Buckingham_Parameters buckingham

        Hertzian_Parameters hertzian

        Gaussian_Parameters gaussian

        DPDParameters dpd_radial
        DPDParameters dpd_trans

        Hat_Parameters hat

        Thole_Parameters thole

    cdef IA_parameters * get_ia_param(int i, int j)
    cdef IA_parameters * get_ia_param_safe(int i, int j)
    cdef string ia_params_get_state()
    cdef void ia_params_set_state(string)
    cdef void reset_ia_params()

cdef extern from "bonded_interactions/bonded_interaction_data.hpp":
    cdef void make_bond_type_exist(int type)

cdef extern from "nonbonded_interactions/lj.hpp":
    cdef int lennard_jones_set_params(int part_type_a, int part_type_b,
                                      double eps, double sig, double cut,
                                      double shift, double offset,
                                      double min)

cdef extern from "nonbonded_interactions/wca.hpp":
    cdef int wca_set_params(int part_type_a, int part_type_b,
                            double eps, double sig)
IF LJCOS:
    cdef extern from "nonbonded_interactions/ljcos.hpp":
        cdef int ljcos_set_params(int part_type_a, int part_type_b,
                                  double eps, double sig,
                                  double cut, double offset)

IF LJCOS2:
    cdef extern from "nonbonded_interactions/ljcos2.hpp":
        cdef int ljcos2_set_params(int part_type_a, int part_type_b,
                                   double eps, double sig, double offset,
                                   double w)

IF GAY_BERNE:
    cdef extern from "nonbonded_interactions/gay_berne.hpp":
        int gay_berne_set_params(int part_type_a, int part_type_b,
                                 double eps, double sig, double cut,
                                 double k1, double k2,
                                 double mu, double nu)

IF THOLE:
    cdef extern from "nonbonded_interactions/thole.hpp":
        int thole_set_params(int part_type_a, int part_type_b, double scaling_coeff, double q1q2)

cdef extern from "nonbonded_interactions/ljgen.hpp":
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
    cdef extern from "nonbonded_interactions/smooth_step.hpp":
        int smooth_step_set_params(int part_type_a, int part_type_b,
                                   double d, int n, double eps,
                                   double k0, double sig,
                                   double cut)
IF BMHTF_NACL:
    cdef extern from "nonbonded_interactions/bmhtf-nacl.hpp":
        int BMHTF_set_params(int part_type_a, int part_type_b,
                             double A, double B, double C,
                             double D, double sig, double cut)

IF MORSE:
    cdef extern from "nonbonded_interactions/morse.hpp":
        int morse_set_params(int part_type_a, int part_type_b,
                             double eps, double alpha,
                             double rmin, double cut)

IF BUCKINGHAM:
    cdef extern from "nonbonded_interactions/buckingham.hpp":
        int buckingham_set_params(int part_type_a, int part_type_b,
                                  double A, double B, double C, double D, double cut,
                                  double discont, double shift)

IF SOFT_SPHERE:
    cdef extern from "nonbonded_interactions/soft_sphere.hpp":
        int soft_sphere_set_params(int part_type_a, int part_type_b,
                                   double a, double n, double cut, double offset)

IF HERTZIAN:
    cdef extern from "nonbonded_interactions/hertzian.hpp":
        int hertzian_set_params(int part_type_a, int part_type_b,
                                double eps, double sig)

IF GAUSSIAN:
    cdef extern from "nonbonded_interactions/gaussian.hpp":
        int gaussian_set_params(int part_type_a, int part_type_b,
                                double eps, double sig, double cut)

IF DPD:
    cdef extern from "dpd.hpp":
        int dpd_set_params(int part_type_a, int part_type_b,
                           double gamma, double r_c, int wf,
                           double tgamma, double tr_c, int twf)

IF HAT:
    cdef extern from "nonbonded_interactions/hat.hpp":
        int hat_set_params(int part_type_a, int part_type_b,
                           double Fmax, double r)

IF MEMBRANE_COLLISION:
    cdef extern from "object-in-fluid/membrane_collision.hpp":
        cdef int membrane_collision_set_params(int part_type_a,
                                               int part_type_b,
                                               double a, double n,
                                               double cut, double offset)

IF SOFT_SPHERE:
    cdef extern from "nonbonded_interactions/soft_sphere.hpp":
        cdef int soft_sphere_set_params(int part_type_a, int part_type_b,
                                        double a, double n,
                                        double cut, double offset)

IF AFFINITY:
    cdef extern from "object-in-fluid/affinity.hpp":
        cdef int affinity_set_params(int part_type_a, int part_type_b,
                                     int afftype, double kappa, double r0,
                                     double Kon, double Koff, double maxBond, double cut)
IF TABULATED:
    cdef extern from "nonbonded_interactions/nonbonded_tab.hpp":
        int tabulated_set_params(int part_type_a, int part_type_b,
                                 double min, double max,
                                 vector[double] energy,
                                 vector[double] force)
IF ROTATION:
    cdef extern from "bonded_interactions/bonded_interaction_data.hpp":
    #* Parameters for the harmonic dumbbell bond potential */
        cdef struct Harmonic_dumbbell_bond_parameters:
            double k1
            double k2
            double r
            double r_cut
ELSE:
    cdef struct Harmonic_dumbbell_bond_parameters:
        double k1
        double k2
        double r
        double r_cut

cdef extern from "bonded_interactions/bonded_interaction_data.hpp":
    #* Parameters for n-body tabulated potential (n=2,3,4). */
    cdef struct Tabulated_bond_parameters:
        int type
        TabulatedPotential * pot

IF ELECTROSTATICS:
    cdef extern from "bonded_interactions/bonded_interaction_data.hpp":
        #* Parameters for Bonded Coulomb p3m sr */
        cdef struct Bonded_coulomb_sr_bond_parameters:
            double q1q2
ELSE:
    cdef struct Bonded_coulomb_sr_bond_parameters:
        double q1q2

cdef extern from "bonded_interactions/bonded_interaction_data.hpp":
    cdef struct Fene_bond_parameters:
        double k
        double drmax
        double r0
        double drmax2
        double drmax2i

    #* Parameters for oif_global_forces */
    cdef struct Oif_global_forces_bond_parameters:
        double A0_g
        double ka_g
        double V0
        double kv

    #* Parameters for oif_local_forces */
    cdef struct Oif_local_forces_bond_parameters:
        double r0
        double ks
        double kslin
        double phi0
        double kb
        double A01
        double A02
        double kal
        double kvisc

    #* Parameters for harmonic bond Potential */
    cdef struct Harmonic_bond_parameters:
        double k
        double r
        double r_cut

    #* Parameters for thermalized  bond */
    cdef struct Thermalized_bond_parameters:
        double temp_com
        double gamma_com
        double temp_distance
        double gamma_distance
        double r_cut

    #* Parameters for Bonded Coulomb */
    cdef struct Bonded_coulomb_bond_parameters:
        double prefactor

    #* Parameters for three body angular potential (bond_angle_harmonic).
    cdef struct Angle_harmonic_bond_parameters:
        double bend
        double phi0

    #* Parameters for three body angular potential (bond_angle_cosine).
    cdef struct Angle_cosine_bond_parameters:
        double bend
        double phi0
        double cos_phi0
        double sin_phi0

    #* Parameters for three body angular potential (bond_angle_cossquare).
    cdef struct Angle_cossquare_bond_parameters:
        double bend
        double phi0
        double cos_phi0

    #* Parameters for four body angular potential (dihedral-angle potentials). */
    cdef struct Dihedral_bond_parameters:
        double mult
        double bend
        double phase

    #* Parameters for n-body overlapped potential (n=2,3,4). */
    cdef struct Overlap_bond_parameters:
        char * filename
        int    type
        double maxval
        int    noverlaps
        double * para_a
        double * para_b
        double * para_c

    #* Parameters for one-directional harmonic potential */
    cdef struct Umbrella_bond_parameters:
        double k
        int    dir
        double r

    #* Parameters for subt-LJ potential */
    cdef struct Subt_lj_bond_parameters:
        double k
        double r
        double r2

    #* Parameters for the rigid_bond/SHAKE/RATTLE ALGORITHM */
    cdef struct Rigid_bond_parameters:
        #*Length of rigid bond/Constrained Bond*/
        # double d
        #*Square of the length of Constrained Bond*/
        double d2
        #*Positional Tolerance/Accuracy value for termination of RATTLE/SHAKE iterations during position corrections*/
        double p_tol
        #*Velocity Tolerance/Accuracy for termination of RATTLE/SHAKE iterations during velocity corrections */
        double v_tol

    #* Parameters for IBM Triel */
    cdef cppclass tElasticLaw:
        pass

    cdef struct IBM_Triel_Parameters:
        double l0
        double lp0
        double sinPhi0
        double cosPhi0
        double area0
        double a1
        double a2
        double b1
        double b2
        double maxDist
        tElasticLaw elasticLaw
        double k1
        double k2

    #* Parameters for IBM Tribend */
    cdef struct IBM_Tribend_Parameters:
        double kb
        double theta0

    #* Parameters for IBM VolCons */
    cdef struct IBM_VolCons_Parameters:
        int softID
        double kappaV
        double volRef

    #* Parameters for Quartic */
    cdef struct Quartic_bond_parameters:
        double k0, k1
        double r
        double r_cut

    #* Union in which to store the parameters of an individual bonded interaction */
    cdef union Bond_parameters:
        Fene_bond_parameters fene
        Oif_global_forces_bond_parameters oif_global_forces
        Oif_local_forces_bond_parameters oif_local_forces
        Thermalized_bond_parameters thermalized_bond
        Bonded_coulomb_bond_parameters bonded_coulomb
        Bonded_coulomb_sr_bond_parameters bonded_coulomb_sr
        Harmonic_bond_parameters harmonic
        Harmonic_dumbbell_bond_parameters harmonic_dumbbell
        Angle_harmonic_bond_parameters angle_harmonic
        Angle_cosine_bond_parameters angle_cosine
        Angle_cossquare_bond_parameters angle_cossquare
        Dihedral_bond_parameters dihedral
        Tabulated_bond_parameters tab
        Overlap_bond_parameters overlap
        Subt_lj_bond_parameters subt_lj
        Rigid_bond_parameters rigid_bond
        IBM_Triel_Parameters ibm_triel
        IBM_Tribend_Parameters ibm_tribend
        IBM_VolCons_Parameters ibmVolConsParameters
        Quartic_bond_parameters quartic

    cdef struct Bonded_ia_parameters:
        int type
        int num
        #* union to store the different bonded interaction parameters. */
        Bond_parameters p

    vector[Bonded_ia_parameters] bonded_ia_params

cdef extern from "bonded_interactions/bonded_interaction_data.hpp" namespace "tElasticLaw":
    cdef tElasticLaw NeoHookean
    cdef tElasticLaw Skalak

cdef extern from "bonded_interactions/fene.hpp":
    int fene_set_params(int bond_type, double k, double drmax, double r0)
cdef extern from "bonded_interactions/harmonic.hpp":
    int harmonic_set_params(int bond_type, double k, double r, double r_cut)
cdef extern from "bonded_interactions/dihedral.hpp":
    int dihedral_set_params(int bond_type, int mult, double bend, double phase)
cdef extern from "bonded_interactions/angle_harmonic.hpp":
    int angle_harmonic_set_params(int bond_type, double bend, double phi0)
cdef extern from "rattle.hpp":
    int rigid_bond_set_params(int bond_type, double d, double p_tol, double v_tol)
cdef extern from "bonded_interactions/angle_cosine.hpp":
    int angle_cosine_set_params(int bond_type, double bend, double phi0)
cdef extern from "bonded_interactions/angle_cossquare.hpp":
    int angle_cossquare_set_params(int bond_type, double bend, double phi0)
cdef extern from "bonded_interactions/subt_lj.hpp":
    int subt_lj_set_params(int bond_type)
cdef extern from "object-in-fluid/oif_global_forces.hpp":
    int oif_global_forces_set_params(int bond_type, double A0_g, double ka_g, double V0, double kv)
cdef extern from "object-in-fluid/oif_local_forces.hpp":
    int oif_local_forces_set_params(int bond_type, double r0, double ks, double kslin, double phi0, double kb, double A01, double A02, double kal, double kvisc)
cdef extern from "object-in-fluid/out_direction.hpp":
    int oif_out_direction_set_params(int bond_type)
cdef extern from "bonded_interactions/thermalized_bond.hpp":
    int thermalized_bond_set_params(int bond_type, double temp_com, double gamma_com, double temp_distance, double gamma_distance, double r_cut)
cdef extern from "bonded_interactions/bonded_coulomb.hpp":
    int bonded_coulomb_set_params(int bond_type, double prefactor)
cdef extern from "bonded_interactions/quartic.hpp":
    int quartic_set_params(int bond_type, double k0, double k1, double r, double r_cut)

cdef extern from "immersed_boundary/ImmersedBoundaries.hpp":
    cppclass ImmersedBoundaries:
        void volume_conservation_set_params(const int bond_type, const int softID, const double kappaV)

cdef extern from "immersed_boundary/ibm_triel.hpp":
    int IBM_Triel_SetParams(const int bond_type, const int ind1, const int ind2, const int ind3, const double max, const tElasticLaw elasticLaw, const double k1, const double k2)
cdef extern from "immersed_boundary/ibm_tribend.hpp":
    int IBM_Tribend_SetParams(const int bond_type, const int ind1, const int ind2, const int ind3, const int ind4, const double kb, const bool flat)

IF ROTATION:
    cdef extern from "bonded_interactions/harmonic_dumbbell.hpp":
        int harmonic_dumbbell_set_params(int bond_type, double k1, double k2, double r, double r_cut)

cdef extern from "bonded_interactions/bonded_interaction_data.hpp":
    cdef enum TabulatedBondedInteraction:
        TAB_UNKNOWN = 0, TAB_BOND_LENGTH, TAB_BOND_ANGLE, TAB_BOND_DIHEDRAL

cdef extern from "bonded_interactions/bonded_tab.hpp":
    int tabulated_bonded_set_params(int bond_type, TabulatedBondedInteraction tab_type, double min, double max, vector[double] energy, vector[double] force)

IF ELECTROSTATICS:
    cdef extern from "bonded_interactions/bonded_coulomb.hpp":
        int bonded_coulomb_set_params(int bond_type, double prefactor)

    cdef extern from "bonded_interactions/bonded_coulomb_sr.hpp":
        int bonded_coulomb_sr_set_params(int bond_type, double q1q2)

cdef extern from "bonded_interactions/bonded_interaction_data.hpp":
    int virtual_set_params(int bond_type)

cdef extern from "bonded_interactions/bonded_interaction_data.hpp":
    cdef enum enum_bonded_interaction "BondedInteraction":
        BONDED_IA_NONE = -1,
        BONDED_IA_FENE,
        BONDED_IA_HARMONIC,
        BONDED_IA_HARMONIC_DUMBBELL,
        BONDED_IA_BONDED_COULOMB,
        BONDED_IA_BONDED_COULOMB_SR,
        BONDED_IA_DIHEDRAL,
        BONDED_IA_TABULATED_DISTANCE,
        BONDED_IA_TABULATED_ANGLE,
        BONDED_IA_TABULATED_DIHEDRAL,
        BONDED_IA_SUBT_LJ,
        BONDED_IA_RIGID_BOND,
        BONDED_IA_VIRTUAL_BOND,
        BONDED_IA_ANGLE_HARMONIC,
        BONDED_IA_ANGLE_COSINE,
        BONDED_IA_ANGLE_COSSQUARE,
        BONDED_IA_OIF_LOCAL_FORCES,
        BONDED_IA_OIF_GLOBAL_FORCES,
        BONDED_IA_OIF_OUT_DIRECTION,
        BONDED_IA_IBM_TRIEL,
        BONDED_IA_IBM_TRIBEND,
        BONDED_IA_IBM_VOLUME_CONSERVATION,
        BONDED_IA_UMBRELLA,
        BONDED_IA_THERMALIZED_DIST
        BONDED_IA_QUARTIC
