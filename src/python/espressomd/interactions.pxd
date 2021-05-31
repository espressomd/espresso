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

# Handling of interactions

from libcpp.memory cimport shared_ptr
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as cbool
from libc cimport stdint

from .thermostat cimport thermalized_bond

include "myconfig.pxi"

# force include of config.hpp
cdef extern from "config.hpp":
    pass

cdef extern from "TabulatedPotential.hpp":
    struct TabulatedPotential:
        double maxval
        double minval
        vector[double] energy_tab
        vector[double] force_tab

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

    cdef struct DPDParameters:
        double gamma
        double k
        double cutoff
        int wf
        double pref

    cdef struct IA_parameters:
        LJ_Parameters lj

        WCA_Parameters wca

        LJcos_Parameters ljcos

        LJcos2_Parameters ljcos2

        LJGen_Parameters ljgen

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
                           double gamma, double k, double r_c, int wf,
                           double tgamma, double tr_c, int twf)

IF HAT:
    cdef extern from "nonbonded_interactions/hat.hpp":
        int hat_set_params(int part_type_a, int part_type_b,
                           double Fmax, double r)

IF SOFT_SPHERE:
    cdef extern from "nonbonded_interactions/soft_sphere.hpp":
        cdef int soft_sphere_set_params(int part_type_a, int part_type_b,
                                        double a, double n,
                                        double cut, double offset)

IF TABULATED:
    cdef extern from "nonbonded_interactions/nonbonded_tab.hpp":
        int tabulated_set_params(int part_type_a, int part_type_b,
                                 double min, double max,
                                 vector[double] energy,
                                 vector[double] force)

    cdef extern from "bonded_interactions/bonded_interaction_data.hpp":
        cdef cppclass CoreTabulatedBond "TabulatedBond":
            CoreTabulatedBond()
            CoreTabulatedBond(double min, double max,
                              const vector[double] & energy,
                              const vector[double] & force) except +
            shared_ptr[TabulatedPotential] pot

        cdef cppclass CoreTabulatedDistanceBond "TabulatedDistanceBond"(CoreTabulatedBond):
            CoreTabulatedDistanceBond()
            CoreTabulatedDistanceBond(double min, double max,
                                      const vector[double] & energy,
                                      const vector[double] & force) except +

        cdef cppclass CoreTabulatedAngleBond "TabulatedAngleBond"(CoreTabulatedBond):
            CoreTabulatedAngleBond()
            CoreTabulatedAngleBond(double min, double max,
                                   const vector[double] & energy,
                                   const vector[double] & force) except +

        cdef cppclass CoreTabulatedDihedralBond "TabulatedDihedralBond"(CoreTabulatedBond):
            CoreTabulatedDihedralBond()
            CoreTabulatedDihedralBond(double min, double max,
                                      const vector[double] & energy,
                                      const vector[double] & force) except +

IF ELECTROSTATICS:
    cdef extern from "bonded_interactions/bonded_coulomb.hpp":
        cdef cppclass CoreBondedCoulomb "BondedCoulomb":
            CoreBondedCoulomb()
            CoreBondedCoulomb(double prefactor)
            double prefactor

    cdef extern from "bonded_interactions/bonded_coulomb_sr.hpp":
        cdef cppclass CoreBondedCoulombSR "BondedCoulombSR":
            CoreBondedCoulombSR()
            CoreBondedCoulombSR(double q1q2)
            double q1q2

ELSE:
    cdef struct CoreBondedCoulomb:
        double prefactor

    cdef struct CoreBondedCoulombSR:
        double q1q2


cdef extern from "bonded_interactions/bonded_interaction_data.hpp":
    cdef struct CoreNoneBond "NoneBond":
        pass

    cdef struct CoreVirtualBond "VirtualBond":
        pass

    # Parameters for FENE bond
    cdef cppclass CoreFeneBond "FeneBond":
        CoreFeneBond()
        CoreFeneBond(double k, double drmax, double r0)
        double k
        double drmax
        double r0
        double drmax2
        double drmax2i

    # Parameters for oif_global_forces
    cdef cppclass CoreOifGlobalForcesBond "OifGlobalForcesBond":
        CoreOifGlobalForcesBond()
        CoreOifGlobalForcesBond(double A0_g, double ka_g,
                                double V0, double kv)
        double A0_g
        double ka_g
        double V0
        double kv

    # Parameters for oif_local_forces
    cdef cppclass CoreOifLocalForcesBond "OifLocalForcesBond":
        CoreOifLocalForcesBond()
        CoreOifLocalForcesBond(double r0, double ks,
                               double kslin, double phi0, double kb,
                               double A01, double A02, double kal,
                               double kvisc)
        double r0
        double ks
        double kslin
        double phi0
        double kb
        double A01
        double A02
        double kal
        double kvisc

    # Parameters for harmonic bond Potential
    cdef cppclass CoreHarmonicBond "HarmonicBond":
        CoreHarmonicBond()
        CoreHarmonicBond(double k, double r, double r_cut)
        double k
        double r
        double r_cut

    # Parameters for thermalized  bond
    cdef cppclass CoreThermalizedBond "ThermalizedBond":
        CoreThermalizedBond()
        CoreThermalizedBond(double temp_com,
                            double gamma_com, double temp_distance,
                            double gamma_distance, double r_cut)
        double temp_com
        double gamma_com
        double temp_distance
        double gamma_distance
        double r_cut

    # Parameters for three body angular potential (bond_angle_harmonic).
    cdef cppclass CoreAngleHarmonicBond "AngleHarmonicBond":
        CoreAngleHarmonicBond()
        CoreAngleHarmonicBond(double bend, double phi0)
        double bend
        double phi0

    # Parameters for three body angular potential (bond_angle_cosine).
    cdef cppclass CoreAngleCosineBond "AngleCosineBond":
        CoreAngleCosineBond()
        CoreAngleCosineBond(double bend, double phi0)
        double bend
        double phi0
        double cos_phi0
        double sin_phi0

    # Parameters for three body angular potential (bond_angle_cossquare).
    cdef cppclass CoreAngleCossquareBond "AngleCossquareBond":
        CoreAngleCossquareBond()
        CoreAngleCossquareBond(double bend, double phi0)
        double bend
        double phi0
        double cos_phi0

    # Parameters for four body angular potential (dihedral-angle potentials).
    cdef cppclass CoreDihedralBond "DihedralBond":
        CoreDihedralBond()
        CoreDihedralBond(int mult, double bend, double phase)
        double mult
        double bend
        double phase

    # Parameters for the rigid_bond/SHAKE/RATTLE ALGORITHM
    cdef cppclass CoreRigidBond "RigidBond":
        CoreRigidBond()
        CoreRigidBond(double d, double p_tol, double v_tol)
        # Length of rigid bond/Constrained Bond
        # double d
        # Square of the length of Constrained Bond
        double d2
        # Positional Tolerance/Accuracy value for termination of RATTLE/SHAKE
        # iterations during position corrections
        double p_tol
        # Velocity Tolerance/Accuracy for termination of RATTLE/SHAKE
        # iterations during velocity corrections
        double v_tol

    # Parameters for IBM Triel
    cdef cppclass tElasticLaw:
        pass

    cdef cppclass CoreIBMTriel "IBMTriel":
        CoreIBMTriel()
        CoreIBMTriel(int ind1, int ind2, int ind3,
                     double maxDist, tElasticLaw elasticLaw, double k1,
                     double k2)
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

    # Parameters for IBM Tribend
    cdef cppclass CoreIBMTribend "IBMTribend":
        CoreIBMTribend()
        CoreIBMTribend(int ind1, int ind2, int ind3, int ind4,
                       double kb, cbool flat)
        double kb
        double theta0

    # Parameters for IBM VolCons
    cdef cppclass CoreIBMVolCons "IBMVolCons":
        CoreIBMVolCons()
        CoreIBMVolCons(int softID, double kappaV)
        int softID
        double kappaV
        double volRef

    # Parameters for Quartic
    cdef cppclass CoreQuarticBond "QuarticBond":
        CoreQuarticBond()
        CoreQuarticBond(double k0, double k1, double r, double r_cut)
        double k0, k1
        double r
        double r_cut

cdef extern from "script_interface/interactions/bonded.hpp":
    T bonded_ia_params_at[T](int bond_id) except +
    cbool bonded_ia_params_is_type[T](int bond_id) except +
    int bonded_ia_params_num_partners(int bond_id) except +
    int bonded_ia_params_zero_based_type(int bond_id)
    int bonded_ia_params_size()

cdef extern from "bonded_interactions/bonded_interaction_utils.hpp":
    void set_bonded_ia_params[T](int bond_id, const T & iaparams)

cdef extern from "bonded_interactions/bonded_interaction_data.hpp" namespace "tElasticLaw":
    cdef tElasticLaw NeoHookean
    cdef tElasticLaw Skalak

cdef extern from "thermostat.hpp":
    void thermalized_bond_set_rng_seed(stdint.uint32_t seed)
    void thermalized_bond_set_rng_counter(stdint.uint64_t counter)

cdef extern from "immersed_boundary/ImmersedBoundaries.hpp":
    cppclass ImmersedBoundaries:
        pass

cdef extern from "immersed_boundaries.hpp":
    extern ImmersedBoundaries immersed_boundaries
