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

cdef extern from "script_interface/interactions/bonded.hpp":
    int bonded_ia_params_zero_based_type(int bond_id) except +
    int bonded_ia_params_size()
    int bonded_ia_params_next_key()

# Map the boost::variant type indices to python type identifiers. These enum
# values must be in the same order as in the definition of the boost::variant.
cdef enum enum_bonded_interaction:
    BONDED_IA_NONE = 0,
    BONDED_IA_FENE,
    BONDED_IA_HARMONIC,
    BONDED_IA_QUARTIC,
    BONDED_IA_BONDED_COULOMB,
    BONDED_IA_BONDED_COULOMB_SR,
    BONDED_IA_ANGLE_HARMONIC,
    BONDED_IA_ANGLE_COSINE,
    BONDED_IA_ANGLE_COSSQUARE,
    BONDED_IA_DIHEDRAL,
    BONDED_IA_TABULATED_DISTANCE,
    BONDED_IA_TABULATED_ANGLE,
    BONDED_IA_TABULATED_DIHEDRAL,
    BONDED_IA_THERMALIZED_DIST,
    BONDED_IA_RIGID_BOND,
    BONDED_IA_IBM_TRIEL,
    BONDED_IA_IBM_VOLUME_CONSERVATION,
    BONDED_IA_IBM_TRIBEND,
    BONDED_IA_OIF_GLOBAL_FORCES,
    BONDED_IA_OIF_LOCAL_FORCES,
    BONDED_IA_VIRTUAL_BOND

cdef extern from "thermostat.hpp":
    void thermalized_bond_set_rng_counter(stdint.uint64_t counter)

cdef extern from "immersed_boundary/ImmersedBoundaries.hpp":
    cppclass ImmersedBoundaries:
        double get_current_volume(int softID)

cdef extern from "immersed_boundaries.hpp":
    extern ImmersedBoundaries immersed_boundaries
