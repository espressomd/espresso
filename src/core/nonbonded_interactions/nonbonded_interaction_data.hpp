/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _INTERACTION_DATA_H
#define _INTERACTION_DATA_H
/** \file
    Various procedures concerning interactions between particles.
*/

#include "particle_data.hpp"
#include "utils.hpp"

#include "TabulatedPotential.hpp"

/** cutoff for deactivated interactions. Below 0, so that even particles on
    top of each other don't interact by chance. */
#define INACTIVE_CUTOFF -1.0

/*@}*/

/** \name Type codes for the type of Coulomb interaction
    Enumeration of implemented methods for the electrostatic
    interaction.
*/
/************************************************************/
/*@{*/

#ifdef ELECTROSTATICS
enum CoulombMethod {
  COULOMB_NONE,      //< Coulomb interaction switched off (NONE)
  COULOMB_DH,        //< Coulomb method is Debye-Hueckel
  COULOMB_P3M,       //< Coulomb method is P3M
  COULOMB_MMM1D,     //< Coulomb method is one-dimensional MMM
  COULOMB_MMM2D,     //< Coulomb method is two-dimensional MMM
  COULOMB_MAGGS,     //< Coulomb method is "Maggs"
  COULOMB_ELC_P3M,   //< Coulomb method is P3M plus ELC
  COULOMB_RF,        //< Coulomb method is Reaction-Field
  COULOMB_P3M_GPU,   //< Coulomb method is P3M with GPU based long range part
                     // calculation
  COULOMB_MMM1D_GPU, //< Coulomb method is one-dimensional MMM running on GPU
  COULOMB_SCAFACOS,  //< Coulomb method is scafacos
};

#endif
/*@}*/

#ifdef DIPOLES
/** \name Type codes for the type of dipolar interaction
  Enumeration of implemented methods for the magnetostatic
  interaction.
 */
/************************************************************/
/*@{*/
enum DipolarInteraction {
  /** dipolar interaction switched off (NONE). */
  DIPOLAR_NONE = 0,
  /** dipolar method is P3M. */
  DIPOLAR_P3M,
  /** Dipolar method is P3M plus DLC. */
  DIPOLAR_MDLC_P3M,
  /** Dipolar method is all with all and no replicas */
  DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA,
  /** Dipolar method is magnetic dipolar direct sum */
  DIPOLAR_DS,
  /** Dipolar method is direct sum plus DLC. */
  DIPOLAR_MDLC_DS,
  /** Direct summation on gpu */
  DIPOLAR_DS_GPU,
#ifdef DIPOLAR_BARNES_HUT
  /** Direct summation on gpu by Barnes-Hut algorithm */
  DIPOLAR_BH_GPU,
#endif
  /** Scafacos library */
  DIPOLAR_SCAFACOS
};
#endif

/* Data Types */
/************************************************************/

/** field containing the interaction parameters for
 *  nonbonded interactions. Access via
 * get_ia_param(i, j), i,j < max_seen_particle_type */
struct IA_parameters {
  /** maximal cutoff for this pair of particle types. This contains
      contributions from the short-ranged interactions, plus any
      cutoffs from global interactions like electrostatics.
  */
  double max_cut = INACTIVE_CUTOFF;

#ifdef LENNARD_JONES
  /** \name Lennard-Jones with shift */
  /*@{*/
  double LJ_eps = 0.0;
  double LJ_sig = 0.0;
  double LJ_cut = 0.0;
  double LJ_shift = 0.0;
  double LJ_offset = 0.0;
  double LJ_min = 0.0;
  /*@}*/

#endif

#ifdef WCA
  double WCA_eps = 0.0;
  double WCA_sig = 0.0;
  double WCA_cut = INACTIVE_CUTOFF;
#endif

  /** flag that tells whether there is any short-ranged interaction,
      i.e. one that contributes to the "nonbonded" section of the
      energy/pressure. Note that even if there is no short-ranged
      interaction present, the \ref max_cut can be non-zero due to
      e.g. electrostatics. */
  int particlesInteract;

#ifdef LENNARD_JONES_GENERIC
  /** \name Generic Lennard-Jones with shift */
  /*@{*/
  double LJGEN_eps = 0.0;
  double LJGEN_sig = 0.0;
  double LJGEN_cut = INACTIVE_CUTOFF;
  double LJGEN_shift = 0.0;
  double LJGEN_offset = 0.0;
  double LJGEN_a1 = 0.0;
  double LJGEN_a2 = 0.0;
  double LJGEN_b1 = 0.0;
  double LJGEN_b2 = 0.0;
  double LJGEN_lambda = 1.0;
  double LJGEN_softrad = 0.0;
/*@}*/
#endif

#ifdef SMOOTH_STEP
  /** \name smooth step potential */
  /*@{*/
  double SmSt_eps = 0.0;
  double SmSt_sig = 0.0;
  double SmSt_cut = INACTIVE_CUTOFF;
  double SmSt_d = 0.0;
  int SmSt_n = 0.0;
  double SmSt_k0 = 0.0;
/*@}*/
#endif

#ifdef HERTZIAN
  /** \name Hertzian potential */
  /*@{*/
  double Hertzian_eps = 0.0;
  double Hertzian_sig = INACTIVE_CUTOFF;
/*@}*/
#endif

#ifdef GAUSSIAN
  /** \name Gaussian potential */
  /*@{*/
  double Gaussian_eps = 0.0;
  double Gaussian_sig = 1.0;
  double Gaussian_cut = INACTIVE_CUTOFF;
/*@}*/
#endif

#ifdef BMHTF_NACL
  /** \name BMHTF NaCl potential */
  /*@{*/
  double BMHTF_A = 0.0;
  double BMHTF_B = 0.0;
  double BMHTF_C = 0.0;
  double BMHTF_D = 0.0;
  double BMHTF_sig = 0.0;
  double BMHTF_cut = INACTIVE_CUTOFF;
  double BMHTF_computed_shift = 0.0;
/*@}*/
#endif

#ifdef MORSE
  /** \name Morse potential */
  /*@{*/
  double MORSE_eps = INACTIVE_CUTOFF;
  double MORSE_alpha = INACTIVE_CUTOFF;
  double MORSE_rmin = INACTIVE_CUTOFF;
  double MORSE_cut = INACTIVE_CUTOFF;
  double MORSE_rest = INACTIVE_CUTOFF;
/*@}*/
#endif

#ifdef BUCKINGHAM
  /** \name Buckingham potential */
  /*@{*/
  double BUCK_A = 0.0;
  double BUCK_B = 0.0;
  double BUCK_C = 0.0;
  double BUCK_D = 0.0;
  double BUCK_cut = INACTIVE_CUTOFF;
  double BUCK_discont = 0.0;
  double BUCK_shift = 0.0;
  double BUCK_F1 = 0.0;
  double BUCK_F2 = 0.0;
/*@}*/
#endif

#ifdef SOFT_SPHERE
  /** \name soft-sphere potential */
  /*@{*/
  double soft_a = 0.0;
  double soft_n = 0.0;
  double soft_cut = INACTIVE_CUTOFF;
  double soft_offset = 0.0;
/*@}*/
#endif

#ifdef AFFINITY
  /** \name affinity potential */
  /*@{*/
  int affinity_type = INACTIVE_CUTOFF;
  double affinity_kappa = INACTIVE_CUTOFF;
  double affinity_r0 = INACTIVE_CUTOFF;
  double affinity_Kon = INACTIVE_CUTOFF;
  double affinity_Koff = INACTIVE_CUTOFF;
  double affinity_maxBond = INACTIVE_CUTOFF;
  double affinity_cut = INACTIVE_CUTOFF;
/*@}*/
#endif

#ifdef MEMBRANE_COLLISION
  /** \name membrane collision potential */
  /*@{*/
  double membrane_a = 0.0;
  double membrane_n = 0.0;
  double membrane_cut = INACTIVE_CUTOFF;
  double membrane_offset = 0.0;
/*@}*/
#endif

#ifdef HAT
  /** \name hat potential */
  /*@{*/
  double HAT_Fmax = 0.0;
  double HAT_r = INACTIVE_CUTOFF;
/*@}*/
#endif

#ifdef LJCOS
  /** \name Lennard-Jones+Cos potential */
  /*@{*/
  double LJCOS_eps = 0.0;
  double LJCOS_sig = 0.0;
  double LJCOS_cut = INACTIVE_CUTOFF;
  double LJCOS_offset = 0.0;
  double LJCOS_alfa = 0.0;
  double LJCOS_beta = 0.0;
  double LJCOS_rmin = 0.0;
/*@}*/
#endif

#ifdef LJCOS2
  /** \name Lennard-Jones with a different Cos potential */
  /*@{*/
  double LJCOS2_eps = 0.0;
  double LJCOS2_sig = 0.0;
  double LJCOS2_cut = INACTIVE_CUTOFF;
  double LJCOS2_offset = 0.0;
  double LJCOS2_w = 0.0;
  double LJCOS2_rchange = 0.0;
/*@}*/
#endif

#ifdef COS2
  /** \name Cos2 potential */
  /*@{*/
  double COS2_eps = INACTIVE_CUTOFF;
  double COS2_cut = INACTIVE_CUTOFF;
  double COS2_offset = INACTIVE_CUTOFF;
  double COS2_w = INACTIVE_CUTOFF;
/*@}*/
#endif

#ifdef GAY_BERNE
  /** \name Gay-Berne potential */
  /*@{*/
  double GB_eps = 0.0;
  double GB_sig = 0.0;
  double GB_cut = INACTIVE_CUTOFF;
  double GB_k1 = 0.0;
  double GB_k2 = 0.0;
  double GB_mu = 0.0;
  double GB_nu = 0.0;
  double GB_chi1 = 0.0;
  double GB_chi2 = 0.0;
/*@}*/
#endif

#ifdef TABULATED
  /** \name Tabulated potential */
  /*@{*/
  TabulatedPotential TAB;
/*@}*/
#endif

#ifdef DPD
  /** \name DPD as interaction */
  /*@{*/
  int dpd_wf = 0;
  int dpd_twf = 0;
  double dpd_gamma = 0.0;
  double dpd_r_cut = INACTIVE_CUTOFF;
  double dpd_pref1 = 0.0;
  double dpd_pref2 = 0.0;
  double dpd_tgamma = 0.0;
  double dpd_tr_cut = INACTIVE_CUTOFF;
  double dpd_pref3 = 0.0;
  double dpd_pref4 = 0.0;
/*@}*/
#endif

#ifdef THOLE
  double THOLE_scaling_coeff;
  double THOLE_q1q2;
#endif

#ifdef SWIMMER_REACTIONS
  double REACTION_range = INACTIVE_CUTOFF;
#endif

#ifdef SHANCHEN
  double affinity[LB_COMPONENTS];
  int affinity_on = 0;
#endif
};

extern std::vector<IA_parameters> ia_params;

/** thermodynamic force parameters */

/** \name Compounds for Coulomb interactions */
/*@{*/

/** field containing the interaction parameters for
 *  the Coulomb  interaction.  */
struct Coulomb_parameters {

#ifdef ELECTROSTATICS
  /** bjerrum length times temperature. */
  double prefactor;

  /** Method to treat Coulomb interaction. */
  CoulombMethod method;
#endif

#ifdef DIPOLES
  double Dprefactor;
  DipolarInteraction Dmethod;
#endif
};

#ifdef ELECTROSTATICS

/** Induced field (for const. potential feature). **/
extern double field_induced;
/** Applied field (for const. potential feature) **/
extern double field_applied;

#endif

/************************************************
 * exported variables
 ************************************************/

/** Maximal particle type seen so far. */
extern int max_seen_particle_type;

/** Structure containing the Coulomb parameters. */
extern Coulomb_parameters coulomb;

/** Maximal interaction cutoff (real space/short range interactions). */
extern double max_cut;
/** Maximal interaction cutoff (real space/short range non-bonded interactions).
 */
extern double max_cut_nonbonded;
/** Cutoff of Coulomb real space part */
extern double coulomb_cutoff;
/** Cutoff of dipolar real space part */
extern double dipolar_cutoff;

/** Minimal global interaction cutoff. Particles with a distance
    smaller than this are guaranteed to be available on the same node
    (through ghosts).  */
extern double min_global_cut;

/************************************************
 * exported functions
 ************************************************/

#ifdef ELECTROSTATICS
/** @brief Set the electrostatics prefactor */
int coulomb_set_prefactor(double prefactor);

/** @brief Deactivates the current Coulomb method
    This was part of coulomb_set_bjerrum()
*/
void deactivate_coulomb_method();
#endif

#ifdef DIPOLES
/** @brief Set the dipolar prefactor */
int dipolar_set_Dprefactor(double prefactor);
#endif

/** get interaction parameters between particle sorts i and j */
inline IA_parameters *get_ia_param(int i, int j) {
  extern std::vector<IA_parameters> ia_params;
  extern int max_seen_particle_type;
  return &ia_params[i * max_seen_particle_type + j];
}

/** get interaction parameters between particle sorts i and j.
    Slower than @ref get_ia_param, but can also be used on not
    yet present particle types*/
IA_parameters *get_ia_param_safe(int i, int j);

/** @brief Get the state of all non bonded interactions.
 */
std::string ia_params_get_state();

/** @brief Set the state of all non bonded interactions.
 */
void ia_params_set_state(std::string const &);

bool is_new_particle_type(int type);
/** Makes sure that ia_params is large enough to cover interactions
    for this particle type. The interactions are initialized with values
    such that no physical interaction occurs. */
void make_particle_type_exist(int type);

void make_particle_type_exist_local(int type);

/** This function increases the LOCAL ia_params field
    to the given size. Better use
    \ref make_particle_type_exist since it takes care of
    the other nodes.  */
void realloc_ia_params(int nsize);

/** calculates the maximal cutoff of all real space
    interactions. these are: bonded, non bonded + real space
    electrostatics. The result is stored in the global variable
    max_cut. The maximal cutoff of the non-bonded + real space
    electrostatic interactions is stored in max_cut_non_bonded. This
    value is used in the Verlet pair list algorithm. */
void recalc_maximal_cutoff();

/**
 * @brief Reset all interaction parameters to their defaults.
 */
void reset_ia_params();

/** check whether all force calculation routines are properly initialized. */
int interactions_sanity_checks();

/**  check if a non bonded interaction is defined */
inline int checkIfInteraction(IA_parameters *data) {
  return data->particlesInteract;
}

/** check if the types of particles i and j have any non bonded
    interaction defined. */
inline int checkIfParticlesInteract(int i, int j) {
  return checkIfInteraction(get_ia_param(i, j));
}

int virtual_set_params(int bond_type);

#ifdef DIPOLES
void set_dipolar_method_local(DipolarInteraction method);
#endif

#include "utils/math/sqr.hpp"

/** Returns true if the particles are to be considered for short range
    interactions */
class VerletCriterion {
  const double m_skin;
  const double m_eff_max_cut2;
  const double m_eff_coulomb_cut2 = 0.;
  const double m_eff_dipolar_cut2 = 0.;
  const double m_collision_cut2 = 0.;

public:
  VerletCriterion(double skin, double max_cut, double coulomb_cut = 0.,
                  double dipolar_cut = 0.,
                  double collision_detection_cutoff = 0.)
      : m_skin(skin), m_eff_max_cut2(Utils::sqr(max_cut + m_skin)),
        m_eff_coulomb_cut2(Utils::sqr(coulomb_cut + m_skin)),
        m_eff_dipolar_cut2(Utils::sqr(dipolar_cut + m_skin)),
        m_collision_cut2(Utils::sqr(collision_detection_cutoff)) {}

  template <typename Distance>
  bool operator()(const Particle &p1, const Particle &p2,
                  Distance const &dist) const {
    auto const &dist2 = dist.dist2;
    if (dist2 > m_eff_max_cut2)
      return false;

// Within real space cutoff of electrostatics and both charged
#ifdef ELECTROSTATICS
    if ((dist2 <= m_eff_coulomb_cut2) && (p1.p.q != 0) && (p2.p.q != 0))
      return true;
#endif

// Within dipolar cutoff and both carry magnetic moments
#ifdef DIPOLES
    if ((dist2 <= m_eff_dipolar_cut2) && (p1.p.dipm != 0) && (p2.p.dipm != 0))
      return true;
#endif

// Collision detection
#ifdef COLLISION_DETECTION
    if (dist2 <= m_collision_cut2)
      return true;
#endif

    // Within short-range distance (incl dpd and the like)
    auto const max_cut = get_ia_param(p1.p.type, p2.p.type)->max_cut;
    if ((max_cut != INACTIVE_CUTOFF) && (dist2 <= Utils::sqr(max_cut + m_skin)))
      return true;

    return false;
  }
};
#endif
