/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file interaction_data.hpp
    Various procedures concerning interactions between particles.
*/

#include "particle_data.hpp"
#include "utils.hpp"

#include "TabulatedPotential.hpp"


/** \name Type codes of bonded interactions
    Enumeration of implemented bonded interactions.
*/
/************************************************************/
/*@{*/

enum BondedInteraction {
  /** This bonded interaction was not set. */
  BONDED_IA_NONE = -1,
  /** Type of bonded interaction is a FENE potential
      (to be combined with Lennard Jones). */
  BONDED_IA_FENE,
  /** Type of bonded interaction is a HARMONIC potential. */
  BONDED_IA_HARMONIC,
  /** Type of bonded interaction is a HARMONIC_DUMBBELL potential. */
  BONDED_IA_HARMONIC_DUMBBELL,
  /** Type of bonded interaction is a QUARTIC potential. */
  BONDED_IA_QUARTIC,
  /** Type of bonded interaction is a BONDED_COULOMB */
  BONDED_IA_BONDED_COULOMB,
  /** Type of bonded interaction is a bond angle potential. */
  BONDED_IA_ANGLE_OLD,
  /** Type of bonded interaction is a dihedral potential. */
  BONDED_IA_DIHEDRAL,
  /** Type of tabulated bonded interaction potential,
      may be of bond length, of bond angle or of dihedral type. */
  BONDED_IA_TABULATED,
  /** Type of bonded interaction is a (-LJ) potential. */
  BONDED_IA_SUBT_LJ,
  /** Type of a Rigid/Constrained bond*/
  BONDED_IA_RIGID_BOND,
  /** Type of a virtual bond*/
  BONDED_IA_VIRTUAL_BOND,
  /** Type of bonded interaction is a bond angle -- constraint distance
     potential. */
  BONDED_IA_ANGLEDIST,
  /** Type of bonded interaction is a bond angle cosine potential. */
  BONDED_IA_ANGLE_HARMONIC,
  /** Type of bonded interaction is a bond angle cosine potential. */
  BONDED_IA_ANGLE_COSINE,
  /** Type of bonded interaction is a bond angle cosine potential. */
  BONDED_IA_ANGLE_COSSQUARE,
  /** Type of bonded interaction: oif local forces. */
  BONDED_IA_OIF_LOCAL_FORCES,
  /** Type of bonded interaction: oif global forces. */
  BONDED_IA_OIF_GLOBAL_FORCES,
  /** Type of bonded interaction: determining outward direction of oif membrane.
     */
  BONDED_IA_OIF_OUT_DIRECTION,
  /** Type of bonded interaction for cg DNA */
  BONDED_IA_CG_DNA_BASEPAIR,
  /** Type of bonded interaction for cg DNA */
  BONDED_IA_CG_DNA_STACKING,
  /** Type of bonded interaction for cg DNA */
  BONDED_IA_CG_DNA_BACKBONE,
  /** Type of bonded interaction is a wall repulsion (immersed boundary). */
  BONDED_IA_IBM_TRIEL,
  /** Type of bonded interaction is volume conservation force (immersed
     boundary). */
  BONDED_IA_IBM_VOLUME_CONSERVATION,
  /** Type of bonded interaction is bending force (immersed boundary). */
  BONDED_IA_IBM_TRIBEND,
  /** Type of bonded interaction is umbrella. */
  BONDED_IA_UMBRELLA,
  /** Type of bonded interaction is thermalized distance bond. */
  BONDED_IA_THERMALIZED_DIST,
  /** Type of bonded interaction is a BONDED_COULOMB_P3M_SR */
  BONDED_IA_BONDED_COULOMB_P3M_SR,
};

/** Specify tabulated bonded interactions  */
enum TabulatedBondedInteraction {
  TAB_UNKNOWN = 0,
  TAB_BOND_LENGTH = 1,
  TAB_BOND_ANGLE = 2,
  TAB_BOND_DIHEDRAL = 3
};

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
  COULOMB_INTER_RF,  //< Coulomb method is Reaction-Field BUT as interaction
  COULOMB_P3M_GPU,   //< Coulomb method is P3M with GPU based long range part
                     // calculation
  COULOMB_MMM1D_GPU, //< Coulomb method is one-dimensional MMM running on GPU
  COULOMB_EK,        //< Coulomb method is electrokinetics
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
  /** dipolar interation switched off (NONE). */
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

#ifdef INTER_RF
  int rf_on = 0;
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
 *  the coulomb  interaction.  */
struct Coulomb_parameters {

#ifdef ELECTROSTATICS
  /** bjerrum length times temperature. */
  double prefactor;

  /** Method to treat coulomb interaction. */
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

/*@}*/
/** Parameters for FENE bond Potential.
k - spring constant.
drmax - maximal bond streching.
r0 - equilibrium bond length.
drmax2 - square of drmax (internal parameter).
*/
struct Fene_bond_parameters {
  double k;
  double drmax;
  double r0;
  double drmax2;
  double drmax2i;
};

#ifdef HYDROGEN_BOND
/** Parameters for the cg_dna potential
    Insert documentation here.
**/
struct Cg_dna_basepair_parameters {
  double r0;
  double alpha;
  double E0;
  double kd;
  double sigma1;
  double sigma2;
  double psi10;
  double psi20;
  /* Parameters for the sugar base interaction */
  double E0sb;
  double r0sb;
  double alphasb;
  double f2;
  double f3;
};
#endif
#ifdef TWIST_STACK
struct Cg_dna_stacking_parameters {
  double rm;
  double epsilon;
  double ref_pot;
  double a[8];
  double b[7];
};
#endif

/** Parameters for oif_global_forces */
struct Oif_global_forces_bond_parameters {
  double A0_g;
  double ka_g;
  double V0;
  double kv;
};

/** Parameters for oif_local_forces */
struct Oif_local_forces_bond_parameters {
    double r0;
    double ks;
    double kslin;
    double phi0;
    double kb;
    double A01;
    double A02;
    double kal;
    double kvisc;
};

/** Parameters for harmonic bond Potential */
struct Harmonic_bond_parameters {
  double k;
  double r;
  double r_cut;
};

/** Parameters for Thermalized bond **/
struct Thermalized_bond_parameters {
    double temp_com;
    double gamma_com;
    double temp_distance;
    double gamma_distance;
    double r_cut;
    double pref1_com;
    double pref2_com;
    double pref1_dist;
    double pref2_dist;
};

#ifdef ROTATION
/** Parameters for harmonic dumbbell bond Potential */
struct Harmonic_dumbbell_bond_parameters {
  double k1;
  double k2;
  double r;
  double r_cut;
};
#endif

/** Parameters for quartic bond Potential */
struct Quartic_bond_parameters {
  double k0, k1;
  double r;
  double r_cut;
};

/** Parameters for coulomb bond Potential */
struct Bonded_coulomb_bond_parameters { double prefactor; };

#ifdef P3M
/** Parameters for coulomb bond p3m shortrange Potential */
struct Bonded_coulomb_p3m_sr_bond_parameters { double q1q2; };
#endif

/** Parameters for three body angular potential (bond-angle potentials).
        ATTENTION: Note that there are different implementations of the bond
   angle
        potential which you may chose with a compiler flag in the file \ref
   config.hpp !
        bend - bending constant.
        phi0 - equilibrium angle (default is 180 degrees / Pi) */
struct Angle_bond_parameters {
  double bend;
  double phi0;
  double cos_phi0;
  double sin_phi0;

};

/** Parameters for three body angular potential (bond_angle_harmonic).
    bend - bending constant.
    phi0 - equilibrium angle (default is 180 degrees / Pi) */
struct Angle_harmonic_bond_parameters {
  double bend;
  double phi0;
};

/** Parameters for three body angular potential (bond_angle_cosine).
    bend - bending constant.
    phi0 - equilibrium angle (default is 180 degrees / Pi) */
struct Angle_cosine_bond_parameters {
  double bend;
  double phi0;
  double cos_phi0;
  double sin_phi0;
};

/** Parameters for three body angular potential (bond_angle_cossquare).
    bend - bending constant.
    phi0 - equilibrium angle (default is 180 degrees / Pi) */
struct Angle_cossquare_bond_parameters {
  double bend;
  double phi0;
  double cos_phi0;
};

/** Parameters for four body angular potential (dihedral-angle potentials). */
struct Dihedral_bond_parameters {
  double mult;
  double bend;
  double phase;
};

/** Parameters for n-body tabulated potential (n=2,3,4). */
struct Tabulated_bond_parameters {
  TabulatedBondedInteraction type;
  TabulatedPotential *pot;
};

#ifdef UMBRELLA
/** Parameters for umbrella potential */
struct Umbrella_bond_parameters {
  double k;
  int dir;
  double r;
};
#endif

/** Dummy parameters for -LJ Potential */
struct Subt_lj_bond_parameters {
};

/**Parameters for the rigid_bond/SHAKE/RATTLE ALGORITHM*/
struct Rigid_bond_parameters {
  /**Square of the length of Constrained Bond*/
  double d2;
  /**Positional Tolerance/Accuracy value for termination of RATTLE/SHAKE
   * iterations during position corrections*/
  double p_tol;
  /**Velocity Tolerance/Accuracy for termination of RATTLE/SHAKE iterations
   * during velocity corrections */
  double v_tol;
};

/** Parameters for three body angular potential (bond-angle potentials) that
    depends on distance to wall constraint.
        ATTENTION: Note that there are different implementations of the bond
   angle
        potential which you may chose with a compiler flag in the file \ref
   config.hpp !
        bend - bending constant.
        phi0 - equilibrium angle (default is 180 degrees / Pi)
        dist0 - equilibrium distance (no default) */
struct Angledist_bond_parameters {
  double bend;
  double phimin;
  double distmin;
  double phimax;
  double distmax;
  double cos_phi0;
  double sin_phi0;
};

enum class tElasticLaw { NeoHookean, Skalak };

/** Parameters for IBM elastic triangle (triel) **/
struct IBM_Triel_Parameters {
  // These values encode the reference state
  double l0;
  double lp0;
  double sinPhi0;
  double cosPhi0;
  double area0;

  // These values are cache values to speed up computation
  double a1;
  double a2;
  double b1;
  double b2;

  // These are interaction parameters
  // k1 is used for Neo-Hookean
  // k1 and k2 are used Skalak
  double maxDist;
  tElasticLaw elasticLaw;
  double k1;
  double k2;

};

/** Parameters for IBM volume conservation bond **/
struct IBM_VolCons_Parameters {
  int softID; // ID of the large soft particle to which this node belongs
  // Reference volume
  double volRef;
  // Spring constant for volume force
  double kappaV;
  // Whether to write out center-of-mass at each time step
  // Actually this is more of an analysis function and does not strictly belong
  // to volume conservation
  //  bool writeCOM;
};

/** Parameters for IBM tribend **/
struct IBM_Tribend_Parameters {
  // Interaction data
  double kb;

  // Reference angle
  double theta0;

};

/** Union in which to store the parameters of an individual bonded interaction
 */
union Bond_parameters {
  Fene_bond_parameters fene;
  Oif_global_forces_bond_parameters oif_global_forces;
  Oif_local_forces_bond_parameters oif_local_forces;
  Harmonic_bond_parameters harmonic;
#ifdef ROTATION
  Harmonic_dumbbell_bond_parameters harmonic_dumbbell;
#endif
  Quartic_bond_parameters quartic;
  Bonded_coulomb_bond_parameters bonded_coulomb;
  Angle_bond_parameters angle;
  Angle_harmonic_bond_parameters angle_harmonic;
  Angle_cosine_bond_parameters angle_cosine;
  Angle_cossquare_bond_parameters angle_cossquare;
  Dihedral_bond_parameters dihedral;
#ifdef TABULATED
  Tabulated_bond_parameters tab;
#endif
#ifdef UMBRELLA
  Umbrella_bond_parameters umbrella;
#endif
  Thermalized_bond_parameters thermalized_bond;
#ifdef P3M
  Bonded_coulomb_p3m_sr_bond_parameters bonded_coulomb_p3m_sr;
#endif
  Subt_lj_bond_parameters subt_lj;
  Rigid_bond_parameters rigid_bond;
  Angledist_bond_parameters angledist;
#if defined(CG_DNA) || defined(HYDROGEN_BOND)
  Cg_dna_basepair_parameters hydrogen_bond;
#endif
#if defined(CG_DNA) || defined(TWIST_STACK)
  Cg_dna_stacking_parameters twist_stack;
#endif
  IBM_Triel_Parameters ibm_triel;
  IBM_VolCons_Parameters ibmVolConsParameters;
  IBM_Tribend_Parameters ibm_tribend;
};

/** Defines parameters for a bonded interaction. */
struct Bonded_ia_parameters {
  /** bonded interaction type. See \ref BONDED_IA_FENE "Type code for bonded" */
  BondedInteraction type;
  /** (Number of particles - 1) interacting for that type */
  int num;
  /** union to store the different bonded interaction parameters. */
  Bond_parameters p;
};

/************************************************
 * exported variables
 ************************************************/

/** Maximal particle type seen so far. */
extern int max_seen_particle_type;

/** Structure containing the coulomb parameters. */
extern Coulomb_parameters coulomb;

/** Field containing the paramters of the bonded ia types */
extern std::vector<Bonded_ia_parameters> bonded_ia_params;

/** Maximal interaction cutoff (real space/short range interactions). */
extern double max_cut;
/** Maximal interaction cutoff (real space/short range non-bonded interactions).
 */
extern double max_cut_nonbonded;
/** Maximal interaction cutoff (real space/short range bonded interactions). */
extern double max_cut_bonded;
/** Cutoff of coulomb real space part */
extern double coulomb_cutoff;
/** Cutoff of dipolar real space part */
extern double dipolar_cutoff;

/** Minimal global interaction cutoff. Particles with a distance
    smaller than this are guaranteed to be available on the same node
    (through ghosts).  */
extern double min_global_cut;

/** Switch for nonbonded interaction exclusion */
extern int ia_excl;

/************************************************
 * exported functions
 ************************************************/

#ifdef ELECTROSTATICS
/** @brief Set the electrostatics prefactor */
int coulomb_set_prefactor(double prefactor);


/** @brief Deactivates the current Coulomb mhthod 
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
void ia_params_set_state(std::string const&);

bool is_new_particle_type(int type);
/** Makes sure that ia_params is large enough to cover interactions
    for this particle type. The interactions are initialized with values
    such that no physical interaction occurs. */
void make_particle_type_exist(int type);

void make_particle_type_exist_local(int type);

/** Makes sure that \ref bonded_ia_params is large enough to cover the
    parameters for the bonded interaction type. Attention: 1: There is
    no initialization done here. 2: Use only in connection with
    creating new or overwriting old bond types*/
void make_bond_type_exist(int type);

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
    value is used in the verlet pair list algorithm. */
void recalc_maximal_cutoff();

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

/** @brief Checks if particle has a pair bond with a given partner  
*  Note that bonds are stored only on one of the two particles in Espresso
* 
* @param P
* @param p          particle on which the bond may be stored
* @param partner    bond partner 
* @param bond_type  numerical bond type */ 
inline bool pair_bond_exists_on(const Particle* const p, const Particle* const partner, int bond_type)
{
  // First check the bonds of p1
  if (p->bl.e) {
    int i = 0;
    while(i < p->bl.n) {
      int size = bonded_ia_params[p->bl.e[i]].num;
      
      if (p->bl.e[i] == bond_type &&
          p->bl.e[i + 1] == partner->p.identity) {
        // There's a bond, already. Nothing to do for these particles
        return true;
      }
      i += size + 1;
    }
  }
  return false;
}

/** @brief Checks both particle for a specific bond. Needs GHOSTS_HAVE_BONDS if particles are ghosts.  
* 
* @param P
* @param p1          particle on which the bond may be stored
* @param p2    	     bond partner
* @param bond        enum bond type */ 
inline bool pair_bond_enum_exists_on(const Particle * const p_bond, const Particle * const p_partner, BondedInteraction bond)
{
    int i = 0;
    while (i < p_bond->bl.n) {
        int type_num = p_bond->bl.e[i];
        Bonded_ia_parameters *iaparams = &bonded_ia_params[type_num];
        if (iaparams->type == (int)bond && p_bond->bl.e[i+1] == p_partner->p.identity) {
            return true;
        } else {
            i+= iaparams->num + 1;
        }
    }
    return false;
}

/** @brief Checks both particle for a specific bond. Needs GHOSTS_HAVE_BONDS if particles are ghosts.  
* 
* @param P
* @param p1          particle on which the bond may be stored
* @param p2    	     particle on which the bond may be stored
* @param bond_type   numerical bond type */ 
inline bool pair_bond_enum_exists_between(const Particle * const p1, const Particle * const p2, BondedInteraction bond)
{
    if (p1==p2)
        return false;
    else {
        //Check if particles have bonds (bl.n > 0) and search for the bond of interest with are_bonded().
        //Could be saved on both sides (and both could have other bonds), so we need to check both.
        return (p1->bl.n > 0 && pair_bond_enum_exists_on(p1, p2, bond)) || (p2->bl.n > 0 && pair_bond_enum_exists_on(p2, p1, bond)); 
    }
}

#include "utils/math/sqr.hpp"

/** Returns true if the particles are to be considered for short range
    interactions */
class VerletCriterion {
  const double m_skin;
  const double m_eff_max_cut2;
  const double m_eff_coulomb_cut2 = 0.;
  const double m_eff_dipolar_cut2 = 0.;
  const double m_collision_cut2 =0.;

public:
  VerletCriterion(double skin, double max_cut, double coulomb_cut = 0.,
                  double dipolar_cut = 0., double collision_detection_cutoff=0.)
      : m_skin(skin), m_eff_max_cut2(Utils::sqr(max_cut + m_skin)),
        m_eff_coulomb_cut2(Utils::sqr(coulomb_cut + m_skin)),
        m_eff_dipolar_cut2(Utils::sqr(dipolar_cut + m_skin)), 
        m_collision_cut2(Utils::sqr(collision_detection_cutoff))
        {}

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

// Within dipolar cutoff and both cary magnetic moments
#ifdef DIPOLES
    if ((dist2 <= m_eff_dipolar_cut2) && (p1.p.dipm != 0) && (p2.p.dipm != 0))
      return true;
#endif


// Collision detectoin
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
