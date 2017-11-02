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
  /** Type of bonded interaction is a bond angle -- chain ends have angle with
     wall constraint */
  BONDED_IA_ENDANGLEDIST,
  /** Type of overlapped bonded interaction potential,
      may be of bond length, of bond angle or of dihedral type. */
  BONDED_IA_OVERLAPPED,
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
  BONDED_IA_UMBRELLA
};

/** Specify tabulated bonded interactions  */
enum TabulatedBondedInteraction {
  TAB_UNKNOWN = 0,
  TAB_BOND_LENGTH = 1,
  TAB_BOND_ANGLE = 2,
  TAB_BOND_DIHEDRAL = 3
};

/** Specify overlapped bonded interactions  */
enum OverlappedBondedInteraction {
  OVERLAP_UNKNOWN = 0,
  OVERLAP_BOND_LENGTH,
  OVERLAP_BOND_ANGLE,
  OVERLAP_BOND_DIHEDRAL
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
  COULOMB_EWALD_GPU, //< Coulomb method is Ewald running on GPU
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
  /** Scafacos library */
  DIPOLAR_SCAFACOS

};
#endif

/* Data Types */
/************************************************************/

/** field containing the interaction parameters for
 *  nonbonded interactions. Access via
 * get_ia_param(i, j), i,j < n_particle_types */
typedef struct {
  /** maximal cutoff for this pair of particle types. This contains
      contributions from the short-ranged interactions, plus any
      cutoffs from global interactions like electrostatics.
  */
  double max_cut;

  /** \name Lennard-Jones with shift */
  /*@{*/
  double LJ_eps;
  double LJ_sig;
  double LJ_cut;
  double LJ_shift;
  double LJ_offset;
  double LJ_min;
  /*@}*/

  /** flag that tells whether there is any short-ranged interaction,
      i.e. one that contributes to the "nonbonded" section of the
      energy/pressure. Note that even if there is no short-ranged
      interaction present, the \ref max_cut can be non-zero due to
      e.g. electrostatics. */
  int particlesInteract;

#ifdef LENNARD_JONES_GENERIC
  /** \name Generic Lennard-Jones with shift */
  /*@{*/
  double LJGEN_eps;
  double LJGEN_sig;
  double LJGEN_cut;
  double LJGEN_shift;
  double LJGEN_offset;
  double LJGEN_a1;
  double LJGEN_a2;
  double LJGEN_b1;
  double LJGEN_b2;
  double LJGEN_lambda;
  double LJGEN_softrad;
/*@}*/
#endif

#ifdef LJ_ANGLE
  /** \name Directional Lennard-Jones */
  /*@{*/
  double LJANGLE_eps;
  double LJANGLE_sig;
  double LJANGLE_cut;
  /* Locate bonded partners */
  int LJANGLE_bonded1type;
  int LJANGLE_bonded1pos;
  int LJANGLE_bonded1neg;
  int LJANGLE_bonded2pos;
  int LJANGLE_bonded2neg;
    /* Optional 2nd environment */
  double LJANGLE_z0;
  double LJANGLE_dz;
  double LJANGLE_kappa;
  double LJANGLE_epsprime;
/*@}*/
#endif

#ifdef SMOOTH_STEP
  /** \name smooth step potential */
  /*@{*/
  double SmSt_eps;
  double SmSt_sig;
  double SmSt_cut;
  double SmSt_d;
  int SmSt_n;
  double SmSt_k0;
/*@}*/
#endif

#ifdef HERTZIAN
  /** \name Hertzian potential */
  /*@{*/
  double Hertzian_eps;
  double Hertzian_sig;
/*@}*/
#endif

#ifdef GAUSSIAN
  /** \name Gaussian potential */
  /*@{*/
  double Gaussian_eps;
  double Gaussian_sig;
  double Gaussian_cut;
/*@}*/
#endif

#ifdef BMHTF_NACL
  /** \name BMHTF NaCl potential */
  /*@{*/
  double BMHTF_A;
  double BMHTF_B;
  double BMHTF_C;
  double BMHTF_D;
  double BMHTF_sig;
  double BMHTF_cut;
  double BMHTF_computed_shift;
/*@}*/
#endif

#ifdef MORSE
  /** \name Morse potential */
  /*@{*/
  double MORSE_eps;
  double MORSE_alpha;
  double MORSE_rmin;
  double MORSE_cut;
  double MORSE_rest;
/*@}*/
#endif

#ifdef BUCKINGHAM
  /** \name Buckingham potential */
  /*@{*/
  double BUCK_A;
  double BUCK_B;
  double BUCK_C;
  double BUCK_D;
  double BUCK_cut;
  double BUCK_discont;
  double BUCK_shift;
  double BUCK_F1;
  double BUCK_F2;
/*@}*/
#endif

#ifdef SOFT_SPHERE
  /** \name soft-sphere potential */
  /*@{*/
  double soft_a;
  double soft_n;
  double soft_cut;
  double soft_offset;
/*@}*/
#endif

#ifdef AFFINITY
  /** \name affinity potential */
  /*@{*/
  int affinity_type;
  double affinity_kappa;
  double affinity_r0;
  double affinity_Kon;
  double affinity_Koff;
  double affinity_maxBond;
  double affinity_cut;
/*@}*/
#endif

#ifdef MEMBRANE_COLLISION
  /** \name membrane collision potential */
  /*@{*/
  double membrane_a;
  double membrane_n;
  double membrane_cut;
  double membrane_offset;
/*@}*/
#endif

#ifdef HAT
  /** \name hat potential */
  /*@{*/
  double HAT_Fmax;
  double HAT_r;
/*@}*/
#endif

#ifdef LJCOS
  /** \name Lennard-Jones+Cos potential */
  /*@{*/
  double LJCOS_eps;
  double LJCOS_sig;
  double LJCOS_cut;
  double LJCOS_offset;
  double LJCOS_alfa;
  double LJCOS_beta;
  double LJCOS_rmin;
/*@}*/
#endif

#ifdef LJCOS2
  /** \name Lennard-Jones with a different Cos potential */
  /*@{*/
  double LJCOS2_eps;
  double LJCOS2_sig;
  double LJCOS2_cut;
  double LJCOS2_offset;
  double LJCOS2_w;
  double LJCOS2_rchange;
/*@}*/
#endif

#ifdef COS2
  /** \name Cos2 potential */
  /*@{*/
  double COS2_eps;
  double COS2_cut;
  double COS2_offset;
  double COS2_w;
/*@}*/
#endif

#ifdef GAY_BERNE
  /** \name Gay-Berne potential */
  /*@{*/
  double GB_eps;
  double GB_sig;
  double GB_cut;
  double GB_k1;
  double GB_k2;
  double GB_mu;
  double GB_nu;
  double GB_chi1;
  double GB_chi2;
/*@}*/
#endif

#ifdef TABULATED
  /** \name Tabulated potential */
  /*@{*/
  int TAB_npoints;
  int TAB_startindex;
  double TAB_minval;
  double TAB_minval2;
  double TAB_maxval;
  double TAB_stepsize;
/** The maximum allowable filename length for a tabulated potential file*/
#define MAXLENGTH_TABFILE_NAME 256
  char TAB_filename[MAXLENGTH_TABFILE_NAME];
/*@}*/
#endif

#ifdef DPD
  /** \name DPD as interaction */
  /*@{*/
  int dpd_wf;
  int dpd_twf;
  double dpd_gamma;
  double dpd_r_cut;
  double dpd_pref1;
  double dpd_pref2;
  double dpd_tgamma;
  double dpd_tr_cut;
  double dpd_pref3;
  double dpd_pref4;
/*@}*/
#endif

#ifdef INTER_RF
  int rf_on;
#endif

#ifdef TUNABLE_SLIP
  double TUNABLE_SLIP_temp;
  double TUNABLE_SLIP_gamma;
  double TUNABLE_SLIP_r_cut;
  double TUNABLE_SLIP_time;
  double TUNABLE_SLIP_vx;
  double TUNABLE_SLIP_vy;
  double TUNABLE_SLIP_vz;
#endif

#ifdef CATALYTIC_REACTIONS
  double REACTION_range;
#endif

#ifdef SHANCHEN
  double affinity[LB_COMPONENTS];
  int affinity_on;
#endif

} IA_parameters;

/** thermodynamic force parameters */

/** \name Compounds for Coulomb interactions */
/*@{*/

/** field containing the interaction parameters for
 *  the coulomb  interaction.  */
typedef struct {

#ifdef ELECTROSTATICS
  /** Bjerrum length. */
  double bjerrum;
  /** bjerrum length times temperature. */
  double prefactor;

  /** Method to treat coulomb interaction. */
  CoulombMethod method;
#endif

#ifdef DIPOLES
  double Dbjerrum;
  double Dprefactor;
  DipolarInteraction Dmethod;
#endif

} Coulomb_parameters;

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
typedef struct {
  double k;
  double drmax;
  double r0;
  double drmax2;
  double drmax2i;
} Fene_bond_parameters;

#ifdef HYDROGEN_BOND
/** Parameters for the cg_dna potential
    Insert documentation here.
**/
typedef struct {
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
} Cg_dna_basepair_parameters;
#endif
#ifdef TWIST_STACK
typedef struct {
  double rm;
  double epsilon;
  double ref_pot;
  double a[8];
  double b[7];
} Cg_dna_stacking_parameters;
#endif

/** Parameters for oif_global_forces */
typedef struct {
  double A0_g;
  double ka_g;
  double V0;
  double kv;
} Oif_global_forces_bond_parameters;

/** Parameters for oif_local_forces */
typedef struct {
  double r0;
  double ks;
  double kslin;
  double phi0;
  double kb;
  double A01;
  double A02;
  double kal;
} Oif_local_forces_bond_parameters;

/** Parameters for oif_out_direction */
typedef struct {

} Oif_out_direction_bond_parameters;

/** Parameters for harmonic bond Potential */
typedef struct {
  double k;
  double r;
  double r_cut;
} Harmonic_bond_parameters;

#ifdef ROTATION
/** Parameters for harmonic dumbbell bond Potential */
typedef struct {
  double k1;
  double k2;
  double r;
  double r_cut;
} Harmonic_dumbbell_bond_parameters;
#endif

/** Parameters for quartic bond Potential */
typedef struct {
  double k0, k1;
  double r;
  double r_cut;
} Quartic_bond_parameters;

/** Parameters for coulomb bond Potential */
typedef struct { double prefactor; } Bonded_coulomb_bond_parameters;

/** Parameters for three body angular potential (bond-angle potentials).
        ATTENTION: Note that there are different implementations of the bond
   angle
        potential which you may chose with a compiler flag in the file \ref
   config.hpp !
        bend - bending constant.
        phi0 - equilibrium angle (default is 180 degrees / Pi) */
typedef struct {
  double bend;
  double phi0;
  double cos_phi0;
  double sin_phi0;

} Angle_bond_parameters;

/** Parameters for three body angular potential (bond_angle_harmonic).
    bend - bending constant.
    phi0 - equilibrium angle (default is 180 degrees / Pi) */
typedef struct {
  double bend;
  double phi0;
} Angle_harmonic_bond_parameters;

/** Parameters for three body angular potential (bond_angle_cosine).
    bend - bending constant.
    phi0 - equilibrium angle (default is 180 degrees / Pi) */
typedef struct {
  double bend;
  double phi0;
  double cos_phi0;
  double sin_phi0;
} Angle_cosine_bond_parameters;

/** Parameters for three body angular potential (bond_angle_cossquare).
    bend - bending constant.
    phi0 - equilibrium angle (default is 180 degrees / Pi) */
typedef struct {
  double bend;
  double phi0;
  double cos_phi0;
} Angle_cossquare_bond_parameters;

/** Parameters for four body angular potential (dihedral-angle potentials). */
typedef struct {
  double mult;
  double bend;
  double phase;
} Dihedral_bond_parameters;

/** Parameters for n-body tabulated potential (n=2,3,4). */
typedef struct {
  char *filename;
  TabulatedBondedInteraction type;
  int npoints;
  double minval;
  double maxval;
  double invstepsize;
  double *f;
  double *e;
} Tabulated_bond_parameters;

/** Parameters for n-body overlapped potential (n=2,3,4). */
typedef struct {
  char *filename;
  OverlappedBondedInteraction type;
  double maxval;
  int noverlaps;
  double *para_a;
  double *para_b;
  double *para_c;
} Overlap_bond_parameters;

#ifdef UMBRELLA
/** Parameters for umbrella potential */
typedef struct {
  double k;
  int dir;
  double r;
} Umbrella_bond_parameters;
#endif

/** Dummy parameters for -LJ Potential */
typedef struct {
} Subt_lj_bond_parameters;

/**Parameters for the rigid_bond/SHAKE/RATTLE ALGORITHM*/
typedef struct {
  /**Length of rigid bond/Constrained Bond*/
  // double d;
  /**Square of the length of Constrained Bond*/
  double d2;
  /**Positional Tolerance/Accuracy value for termination of RATTLE/SHAKE
   * iterations during position corrections*/
  double p_tol;
  /**Velocity Tolerance/Accuracy for termination of RATTLE/SHAKE iterations
   * during velocity corrections */
  double v_tol;
} Rigid_bond_parameters;

/** Parameters for three body angular potential (bond-angle potentials) that
    depends on distance to wall constraint.
        ATTENTION: Note that there are different implementations of the bond
   angle
        potential which you may chose with a compiler flag in the file \ref
   config.hpp !
        bend - bending constant.
        phi0 - equilibrium angle (default is 180 degrees / Pi)
        dist0 - equilibrium distance (no default) */
typedef struct {
  double bend;
  double phimin;
  double distmin;
  double phimax;
  double distmax;
  double cos_phi0;
  double sin_phi0;
} Angledist_bond_parameters;

/** Parameters for chainend angular potential with wall  */
typedef struct {
  double bend;
  double phi0;
  double distmin;
  double distmax;
} Endangledist_bond_parameters;

typedef enum { NeoHookean, Skalak } tElasticLaw;

/** Parameters for IBM elastic triangle (triel) **/
typedef struct {
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
  double maxdist;
  tElasticLaw elasticLaw;
  double k1;
  double k2;

} IBM_Triel_Parameters;

/** Parameters for IBM volume conservation bond **/
typedef struct {
  int softID; // ID of the large soft particle to which this node belongs
  // Reference volume
  double volRef;
  // Spring constant for volume force
  double kappaV;
  // Whether to write out center-of-mass at each time step
  // Actually this is more of an analysis function and does not strictly belong
  // to volume conservation
  //  bool writeCOM;
} IBM_VolCons_Parameters;

typedef enum { TriangleNormals, NodeNeighbors } tBendingMethod;

/** Parameters for IBM tribend **/
typedef struct {
  // Interaction data
  double kb;
  tBendingMethod method;

  // Reference angle
  double theta0;

} IBM_Tribend_Parameters;

/** Union in which to store the parameters of an individual bonded interaction
 */
typedef union {
  Fene_bond_parameters fene;
  Oif_global_forces_bond_parameters oif_global_forces;
  Oif_local_forces_bond_parameters oif_local_forces;
  Oif_out_direction_bond_parameters oif_out_direction;
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
  Tabulated_bond_parameters tab;
  Overlap_bond_parameters overlap;
#ifdef UMBRELLA
  Umbrella_bond_parameters umbrella;
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
  Endangledist_bond_parameters endangledist;
  IBM_Triel_Parameters ibm_triel;
  IBM_VolCons_Parameters ibmVolConsParameters;
  IBM_Tribend_Parameters ibm_tribend;
} Bond_parameters;

/** Defines parameters for a bonded interaction. */
typedef struct {
  /** bonded interaction type. See \ref BONDED_IA_FENE "Type code for bonded" */
  BondedInteraction type;
  /** (Number of particles - 1) interacting for that type */
  int num;
  /** union to store the different bonded interaction parameters. */
  Bond_parameters p;
} Bonded_ia_parameters;

/************************************************
 * exported variables
 ************************************************/

/** Maximal particle type seen so far. */
extern int n_particle_types;
/* Number of nonbonded (short range) interactions. Not used so far.*/
extern int n_interaction_types;

/** Structure containing the coulomb parameters. */
extern Coulomb_parameters coulomb;

/** number of bonded interactions. Not used so far. */
extern int n_bonded_ia;
/** Field containing the paramters of the bonded ia types */
extern Bonded_ia_parameters *bonded_ia_params;

/** Array containing all tabulated forces*/
extern DoubleList tabulated_forces;
/** Array containing all tabulated energies*/
extern DoubleList tabulated_energies;

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
/** Function for initializing force and energy tables */
void force_and_energy_tables_init();

#ifdef ELECTROSTATICS
int coulomb_set_bjerrum(double bjerrum);
#endif

#ifdef DIPOLES
int dipolar_set_Dbjerrum(double bjerrum);
#endif

/** copy a set of interaction parameters. */
void copy_ia_params(IA_parameters *dst, IA_parameters *src);

/** get interaction parameters between particle sorts i and j */
inline IA_parameters *get_ia_param(int i, int j) {
  extern IA_parameters *ia_params;
  extern int n_particle_types;
  return &ia_params[i * n_particle_types + j];
}

/** get interaction parameters between particle sorts i and j.
    Slower than @ref get_ia_param, but can also be used on not
    yet present particle types*/
IA_parameters *get_ia_param_safe(int i, int j);

/** Makes sure that ia_params is large enough to cover interactions
    for this particle type. The interactions are initialized with values
    such that no physical interaction occurs. */
void make_particle_type_exist(int type);

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

/** call when the temperature changes, for Bjerrum length adjusting. */
void recalc_coulomb_prefactor();

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

///
const char *get_name_of_bonded_ia(BondedInteraction type);

#ifdef BOND_VIRTUAL
int virtual_set_params(int bond_type);
#endif

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

public:
  VerletCriterion(double skin, double max_cut, double coulomb_cut = 0.,
                  double dipolar_cut = 0.)
      : m_skin(skin), m_eff_max_cut2(Utils::sqr(max_cut + m_skin)),
        m_eff_coulomb_cut2(Utils::sqr(coulomb_cut + m_skin)),
        m_eff_dipolar_cut2(Utils::sqr(dipolar_cut + m_skin)) {}

  template<typename Distance>
  bool operator()(const Particle &p1, const Particle &p2, Distance const& dist) const {
    auto const& dist2 = dist.dist2;
    if (dist2 > m_eff_max_cut2)
      return false;

#ifdef EXCLUSIONS
    if (!do_nonbonded(&p1, &p2))
      return false;
#endif

    // Within short-range distance (incl dpd and the like)
    if (dist2 <= SQR(get_ia_param(p1.p.type, p2.p.type)->max_cut + m_skin))
      return true;

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

    return false;
  }
};
#endif
