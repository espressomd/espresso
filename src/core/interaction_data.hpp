/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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

#include "utils.hpp"
#include "particle_data.hpp" /* needed for constraints */

/** \name Type codes of bonded interactions
    Enumeration of implemented bonded interactions.
*/
/************************************************************/
/*@{*/

enum BondedInteraction{
    /** This bonded interaction was not set. */
    BONDED_IA_NONE = -1,
    /** Type of bonded interaction is a FENE potential
        (to be combined with Lennard Jones). */
    BONDED_IA_FENE,
    /** Type of bonded interaction is a HARMONIC potential. */
    BONDED_IA_HARMONIC,
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
    /** Type of bonded interaction is a bond angle -- constraint distance potential. */
    BONDED_IA_ANGLEDIST,
    /** Type of bonded interaction is a bond angle -- chain ends have angle with wall constraint */
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
    /** Type of bonded interaction is a hyperelastic stretching force. */
    BONDED_IA_STRETCHING_FORCE,
    /** Type of bonded interaction is a local area force. */
    BONDED_IA_AREA_FORCE_LOCAL,
    /** Type of bonded interaction is a bending force. */
    BONDED_IA_BENDING_FORCE,
    /** Type of bonded interaction is a bending force. */
    BONDED_IA_VOLUME_FORCE,
    /** Type of bonded interaction is a global area force. */
    BONDED_IA_AREA_FORCE_GLOBAL,
    /** Type of bonded interaction is a linear stretching force. */
    BONDED_IA_STRETCHLIN_FORCE
};

/** Specify tabulated bonded interactions  */
enum TabulatedBondedInteraction{
    TAB_UNKNOWN = 0,
    TAB_BOND_LENGTH,
    TAB_BOND_ANGLE,
    TAB_BOND_DIHEDRAL
};

/** Specify overlapped bonded interactions  */
enum OverlappedBondedInteraction{
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
		COULOMB_NONE, //< Coulomb interaction switched off (NONE)
		COULOMB_DH, //< Coulomb method is Debye-Hueckel
		COULOMB_P3M, //< Coulomb method is P3M
		COULOMB_MMM1D, //< Coulomb method is one-dimensional MMM
		COULOMB_MMM2D, //< Coulomb method is two-dimensional MMM
		COULOMB_MAGGS, //< Coulomb method is "Maggs"
		COULOMB_ELC_P3M, //< Coulomb method is P3M plus ELC
		COULOMB_RF, //< Coulomb method is Reaction-Field
		COULOMB_INTER_RF, //< Coulomb method is Reaction-Field BUT as interaction
		COULOMB_P3M_GPU, //< Coulomb method is P3M with GPU based long range part calculation
		COULOMB_MMM1D_GPU, //< Coulomb method is one-dimensional MMM running on GPU
		COULOMB_EWALD_GPU, //< Coulomb method is Ewald running on GPU
	};

#endif
/*@}*/


#ifdef  DIPOLES
  /** \name Type codes for the type of dipolar interaction
    Enumeration of implemented methods for the magnetostatic
    interaction.
   */
  /************************************************************/
  /*@{*/
enum DipolarInteration{
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
    DIPOLAR_MDLC_DS
};
   /*@}*/
#endif 



/** \name Type codes for constraints
    Enumeration of implemented constraint types.
*/
/************************************************************/
/*@{*/
enum ConstraintApplied{
/** No constraint applied */
    CONSTRAINT_NONE = 0,
/** wall constraint applied */
    CONSTRAINT_WAL,
/** spherical constraint applied */
    CONSTRAINT_SPH,
/** (finite) cylinder shaped constraint applied */
    CONSTRAINT_CYL,
/** Rod-like charge. */
    CONSTRAINT_ROD,
/** Plate-like charge. */
    CONSTRAINT_PLATE,
/** maze-like constraint applied */
    CONSTRAINT_MAZE,
/** pore constraint applied */
    CONSTRAINT_PORE,
//ER
/** External magnetic field constraint applied */
    CONSTRAINT_EXT_MAGN_FIELD,
//end ER
/** Constraint for tunable-lsip boundary conditions */
    CONSTRAINT_PLANE,
/** Constraint for tunable-lsip boundary conditions */
    CONSTRAINT_RHOMBOID,
/** Constraint for a stomatocyte boundary */
    CONSTRAINT_STOMATOCYTE,
/** slitpore constraint applied */
    CONSTRAINT_SLITPORE,
/** Constraint for a hollow cone boundary */
    CONSTRAINT_HOLLOW_CONE
};
/*@}*/

/* Data Types */
/************************************************************/

/** field containing the interaction parameters for
 *  nonbonded interactions. Access via
 * get_ia_param(i, j), i,j < n_particle_types */
typedef struct {

  /** flag that tells whether there is any short-ranged interaction,
   i.e. one that contributes to the "nonbonded" section of the
   energy/pressure. Note that even if there is no short-ranged
   interaction present, the \ref max_cut can be non-zero due to
   e.g. electrostatics. */
  int particlesInteract;

  /** maximal cutoff for this pair of particle types. This contains
      contributions from the short-ranged interactions, plus any
      cutoffs from global interactions like electrostatics.
  */
  double max_cut;

#ifdef LENNARD_JONES
  /** \name Lennard-Jones with shift */
  /*@{*/
  double LJ_eps;
  double LJ_sig;
  double LJ_cut;
  double LJ_shift;
  double LJ_offset;
  double LJ_capradius;
  double LJ_min;
  /*@}*/
#endif

#ifdef LENNARD_JONES_GENERIC
  /** \name Generic Lennard-Jones with shift */
  /*@{*/
  double LJGEN_eps;
  double LJGEN_sig;
  double LJGEN_cut;
  double LJGEN_shift;
  double LJGEN_offset;
  double LJGEN_capradius;
  int LJGEN_a1;
  int LJGEN_a2;
  double LJGEN_b1;
  double LJGEN_b2;
#ifdef LJGEN_SOFTCORE
  double LJGEN_lambda;
  double LJGEN_softrad;
#endif
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
  /* Cap */
  double LJANGLE_capradius;
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
  int    SmSt_n;
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
  double MORSE_capradius;
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
  double BUCK_capradius;
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
  double LJCOS2_capradius;
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

#ifdef COMFORCE
  /** \name center of mass directed force */
  /*@{*/
  int COMFORCE_flag;
  int COMFORCE_dir;
  double COMFORCE_force;
  double COMFORCE_fratio;
  /*@}*/
#endif
  
#ifdef COMFIXED
  /** \name center of mass directed force */
  /*@{*/
  int COMFIXED_flag;
  /*@}*/
#endif

#ifdef INTER_DPD
  /** \name DPD as interaction */
  /*@{*/
  double dpd_gamma;
  double dpd_r_cut;
  int dpd_wf;
  double dpd_pref1;
  double dpd_pref2;
  double dpd_tgamma;
  double dpd_tr_cut;
  int dpd_twf;
  double dpd_pref3;
  double dpd_pref4;
  /*@}*/  
#endif

#ifdef INTER_RF
  int rf_on;
#endif

#ifdef MOL_CUT
  int mol_cut_type;
  double mol_cut_cutoff;
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
  DipolarInteration    Dmethod;
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


/** Parameters for hyperelastic stretching_force */
typedef struct {
  double r0;
  double ks;
} Stretching_force_bond_parameters;


/** Parameters for linear stretching_force */
typedef struct {
  double r0;
  double kslin;
} Stretchlin_force_bond_parameters;

/** Parameters for area_force_local */
typedef struct {
  double A0_l;
  double ka_l;
} Area_force_local_bond_parameters;
/** Parameters for area_force_global */
typedef struct {
  double A0_g;
  double ka_g;
} Area_force_global_bond_parameters;

/** Parameters for bending_force */
typedef struct {
   double phi0;
   double kb;
} Bending_force_bond_parameters;

/** Parameters for volume_force */
typedef struct {
  double V0;
  double kv;
} Volume_force_bond_parameters;
    
    
/** Parameters for harmonic bond Potential */
typedef struct {
      double k;
      double r;
      double r_cut;
} Harmonic_bond_parameters;

/** Parameters for quartic bond Potential */
typedef struct {
      double k0, k1;
      double r;
      double r_cut;
} Quartic_bond_parameters;

/** Parameters for coulomb bond Potential */
typedef struct {
      double prefactor;
} Bonded_coulomb_bond_parameters;

/** Parameters for three body angular potential (bond-angle potentials). 
	ATTENTION: Note that there are different implementations of the bond angle
	potential which you may chose with a compiler flag in the file \ref config.hpp !
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
      char   *filename;
      TabulatedBondedInteraction    type;
      int    npoints;
      double minval;
      double maxval;
      double invstepsize;
      double *f;
      double *e;
} Tabulated_bond_parameters;

/** Parameters for n-body overlapped potential (n=2,3,4). */
typedef struct {
      char   *filename;
      OverlappedBondedInteraction    type;
      double maxval;
      int    noverlaps;
      double *para_a;
      double *para_b;
      double *para_c;
} Overlap_bond_parameters;


/** Dummy parameters for -LJ Potential */
typedef struct {
      double k;
      double r;
      double r2;
} Subt_lj_bond_parameters;
    
    
/**Parameters for the rigid_bond/SHAKE/RATTLE ALGORITHM*/
typedef struct {
      /**Length of rigid bond/Constrained Bond*/
      //double d;
      /**Square of the length of Constrained Bond*/
      double d2;
      /**Positional Tolerance/Accuracy value for termination of RATTLE/SHAKE iterations during position corrections*/
      double p_tol;
      /**Velocity Tolerance/Accuracy for termination of RATTLE/SHAKE iterations during velocity corrections */
      double v_tol;
} Rigid_bond_parameters;


/** Parameters for three body angular potential (bond-angle potentials) that 
    depends on distance to wall constraint.
	ATTENTION: Note that there are different implementations of the bond angle
	potential which you may chose with a compiler flag in the file \ref config.hpp !
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

/** Union in which to store the parameters of an individual bonded interaction */
typedef union {
    Fene_bond_parameters fene;
    Stretchlin_force_bond_parameters stretchlin_force;
    Stretching_force_bond_parameters stretching_force;
    Area_force_local_bond_parameters area_force_local;
    Area_force_global_bond_parameters area_force_global;
    Bending_force_bond_parameters bending_force;
    Volume_force_bond_parameters volume_force;
    Harmonic_bond_parameters harmonic;
    Quartic_bond_parameters quartic;
    Bonded_coulomb_bond_parameters bonded_coulomb;
    Angle_bond_parameters angle;
    Angle_harmonic_bond_parameters angle_harmonic;
    Angle_cosine_bond_parameters angle_cosine;
    Angle_cossquare_bond_parameters angle_cossquare;
    Dihedral_bond_parameters dihedral;
    Tabulated_bond_parameters tab;
    Overlap_bond_parameters overlap;
    Subt_lj_bond_parameters subt_lj;
    Rigid_bond_parameters rigid_bond;
    Angledist_bond_parameters angledist;
    Endangledist_bond_parameters endangledist;
    
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

#ifdef CONSTRAINTS
/** \name Compounds for constraints */
/*@{*/

/** Parameters for a WALL constraint (or a plane if you like that more). */
typedef struct {
  /** normal vector on the plane. */
  double n[3];
  /** distance of the wall from the origin. */
  double d;
  /** whether the constraint is penetrable 1 or not 0*/
  int penetrable; 
  int reflecting;
  int only_positive;
  /** whether to calculate tunable slip forces 1 or not 0 */
  int tunable_slip;
} Constraint_wall;

/** Parameters for a SPHERE constraint. */
typedef struct {
  /** sphere center. */
  double pos[3];
  /** sphere radius. */
  double rad;  
  /** sphere direction. (+1 outside -1 inside interaction direction)*/
  double direction;
  /** whether the constraint is penetrable 1 or not 0*/
  int penetrable; 
  int reflecting;
} Constraint_sphere;

/** Parameters for a CYLINDER constraint. */
typedef struct {
  /** center of the cylinder. */
  double pos[3];
  /** Axis of the cylinder .*/
  double axis[3];
  /** cylinder radius. */
  double rad;
  /** cylinder length. (!!!NOTE this is only the half length of the cylinder.)*/
  double length;
  /** cylinder direction. (+1 outside -1 inside interaction direction)*/
  double direction;
  /** whether the constraint is penetrable 1 or not 0*/
  int penetrable; 
  int reflecting;
} Constraint_cylinder;

/** Parameters for a RHOMBOID constraint. */
typedef struct {
  /** corner of the rhomboid */
  double pos[3];
  /** edges adjacent to the corner */
  double a[3];
  double b[3];
  double c[3];
  /** rhomboid direction. (+1 outside -1 inside interaction direction)*/
  double direction;
  /** whether the constraint is penetrable 1 or not 0*/
  int penetrable; 
  int reflecting;
} Constraint_rhomboid;

/** Parameters for a PORE constraint. */
typedef struct {
  /** center of the cylinder. */
  double pos[3];
  /** Axis of the cylinder .*/
  double axis[3];
  /** cylinder radius. */
  double rad_left;
  double rad_right;
  double smoothing_radius;
  /** cylinder length. (!!!NOTE this is only the half length of the cylinder.)*/
  double length;
  int reflecting;
  double outer_rad_left;
  double outer_rad_right;
} Constraint_pore;


/** Parameters for a SLITPORE constraint. */
typedef struct {
  /** center of the cylinder. */
  double pore_mouth;
  /** Axis of the cylinder .*/
  double upper_smoothing_radius;
  double lower_smoothing_radius;
  /** cylinder length. (!!!NOTE this is only the half length of the cylinder.)*/
  double channel_width;
  double pore_width;
  double pore_length;
  int reflecting;
} Constraint_slitpore;

/** Parameters for a ROD constraint. */
typedef struct {
  /** center of the cylinder in the x-y plane. */
  double pos[2];
  /** line charge density. Only makes sense if the axis along the rod is
      periodically replicated and only with MMM1D. */
  double lambda;
} Constraint_rod;

/** Parameters for a PLATE constraint. */
typedef struct {
  /** height of plane in z-axis. */
  double pos;
  /** charge density. Only makes sense if the axis along the rod is
      periodically replicated and only with MMM2D. */
  double sigma;
} Constraint_plate;

/** Parameters for a MAZE constraint. */
typedef struct {
  /** number of spheres. */
  double nsphere;
  /** dimension of the maze. */
  double dim;
  /** sphere radius. */
  double sphrad;
  /** cylinder (connecting the spheres) radius*/
  double cylrad;
  /** whether the constraint is penetrable 1 or not 0*/
  int penetrable; 
} Constraint_maze;

/** Parameters for a STOMATOCYTE constraint. */
typedef struct {

  /** Stomatocyte position. */

  double position_x;
  double position_y;
  double position_z;

  /** Stomatocyte orientation. */

  double orientation_x;
  double orientation_y;
  double orientation_z;

  /** Stomatocyte dimensions. */

  double outer_radius;
  double inner_radius;
  double layer_width;

  /** Inside/Outside (+1 outside -1 inside interaction direction)*/

  double direction;

  /** whether the constraint is penetrable 1 or not 0*/

  int penetrable; 
  int reflecting;

} Constraint_stomatocyte;

/** Parameters for a HOLLOW_CONE constraint. */
typedef struct {

  /** Hollow cone position. */

  double position_x;
  double position_y;
  double position_z;

  /** Hollow cone orientation. */

  double orientation_x;
  double orientation_y;
  double orientation_z;

  /** Hollow cone dimensions. */

  double inner_radius;
  double outer_radius;
  double width;
  double opening_angle;

  /** Inside/Outside (+1 outside -1 inside interaction direction)*/

  double direction;

  /** whether the constraint is penetrable 1 or not 0*/

  int penetrable; 
  int reflecting;

} Constraint_hollow_cone;

/** Parameters for a BOX constraint. */
typedef struct {
  int value;
} Constraint_box;

//ER
/** Parameters for a EXTERNAL MAGNETIC FIELD constraint */
typedef struct{
  /** vector (direction and magnitude) of the external magnetic field */
  double ext_magn_field[3];
} Constraint_ext_magn_field;
//end ER

/** Parameters for a plane constraint which is needed for tunable-slip boundary conditions. */
typedef struct {
  /** Position of the plain. Negative values mean non-existing in that direction. */
  double pos[3];
} Constraint_plane;

typedef struct {
  double omega;
  double Prefactor;
} SinusoidalField;

/** Structure to specify a constraint. */
typedef struct {
  /** type of the constraint. */
  ConstraintApplied type;

  union {
    Constraint_wall wal;
    Constraint_sphere sph;
    Constraint_cylinder cyl;
    Constraint_rhomboid rhomboid;
    Constraint_rod rod;
    Constraint_plate plate;
    Constraint_maze maze;
    Constraint_pore pore;
    Constraint_slitpore slitpore;
    Constraint_stomatocyte stomatocyte;
    Constraint_hollow_cone hollow_cone;
    //ER
    Constraint_ext_magn_field emfield;
    //end ER
    Constraint_plane plane;
  } c;

  /** particle representation of this constraint. Actually needed are only the identity,
      the type and the force. */
  Particle part_rep;
} Constraint;
/*@}*/
#endif


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
/** Maximal interaction cutoff (real space/short range non-bonded interactions). */
extern double max_cut_nonbonded;
/** Maximal interaction cutoff (real space/short range bonded interactions). */
extern double max_cut_bonded;
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
  return &ia_params[i*n_particle_types + j];
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
    value is used in the verlet pair list algorithm (see \ref
    verlet.hpp). */
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

#endif
