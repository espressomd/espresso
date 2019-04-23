#ifndef _BONDED_INTERACTION_DATA_HPP
#define _BONDED_INTERACTION_DATA_HPP

#include "TabulatedPotential.hpp"
#include "particle_data.hpp"
/** \name Type codes of bonded interactions
 *  Enumeration of implemented bonded interactions.
 */
/************************************************************/
/*@{*/

enum BondedInteraction {
  /** This bonded interaction was not set. */
  BONDED_IA_NONE = -1,
  /** Type of bonded interaction is a FENE potential
      (to be combined with Lennard-Jones). */
  BONDED_IA_FENE,
  /** Type of bonded interaction is a HARMONIC potential. */
  BONDED_IA_HARMONIC,
  /** Type of bonded interaction is a HARMONIC_DUMBBELL potential. */
  BONDED_IA_HARMONIC_DUMBBELL,
  /** Type of bonded interaction is a QUARTIC potential. */
  BONDED_IA_QUARTIC,
  /** Type of bonded interaction is a BONDED_COULOMB */
  BONDED_IA_BONDED_COULOMB,
  /** Type of bonded interaction is a BONDED_COULOMB_SR */
  BONDED_IA_BONDED_COULOMB_SR,
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
};

/** Specify tabulated bonded interactions  */
enum TabulatedBondedInteraction {
  TAB_UNKNOWN = 0,
  TAB_BOND_LENGTH = 1,
  TAB_BOND_ANGLE = 2,
  TAB_BOND_DIHEDRAL = 3
};

/*@}*/
/** Parameters for FENE bond Potential. */
struct Fene_bond_parameters {
  /** spring constant */
  double k;
  /** maximal bond stretching */
  double drmax;
  /** equilibrium bond length */
  double r0;
  /** square of @p drmax (internal parameter) */
  double drmax2;
  /** inverse square of @p drmax (internal parameter) */
  double drmax2i;
};

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
  /** spring constant */
  double k;
  /** equilibrium bond length */
  double r;
  /** cutoff bond length */
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
  /** spring constant */
  double k1;
  /** rotation constant */
  double k2;
  /** equilibrium bond length */
  double r;
  /** cutoff bond length */
  double r_cut;
};
#endif

/** Parameters for quartic bond Potential */
struct Quartic_bond_parameters {
  double k0, k1;
  double r;
  double r_cut;
};

/** Parameters for Coulomb bond Potential */
struct Bonded_coulomb_bond_parameters {
  /** Coulomb prefactor */
  double prefactor;
};

/** Parameters for Coulomb bond short-range Potential */
struct Bonded_coulomb_sr_bond_parameters {
  /** charge factor */
  double q1q2;
};

/** Parameters for three-body angular potential.
 *  @note
 *      ATTENTION: there are different implementations of the bond angle
 *      potential which you may choose with a compiler flag in the file
 *      \ref config.hpp !
 */
struct Angle_bond_parameters {
  /** bending constant */
  double bend;
  /** equilibrium angle (default is 180 degrees) */
  double phi0;
  /** cosine of @p phi0 (internal parameter) */
  double cos_phi0;
  /** sine of @p phi0 (internal parameter) */
  double sin_phi0;
};

/** Parameters for three-body angular potential (harmonic). */
struct Angle_harmonic_bond_parameters {
  /** bending constant */
  double bend;
  /** equilibrium angle (default is 180 degrees) */
  double phi0;
};

/** Parameters for three-body angular potential (cosine). */
struct Angle_cosine_bond_parameters {
  /** bending constant */
  double bend;
  /** equilibrium angle (default is 180 degrees) */
  double phi0;
  /** cosine of @p phi0 (internal parameter) */
  double cos_phi0;
  /** sine of @p phi0 (internal parameter) */
  double sin_phi0;
};

/** Parameters for three-body angular potential (cossquare). */
struct Angle_cossquare_bond_parameters {
  /** bending constant */
  double bend;
  /** equilibrium angle (default is 180 degrees) */
  double phi0;
  /** cosine of @p phi0 (internal parameter) */
  double cos_phi0;
};

/** Parameters for four-body angular potential (dihedral-angle potentials). */
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
struct Subt_lj_bond_parameters {};

/** Parameters for the rigid_bond/SHAKE/RATTLE ALGORITHM */
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
  /** ID of the large soft particle to which this node belongs */
  int softID;
  /** Reference volume */
  double volRef;
  /** Spring constant for volume force */
  double kappaV;
  // Whether to write out center-of-mass at each time step
  // Actually this is more of an analysis function and does not strictly belong
  // to volume conservation
  //  bool writeCOM;
};

/** Parameters for IBM tribend **/
struct IBM_Tribend_Parameters {
  /** Interaction data */
  double kb;

  /** Reference angle */
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
  Bonded_coulomb_sr_bond_parameters bonded_coulomb_sr;
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
  Subt_lj_bond_parameters subt_lj;
  Rigid_bond_parameters rigid_bond;
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

/** Field containing the parameters of the bonded ia types */
extern std::vector<Bonded_ia_parameters> bonded_ia_params;

/** Maximal interaction cutoff (real space/short range bonded interactions). */
extern double max_cut_bonded;

/** Makes sure that \ref bonded_ia_params is large enough to cover the
    parameters for the bonded interaction type. Attention: 1: There is
    no initialization done here. 2: Use only in connection with
    creating new or overwriting old bond types*/
void make_bond_type_exist(int type);

/** @brief Checks if particle has a pair bond with a given partner
 *  Note that bonds are stored only on one of the two particles in Espresso
 *
 * @param p          particle on which the bond may be stored
 * @param partner    bond partner
 * @param bond_type  numerical bond type */
inline bool pair_bond_exists_on(const Particle *const p,
                                const Particle *const partner, int bond_type) {
  // First check the bonds of p1
  if (p->bl.e) {
    int i = 0;
    while (i < p->bl.n) {
      int size = bonded_ia_params[p->bl.e[i]].num;

      if (p->bl.e[i] == bond_type && p->bl.e[i + 1] == partner->p.identity) {
        // There's a bond, already. Nothing to do for these particles
        return true;
      }
      i += size + 1;
    }
  }
  return false;
}

/** @brief Checks both particle for a specific bond. Needs GHOSTS_HAVE_BONDS if
 *  particles are ghosts.
 *
 *  @param p_bond      particle on which the bond may be stored
 *  @param p_partner   bond partner
 *  @param bond        enum bond type
 */
inline bool pair_bond_enum_exists_on(const Particle *const p_bond,
                                     const Particle *const p_partner,
                                     BondedInteraction bond) {
#ifdef ADDITIONAL_CHECKS
  extern int ghosts_have_bonds;
  assert(ghosts_have_bonds);
#endif

  int i = 0;
  while (i < p_bond->bl.n) {
    int type_num = p_bond->bl.e[i];
    Bonded_ia_parameters *iaparams = &bonded_ia_params[type_num];
    if (iaparams->type == (int)bond &&
        p_bond->bl.e[i + 1] == p_partner->p.identity) {
      return true;
    }
    i += iaparams->num + 1;
  }
  return false;
}

/** @brief Checks both particle for a specific bond. Needs GHOSTS_HAVE_BONDS if
 *  particles are ghosts.
 *
 *  @param p1     particle on which the bond may be stored
 *  @param p2     particle on which the bond may be stored
 *  @param bond   numerical bond type
 */
inline bool pair_bond_enum_exists_between(const Particle *const p1,
                                          const Particle *const p2,
                                          BondedInteraction bond) {
  if (p1 == p2)
    return false;

  // Check if particles have bonds (bl.n > 0) and search for the bond of
  // interest with are_bonded(). Could be saved on both sides (and both could
  // have other bonds), so we need to check both.
  return (p1->bl.n > 0 && pair_bond_enum_exists_on(p1, p2, bond)) ||
         (p2->bl.n > 0 && pair_bond_enum_exists_on(p2, p1, bond));
}

void recalc_maximal_cutoff_bonded();
int virtual_set_params(int bond_type);
#endif
