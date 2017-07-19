#ifndef BOND_TYPE_H
#define BOND_TYPE_H

enum class BondType {

  FENE_BOND,
  /** Type of bonded interaction is a HARMONIC potential. */
  HARMONIC_BOND,
  /** Type of bonded interaction is a HARMONIC_DUMBBELL potential. */
  HARMONIC_DUMBBELL_BOND,
  /** Type of bonded interaction is a QUARTIC potential. */
  QUARTIC_BOND,
  /** Type of bonded interaction is a BONDED_COULOMB */
  BONDED_COULOMB_BOND,
  /** Type of bonded interaction is a bond angle potential. */
  ANGLE_OLD_BOND,
  /** Type of bonded interaction is a dihedral potential. */
  DIHEDRAL_BOND,
  /** Type of tabulated bonded interaction potential,
      may be of bond length, of bond angle or of dihedral type. */
  TABULATED_BOND,
  /** Type of bonded interaction is a (-LJ) potential. */
  SUBT_LJ_BOND,
  /** Type of a Rigid/Constrained bond*/
  RIGID_BOND,
  /** Type of a virtual bond*/
  VIRTUAL_BOND,
  /** Type of bonded interaction is a bond angle -- constraint distance
     potential. */
  ANGLEDIST_BOND,
  /** Type of bonded interaction is a bond angle -- chain ends have angle with
     wall constraint */
  ENDANGLEDIST_BOND,
  /** Type of overlapped bonded interaction potential,
      may be of bond length, of bond angle or of dihedral type. */
  OVERLAPPED_BOND,
  /** Type of bonded interaction is a bond angle cosine potential. */
  ANGLE_HARMONIC_BOND,
  /** Type of bonded interaction is a bond angle cosine potential. */
  ANGLE_COSINE_BOND,
  /** Type of bonded interaction is a bond angle cosine potential. */
  ANGLE_COSSQUARE_BOND,
  /** Type of bonded interaction: oif local forces. */
  OIF_LOCAL_FORCES_BOND,
  /** Type of bonded interaction: oif global forces. */
  OIF_GLOBAL_FORCES_BOND,
  /** Type of bonded interaction: determining outward direction of oif membrane.
     */
  OIF_OUT_DIRECTION_BOND,
  /** Type of bonded interaction for cg DNA */
  CG_DNA_BASEPAIR_BOND,
  /** Type of bonded interaction for cg DNA */
  CG_DNA_STACKING_BOND,
  /** Type of bonded interaction for cg DNA */
  CG_DNA_BACKBONE_BOND,
  /** Type of bonded interaction is a wall repulsion (immersed boundary). */
  IBM_TRIEL_BOND,
  /** Type of bonded interaction is volume conservation force (immersed
     boundary). */
  IBM_VOLUME_CONSERVATION_BOND,
  /** Type of bonded interaction is bending force (immersed boundary). */
  IBM_TRIBEND_BOND,
  /** Type of bonded interaction is umbrella. */
  UMBRELLA_BOND
};
#endif
