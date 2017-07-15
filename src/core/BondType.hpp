#ifndef BOND_TYPE_H
#define BOND_TYPE_H
enum class BondType {
  /** This bonded interaction was not set. */
  NONE,
  /** Type of bonded interaction is a FENE potential
      (to be combined with Lennard Jones). */
  FENE,
  /** Type of bonded interaction is a HARMONIC potential. */
  HARMONIC,
  /** Type of bonded interaction is a HARMONIC_DUMBBELL potential. */
  HARMONIC_DUMBBELL,
  /** Type of bonded interaction is a QUARTIC potential. */
  QUARTIC,
  /** Type of bonded interaction is a BONDED_COULOMB */
  BONDED_COULOMB,
  /** Type of bonded interaction is a bond angle potential. */
  ANGLE_OLD,
  /** Type of bonded interaction is a dihedral potential. */
  DIHEDRAL,
  /** Type of tabulated bonded interaction potential,
      may be of bond length, of bond angle or of dihedral type. */
  TABULATED,
  /** Type of bonded interaction is a (-LJ) potential. */
  SUBT_LJ,
  /** Type of a Rigid/Constrained bond*/
  RIGID_BOND,
  /** Type of a virtual bond*/
  VIRTUAL_BOND,
  /** Type of bonded interaction is a bond angle -- constraint distance
     potential. */
  ANGLEDIST,
  /** Type of bonded interaction is a bond angle -- chain ends have angle with
     wall constraint */
  ENDANGLEDIST,
  /** Type of overlapped bonded interaction potential,
      may be of bond length, of bond angle or of dihedral type. */
  OVERLAPPED,
  /** Type of bonded interaction is a bond angle cosine potential. */
  ANGLE_HARMONIC,
  /** Type of bonded interaction is a bond angle cosine potential. */
  ANGLE_COSINE,
  /** Type of bonded interaction is a bond angle cosine potential. */
  ANGLE_COSSQUARE,
  /** Type of bonded interaction: oif local forces. */
  OIF_LOCAL_FORCES,
  /** Type of bonded interaction: oif global forces. */
  OIF_GLOBAL_FORCES,
  /** Type of bonded interaction: determining outward direction of oif membrane.
     */
  OIF_OUT_DIRECTION,
  /** Type of bonded interaction for cg DNA */
  CG_DNA_BASEPAIR,
  /** Type of bonded interaction for cg DNA */
  CG_DNA_STACKING,
  /** Type of bonded interaction for cg DNA */
  CG_DNA_BACKBONE,
  /** Type of bonded interaction is a wall repulsion (immersed boundary). */
  IBM_TRIEL,
  /** Type of bonded interaction is volume conservation force (immersed
     boundary). */
  IBM_VOLUME_CONSERVATION,
  /** Type of bonded interaction is bending force (immersed boundary). */
  IBM_TRIBEND,
  /** Type of bonded interaction is umbrella. */
  UMBRELLA
};
#endif
