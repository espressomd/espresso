#ifndef BOND_TYPE_CLASS_H
#define BOND_TYPE_CLASS_H

namespace Bond {

  enum class BondType {
   /** This bonded interaction was not set. */
   BONDED_IA_NONE=-1,
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
     BONDED_IA_UMBRELLA,
  /** Type of bonded interaction is thermalized distance bond. */
  BONDED_IA_THERMALIZED_DIST,
  /** Type of bonded interaction is a BONDED_COULOMB_P3M_SR */
  BONDED_IA_BONDED_COULOMB_P3M_SR
     
  };

}

#endif
