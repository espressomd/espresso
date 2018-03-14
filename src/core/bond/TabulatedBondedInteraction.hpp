#ifndef TABULATED_BONDED_INTERACTION_CLASS_H
#define TABULATED_BONDED_INTERACTION_CLASS_H

namespace Bond {

  enum class TabulatedBondedInteraction {
    
    TAB_UNKNOWN = 0,
    TAB_BOND_LENGTH = 1,
    TAB_BOND_ANGLE = 2,
    TAB_BOND_DIHEDRAL = 3
      
  };

}

#endif
