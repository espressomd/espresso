#ifndef ABSTRACT_BOND_CLASS_H
#define ABSTRACT_BOND_CLASS_H
#include "particle_data.hpp"
#include <vector>
#include "BondType.hpp"

/*
This is the basis class for bonded interactions.
It is abstract and all classes inheriting from it
are concrete classes.
 */
namespace Bond {
class Bond {
public:
  // virtual desctructor
  virtual ~Bond() = default;
  //no definition of these functions because this is an abstract class,
  // which is used as an interface
  //force calculation
  virtual int add_bonded_force(Particle *p1, int bl_id) const=0;
  //energy calculation
  virtual int add_bonded_energy(Particle *p1, int bl_id, double *_energy) const=0;
  
  // get the total number of bond partners
  virtual int get_number_of_bond_partners() const=0;

  // get bond type
  BondType get_Bond_Type() {return m_bondtype;};

  // Variables
  BondType m_bondtype = BondType::BONDED_IA_NONE;

};

}
#endif
