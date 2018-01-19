#ifndef ABSTRACT_BOND_CLASS_H
#define ABSTRACT_BOND_CLASS_H
#include "BondType.hpp"
#include "particle_data.hpp"
#include "errorhandling.hpp"

#include <vector> //for std::vector
#include <array> //for std::array
#include <boost/optional.hpp> // for boost::optional

/*
This is the base class for bonded interactions.
It is abstract and all classes inheriting from it
are concrete classes.
 */
namespace Bond {
class Bond {
public:

  // constructor
  Bond(int n_partners){m_npartners = n_partners;}
  // virtual desctructor
  virtual ~Bond()=default;

  //### FORCE AND ENERGY ###
  //no definition of these functions because this is an abstract class,
  // which is used as an interface
  //force calculation
  // return value: 0: ok, 1: bond broken, 2: return from "add_bonded_force" in forces_inline.cpp
  virtual int add_bonded_force(Particle *p1, int bl_id) = 0;
  //energy calculation
  // return value: 0: ok, 1: bond broken, 2: return from "add_bonded_force" in forces_inline.cpp
  virtual int add_bonded_energy(Particle *p1, int bl_id) = 0;

  //### BOND PARTNER ###
  // get the total number of bond partners
  int get_number_of_bond_partners() const {return m_npartners;};

  //!!!Definition of template function has to be in header file Bond.hpp!!!
  // get n bond partners
  // p1: Particle, which bond partners are determined
  // bl_id: index in p1->bl.e where the bond is
  template<size_t N>
  boost::optional<std::array<Particle *, N> > get_n_bond_partners(Particle *p1, int bl_id) const
  {

    //definition and initialization of return value
    //if there are null pointers, the function returns nothing
    std::array<Particle *, N> bond_partners;
    // fills entries with NULL pointer
    bond_partners.fill(NULL);

    for(int i=0;i<bond_partners.size();i++){
      // get ith bond partner
      bond_partners[i] = local_particles[p1->bl.e[bl_id+(i+1)]];
      // error message if one partner doesnt exist
      if(!bond_partners[i]){
	runtimeErrorMsg() << "bond broken between particles " << p1->p.identity;
	for(int j=1; j<i+1;j++){
	  runtimeErrorMsg() << ", " << p1->bl.e[bl_id+j] << ", ";
	};
	runtimeErrorMsg() << " (particles not stored on the same node)";
	return boost::none;
      };

    };
    return bond_partners;
  }

  //### BOND TYPE ###
  // get bond type
  BondType get_Bond_Type() const {return m_bondtype;};

  //### VARIABLES ###
  //bond type
  BondType m_bondtype{BondType::BONDED_IA_NONE};
  // number of bond partners
  size_t m_npartners;

};

}
#endif
