#ifndef ABSTRACT_BOND_CLASS_H
#define ABSTRACT_BOND_CLASS_H
#include "particle_data.hpp"
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
  //no definition of these functions because this is an abstract class
  //force calculation
  virtual int add_bonded_force(Particle *p1, Particle *p2, double dx[3], 
			       double force[3]) const=0;
  //energy calculation
  virtual int add_bonded_energy(Particle *p1, Particle *p2, double dx[3], 
				double *_energy) const=0;
};

}
#endif
