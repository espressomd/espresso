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

  // constructor
  Bond(int n_partners);
  // virtual desctructor
  virtual ~Bond();

  //### FORCE AND ENERGY ###
  //no definition of these functions because this is an abstract class,
  // which is used as an interface
  //force calculation
  // return value: 0: ok, 1: bond broken, 2: return from "add_bonded_force" in forces_inline.cpp
  virtual int add_bonded_force(Particle *p1, int bl_id) = 0;
  //energy calculation
  // return value: 0: ok, 1: bond broken, 2: return from "add_bonded_force" in forces_inline.cpp
  virtual int add_bonded_energy(Particle *p1, int bl_id, double *_energy) = 0;

  //### PRESSURE CALCULATION ###
  // for pressure calculation in 
  // add_bonded_virials 
  // in pressure.hpp
  // return value: 0: ok, 1: bond broken, 2: return from "add_bonded_virials" in pressure.hpp
  virtual int add_virial(Particle *p1, int bl_id);
  //for add_three_body_bonded_stress
  virtual int add_three_body_pressure(Particle *p1, int bl_id);
  // for calculating stress tensor in pressure.cpp: int local_stress_tensor_calc()
  virtual int calc_pair_force(Particle *p1, Particle *p2, int bl_id, double force[3]);

  //### BOND PARTNER ###
  // get the total number of bond partners
  int get_number_of_bond_partners() const {return m_npartners;};

  // get n bond partners
  // p1: Particle, which bond partners are determined
  // bl_id: index in p1->bl.e where the bond is
  int get_n_bond_partners(Particle *p1, int bl_id);

  //### BOND TYPE ###
  // get bond type
  BondType get_Bond_Type() const {return m_bondtype;};

  //### VARIABLES ###
  //bond type
  BondType m_bondtype{BondType::BONDED_IA_NONE};
  // number of bond partners
  int m_npartners;
  // bond partners
  Particle** m_bond_partners;

private:
  // asign value of bond partners to NULL
  void reset_bond_partners();

};

}
#endif
