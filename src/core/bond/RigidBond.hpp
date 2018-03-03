#ifndef RIGID_BOND_CLASS_H
#define RIGID_BOND_CLASS_H
#include "Bond.hpp"
#include "CutoffBond.hpp"

namespace Bond{

  class RigidBond : public Bond, public CutoffBond{

  public:
    
    RigidBond(double d, double p_tol, double v_tol, double d2) :
      Bond(1), CutoffBond(sqrt(d2)), m_d{d}, m_p_tol{2*p_tol}, m_v_tol{v_tol}, m_d2{d2}
    {m_bondtype = BondType::BONDED_IA_RIGID_BOND;};

    int pos_corr(Particle *p1, int bl_id, int* repeat, int &cnt);
    int vel_corr(Particle *p1, int bl_id, int* repeat);
    
    //virtual functions
    // return value: 0: ok, 1: bond broken, 2: return from "add_bonded_force" in forces_inline.cpp
    virtual int add_bonded_force(Particle *p1, int bl_id) override;
    //energy calculation
    // return value: 0: ok, 1: bond broken, 2: return from "add_bonded_force" in forces_inline.cpp
    virtual int add_bonded_energy(Particle *p1, int bl_id) override;

  private:
    //member variables
    double m_d;
    double m_p_tol;
    double m_v_tol;
    double m_d2;
    
  };
  
}

#endif
