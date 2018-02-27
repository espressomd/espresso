#ifndef IBM_TRIBEND_BOND_CLASS_H
#define IBM_TRIBEND_BOND_CLASS_H
#include "Bond.hpp"
#include "tBendingMethod.hpp"

namespace Bond {

  class IbmTribend : public Bond {
  public :
    //constructor
    IbmTribend(double kb, tBendingMethod method, double theta0, int numNeighbors) : 
      Bond(numNeighbors), m_kb{kb}, m_method{method}, m_theta0{theta0}, m_numNeighbors{numNeighbors}
    {m_bondtype = BondType::BONDED_IA_IBM_TRIBEND; m_npartners = numNeighbors;}

    //functions from bond
    int add_bonded_force(Particle *p1, int bl_id) override;
    int add_bonded_energy(Particle *p1, int bl_id) override;
    //reset function
    int ResetParams(const double kb);

  private:
    //internal functions
    void CalcForceGompper(Particle *p1, Particle **const partners) const;
    void IBM_Tribend_CalcForce(Particle *p1, Particle **const partners) const;

    //variables
    double m_kb;
    tBendingMethod m_method;
    double m_theta0;
    int m_numNeighbors;
  };

}

#endif
