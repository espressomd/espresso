#ifndef TWIST_STACK_BOND_CLASS_H
#define TWIST_STACK_BOND_CLASS_H
#include "EightParticleBond.hpp"
#include <stdlib.h> //free

namespace Bond {
  
  class TwistStack : public EightParticleBond {
  public:
    TwistStack(double rm, double epsilon, double ref_pot, double a[8], double b[7]) : 
      m_rm{rm}, m_epsilon{epsilon}, m_ref_pot{ref_pot}, m_a{a}, m_b{b} 
    {m_bondtype = BondType::BONDED_IA_CG_DNA_STACKING;}

    ~TwistStack(){
      free(m_a);
      free(m_b);
    }

    int calc_bonded_eight_particle_force(Particle *p1, Particle *p2, Particle *p3, 
					Particle *p4, Particle *p5, Particle *p6, 
					Particle *p7, Particle *p8, double force[3], 
					double force2[3], double force3[3], 
					double force4[3], double force5[3],
					double force6[3], double force7[3],
					double force8[3]) const override;

    int calc_bonded_eight_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
					 Particle *p4, Particle *p5, Particle *p6, 
					 Particle *p7, Particle *p8, 
					 double *_energy) const override;
    
  private:
    double m_rm;
    double m_epsilon;
    double m_ref_pot;
    double* m_a;
    double* m_b;
  };

}

#endif
