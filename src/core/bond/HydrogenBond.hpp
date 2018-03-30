#ifndef HYDROGEN_BOND_BOND_CLASS_H
#define HYDROGEN_BOND_BOND_CLASS_H
#include "FourParticleBond.hpp"

namespace Bond {

  class HydrogenBond : public FourParticleBond {
    
  public:
    //constructor
    HydrogenBond(double r0, double alpha, double E0, double kd, double sigma1, 
		 double sigma2, double psi10, double psi20, double E0sb, 
		 double r0sb, double alphasb, double f2, double f3) : 
      m_r0{r0}, m_alpha{alpha}, m_E0{E0}, m_kd{kd}, m_sigma1{sigma1}, 
      m_sigma2{sigma2}, m_psi10{psi10}, m_psi20{psi20}, m_E0sb{E0sb}, 
      m_r0sb{r0sb}, m_alphasb{alphasb}, m_f2{f2}, m_f3{f3} {m_bondtype = BondType::BONDED_IA_CG_DNA_BASEPAIR;}
    
    //force calculation
    int calc_bonded_four_particle_force(Particle *p1, Particle *p2, Particle *p3, 
					       Particle *p4, double force[3], double force2[3], 
					       double force3[3], double force4[3]) const override;
    //energy calculation
    int calc_bonded_four_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
						Particle *p4, double *_energy) const override;
    double &r0(){return m_r0;}
    double &alpha(){return m_alpha;}
    double &E0(){return m_E0;}
    double &kd(){return m_kd;}
    double &sigma1(){return m_sigma1;}
    double &sigma2(){return m_sigma2;}
    double &psi10(){return m_psi10;}
    double &psi20(){return m_psi20;}
    double &E0sb(){return m_E0sb;}
    double &r0sb(){return m_r0sb;}
    double &alphasb(){return m_alphasb;}
    double &f2(){return m_f2;}
    double &f3(){return m_f3;}  
    
  private:  
    double m_r0;
    double m_alpha;
    double m_E0;
    double m_kd;
    double m_sigma1;
    double m_sigma2;
    double m_psi10;
    double m_psi20;
    double m_E0sb;
    double m_r0sb;
    double m_alphasb;
    double m_f2;
    double m_f3;  

  };

}

#endif
