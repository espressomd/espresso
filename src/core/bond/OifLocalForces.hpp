#ifndef OIF_LOCAL_FORCES_BOND_CLASS_H
#define OIF_LOCAL_FORCES_BOND_CLASS_H
#include "FourParticleBond.hpp"

namespace Bond {

  class OifLocalForces : public FourParticleBond {
    
  public:
    //constructor
    OifLocalForces(double phi0, double kb, double r0, 
		   double ks, double kslin, double A01, 
		   double A02, double kal) : m_phi0{phi0}, 
      m_kb{kb}, m_r0{r0}, m_ks{ks}, m_kslin{kslin}, 
      m_A01{A01}, m_A02{A02}, m_kal{kal} {m_bondtype = BondType::BONDED_IA_OIF_LOCAL_FORCES;}
    
    //force calculation
    int calc_bonded_four_particle_force(Particle *p1, Particle *p2, Particle *p3, 
					       Particle *p4, double force[3], double force2[3], 
					       double force3[3], double force4[3]) const override;
    //energy calculation
    int calc_bonded_four_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
						Particle *p4, double *_energy) const override;
    
    // is not default -> has to be modified
    void write_force_to_particle(Particle *p1, Particle *p2, Particle *p3, 
				 Particle *p4, double force[3], double force2[3],
				 double force3[3], double force4[3]) const override;


    boost::any get_bond_parameters_from_bond() const override;
    
    double &phi0(){return m_phi0;}
    double &kb(){return m_kb;}
    double &r0(){return m_r0;}
    double &ks(){return m_ks;}
    double &kslin(){return m_kslin;}
    double &A01(){return m_A01;}
    double &A02(){return m_A02;}
    double &kal(){return m_kal;}
    
  private:  
    double m_phi0;
    double m_kb;
    double m_r0;
    double m_ks;
    double m_kslin;
    double m_A01;
    double m_A02;
    double m_kal;

  };

}

#endif
