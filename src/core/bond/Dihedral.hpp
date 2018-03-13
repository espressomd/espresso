#ifndef DIHEDRAL_BOND_CLASS_H
#define DIHEDRAL_BOND_CLASS_H
#include "FourParticleBond.hpp"
#include "CutoffBond.hpp"

namespace Bond {

  class Dihedral : public FourParticleBond, public CutoffBond {
    
  public:
    //constructor
    Dihedral(double mult, double bend, double phase) :
      CutoffBond(0.0),
      m_mult{mult}, m_bend{bend}, m_phase{phase} 
    {m_bondtype = BondType::BONDED_IA_DIHEDRAL;}
    
    //force calculation
    int calc_bonded_four_particle_force(Particle *p2, Particle *p1, Particle *p3, 
					Particle *p4, double force2[3], double force[3], 
					double force3[3], double force4[3]) const override;
    //energy calculation
    int calc_bonded_four_particle_energy(Particle *p2, Particle *p1, Particle *p3, 
					 Particle *p4, double *_energy) const override;

    void write_force_to_particle(Particle *p1, Particle *p2, Particle *p3, 
				 Particle *p4, double force[3], double force2[3],
				 double force3[3], double force4[3]) const override;

    boost::any get_bond_parameters_from_bond() const override;
    
  private:  
    //internal function
    void calc_dihedral_angle(Particle *p1, Particle *p2, Particle *p3, Particle *p4, 
			     double a[3], double b[3], double c[3], 
			     double aXb[3], double *l_aXb, double bXc[3], double *l_bXc, 
			     double *cosphi, double *phi) const;

    // variables
    double m_mult;
    double m_bend;
    double m_phase;

  };

}

#endif
