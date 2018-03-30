#ifndef ANGLE_DIST_BOND_CLASS_H
#define ANGLE_DIST_BOND_CLASS_H
#include "ThreeParticleBond.hpp"

#ifdef BOND_ANGLEDIST

namespace Bond {

  class AngleDist : public ThreeParticleBond {

  public:
    AngleDist(double bend, double phimin, double distmin, double phimax, double distmax) :
      m_bend{bend}, m_phimin{phimin}, m_distmin{distmin}, m_phimax{phimax}, m_distmax{distmax}
    {m_bondtype = BondType::BONDED_IA_ANGLEDIST;}

    //force *
    int calc_bonded_three_particle_force(Particle *p1, Particle *p2, Particle *p3, double force[3],
					 double force2[3]) const override;

    //energy *
    int calc_bonded_three_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
					  double *_energy) const override;

    double &bend(){return m_bend;}
    double &phimin(){return m_phimin;}
    double &distmin(){return m_distmin;}
    double &phimax(){return m_phimax;}
    double &distmax(){return m_distmax;}
    
  private:
    /** Function to calculate wall distance and phi0(dist).
	Called by \ref calc_angledist_force */
    double calc_angledist_param(Particle *p1, Particle *p2,
				Particle *p3) const;
    //variables
    double m_bend;
    double m_phimin;
    double m_distmin;
    double m_phimax;
    double m_distmax;

  };

}

#endif //BOND_ANGLE_DIST

#endif
