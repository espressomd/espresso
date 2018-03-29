#ifndef TABULATED_BOND_DIHEDRAL_BOND_CLASS_H
#define TABULATED_BOND_DIHEDRAL_BOND_CLASS_H
#include "FourParticleBond.hpp"
#include "Tabulated.hpp"
#include "CutoffBond.hpp"

namespace Bond {
  class TabulatedBondDihedral : public FourParticleBond, public Tabulated, public CutoffBond {
  public:

    //constructor
    TabulatedBondDihedral(double min, double max, std::vector<double> const &energy,
			  std::vector<double> const &force) : 
      Tabulated{min, max, energy, force, TabulatedBondedInteraction::TAB_BOND_DIHEDRAL},
      CutoffBond(0.0)
    {m_bondtype = BondType::BONDED_IA_TABULATED;}

    //force calculation
    int calc_bonded_four_particle_force(Particle *p1, Particle *p2, Particle *p3, 
				       Particle *p4, double force[3], double force2[3], 
				       double force3[3], double force4[3]) const override;
    //energy calculation
    int calc_bonded_four_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
					Particle *p4, double *_energy) const override;

    boost::any get_bond_parameters_from_bond() const override;

    
    
  private:  
    //internal function
    void calc_dihedral_angle(Particle *p1, Particle *p2, Particle *p3, Particle *p4, 
			     double a[3], double b[3], double c[3], 
			     double aXb[3], double *l_aXb, double bXc[3], double *l_bXc, 
			     double *cosphi, double *phi) const;
  };
}
#endif
