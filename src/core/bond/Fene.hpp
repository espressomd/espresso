#ifndef FENE_BOND_CLASS_H
#define FENE_BOND_CLASS_H
#include"PairBond.hpp"
#include"CutoffBond.hpp"

/*
The definition of the concrete classes.
only possible to inherit public from abstact class!
*/

//FENE Bond
namespace Bond {
  class Fene : public PairBond, public CutoffBond{
  public:
    //constructor
    // initializer list -> faster
    Fene(double r0, double drmax, double k) : CutoffBond(r0 + drmax), 
			   m_r0{r0}, m_drmax{drmax},
		     m_drmax2{drmax*drmax}, m_drmax2i{1./(drmax*drmax)}, m_k{k}
    {m_bondtype = BondType::BONDED_IA_FENE;}
    // Member function 
    int calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3],
			       double force[3]) const override;
    int calc_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3],
				double *_energy) const override;
    double &r0() { return m_r0; }
    
    const double &dr_max() const { return m_drmax; }
    
    void set_dr_max(double dr_max) {
      m_drmax  = dr_max;
      m_drmax2 = dr_max*dr_max;
      m_drmax2i = 1. / m_drmax2;
    }
    double &k() { return m_k; }

    //get parameters
    boost::any get_bond_parameters_from_bond() const override;
    
  private:
    // specific bond parameters
    double m_r0;
    double m_drmax;
    double m_drmax2;
    double m_drmax2i;
    double m_k;
    
  };
  
}
#endif
