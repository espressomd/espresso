#ifndef BOND_H
#define BOND_H
#include "particle_data.hpp"
#include "BondType.hpp"

/*
This is the basis class for bonded interactions.
It is abstract and all classes inheriting from it
are concrete classes.
 */
class Bond {
public:

  // virtual desctructor
  virtual ~Bond() = default;
  //no definition of these functions because this is an abstract class
  //force calculation
  virtual int add_bonded_force(Particle *p1, Particle *p2, double dx[3], 
			       double force[3])=0;
  //energy calculation
  virtual int add_bonded_energy(Particle *p1, Particle *p2, double dx[3], 
				double *_energy)=0;
  
  //get interaction type of bond
  virtual BondType get_interaction_code();
};

/*
The definition of the concrete classes.
only possible to inherit public from abstact class!
*/

//FENE Bond
class Fene : public Bond{
public:
  //constructor
  Fene(double r0_input, double drmax_input, double drmax2_input, double drmax2i_input, double k_input);
  // Member function
  int add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]) override;
  int add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) override;
  BondType get_interaction_code() override;
private:
  // specific bond parameters
  double m_r0;
  double m_drmax;
  double m_drmax2;
  double m_drmax2i;
  double m_k;
  BondType m_interaction_type_code;

};

// Harmonic dumbbell bond
class HarmonicDumbbell : public Bond{
public: 
  HarmonicDumbbell(double k1_i, double k_2_i, double r_i, double r_cut_i);
    // Member function
  int add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]) override;
  int add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) override;
  BondType get_interaction_code() override;
  //bond parameters
private:
  double m_k1;
  double m_k2;
  double m_r;
  double m_r_cut;
  BondType m_interaction_type_code;
};

// Harmonic bond
class Harmonic : public Bond{
public: 
  Harmonic(double k_i, double r_i, double r_cut_i);
    // Member function
  int add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]) override;
  int add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) override;
  BondType get_interaction_code() override;
  //bond parameters
private:
  double m_k;
  double m_r;
  double m_r_cut;
  BondType m_interaction_type_code;
};

// quartic bond
class Quartic : public Bond{
public: 
  Quartic(double k0_i, double k_1_i, double r_i, double r_cut_i);
    // Member function
  int add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]) override;
  int add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) override;
  BondType get_interaction_code() override;
  //bond parameters
private:
  double m_k0;
  double m_k1;
  double m_r;
  double m_r_cut;
  BondType m_interaction_type_code;
};

class BondedCoulomb : public Bond{
public:
  BondedCoulomb(double prefactor_i);
    // Member function
  int add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]) override;
  int add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) override;
  BondType get_interaction_code() override;
  //bond parameters
private:
  double m_prefactor;
  BondType m_interaction_type_code;
};
#endif
