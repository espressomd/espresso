#ifndef BOND_H
#define BOND_H
//headers from fene.hpp
//#include "utils.hpp"
//#include "interaction_data.hpp"
#include "particle_data.hpp"
//#include "random.hpp"
//#include "errorhandling.hpp"
//#include "debug.hpp"
#include "bonded_interaction.hpp"

/*
This is the basis class for bonded interactions.
It is abstract and all classes inheriting from it
are concrete classes.
 */
class Bond{
public:
  //no definition of these functions because this is an abstract class
  //force calculation
  virtual int add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3])=0;
  //energy calculation
  virtual int add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy)=0;
  
  //get interaction type of bond
  BondedInteraction get_interaction_code();

  // only subclasses have acces to that
protected:
  BondedInteraction interaction_type_code {BONDED_IA_NONE};
  
};

/*
The definition of the concrete classes.
only possible to inherit public from abstact class!
*/

/*NON-BOND - just a placeholder class in order to copy and paste the old algorithm
for the new bond vector instead of bond_ia_params. For the future this should be replaced.
 */
class NO_BOND : public Bond{
public:
  NO_BOND();
  // member functions will be empty, they only exist for preventing runtime errors(if they are accidently called)
  virtual int add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]) override;
  virtual int add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) override;  
};

//FENE Bond
class FENE : public Bond{
public:
  //constructor
  FENE(double r0_input, double drmax_input, double drmax2_input, double drmax2i_input, double k_input);
  // Member function
  virtual int add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]) override;
  virtual int add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) override;
private:
  // specific bond parameters
  double r0;
  double drmax;
  double drmax2;
  double drmax2i;
  double k;

};

// Harmonic dumbbell bond
class HARMONIC_DUMBBELL : public Bond{
public: 
  HARMONIC_DUMBBELL(double k1_i, double k_2_i, double r_i, double r_cut_i);
    // Member function
  virtual int add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]) override;
  virtual int add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) override;
  //bond parameters
private:
  double k1;
  double k2;
  double r;
  double r_cut;

};

// Harmonic bond
class HARMONIC : public Bond{
public: 
  HARMONIC(double k_i, double r_i, double r_cut_i);
    // Member function
  virtual int add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]) override;
  virtual int add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) override;
  //bond parameters
private:
  double k;
  double r;
  double r_cut;

};

// quartic bond
class QUARTIC : public Bond{
public: 
  QUARTIC(double k0_i, double k_1_i, double r_i, double r_cut_i);
    // Member function
  virtual int add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]) override;
  virtual int add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) override;
  //bond parameters
private:
  double k0;
  double k1;
  double r;
  double r_cut;

};

class BONDED_COULOMB : public Bond{
public:
  BONDED_COULOMB(double prefactor_i);
    // Member function
  virtual int add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]) override;
  virtual int add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) override;
  //bond parameters
private:
  double prefactor;
};
#endif
