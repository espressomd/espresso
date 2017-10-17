#include "Bond.hpp"

Bond::Bond::Bond(int n_partners)
{
  m_npartners = n_partners;
  //zero initialize
  m_bond_partners = new Particle* [m_npartners]();

}

Bond::Bond::~Bond()
{

  delete[] m_bond_partners;

}

int Bond::Bond::get_n_bond_partners(Particle *p1, int bl_id)
{

  //reset partners to get null pointers
  //not necessary but then we know the pointers
  reset_bond_partners();

  for(int i=0;i<m_npartners;i++){
    // get ith bond partner
    m_bond_partners[i] = local_particles[p1->bl.e[bl_id+(i+1)]];
    // error message if one partner doesnt exist
    if(!m_bond_partners[i]){
      runtimeErrorMsg() << "bond broken between particles " << p1->p.identity;
	for(int j=1; j<i+1;j++){
	  runtimeErrorMsg() << ", " << p1->bl.e[bl_id+j] << ", ";
	}
	runtimeErrorMsg() << " (particles not stored on the same node)";
      return 1;
    };

  };

  return 0;
}
void Bond::Bond::reset_bond_partners()
{
  for(int i=0;i<m_npartners;i++){
    m_bond_partners[i]=NULL;
  };
}

// as default: no contribution to virial pressure
// exception:pair bonds
int Bond::Bond::add_virial(Particle *p1, int bl_id){

  return 0;

}

//default no contribution to three body pressure
//exception: angle bonds and tabulated angle
int Bond::Bond::add_three_body_pressure(Particle *p1, int bl_id){

  return 0;

}

// default
// only for pair bonds there is a modification
int Bond::Bond::calc_pair_force(Particle *p1, Particle *p2,  int bl_id, double force[3])
{

  p2 = NULL;
  for(int i =0;i<3;i++)
    {
      force[i] = 0.0;
    };
  return 0;

}


