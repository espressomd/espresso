#include "RigidBond.hpp"
#include "debug.hpp"
#include "grid.hpp" //get_mi_vector
#include "interaction_data.hpp"

//there is neither a direct force or energy contribution
int Bond::RigidBond::add_bonded_energy(Particle *p1, int bl_id)
{

  return 0;
  
}

int Bond::RigidBond::add_bonded_force(Particle *p1, int bl_id)
{
  
  return 0;
  
}

//functions which called by BondContainer
int Bond::RigidBond::pos_corr(Particle *p1, int bl_id, int* repeat, int &cnt)
{
#ifdef BOND_CONSTRAINT  
  double r_ij_t[3], r_ij[3], r_ij_dot, G, pos_corr, r_ij2;
  
  cnt++;
  //get bond partner
  if(auto bond_partners = get_n_bond_partners<1>(p1, bl_id)){

    // get bond partners from bond partners list
    Particle* p2 = (*bond_partners)[0];
    
    get_mi_vector(r_ij, p1->r.p, p2->r.p);
    r_ij2 = sqrlen(r_ij);
    if (fabs(1.0 - r_ij2 / m_d2) > m_p_tol){
      get_mi_vector(r_ij_t, p1->r.p_old, p2->r.p_old);
      r_ij_dot = scalar(r_ij_t, r_ij);
      G = 0.50 * (m_d2 - r_ij2) / r_ij_dot;
#ifdef MASS
      G /= ((*p1).p.mass + (*p2).p.mass);
#else
      G /= 2;
#endif
      for (int j = 0; j < 3; j++) {
	pos_corr = G * r_ij_t[j];
	p1->f.f[j] += pos_corr * (*p2).p.mass;
	p2->f.f[j] -= pos_corr * (*p1).p.mass;
      }
      /*Increase the 'repeat' flag by one */
      *repeat = *repeat + 1;
    };
  }
  else{
    return 2;
  };
#endif
  return 0;
}

int Bond::RigidBond::vel_corr(Particle *p1, int bl_id, int* repeat)
{
#ifdef BOND_CONSTRAINT
  double v_ij[3], r_ij[3], K, vel_corr;
  
  //get bond partner
  if(auto bond_partners = get_n_bond_partners<1>(p1, bl_id)){
    
    // get bond partners from bond partners list
    Particle* p2 = (*bond_partners)[0];
    
    vecsub(p1->m.v, p2->m.v, v_ij);
    get_mi_vector(r_ij, p1->r.p, p2->r.p);
    if (fabs(scalar(v_ij, r_ij)) > m_v_tol) {
      K = scalar(v_ij, r_ij) / m_d2;
#ifdef MASS
      K /= ((*p1).p.mass + (*p2).p.mass);
#else
      K /= 2.0;
#endif
      for (int j = 0; j < 3; j++) {
	vel_corr = K * r_ij[j];
	p1->f.f[j] -= vel_corr * (*p2).p.mass;
	p2->f.f[j] += vel_corr * (*p1).p.mass;
      }
      *repeat = *repeat + 1;
    };
  }
  else{
    return 2;
  };
#endif
  return 0;
}

boost::any Bond::RigidBond::get_bond_parameters_from_bond() const
{

  Rigid_bond_parameters params = {m_d2, m_p_tol, m_v_tol};
  return boost::any(params);
  
}

