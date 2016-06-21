#include "ChargedRod.hpp"
#include "interaction_data.hpp"
#include "grid.hpp"

#define C_GAMMA 0.57721566490153286060651209008

namespace Constraints {

  void ChargedRod::add_energy(Particle *p, const double *folded_pos, Observable_stat &energy) {
#ifdef ELECTROSTATICS
    int i;
    double vec[2], c_dist_2;

    c_dist_2 = 0.0;
    for(i=0;i<2;i++) {
      vec[i] = folded_pos[i] - center[i];
      c_dist_2 += SQR(vec[i]);
    }

    if ((coulomb.prefactor != 0.0) && (p->p.q != 0.0) && (lambda != 0.0)) {
      energy.coulomb[0] += coulomb.prefactor*p->p.q*lambda*(-log(c_dist_2*SQR(box_l_i[2])) + 2*(M_LN2 - C_GAMMA));
    }
#endif
  }

  void ChargedRod::add_force(Particle *p, const double *folded_pos) {
#ifdef ELECTROSTATICS
    int i;
    double fac, vec[2], c_dist_2;

    c_dist_2 = 0.0;
    for(i=0;i<2;i++) {
      vec[i] = folded_pos[i] - center[i];
      c_dist_2 += SQR(vec[i]);
    }

    if ((coulomb.prefactor != 0.0) && (p->p.q != 0.0) && (lambda != 0.0)) {
      fac = 2*coulomb.prefactor*lambda*p->p.q/c_dist_2;
      p->f.f[0]  += fac*vec[0];
      p->f.f[1]  += fac*vec[1];
      total_force[0] -= fac*vec[0];
      total_force[1] -= fac*vec[1];
    }
#endif
  };

std::map<std::string, Variant> ChargedRod::get_parameters() {
  std::map<std::string, Variant> p;

  p["center"] = center;    
  p["lambda"] = lambda;
      
  return p;
}

void ChargedRod::set_parameter(const std::string &name, const Variant &value) {
  SET_PARAMETER_HELPER("center", center);
  SET_PARAMETER_HELPER("lambda", lambda);
}

Parameters ChargedRod::all_parameters() const {
  static bool init = false;
  static Parameters p;
  if(!init) {
    p["center"] = Parameter(Parameter::Type::DOUBLE_VECTOR, 3, true);
    p["lambda"] = Parameter(Parameter::Type::DOUBLE, true);
    init = true;
  }  

  return p;
}

}
