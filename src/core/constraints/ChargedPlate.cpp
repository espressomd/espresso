#include "ChargedPlate.hpp"
#include "interaction_data.hpp"
#include "grid.hpp"

namespace Constraints {
  void ChargedPlate::add_energy(Particle *p, const double *folded_pos, Observable_stat &energy) {
#ifdef ELECTROSTATICS
    if (coulomb.prefactor != 0.0 && p->p.q != 0.0 && sigma != 0.0)
      energy.coulomb[0] += -2*M_PI*coulomb.prefactor*sigma*p->p.q*fabs(folded_pos[2] - pos);
#endif
  }

  void ChargedPlate::add_force(Particle *p, const double *folded_pos) {
#ifdef ELECTROSTATICS
    double f;

    if (coulomb.prefactor != 0.0 && p->p.q != 0.0 && sigma != 0.0) {
      f = 2*M_PI*coulomb.prefactor*sigma*p->p.q;
      if (folded_pos[2] < pos)
        f = -f;
      p->f.f[2]  += f;
      total_force[2] -= f;
    }
#endif
  }

std::map<std::string, Variant> ChargedPlate::get_parameters() {
  std::map<std::string, Variant> p;

  p["height"] = pos;    
  p["sigma"] = sigma;
      
  return p;
}

void ChargedPlate::set_parameter(const std::string &name, const Variant &value) {
  SET_PARAMETER_HELPER("height", pos);
  SET_PARAMETER_HELPER("sigma", sigma);
}

Parameters ChargedPlate::all_parameters() const {
  static bool init = false;
  static Parameters p;
  if(!init) {
    p["height"] = Parameter(Parameter::Type::DOUBLE, true);
    p["sigma"] = Parameter(Parameter::Type::DOUBLE, true);
    init = true;
  }  

  return p;
}

}
