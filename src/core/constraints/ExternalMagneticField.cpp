#include "ExternalMagneticField.hpp"
#include "interaction_data.hpp"
#include "grid.hpp"

namespace Constraints {

  void ExternalMagneticField::add_energy(Particle *p, const double *folded_pos, Observable_stat &energy) {
#ifdef DIPOLES
    energy.dipolar[0] += -1.0 * ext_magn_field.dot(Vector3d(p->r.dip));
#endif
  }

void ExternalMagneticField::add_force(Particle *p, const double *folded_pos) {
#ifdef ROTATION
#ifdef DIPOLES
    p->f.torque[0] += p->r.dip[1]*ext_magn_field[2]-p->r.dip[2]*ext_magn_field[1];
    p->f.torque[1] += p->r.dip[2]*ext_magn_field[0]-p->r.dip[0]*ext_magn_field[2];
    p->f.torque[2] += p->r.dip[0]*ext_magn_field[1]-p->r.dip[1]*ext_magn_field[0];
#endif
#endif
  }

std::map<std::string, Variant> ExternalMagneticField::get_parameters() {
  std::map<std::string, Variant> p;
  p["field"] = ext_magn_field;
  
  return p;
}

void ExternalMagneticField::set_parameter(const std::string &name, const Variant &value) {
  SET_PARAMETER_HELPER("field", ext_magn_field);
}

Parameters ExternalMagneticField::all_parameters() const {
  static bool init = false;
  static Parameters p;
  if(!init) {
    p["field"] = Parameter(Parameter::Type::DOUBLE_VECTOR, 3, true);
    init = true;
  }  

  return p;
}


}
