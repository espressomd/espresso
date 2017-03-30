#include "HomogeneousMagneticField.hpp"
#include "energy_inline.hpp"

namespace Constraints {

void HomogeneousMagneticField::add_force(Particle *p, double *folded_pos) {
#ifdef ROTATION
#ifdef DIPOLES
    p->f.torque[0] += p->r.dip[1]*m_field[2]-p->r.dip[2]*m_field[1];
    p->f.torque[1] += p->r.dip[2]*m_field[0]-p->r.dip[0]*m_field[2];
    p->f.torque[2] += p->r.dip[0]*m_field[1]-p->r.dip[1]*m_field[0];
#endif
#endif
}

void HomogeneousMagneticField::add_energy(Particle *p, double *folded_pos, Observable_stat &energy) const {
#ifdef ROTATION
#ifdef DIPOLES
  /* TODO:
   * see TCL:
   * double ext_magn_field_energy(Particle *p1, Constraint_ext_magn_field *c)
   * {
   * #ifdef DIPOLES
   *   return -1.0 * scalar(c->ext_magn_field,p1->r.dip);
   * #endif
   *   return 0;
   * }
   */
#endif
#endif    
}

}
