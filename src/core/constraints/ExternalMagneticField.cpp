#include "ExternalMagneticField.hpp"
#include "interaction_data.hpp"
#include "grid.hpp"

namespace Constraints {

  void ExternalMagneticField::add_energy(const Particle *p, const double *folded_pos, Observable_stat &energy) {
#ifdef DIPOLES
    energy.dipolar[0] += -1.0 * scalar(ext_magn_field,p->r.dip);
#endif
  }

    void ExternalMagneticField::add_force(Particle *p, const double *folded_pos) {
    ;
#ifdef ROTATION
#ifdef DIPOLES
    // p1->f.torque[0] += p1->r.dip[1]*c->ext_magn_field[2]-p1->r.dip[2]*c->ext_magn_field[1];
    // p1->f.torque[1] += p1->r.dip[2]*c->ext_magn_field[0]-p1->r.dip[0]*c->ext_magn_field[2];
    // p1->f.torque[2] += p1->r.dip[0]*c->ext_magn_field[1]-p1->r.dip[1]*c->ext_magn_field[0];
#endif
#endif

  }
  }
