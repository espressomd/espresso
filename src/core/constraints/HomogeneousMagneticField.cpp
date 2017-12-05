#include "HomogeneousMagneticField.hpp"
#include "energy.hpp"

namespace Constraints {

void HomogeneousMagneticField::add_force(Particle *p, double *folded_pos) {
#ifdef ROTATION
#ifdef DIPOLES
    double c[3];
    Utils::cross_product(p->r.dip, &m_field.front(), c);
    for (int i=0; i<3; ++i) {
        p->f.torque[i] += c[i];
    }
#endif
#endif
}

void HomogeneousMagneticField::add_energy(Particle *p, double *folded_pos, Observable_stat &energy) const {
#ifdef DIPOLES
    energy.dipolar[0] += -1.0 * Utils::dot_product(&m_field.front(), p->r.dip);
#endif
}

}
