#include "HomogeneousMagneticField.hpp"
#include "energy.hpp"

namespace Constraints {

ParticleForce HomogeneousMagneticField::force(const Particle &p, const Vector3d &folded_pos) {
#ifdef ROTATION
#ifdef DIPOLES
    return {Vector3d{}, Vector3d::cross(p.r.dip, m_field)};
#endif
#else
    return {Vector3d{}};
#endif
}

void HomogeneousMagneticField::add_energy(const Particle &p, const Vector3d &folded_pos, Observable_stat &energy) const {
#ifdef DIPOLES
    energy.dipolar[0] += -1.0 * m_field * p.r.dip;
#endif
}

}
