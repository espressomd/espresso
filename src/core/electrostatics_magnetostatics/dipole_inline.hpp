#ifndef ESPRESSO_DIPOLE_INLINE_HPP
#define ESPRESSO_DIPOLE_INLINE_HPP

#include "electrostatics_magnetostatics/dipole.hpp"
#include "electrostatics_magnetostatics/p3m-dipolar.hpp"
#include "integrate.hpp"
#include "npt.hpp"

namespace Dipole {
// forces_inline
inline void calc_pair_force(Particle *p1, Particle *p2, double *d, double dist,
                            double dist2, Utils::Vector3d &force) {
  switch (dipole.method) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    // fall trough
  case DIPOLAR_P3M: {
#ifdef NPT
    double eng = dp3m_add_pair_force(p1, p2, d, dist2, dist, force.data());
    if (integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir[0] += eng;
#else
    dp3m_add_pair_force(p1, p2, d, dist2, dist, force.data());
#endif
    break;
  }
#endif /*ifdef DP3M */
  default:
    break;
  }
}

// energy_inline
inline void add_pair_energy(Particle *p1, Particle *p2, double *d, double dist,
                            double dist2, Observable_stat &energy) {
  double ret = 0;
  if (dipole.method != DIPOLAR_NONE) {
    // ret=0;
    switch (dipole.method) {
#ifdef DP3M
    case DIPOLAR_MDLC_P3M:
      // fall trough
    case DIPOLAR_P3M:
      ret = dp3m_pair_energy(p1, p2, d, dist2, dist);
      break;
#endif
    default:
      ret = 0;
    }
    energy.dipolar[0] += ret;
  }
}

} // namespace Dipole

#endif // ESPRESSO_DIPOLE_INLINE_HPP
