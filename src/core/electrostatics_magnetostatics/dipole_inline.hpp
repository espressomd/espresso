#ifndef ESPRESSO_DIPOLE_INLINE_HPP
#define ESPRESSO_DIPOLE_INLINE_HPP

#include "electrostatics_magnetostatics/dipole.hpp"
#include "electrostatics_magnetostatics/p3m-dipolar.hpp"
#include "integrate.hpp"
#include "npt.hpp"

namespace Dipole {
// forces_inline
inline std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
pair_force(Particle const &p1, Particle const &p2, Utils::Vector3d const &d,
           double dist, double dist2) {
  Utils::Vector3d force{}, torque1{}, torque2{};
  switch (dipole.method) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    // fall trough
  case DIPOLAR_P3M: {
    double eng;
    std::tie(eng, force, torque1, torque2) =
        dp3m_pair_force(p1, p2, d, dist2, dist);
#ifdef NPT
    if (integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir[0] += eng;
#endif
    break;
  }
#endif /*ifdef DP3M */
  default:
    break;
  }
  return std::make_tuple(force, torque1, torque2);
}

// energy_inline
inline double pair_energy(Particle const &p1, Particle const &p2,
                          Utils::Vector3d const &d, double dist, double dist2) {
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
  }
  return ret;
}

} // namespace Dipole

#endif // ESPRESSO_DIPOLE_INLINE_HPP
