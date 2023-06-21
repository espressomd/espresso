/*
 * Copyright (C) 2010-2022 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config/config.hpp"

#ifdef DIPOLES

#include "magnetostatics/dipoles.hpp"

#include "ParticleRange.hpp"
#include "actor/traits.hpp"
#include "actor/visit_try_catch.hpp"
#include "actor/visitors.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "npt.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/optional.hpp>

#include <cassert>
#include <cstdio>
#include <stdexcept>

boost::optional<MagnetostaticsActor> magnetostatics_actor;

namespace Dipoles {

void sanity_checks() {
  if (magnetostatics_actor) {
    boost::apply_visitor([](auto &actor) { actor->sanity_checks(); },
                         *magnetostatics_actor);
  }
}

void on_dipoles_change() {
  visit_active_actor_try_catch([](auto &actor) { actor->init(); },
                               magnetostatics_actor);
}

void on_boxl_change() {
  visit_active_actor_try_catch([](auto &actor) { actor->on_boxl_change(); },
                               magnetostatics_actor);
}

void on_node_grid_change() {
  if (magnetostatics_actor) {
    boost::apply_visitor([](auto &actor) { actor->on_node_grid_change(); },
                         *magnetostatics_actor);
  }
}

void on_periodicity_change() {
  visit_active_actor_try_catch(
      [](auto &actor) { actor->on_periodicity_change(); },
      magnetostatics_actor);
}

void on_cell_structure_change() {
  visit_active_actor_try_catch(
      [](auto &actor) { actor->on_cell_structure_change(); },
      magnetostatics_actor);
}

void calc_pressure_long_range() {
  if (magnetostatics_actor) {
    runtimeWarningMsg() << "pressure calculated, but pressure not implemented.";
  }
}

double cutoff() {
#ifdef DP3M
  if (auto dp3m = get_actor_by_type<DipolarP3M>(magnetostatics_actor)) {
    return dp3m->dp3m.params.r_cut;
  }
#endif
  return -1.;
}

void on_observable_calc() {
#ifdef DP3M
  if (auto dp3m = get_actor_by_type<DipolarP3M>(magnetostatics_actor)) {
    dp3m->count_magnetic_particles();
  }
#endif
}

struct LongRangeForce : public boost::static_visitor<void> {
  ParticleRange const &m_particles;
  explicit LongRangeForce(ParticleRange const &particles)
      : m_particles(particles) {}

#ifdef DP3M
  void operator()(std::shared_ptr<DipolarP3M> const &actor) const {
    actor->dipole_assign(m_particles);
#ifdef NPT
    if (integ_switch == INTEG_METHOD_NPT_ISO) {
      auto const energy = actor->kernel(true, true, m_particles);
      npt_add_virial_contribution(energy);
      fprintf(stderr, "dipolar_P3M at this moment is added to p_vir[0]\n");
    } else
#endif // NPT
      actor->kernel(true, false, m_particles);
  }
#endif // DP3M
  void operator()(std::shared_ptr<DipolarLayerCorrection> const &actor) const {
    actor->add_force_corrections(m_particles);
    boost::apply_visitor(*this, actor->base_solver);
  }
  void operator()(std::shared_ptr<DipolarDirectSum> const &actor) const {
    actor->add_long_range_forces(m_particles);
  }
#ifdef DIPOLAR_DIRECT_SUM
  void operator()(std::shared_ptr<DipolarDirectSumGpu> const &actor) const {
    actor->add_long_range_forces();
  }
#endif
#ifdef DIPOLAR_BARNES_HUT
  void operator()(std::shared_ptr<DipolarBarnesHutGpu> const &actor) const {
    actor->add_long_range_forces();
  }
#endif
#ifdef SCAFACOS_DIPOLES
  void operator()(std::shared_ptr<DipolarScafacos> const &actor) const {
    actor->add_long_range_forces();
  }
#endif
};

struct LongRangeEnergy : public boost::static_visitor<double> {
  ParticleRange const &m_particles;
  explicit LongRangeEnergy(ParticleRange const &particles)
      : m_particles(particles) {}

#ifdef DP3M
  double operator()(std::shared_ptr<DipolarP3M> const &actor) const {
    actor->dipole_assign(m_particles);
    return actor->kernel(false, true, m_particles);
  }
#endif // DP3M
  double
  operator()(std::shared_ptr<DipolarLayerCorrection> const &actor) const {
    auto energy = boost::apply_visitor(*this, actor->base_solver);
    return energy + actor->energy_correction(m_particles);
  }
  double operator()(std::shared_ptr<DipolarDirectSum> const &actor) const {
    return actor->long_range_energy(m_particles);
  }
#ifdef DIPOLAR_DIRECT_SUM
  double operator()(std::shared_ptr<DipolarDirectSumGpu> const &actor) const {
    actor->long_range_energy();
    return 0.;
  }
#endif
#ifdef DIPOLAR_BARNES_HUT
  double operator()(std::shared_ptr<DipolarBarnesHutGpu> const &actor) const {
    actor->long_range_energy();
    return 0.;
  }
#endif
#ifdef SCAFACOS_DIPOLES
  double operator()(std::shared_ptr<DipolarScafacos> const &actor) const {
    return actor->long_range_energy();
  }
#endif
};

#ifdef DIPOLE_FIELD_TRACKING
struct LongRangeField : public boost::static_visitor<void> {
  ParticleRange const &m_particles;
  explicit LongRangeField(ParticleRange const &particles)
      : m_particles(particles) {}

  void operator()(std::shared_ptr<DipolarDirectSum> const &actor) const {
    actor->dipole_field_at_part(m_particles);
  }

  template <typename T,
            std::enable_if_t<!traits::has_dipole_fields<T>::value> * = nullptr>
  void operator()(std::shared_ptr<T> const &) const {
    runtimeErrorMsg() << "Dipoles field calculation not implemented by "
                      << "dipolar method " << Utils::demangle<T>();
  }
};
#endif

void calc_long_range_force(ParticleRange const &particles) {
  if (magnetostatics_actor) {
    boost::apply_visitor(LongRangeForce(particles), *magnetostatics_actor);
  }
}

double calc_energy_long_range(ParticleRange const &particles) {
  if (magnetostatics_actor) {
    return boost::apply_visitor(LongRangeEnergy(particles),
                                *magnetostatics_actor);
  }
  return 0.;
}

#ifdef DIPOLE_FIELD_TRACKING
void calc_long_range_field(ParticleRange const &particles) {
  if (magnetostatics_actor) {
    boost::apply_visitor(LongRangeField(particles), *magnetostatics_actor);
  }
}
#endif

namespace detail {
bool flag_all_reduce(bool flag) {
  return boost::mpi::all_reduce(comm_cart, flag, std::logical_or<>());
}
} // namespace detail

} // namespace Dipoles
#endif // DIPOLES
