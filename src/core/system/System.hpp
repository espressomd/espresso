/*
 * Copyright (C) 2014-2022 The ESPResSo project
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

#pragma once

#include "config/config.hpp"

#include "GpuParticleData.hpp"
#include "ResourceCleanup.hpp"

#include "electrostatics/solver.hpp"
#include "magnetostatics/solver.hpp"

#include "ek/Solver.hpp"
#include "lb/Solver.hpp"

#include <utils/Vector.hpp>

#include <memory>

class BoxGeometry;
class LocalBox;
struct CellStructure;
class InteractionsNonBonded;

namespace System {

class System {
public:
  System();
#ifdef CUDA
  GpuParticleData gpu;
#endif
  ResourceCleanup cleanup_queue;

  Utils::Vector3d box() const;
  void init();

  /** @brief Get @ref force_cap. */
  auto get_force_cap() { return force_cap; }

  /** @brief Set @ref force_cap. */
  void set_force_cap(double value);

  /** @brief Get @ref min_global_cut. */
  double get_min_global_cut() const { return min_global_cut; }

  /** @brief Set @ref min_global_cut. */
  void set_min_global_cut(double new_value);

  Coulomb::Solver coulomb;
  Dipoles::Solver dipoles;
  LB::Solver lb;
  EK::Solver ek;
  std::shared_ptr<BoxGeometry> box_geo;
  std::shared_ptr<LocalBox> local_geo;
  std::shared_ptr<CellStructure> cell_structure;
  std::shared_ptr<InteractionsNonBonded> nonbonded_ias;

protected:
  /** @brief Molecular dynamics integrator force capping. */
  double force_cap;
  /**
   * @brief Minimal global interaction cutoff.
   * Particles with a distance smaller than this are guaranteed
   * to be available on the same node (through ghosts).
   */
  double min_global_cut;
};

System &get_system();
void set_system(std::shared_ptr<System> new_instance);
void reset_system();
bool is_system_set();

} // namespace System
