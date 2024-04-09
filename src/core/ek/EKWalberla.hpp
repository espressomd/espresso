/*
 * Copyright (C) 2023 The ESPResSo project
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

#ifdef WALBERLA

#include <utils/Vector.hpp>

#include <memory>
#include <stdexcept>
#include <utility>

// forward declarations
class EKinWalberlaBase;
template <class Base> class EKContainer;
template <class Base> class EKReactions;
namespace walberla {
class EKReactionBase;
}
namespace System {
class System;
}

namespace EK {

struct EKWalberla {
  using ek_container_type = EKContainer<EKinWalberlaBase>;
  using ek_reactions_type = EKReactions<walberla::EKReactionBase>;
  std::shared_ptr<ek_container_type> ek_container;
  std::shared_ptr<ek_reactions_type> ek_reactions;

  EKWalberla(std::shared_ptr<ek_container_type> ek_container_instance,
             std::shared_ptr<ek_reactions_type> ek_reactions_instance)
      : ek_container{std::move(ek_container_instance)},
        ek_reactions{std::move(ek_reactions_instance)} {}

  double get_tau() const;
  void veto_time_step(double time_step) const;
  void veto_kT(double kT) const;
  void sanity_checks(System::System const &system) const;
  bool is_ready_for_propagation() const noexcept;
  void propagate();
  void perform_reactions();

  void on_cell_structure_change() const {}
  void veto_boxl_change() const {
    throw std::runtime_error("MD cell geometry change not supported by EK");
  }
  void on_boxl_change() const { veto_boxl_change(); }
  void on_node_grid_change() const {
    throw std::runtime_error("MPI topology change not supported by EK");
  }
  void on_timestep_change() const {}
  void on_temperature_change() const {}
};

} // namespace EK

#endif // WALBERLA
