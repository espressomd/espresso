/*
 * Copyright (C) 2022 The ESPResSo project
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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_WALBERLA_EK_REACTION_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_WALBERLA_EK_REACTION_HPP

#include "config/config.hpp"

#ifdef WALBERLA

#include "EKReactant.hpp"
#include "LatticeIndices.hpp"
#include "LatticeWalberla.hpp"

#include "walberla_bridge/electrokinetics/reactions/EKReactionBase.hpp"
#include "walberla_bridge/src/electrokinetics/reactions/EKReactionImplBulk.hpp"
#include "walberla_bridge/src/electrokinetics/reactions/EKReactionImplIndexed.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameter.hpp"
#include "script_interface/communication.hpp"

#include <algorithm>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface::walberla {

class EKReaction : public AutoParameters<EKReaction, LatticeIndices> {
public:
  [[nodiscard]] std::shared_ptr<::walberla::EKReactionBase>
  get_instance() const {
    return m_ekreaction;
  }

protected:
  template <typename T>
  std::shared_ptr<T> make_instance(VariantMap const &args) const {
    auto lattice =
        get_value<std::shared_ptr<LatticeWalberla>>(args, "lattice")->lattice();

    auto reactant = get_value<std::vector<Variant>>(args, "reactants");
    auto output =
        std::vector<std::shared_ptr<::walberla::EKReactant>>(reactant.size());
    auto get_instance = [](Variant const &v) {
      return get_value<std::shared_ptr<EKReactant>>(v)->get_instance();
    };
    std::transform(reactant.begin(), reactant.end(), output.begin(),
                   get_instance);

    return std::make_shared<T>(lattice, output,
                               get_value<double>(args, "coefficient"));
  }

  std::shared_ptr<::walberla::EKReactionBase> m_ekreaction;
};

class EKBulkReaction : public EKReaction {
public:
  EKBulkReaction() {
    add_parameters({{"coefficient",
                     [this](Variant const &v) {
                       get_instance()->set_coefficient(get_value<double>(v));
                     },
                     [this]() { return get_instance()->get_coefficient(); }}});
  }

  void do_construct(VariantMap const &args) override {
    m_ekreaction = make_instance<::walberla::EKReactionImplBulk>(args);
  }
};

class EKIndexedReaction : public EKReaction {
public:
  EKIndexedReaction() {
    add_parameters(
        {{"coefficient",
          [this](Variant const &v) {
            get_instance()->set_coefficient(get_value<double>(v));
          },
          [this]() { return get_instance()->get_coefficient(); }},
         {"shape", AutoParameter::read_only, [this]() {
            return get_instance()->get_lattice()->get_grid_dimensions();
          }}});
  }

  void do_construct(VariantMap const &args) override {
    m_ekreaction = make_instance<::walberla::EKReactionImplIndexed>(args);
    m_ekreaction_impl =
        std::dynamic_pointer_cast<::walberla::EKReactionImplIndexed>(
            get_instance());
  }

  [[nodiscard]] Variant do_call_method(std::string const &method,
                                       VariantMap const &parameters) override {
    if (method == "set_node_is_boundary") {
      auto const index = get_mapped_index(
          get_value<Utils::Vector3i>(parameters, "node"),
          get_instance()->get_lattice()->get_grid_dimensions());
      m_ekreaction_impl->set_node_is_boundary(
          index, get_value<bool>(parameters, "is_boundary"));
      return none;
    }
    if (method == "get_node_is_boundary") {
      auto const index = get_mapped_index(
          get_value<Utils::Vector3i>(parameters, "node"),
          get_instance()->get_lattice()->get_grid_dimensions());
      auto const result = m_ekreaction_impl->get_node_is_boundary(index);
      return mpi_reduce_optional(context()->get_comm(), result);
    }
    return {};
  }

private:
  std::shared_ptr<::walberla::EKReactionImplIndexed> m_ekreaction_impl;
};

} // namespace ScriptInterface::walberla

#endif // WALBERLA
#endif
