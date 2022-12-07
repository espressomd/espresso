/*
 * Copyright (C) 2021-2022 The ESPResSo project
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
#ifndef SCRIPT_INTERFACE_WALBERLA_EKSPECIESNODE_HPP
#define SCRIPT_INTERFACE_WALBERLA_EKSPECIESNODE_HPP

#include "config/config.hpp"

#ifdef WALBERLA

#include <walberla_bridge/electrokinetics/EKinWalberlaBase.hpp>

#include "LatticeIndices.hpp"

#include "core/errorhandling.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/communication.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/optional.hpp>
#include <boost/serialization/vector.hpp>

#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>

namespace ScriptInterface::walberla {

class EKSpeciesNode : public AutoParameters<EKSpeciesNode, LatticeIndices> {
  std::shared_ptr<::EKinWalberlaBase> m_ek_species;
  Utils::Vector3i m_index;
  Utils::Vector3i m_grid_size;
  double m_conv_dens;
  double m_conv_flux;
  double m_conv_velocity;

public:
  EKSpeciesNode() {
    add_parameters(
        {{"_index", AutoParameter::read_only, [this]() { return m_index; }}});
  }

  void do_construct(VariantMap const &params) override {
    try {
      auto const ek_sip =
          get_value<std::shared_ptr<EKSpecies>>(params, "parent_sip");
      m_ek_species = ek_sip->get_ekinstance();
      assert(m_ek_species);
      m_conv_dens = 1.;
      m_conv_flux = 1.;
      m_conv_velocity = 1.;
    } catch (std::exception const &e) {
      if (context()->is_head_node()) {
        runtimeErrorMsg() << "EKSpeciesNode failed: " << e.what();
      }
      m_ek_species.reset();
      return;
    }
    m_grid_size = m_ek_species->get_lattice().get_grid_dimensions();
    m_index = get_mapped_index(get_value<Utils::Vector3i>(params, "index"),
                               m_grid_size);
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "override_index") {
      // this hidden feature is used to iterate an EK slice without
      // rebuilding a EKSpeciesNode for each node in the slice
      auto const index = get_value<Utils::Vector3i>(params, "index");
      if (not is_index_valid(index, m_grid_size)) {
        return ES_ERROR;
      }
      m_index = index;
      return ES_OK;
    }
    if (name == "set_density") {
      auto const dens = get_value<double>(params, "value");
      m_ek_species->set_node_density(m_index, dens * m_conv_dens);
      m_ek_species->ghost_communication();
      return {};
    }
    if (name == "get_density") {
      auto const result = m_ek_species->get_node_density(m_index);
      return mpi_reduce_optional(context()->get_comm(), result) / m_conv_dens;
    }
    if (name == "get_is_boundary") {
      auto const result = m_ek_species->get_node_is_boundary(m_index);
      return mpi_reduce_optional(context()->get_comm(), result);
    }
    if (name == "get_node_density_at_boundary") {
      if (is_boundary_all_reduce(
              m_ek_species->get_node_is_density_boundary(m_index))) {
        auto const result = m_ek_species->get_node_density_at_boundary(m_index);
        return mpi_reduce_optional(context()->get_comm(), result) / m_conv_dens;
      }
      return Variant{None{}};
    }
    if (name == "set_node_density_at_boundary") {
      if (is_none(params.at("value"))) {
        m_ek_species->remove_node_from_density_boundary(m_index);
      } else {
        auto const dens = get_value<double>(params, "value") * m_conv_dens;
        m_ek_species->set_node_density_boundary(m_index, dens);
      }
      return {};
    }
    if (name == "get_node_flux_at_boundary") {
      if (is_boundary_all_reduce(
              m_ek_species->get_node_is_flux_boundary(m_index))) {
        auto const result = m_ek_species->get_node_flux_at_boundary(m_index);
        return mpi_reduce_optional(context()->get_comm(), result) / m_conv_flux;
      }
      return Variant{None{}};
    }
    if (name == "set_node_flux_at_boundary") {
      if (is_none(params.at("value"))) {
        m_ek_species->remove_node_from_flux_boundary(m_index);
      } else {
        auto const flux =
            get_value<Utils::Vector3d>(params, "value") * m_conv_flux;
        m_ek_species->set_node_flux_boundary(m_index, flux);
      }
      return {};
    }
    if (name == "get_lattice_speed") {
      return 1. / m_conv_velocity;
    }

    return {};
  }

private:
  bool is_boundary_all_reduce(boost::optional<bool> const &is_boundary) const {
    return boost::mpi::all_reduce(context()->get_comm(),
                                  is_boundary ? *is_boundary : false,
                                  std::logical_or<>());
  }
};
} // namespace ScriptInterface::walberla

#endif // WALBERLA
#endif
