/*
 * Copyright (C) 2026 The ESPResSo project
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

#ifdef ESPRESSO_H5MD

#include "Shape.hpp"

#include "core/BoxGeometry.hpp"
#include "core/LocalBox.hpp"
#include "core/io/voxel_mask/VoxelMaskFile.hpp"
#include "core/system/System.hpp"

#include <shapes/VoxelMask.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/collectives/gather.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/operations.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ScriptInterface {
namespace Shapes {

namespace detail {

/**
 * @brief Throw on all ranks if any rank has an error.
 * All ranks throw the message of the lowest failing rank.
 */
inline void throw_on_any_rank(boost::mpi::communicator const &comm,
                              std::optional<std::string> const &error) {
  auto const failing = boost::mpi::all_reduce(
      comm, error ? comm.rank() : comm.size(), boost::mpi::minimum<int>());
  if (failing == comm.size()) {
    return;
  }
  auto message = error.value_or(std::string{});
  boost::mpi::broadcast(comm, message, failing);
  throw std::runtime_error(message);
}

/** Local geometry of the last load. */
struct VoxelMaskLoadState {
  Utils::Vector3d my_left;
  Utils::Vector3d my_right;
  Utils::Vector3d box_l;
  std::array<bool, 3> periodic;
  bool operator==(VoxelMaskLoadState const &) const = default;
};

/**
 * @brief Read the stored region of this rank. Collective.
 * Sets the dataset shape on the first load, and checks it afterwards.
 */
inline void load_voxel_mask(::Shapes::VoxelMask &mask, std::string const &path,
                            std::string const &dataset,
                            boost::mpi::communicator const &comm,
                            VoxelMaskLoadState const &state) {
  IO::VoxelMaskFile const file(path, dataset, comm);
  auto g = mask.geometry();
  if (g.mask_shape == Utils::Vector3i{}) {
    g.mask_shape = file.shape();
    mask.set_geometry(g);
  } else if (file.shape() != g.mask_shape) {
    throw std::runtime_error("The shape of dataset '" + dataset +
                             "' in file '" + path + "' has changed");
  }
  auto const region =
      ::Shapes::VoxelMask::stored_region(g, state.my_left, state.my_right);
  IO::read_voxel_mask(file, mask, region, comm);
}

/**
 * @brief Make a loader that reads the stored region of this rank.
 * The loader is collective. It re-reads the file only if the local box
 * or the box geometry changed on any rank since the last load.
 */
inline auto make_voxel_mask_loader(std::string path, std::string dataset,
                                   boost::mpi::communicator comm) {
  auto last = std::make_shared<std::optional<VoxelMaskLoadState>>();
  return [path = std::move(path), dataset = std::move(dataset),
          comm = std::move(comm), last](::Shapes::VoxelMask &mask) {
    auto const &system = ::System::get_system();
    auto const &box_geo = *system.box_geo;
    auto const &local_geo = *system.local_geo;
    VoxelMaskLoadState const state{
        local_geo.my_left(),
        local_geo.my_right(),
        box_geo.length(),
        {{box_geo.periodic(0u), box_geo.periodic(1u), box_geo.periodic(2u)}}};
    auto const changed = not last->has_value() or * *last != state;
    if (not boost::mpi::all_reduce(comm, changed, std::logical_or<>())) {
      return;
    }
    auto g = mask.geometry();
    g.box_l = state.box_l;
    g.periodic = state.periodic;
    // discard the stale data; after a failure, the next call retries
    mask.set_geometry(g);
    last->reset();
    std::optional<std::string> error;
    try {
      load_voxel_mask(mask, path, dataset, comm, state);
    } catch (std::exception const &err) {
      error = err.what();
    }
    if (error) {
      mask.set_geometry(g);
    }
    // all ranks agree on the outcome
    throw_on_any_rank(comm, error);
    *last = state;
  };
}

} // namespace detail

class VoxelMask : public Shape {
public:
  VoxelMask() : m_shape(std::make_shared<::Shapes::VoxelMask>()) {
    add_parameters(
        {{"file", AutoParameter::read_only, [this]() { return m_file; }},
         {"dataset", AutoParameter::read_only, [this]() { return m_dataset; }},
         {"voxel_size", AutoParameter::read_only,
          [this]() { return m_shape->geometry().voxel_size; }},
         {"origin", AutoParameter::read_only,
          [this]() { return m_shape->geometry().origin; }},
         {"halo", AutoParameter::read_only,
          [this]() { return m_shape->geometry().halo; }},
         {"mask_shape", AutoParameter::read_only,
          [this]() { return m_shape->geometry().mask_shape; }}});
  }

  std::shared_ptr<::Shapes::Shape> shape() const override { return m_shape; }

  void do_construct(VariantMap const &params) override {
    context()->parallel_try_catch([&]() {
      m_file = get_value<std::string>(params, "file");
      m_dataset = get_value<std::string>(params, "dataset");
      ::Shapes::VoxelMask::Geometry g;
      g.voxel_size = get_value<double>(params, "voxel_size");
      g.origin = get_value_or<Utils::Vector3d>(params, "origin", {});
      g.halo = get_value_or<int>(params, "halo", 2);
      if (not(g.voxel_size > 0.)) {
        throw std::domain_error("Parameter 'voxel_size' must be > 0");
      }
      if (g.halo < 1) {
        throw std::domain_error("Parameter 'halo' must be >= 1");
      }
      auto const &comm = context()->get_comm();
      // the loader opens the file and sets the dataset shape
      auto const &box_geo = *::System::get_system().box_geo;
      g.box_l = box_geo.length();
      for (unsigned int a = 0u; a < 3u; ++a) {
        g.periodic[a] = box_geo.periodic(a);
      }
      m_shape->set_geometry(g);
      m_shape->set_loader(
          detail::make_voxel_mask_loader(m_file, m_dataset, comm));
      m_shape->on_local_geometry_change();
    });
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    context()->parallel_try_catch(
        [this]() { m_shape->on_local_geometry_change(); });
    auto const &comm = context()->get_comm();

    if (name == "calc_distance" or name == "is_inside") {
      context()->parallel_try_catch([&]() {
        if (not params.contains("position")) {
          throw std::runtime_error("Parameter 'position' is missing");
        }
      });
      auto const &box_geo = *::System::get_system().box_geo;
      auto const pos = box_geo.folded_position(
          get_value<Utils::Vector3d>(params.at("position")));
      auto const owner = owner_rank(pos);
      // result: distance, vector, inside
      std::array<double, 5> result{};
      context()->parallel_try_catch([&]() {
        if (comm.rank() != owner) {
          return;
        }
        if (name == "calc_distance") {
          double dist;
          Utils::Vector3d vec;
          m_shape->calculate_dist(pos, dist, vec);
          result = {dist, vec[0], vec[1], vec[2], 0.};
        } else {
          result[4] = m_shape->is_inside(pos) ? 1. : 0.;
        }
      });
      boost::mpi::broadcast(comm, result.data(),
                            static_cast<int>(result.size()), owner);
      if (comm.rank() != 0) {
        return {};
      }
      if (name == "calc_distance") {
        return std::vector<Variant>{
            result[0], Utils::Vector3d{result[1], result[2], result[3]}};
      }
      return {result[4] != 0.};
    }

    if (name == "rasterize") {
      context()->parallel_try_catch([]() {
        throw std::runtime_error("VoxelMask cannot be rasterized globally");
      });
      return {};
    }

    if (name == "get_local_voxel_counts") {
      return gather_on_head(m_shape->local_voxel_count());
    }

    if (name == "get_local_memory_bytes") {
      return gather_on_head(m_shape->local_memory_bytes());
    }

    return Shape::do_call_method(name, params);
  }

private:
  /** @brief Rank whose local box contains a folded position. */
  int owner_rank(Utils::Vector3d const &folded_pos) const {
    auto const &comm = context()->get_comm();
    auto const &system = ::System::get_system();
    auto const &box_geo = *system.box_geo;
    auto const &local_geo = *system.local_geo;
    auto folded = folded_pos;
    auto is_local = true;
    for (unsigned int a = 0u; a < 3u; ++a) {
      if (not box_geo.periodic(a)) {
        // positions beyond a non-periodic boundary go to the boundary rank
        folded[a] = std::clamp(folded[a], 0., box_geo.length()[a]);
      }
      // boundaries are inclusive to close floating-point gaps between ranks
      is_local &= local_geo.my_left()[a] <= folded[a] and
                  folded[a] <= local_geo.my_right()[a];
    }
    auto const owner = boost::mpi::all_reduce(
        comm, is_local ? comm.rank() : comm.size(), boost::mpi::minimum<int>());
    context()->parallel_try_catch([&]() {
      if (owner == comm.size()) {
        throw std::runtime_error("No MPI rank owns the position");
      }
    });
    return owner;
  }

  Variant gather_on_head(std::size_t value) const {
    auto const &comm = context()->get_comm();
    std::vector<std::size_t> values;
    boost::mpi::gather(comm, value, values, 0);
    if (comm.rank() != 0) {
      return {};
    }
    return make_vector_of_variants(values);
  }

  std::shared_ptr<::Shapes::VoxelMask> m_shape;
  std::string m_file;
  std::string m_dataset;
};

} /* namespace Shapes */
} /* namespace ScriptInterface */

#endif // ESPRESSO_H5MD
