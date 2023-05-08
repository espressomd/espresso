/*
 * Copyright (C) 2021-2023 The ESPResSo project
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

#ifdef WALBERLA

#include "EKSpecies.hpp"
#include "WalberlaCheckpoint.hpp"

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/EKWalberlaNodeState.hpp>
#include <walberla_bridge/electrokinetics/ek_walberla_init.hpp>

#include <boost/mpi.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/optional.hpp>

#include <algorithm>
#include <cassert>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface::walberla {

std::unordered_map<std::string, int> const EKVTKHandle::obs_map = {
    {"density", static_cast<int>(EKOutputVTK::density)},
};

Variant EKSpecies::do_call_method(std::string const &method,
                                  VariantMap const &parameters) {
  if (method == "update_flux_boundary_from_shape") {
    auto values = get_value<std::vector<double>>(parameters, "values");
    std::transform(values.begin(), values.end(), values.begin(),
                   [this](double v) { return v * m_conv_flux; });

    m_instance->update_flux_boundary_from_shape(
        get_value<std::vector<int>>(parameters, "raster"), values);
    return {};
  }
  if (method == "update_density_boundary_from_shape") {
    auto values = get_value<std::vector<double>>(parameters, "values");
    std::transform(values.begin(), values.end(), values.begin(),
                   [this](double v) { return v * m_conv_density; });
    m_instance->update_density_boundary_from_shape(
        get_value<std::vector<int>>(parameters, "raster"), values);
    return {};
  }
  if (method == "clear_flux_boundaries") {
    m_instance->clear_flux_boundaries();
    return {};
  }
  if (method == "clear_density_boundaries") {
    m_instance->clear_density_boundaries();
    return {};
  }
  if (method == "save_checkpoint") {
    auto const path = get_value<std::string>(parameters, "path");
    auto const mode = get_value<int>(parameters, "mode");
    save_checkpoint(path, mode);
    return {};
  }
  if (method == "load_checkpoint") {
    auto const path = get_value<std::string>(parameters, "path");
    auto const mode = get_value<int>(parameters, "mode");
    load_checkpoint(path, mode);
    return {};
  }
  return Base::do_call_method(method, parameters);
}

void EKSpecies::do_construct(VariantMap const &args) {
  m_lattice = get_value<std::shared_ptr<LatticeWalberla>>(args, "lattice");
  m_vtk_writers =
      get_value_or<decltype(m_vtk_writers)>(args, "vtk_writers", {});
  auto const single_precision = get_value<bool>(args, "single_precision");
  auto const agrid = get_value<double>(m_lattice->get_parameter("agrid"));
  auto const diffusion = get_value<double>(args, "diffusion");
  auto const ext_efield = get_value<Utils::Vector3d>(args, "ext_efield");
  auto const density = get_value<double>(args, "density");
  auto const kT = get_value<double>(args, "kT");
  auto const tau = m_tau = get_value<double>(args, "tau");
  context()->parallel_try_catch([&]() {
    if (tau <= 0.) {
      throw std::domain_error("Parameter 'tau' must be > 0");
    }
    if (kT < 0.) {
      throw std::domain_error("Parameter 'kT' must be >= 0");
    }
    if (density < 0.) {
      throw std::domain_error("Parameter 'density' must be >= 0");
    }
    m_conv_energy = Utils::int_pow<2>(tau) / Utils::int_pow<2>(agrid);
    m_conv_diffusion = tau / Utils::int_pow<2>(agrid);
    m_conv_ext_efield = Utils::int_pow<2>(tau) / agrid;
    m_conv_density = Utils::int_pow<3>(agrid);
    m_conv_flux = tau * Utils::int_pow<2>(agrid);
    auto const ek_diffusion = diffusion * m_conv_diffusion;
    auto const ek_ext_efield = ext_efield * m_conv_ext_efield;
    auto const ek_density = m_density = density * m_conv_density;
    auto const ek_kT = kT * m_conv_energy;
    m_instance = new_ek_walberla(
        m_lattice->lattice(), ek_diffusion, ek_kT,
        get_value<double>(args, "valency"), ek_ext_efield, ek_density,
        get_value<bool>(args, "advection"),
        get_value<bool>(args, "friction_coupling"), single_precision);
    for (auto &vtk : m_vtk_writers) {
      vtk->attach_to_lattice(m_instance, get_latice_to_md_units_conversion());
    }
  });
}

void EKSpecies::load_checkpoint(std::string const &filename, int mode) {
  auto &ek_obj = *m_instance;

  auto const read_metadata = [&ek_obj](CheckpointFile &cpfile) {
    auto const expected_grid_size = ek_obj.get_lattice().get_grid_dimensions();
    Utils::Vector3i read_grid_size;
    cpfile.read(read_grid_size);
    if (read_grid_size != expected_grid_size) {
      std::stringstream message;
      message << "grid dimensions mismatch, "
              << "read [" << read_grid_size << "], "
              << "expected [" << expected_grid_size << "].";
      throw std::runtime_error(message.str());
    }
  };

  auto const read_data = [&ek_obj](CheckpointFile &cpfile) {
    auto const grid_size = ek_obj.get_lattice().get_grid_dimensions();
    auto const i_max = grid_size[0];
    auto const j_max = grid_size[1];
    auto const k_max = grid_size[2];
    EKWalberlaNodeState cpnode;
    for (int i = 0; i < i_max; i++) {
      for (int j = 0; j < j_max; j++) {
        for (int k = 0; k < k_max; k++) {
          auto const ind = Utils::Vector3i{{i, j, k}};
          cpfile.read(cpnode.density);
          cpfile.read(cpnode.is_boundary_density);
          if (cpnode.is_boundary_density) {
            cpfile.read(cpnode.density_boundary);
          }
          cpfile.read(cpnode.is_boundary_flux);
          if (cpnode.is_boundary_flux) {
            cpfile.read(cpnode.flux_boundary);
          }
          ek_obj.set_node_density(ind, cpnode.density);
          if (cpnode.is_boundary_density) {
            ek_obj.set_node_density_boundary(ind, cpnode.density_boundary);
          }
          if (cpnode.is_boundary_flux) {
            ek_obj.set_node_flux_boundary(ind, cpnode.flux_boundary);
          }
        }
      }
    }
  };

  auto const on_success = [&ek_obj]() { ek_obj.ghost_communication(); };

  load_checkpoint_common(*context(), "EK", filename, mode, read_metadata,
                         read_data, on_success);
}

void EKSpecies::save_checkpoint(std::string const &filename, int mode) {
  auto &ek_obj = *m_instance;

  auto const write_metadata = [&ek_obj,
                               mode](std::shared_ptr<CheckpointFile> cpfile_ptr,
                                     Context const &context) {
    auto const grid_size = ek_obj.get_lattice().get_grid_dimensions();
    if (context.is_head_node()) {
      cpfile_ptr->write(grid_size);
      unit_test_handle(mode);
    }
  };

  auto const on_failure = [](std::shared_ptr<CheckpointFile> const &,
                             Context const &context) {
    if (context.is_head_node()) {
      auto failure = true;
      boost::mpi::broadcast(context.get_comm(), failure, 0);
    }
  };

  auto const write_data = [&ek_obj,
                           mode](std::shared_ptr<CheckpointFile> cpfile_ptr,
                                 Context const &context) {
    auto const get_node_checkpoint = [&](Utils::Vector3i const &ind)
        -> boost::optional<EKWalberlaNodeState> {
      auto const density = ek_obj.get_node_density(ind);
      auto const is_b_d = ek_obj.get_node_is_density_boundary(ind);
      auto const dens_b = ek_obj.get_node_density_at_boundary(ind);
      auto const is_b_f = ek_obj.get_node_is_flux_boundary(ind);
      auto const flux_b = ek_obj.get_node_flux_at_boundary(ind);
      if (density and is_b_d and is_b_f and
          ((*is_b_d) ? dens_b.has_value() : true) and
          ((*is_b_f) ? flux_b.has_value() : true)) {
        EKWalberlaNodeState cpnode;
        cpnode.density = *density;
        cpnode.is_boundary_density = *is_b_d;
        if (*is_b_d) {
          cpnode.density_boundary = *dens_b;
        }
        cpnode.is_boundary_flux = *is_b_f;
        if (*is_b_f) {
          cpnode.flux_boundary = *flux_b;
        }
        return cpnode;
      }
      return {boost::none};
    };

    auto failure = false;
    auto const &comm = context.get_comm();
    auto const is_head_node = context.is_head_node();
    auto const unit_test_mode = (mode != static_cast<int>(CptMode::ascii)) and
                                (mode != static_cast<int>(CptMode::binary));
    auto const grid_size = ek_obj.get_lattice().get_grid_dimensions();
    auto const i_max = grid_size[0];
    auto const j_max = grid_size[1];
    auto const k_max = grid_size[2];
    EKWalberlaNodeState cpnode;
    for (int i = 0; i < i_max; i++) {
      for (int j = 0; j < j_max; j++) {
        for (int k = 0; k < k_max; k++) {
          auto const ind = Utils::Vector3i{{i, j, k}};
          auto const result = get_node_checkpoint(ind);
          if (!unit_test_mode) {
            assert(1 == boost::mpi::all_reduce(comm, static_cast<int>(!!result),
                                               std::plus<>()) &&
                   "Incorrect number of return values");
          }
          if (is_head_node) {
            if (result) {
              cpnode = *result;
            } else {
              comm.recv(boost::mpi::any_source, 42, cpnode);
            }
            auto &cpfile = *cpfile_ptr;
            cpfile.write(cpnode.density);
            cpfile.write(cpnode.is_boundary_density);
            if (cpnode.is_boundary_density) {
              cpfile.write(cpnode.density_boundary);
            }
            cpfile.write(cpnode.is_boundary_flux);
            if (cpnode.is_boundary_flux) {
              cpfile.write(cpnode.flux_boundary);
            }
            boost::mpi::broadcast(comm, failure, 0);
          } else {
            if (result) {
              comm.send(0, 42, *result);
            }
            boost::mpi::broadcast(comm, failure, 0);
            if (failure) {
              return;
            }
          }
        }
      }
    }
  };

  save_checkpoint_common(*context(), "EK", filename, mode, write_metadata,
                         write_data, on_failure);
}

} // namespace ScriptInterface::walberla

#endif // WALBERLA
