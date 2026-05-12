/*
 * Copyright (C) 2022-2026 The ESPResSo project
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

#include <config/config.hpp>

#ifdef ESPRESSO_WALBERLA

#include "EKPoissonSolver.hpp"
#include "WalberlaCheckpoint.hpp"

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/ek_walberla_init.hpp>

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameters.hpp>
#include <script_interface/code_info/CodeInfo.hpp>

#include <boost/mpi.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>

#include <cassert>
#include <filesystem>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

namespace ScriptInterface::walberla {

class EKNone : public EKPoissonSolver {
  std::shared_ptr<::walberla::PoissonSolver> m_instance;
  std::shared_ptr<LatticeWalberla> m_lattice;
  bool m_gpu;
  bool m_single_precision;

protected:
  void make_instance(VariantMap const &args) override {

    auto *make_new_instance = &::walberla::new_ek_poisson_none;
    if (m_gpu) {
      std::vector<std::string> required_features;
      required_features.emplace_back("CUDA");
      CodeInfo::check_features(required_features);
#ifdef ESPRESSO_CUDA
      make_new_instance = &::walberla::new_ek_poisson_none_cuda;
#endif
    }
    m_instance = make_new_instance(m_lattice->lattice(), m_single_precision);
  }

public:
  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "save_checkpoint") {
      auto const path = get_value<std::filesystem::path>(parameters, "path");
      auto const mode = get_value<int>(parameters, "mode");
      save_checkpoint(path, mode);
      return {};
    }
    if (method == "load_checkpoint") {
      auto const path = get_value<std::filesystem::path>(parameters, "path");
      auto const mode = get_value<int>(parameters, "mode");
      load_checkpoint(path, mode);
      return {};
    }
    return EKPoissonSolver::do_call_method(method, parameters);
  }

  void do_construct(VariantMap const &args) override {
    m_gpu = get_value_or<bool>(args, "gpu", false);
    m_single_precision = get_value_or<bool>(args, "single_precision", m_gpu);
    m_lattice = get_value<decltype(m_lattice)>(args, "lattice");
    auto const agrid = get_value<double>(m_lattice->get_parameter("agrid"));
    m_tau = get_value<double>(args, "tau");
    set_potential_conversion(agrid, m_tau);

    make_instance(args);
    add_parameters({
        {"tau", AutoParameter::read_only, [this]() { return m_tau; }},
        {"single_precision", AutoParameter::read_only,
         [this]() { return m_single_precision; }},
        {"gpu", AutoParameter::read_only, [this]() { return m_gpu; }},
        {"lattice", AutoParameter::read_only, [this]() { return m_lattice; }},
        {"shape", AutoParameter::read_only,
         [this]() { return m_instance->get_lattice().get_grid_dimensions(); }},
    });
  }

  [[nodiscard]] std::shared_ptr<::walberla::PoissonSolver>
  get_instance() const noexcept override {
    return m_instance;
  }

private:
  void load_checkpoint(std::filesystem::path const &path, int mode) {
    auto &solver = *m_instance;

    auto const read_metadata = [&solver](CheckpointFile &cpfile) {
      auto const expected_grid_size =
          solver.get_lattice().get_grid_dimensions();
      Utils::Vector3i read_grid_size;
      cpfile.read(read_grid_size);
      if (read_grid_size != expected_grid_size) {
        std::stringstream message;
        message << "grid dimensions mismatch, read [" << read_grid_size << "], "
                << "expected [" << expected_grid_size << "].";
        throw std::runtime_error(message.str());
      }
    };

    auto const read_data = [&solver](CheckpointFile &cpfile) {
      auto const grid_size = solver.get_lattice().get_grid_dimensions();
      auto const i_max = grid_size[0];
      auto const j_max = grid_size[1];
      auto const k_max = grid_size[2];
      for (int i = 0; i < i_max; i++) {
        for (int j = 0; j < j_max; j++) {
          for (int k = 0; k < k_max; k++) {
            auto const ind = Utils::Vector3i{{i, j, k}};
            double potential{};
            cpfile.read(potential);
            static_cast<void>(solver.set_node_potential(ind, potential));
          }
        }
      }
    };

    auto const on_success = []() {};
    load_checkpoint_common(*context(), "EKNone", path, mode, read_metadata,
                           read_data, on_success);
  }

  void save_checkpoint(std::filesystem::path const &path, int mode) {
    auto &solver = *m_instance;

    auto const write_metadata =
        [&solver, mode](std::shared_ptr<CheckpointFile> cpfile_ptr,
                        Context const &context) {
          auto const grid_size = solver.get_lattice().get_grid_dimensions();
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

    auto const write_data = [&solver,
                             mode](std::shared_ptr<CheckpointFile> cpfile_ptr,
                                   Context const &context) {
      auto const get_node_checkpoint =
          [&solver](Utils::Vector3i const &ind) -> std::optional<double> {
        return solver.get_node_potential(ind);
      };

      auto failure = false;
      auto const &comm = context.get_comm();
      auto const is_head_node = context.is_head_node();
      auto const unit_test_mode = (mode != static_cast<int>(CptMode::ascii)) and
                                  (mode != static_cast<int>(CptMode::binary));
      auto const grid_size = solver.get_lattice().get_grid_dimensions();
      auto const i_max = grid_size[0];
      auto const j_max = grid_size[1];
      auto const k_max = grid_size[2];
      double cpnode = 0.;
      for (int i = 0; i < i_max; i++) {
        for (int j = 0; j < j_max; j++) {
          for (int k = 0; k < k_max; k++) {
            auto const ind = Utils::Vector3i{{i, j, k}};
            auto const result = get_node_checkpoint(ind);
            if (!unit_test_mode) {
              assert(1 == boost::mpi::all_reduce(comm,
                                                 static_cast<int>(!!result),
                                                 std::plus<>()) &&
                     "Incorrect number of return values");
            }
            if (is_head_node) {
              if (result) {
                cpnode = *result;
              } else {
                comm.recv(boost::mpi::any_source, 42, cpnode);
              }
              cpfile_ptr->write(cpnode);
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

    save_checkpoint_common(*context(), "EKNone", path, mode, write_metadata,
                           write_data, on_failure);
  }
};

} // namespace ScriptInterface::walberla

#endif // ESPRESSO_WALBERLA
