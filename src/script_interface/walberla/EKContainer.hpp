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

#include "EKFFT.hpp"
#include "EKNone.hpp"
#include "EKReactions.hpp"
#include "EKSpecies.hpp"

#include <walberla_bridge/electrokinetics/EKContainer.hpp>
#include <walberla_bridge/electrokinetics/EKinWalberlaBase.hpp>

#include "core/communication.hpp"
#include "core/ek/EKWalberla.hpp"
#include "core/ek/Solver.hpp"
#include "core/system/System.hpp"

#include <script_interface/ObjectList.hpp>
#include <script_interface/ScriptInterface.hpp>
#include <script_interface/ek/Container.hpp>

#include <cassert>
#include <memory>
#include <optional>
#include <string>
#include <variant>

namespace ScriptInterface::walberla {

class EKContainer : public ObjectList<EKSpecies, EK::Container> {
  using Base = ObjectList<EKSpecies, EK::Container>;
  using Base::value_type;

  std::variant<
#ifdef ESPRESSO_WALBERLA_FFT
      std::shared_ptr<EKFFT>,
#endif
      std::shared_ptr<EKNone>>
      m_poisson_solver;

  std::shared_ptr<EKReactions> m_ek_reactions;
  std::shared_ptr<::EK::EKWalberla> m_ek_instance;
  std::shared_ptr<::EK::EKWalberla::ek_container_type> m_ek_container;
  bool m_is_active;

  auto get_precision(decltype(m_poisson_solver) const &solver) const {
    std::optional<bool> result = std::nullopt;
    std::visit(
        [&](auto const &ptr) {
          using SolverType = std::remove_cvref_t<decltype(ptr)>::element_type;
          if (ptr and not std::is_same_v<SolverType, EKNone>) {
            result = get_value<bool>(ptr->get_parameter("single_precision"));
          }
        },
        solver);
    return result;
  }

  auto get_precision(value_type const &species) const {
    return get_value<bool>(species->get_parameter("single_precision"));
  }

  auto get_precision(std::vector<value_type> const &species_list) const {
    std::optional<bool> result = std::nullopt;
    for (auto const &species : species_list) {
      result = get_precision(species);
    }
    return result;
  }

  bool has_in_core(value_type const &obj_ptr) const override {
    return m_ek_container->contains(obj_ptr->get_ekinstance());
  }
  void add_in_core(value_type const &obj_ptr) override {
    context()->parallel_try_catch([this, &obj_ptr]() {
      auto const prec_new_species = get_precision(obj_ptr);
      auto const prec_old_species = get_precision(elements());
      auto const prec_solver = get_precision(m_poisson_solver);
      if ((prec_solver and *prec_solver != prec_new_species) or
          (prec_old_species and *prec_old_species != prec_new_species)) {
        throw std::runtime_error(
            "Cannot mix single and double precision kernels");
      }
      ek_throw_if_expired(obj_ptr->get_mpi_cart_comm_observer());
      m_ek_container->add(obj_ptr->get_ekinstance());
    });
  }
  void remove_in_core(value_type const &obj_ptr) final {
    m_ek_container->remove(obj_ptr->get_ekinstance());
  }

  struct GetPoissonSolverAsVariant {
    template <typename T>
    auto operator()(std::shared_ptr<T> const &solver) const {
      return (solver) ? Variant{solver} : Variant{none};
    }
  };

  Variant get_solver() const {
    return std::visit(GetPoissonSolverAsVariant(), m_poisson_solver);
  }

  struct GetPoissonSolverCoreInstance {
    template <typename T>
    std::shared_ptr<::walberla::PoissonSolver>
    operator()(std::shared_ptr<T> const &solver) const {
      return solver->get_instance();
    }
  };

  auto extract_solver(Variant const &v) {
    std::optional<decltype(m_poisson_solver)> solver;
    auto so_ptr = get_value<ObjectRef>(v);
    if (auto ptr = std::dynamic_pointer_cast<EKNone>(so_ptr)) {
      solver = std::move(ptr);
    }
#ifdef ESPRESSO_WALBERLA_FFT
    else if (auto ptr = std::dynamic_pointer_cast<EKFFT>(so_ptr)) {
      solver = std::move(ptr);
    }
#endif // ESPRESSO_WALBERLA_FFT
    assert(solver.has_value());
    return *solver;
  }

  void set_solver(Variant const &v) {
    auto new_solver = extract_solver(v);
    auto handle = std::visit(GetPoissonSolverCoreInstance{}, new_solver);
    context()->parallel_try_catch([&]() {
      auto const prec_solver = get_precision(new_solver);
      auto const prec_species = get_precision(elements());
      if (prec_solver and prec_species and *prec_solver != *prec_species) {
        throw std::runtime_error(
            "Cannot mix single and double precision kernels");
      }
      m_ek_container->set_poisson_solver(handle);
    });
    m_poisson_solver = new_solver;
  }

public:
  EKContainer() : Base::ObjectList() {
    add_parameters({
        {"tau", AutoParameter::read_only,
         [this]() { return m_ek_container->get_tau(); }},
        {"solver", [this](Variant const &v) { set_solver(v); },
         [this]() { return get_solver(); }},
        {"reactions", AutoParameter::read_only,
         [this]() { return m_ek_reactions; }},
        {"is_active", AutoParameter::read_only,
         [this]() { return m_is_active; }},
    });
  }

  ~EKContainer() override { do_destruct(); }

  void do_construct(VariantMap const &params) override {
    m_is_active = false;
    auto const tau = get_value<double>(params, "tau");
    context()->parallel_try_catch([tau]() {
      if (tau <= 0.) {
        throw std::domain_error("Parameter 'tau' must be > 0");
      }
    });
    m_poisson_solver = extract_solver(
        params.contains("solver") ? params.at("solver") : Variant{none});
    m_ek_container = std::make_shared<::EK::EKWalberla::ek_container_type>(
        tau, std::visit(GetPoissonSolverCoreInstance{}, m_poisson_solver));
    m_ek_reactions = get_value<decltype(m_ek_reactions)>(params, "reactions");
    m_ek_instance = std::make_shared<::EK::EKWalberla>(
        m_ek_container, m_ek_reactions->get_handle());
    // EK species must be added after tau
    Base::do_construct(params);
  }

protected:
  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "activate") {
      get_system().ek.set<::EK::EKWalberla>(m_ek_instance);
      m_is_active = true;
      return {};
    }
    if (method == "deactivate") {
      get_system().ek.reset();
      m_is_active = false;
      return {};
    }

    return Base::do_call_method(method, parameters);
  }
};

} // namespace ScriptInterface::walberla

#endif // ESPRESSO_WALBERLA
