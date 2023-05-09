/*
 * Copyright (C) 2022-2023 The ESPResSo project
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

#include "EKFFT.hpp"
#include "EKNone.hpp"
#include "EKSpecies.hpp"

#include "core/grid_based_algorithms/ek_container.hpp"

#include <script_interface/ObjectList.hpp>
#include <script_interface/ScriptInterface.hpp>

#include <memory>
#include <string>

namespace ScriptInterface::walberla {

class EKContainer : public ObjectList<EKSpecies> {
  using Base = ObjectList<EKSpecies>;

  boost::variant<
#ifdef WALBERLA_FFT
      std::shared_ptr<EKFFT>,
#endif
      std::shared_ptr<EKNone>>
      m_poisson_solver{std::shared_ptr<EKNone>()};

  void add_in_core(std::shared_ptr<EKSpecies> const &obj_ptr) override {
    context()->parallel_try_catch(
        [&obj_ptr]() { EK::ek_container.add(obj_ptr->get_ekinstance()); });
  }
  void remove_in_core(std::shared_ptr<EKSpecies> const &obj_ptr) override {
    EK::ek_container.remove(obj_ptr->get_ekinstance());
  }

public:
  EKContainer() : Base::ObjectList() {
    add_parameters({
        {"tau",
         [this](Variant const &v) {
           if (is_none(v)) {
             if (get_value<int>(do_call_method("size", {})) == 0) {
               EK::ek_container.set_tau(0.);
               return;
             }
             context()->parallel_try_catch([]() {
               throw std::domain_error(
                   "Parameter 'tau' is required when container isn't empty");
             });
           }
           auto const tau = get_value<double>(v);
           context()->parallel_try_catch([tau]() {
             if (tau <= 0.) {
               throw std::domain_error("Parameter 'tau' must be > 0");
             }
           });
           EK::ek_container.set_tau(get_value<double>(v));
         },
         []() {
           auto const tau = EK::ek_container.get_tau();
           return (tau == 0.) ? Variant{none} : Variant{tau};
         }},
        {"solver", [this](Variant const &v) { set_solver(v); },
         [this]() { return get_solver(); }},
    });
  }

  void do_construct(VariantMap const &params) override {
    if (params.count("solver")) {
      set_solver(params.at("solver"));
    }
    if (params.count("tau")) {
      do_set_parameter("tau", params.at("tau"));
    }
    // EK species must be added after tau
    Base::do_construct(params);
  }

protected:
  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "is_poisson_solver_set") {
      return EK::ek_container.is_poisson_solver_set();
    }

    return Base::do_call_method(method, parameters);
  }

private:
  struct GetPoissonSolverVariant : public boost::static_visitor<Variant> {
    template <typename T>
    auto operator()(std::shared_ptr<T> const &solver) const {
      return (solver) ? Variant{solver} : Variant{none};
    }
  };

  struct GetPoissonSolverInstance
      : public boost::static_visitor<
            std::shared_ptr<::walberla::PoissonSolver>> {
    template <typename T>
    auto operator()(std::shared_ptr<T> const &solver) const {
      return (solver) ? solver->get_instance()
                      : std::shared_ptr<::walberla::PoissonSolver>();
    }
  };

  Variant get_solver() const {
    auto const visitor = GetPoissonSolverVariant();
    return boost::apply_visitor(visitor, m_poisson_solver);
  }

  void set_solver(Variant const &solver_variant) {
    boost::optional<decltype(m_poisson_solver)> solver;
    if (is_none(solver_variant)) {
      solver = std::shared_ptr<EKNone>();
    } else {
#ifdef WALBERLA_FFT
      try {
        solver = get_value<std::shared_ptr<EKFFT>>(solver_variant);
      } catch (...) {
      }
#endif
      if (not solver) {
        solver = get_value<std::shared_ptr<EKNone>>(solver_variant);
      }
    }
    m_poisson_solver = *solver;
    auto const visitor = GetPoissonSolverInstance();
    auto const instance = boost::apply_visitor(visitor, m_poisson_solver);
    context()->parallel_try_catch(
        [&instance]() { EK::ek_container.set_poisson_solver(instance); });
  }
};

} // namespace ScriptInterface::walberla

#endif // WALBERLA
