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

/** @file
 *  The ScriptInterface counterparts of the non-bonded interactions parameters
 *  structs from the core are defined here.
 */

#pragma once

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/get_value.hpp"

#include "core/nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "core/system/System.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

namespace ScriptInterface {
namespace Interactions {

class NonBondedInteractionHandle;

template <class CoreIA>
class InteractionPotentialInterface
    : public AutoParameters<InteractionPotentialInterface<CoreIA>> {
public:
  using CoreInteraction = CoreIA;
  using CoreInteractionIAParamsPtr = CoreInteraction IA_parameters::*;

protected:
  using BaseClass = AutoParameters<InteractionPotentialInterface<CoreIA>>;
  using BaseClass::context;
  using BaseClass::valid_parameters;

  /** @brief Managed object. */
  std::shared_ptr<CoreInteraction> m_handle;
  /** @brief Handle to the container whose members have to be synchronized. */
  std::weak_ptr<::IA_parameters> m_core_struct;
  /** @brief Handle to the interface used to synchronize data members. */
  std::weak_ptr<NonBondedInteractionHandle *> m_si_struct;
  /** @brief Callback to notify changes to the interaction range. */
  std::weak_ptr<std::function<void()>> m_notify_non_bonded_ia_change;
  /** @brief Pointer to the corresponding member in a handle. */
  virtual CoreInteractionIAParamsPtr get_ptr_offset() const = 0;
  /** @brief Create a new instance using the constructor with range checks. */
  virtual void make_new_instance(VariantMap const &params) = 0;
  /** @brief Which parameter indicates whether the potential is inactive. */
  virtual std::string inactive_parameter() const { return "cutoff"; }
  /** @brief Which magic value indicates the potential is inactive. */
  virtual double get_inactive_cutoff() const { return inactive_cutoff; }

  template <typename T>
  auto make_autoparameter(T CoreInteraction::*ptr, char const *name) {
    return AutoParameter{name, AutoParameter::read_only,
                         [this, ptr]() { return m_handle.get()->*ptr; }};
  }

  std::set<std::string> get_valid_parameters() const {
    auto const vec = valid_parameters();
    return {vec.begin(), vec.end()};
  }

private:
  void check_valid_parameters(VariantMap const &params) const {
    auto const valid_keys = get_valid_parameters();
    for (auto const &key : valid_keys) {
      if (not params.contains(key)) {
        throw std::runtime_error("Parameter '" + key + "' is missing");
      }
    }
    for (auto const &key : std::views::elements<0>(params)) {
      if (not valid_keys.contains(key)) {
        throw std::runtime_error("Parameter '" + key + "' is not recognized");
      }
    }
  }

public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "set_params") {
      context()->parallel_try_catch([this, &params]() {
        check_valid_parameters(params);
        make_new_instance(params);
      });
      // copying the new value to the core may queue a runtime error message,
      // but we can't detect it to roll back to the last valid state
      update_core();
      return {};
    }
    if (name == "deactivate") {
      m_handle = std::make_shared<CoreInteraction>();
      update_core(get_value_or<bool>(params, "notify", true));
      return {};
    }
    return {};
  }

  void do_construct(VariantMap const &params) final {
    if (params.empty()) {
      m_handle = std::make_shared<CoreInteraction>();
    } else {
      if (std::abs(get_value<double>(params, inactive_parameter()) -
                   get_inactive_cutoff()) < 1e-9) {
        m_handle = std::make_shared<CoreInteraction>();
      } else {
        context()->parallel_try_catch([this, &params]() {
          check_valid_parameters(params);
          make_new_instance(params);
        });
      }
    }
  }

  void attach(std::weak_ptr<NonBondedInteractionHandle *> si_struct,
              std::weak_ptr<::IA_parameters> core_struct) {
    m_si_struct = si_struct;
    m_core_struct = core_struct;
  }

  void update_core(bool notify = true);
};

#ifdef ESPRESSO_WCA
class InteractionWCA : public InteractionPotentialInterface<::WCA_Parameters> {
protected:
  CoreInteractionIAParamsPtr get_ptr_offset() const override {
    return &::IA_parameters::wca;
  }

public:
  InteractionWCA() {
    add_parameters({
        make_autoparameter(&CoreInteraction::eps, "epsilon"),
        make_autoparameter(&CoreInteraction::sig, "sigma"),
    });
  }

private:
  std::string inactive_parameter() const override { return "sigma"; }
  double get_inactive_cutoff() const override { return 0.; }

  void make_new_instance(VariantMap const &params) override {
    m_handle = make_shared_from_args<CoreInteraction, double, double>(
        params, "epsilon", "sigma");
  }

public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "get_cutoff") {
      return m_handle.get()->cut;
    }
    return InteractionPotentialInterface<CoreInteraction>::do_call_method(
        name, params);
  }
};
#endif // ESPRESSO_WCA

#ifdef ESPRESSO_LENNARD_JONES
class InteractionLJ : public InteractionPotentialInterface<::LJ_Parameters> {
protected:
  CoreInteractionIAParamsPtr get_ptr_offset() const override {
    return &::IA_parameters::lj;
  }

public:
  InteractionLJ() {
    add_parameters({
        make_autoparameter(&CoreInteraction::eps, "epsilon"),
        make_autoparameter(&CoreInteraction::sig, "sigma"),
        make_autoparameter(&CoreInteraction::cut, "cutoff"),
        make_autoparameter(&CoreInteraction::shift, "shift"),
        make_autoparameter(&CoreInteraction::offset, "offset"),
        make_autoparameter(&CoreInteraction::min, "min"),
    });
  }

private:
  void make_new_instance(VariantMap const &params) override {
    auto new_params = params;
    auto const *shift_string = std::get_if<std::string>(&params.at("shift"));
    if (shift_string != nullptr) {
      if (*shift_string != "auto") {
        throw std::invalid_argument(
            "LJ parameter 'shift' has to be 'auto' or a float");
      }
      new_params["shift"] = 0.;
    }
    m_handle = make_shared_from_args<CoreInteraction, double, double, double,
                                     double, double, double>(
        new_params, "epsilon", "sigma", "cutoff", "offset", "min", "shift");
    if (shift_string != nullptr) {
      m_handle->shift = m_handle->get_auto_shift();
    }
  }
};
#endif // ESPRESSO_LENNARD_JONES

#ifdef ESPRESSO_LENNARD_JONES_GENERIC
class InteractionLJGen
    : public InteractionPotentialInterface<::LJGen_Parameters> {
protected:
  CoreInteractionIAParamsPtr get_ptr_offset() const override {
    return &::IA_parameters::ljgen;
  }

public:
  InteractionLJGen() {
    add_parameters({
        make_autoparameter(&CoreInteraction::eps, "epsilon"),
        make_autoparameter(&CoreInteraction::sig, "sigma"),
        make_autoparameter(&CoreInteraction::cut, "cutoff"),
        make_autoparameter(&CoreInteraction::shift, "shift"),
        make_autoparameter(&CoreInteraction::offset, "offset"),
#ifdef ESPRESSO_LJGEN_SOFTCORE
        make_autoparameter(&CoreInteraction::lambda, "lam"),
        make_autoparameter(&CoreInteraction::softrad, "delta"),
#endif
        make_autoparameter(&CoreInteraction::a1, "e1"),
        make_autoparameter(&CoreInteraction::a2, "e2"),
        make_autoparameter(&CoreInteraction::b1, "b1"),
        make_autoparameter(&CoreInteraction::b2, "b2"),
    });
  }

private:
  void make_new_instance(VariantMap const &params) override {
    auto new_params = params;
    auto const *shift_string = std::get_if<std::string>(&params.at("shift"));
    if (shift_string != nullptr) {
      if (*shift_string != "auto") {
        throw std::invalid_argument(
            "Generic LJ parameter 'shift' has to be 'auto' or a float");
      }
      new_params["shift"] = 0.;
    }
    m_handle = make_shared_from_args<CoreInteraction, double, double, double,
                                     double, double,
#ifdef ESPRESSO_LJGEN_SOFTCORE
                                     double, double,
#endif
                                     double, double, double, double>(
        new_params, "epsilon", "sigma", "cutoff", "shift", "offset",
#ifdef ESPRESSO_LJGEN_SOFTCORE
        "lam", "delta",
#endif
        "e1", "e2", "b1", "b2");
    if (shift_string != nullptr) {
      m_handle->shift = m_handle->get_auto_shift();
    }
  }
};
#endif // ESPRESSO_LENNARD_JONES_GENERIC

#ifdef ESPRESSO_LJCOS
class InteractionLJcos
    : public InteractionPotentialInterface<::LJcos_Parameters> {
protected:
  CoreInteractionIAParamsPtr get_ptr_offset() const override {
    return &::IA_parameters::ljcos;
  }

public:
  InteractionLJcos() {
    add_parameters({
        make_autoparameter(&CoreInteraction::eps, "epsilon"),
        make_autoparameter(&CoreInteraction::sig, "sigma"),
        make_autoparameter(&CoreInteraction::cut, "cutoff"),
        make_autoparameter(&CoreInteraction::offset, "offset"),
    });
  }

private:
  void make_new_instance(VariantMap const &params) override {
    m_handle =
        make_shared_from_args<CoreInteraction, double, double, double, double>(
            params, "epsilon", "sigma", "cutoff", "offset");
  }
};
#endif // ESPRESSO_LJCOS

#ifdef ESPRESSO_LJCOS2
class InteractionLJcos2
    : public InteractionPotentialInterface<::LJcos2_Parameters> {
protected:
  CoreInteractionIAParamsPtr get_ptr_offset() const override {
    return &::IA_parameters::ljcos2;
  }

public:
  InteractionLJcos2() {
    add_parameters({
        make_autoparameter(&CoreInteraction::eps, "epsilon"),
        make_autoparameter(&CoreInteraction::sig, "sigma"),
        make_autoparameter(&CoreInteraction::offset, "offset"),
        make_autoparameter(&CoreInteraction::w, "width"),
    });
  }

private:
  std::string inactive_parameter() const override { return "sigma"; }
  double get_inactive_cutoff() const override { return 0.; }

  void make_new_instance(VariantMap const &params) override {
    m_handle =
        make_shared_from_args<CoreInteraction, double, double, double, double>(
            params, "epsilon", "sigma", "offset", "width");
  }

public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "get_cutoff") {
      return m_handle.get()->cut;
    }
    return InteractionPotentialInterface<CoreInteraction>::do_call_method(
        name, params);
  }
};
#endif // ESPRESSO_LJCOS2

#ifdef ESPRESSO_HERTZIAN
class InteractionHertzian
    : public InteractionPotentialInterface<::Hertzian_Parameters> {
protected:
  CoreInteractionIAParamsPtr get_ptr_offset() const override {
    return &::IA_parameters::hertzian;
  }

public:
  InteractionHertzian() {
    add_parameters({
        make_autoparameter(&CoreInteraction::eps, "eps"),
        make_autoparameter(&CoreInteraction::sig, "sig"),
    });
  }

private:
  std::string inactive_parameter() const override { return "sig"; }

  void make_new_instance(VariantMap const &params) override {
    m_handle = make_shared_from_args<CoreInteraction, double, double>(
        params, "eps", "sig");
  }
};
#endif // ESPRESSO_HERTZIAN

#ifdef ESPRESSO_GAUSSIAN
class InteractionGaussian
    : public InteractionPotentialInterface<::Gaussian_Parameters> {
protected:
  CoreInteractionIAParamsPtr get_ptr_offset() const override {
    return &::IA_parameters::gaussian;
  }

public:
  InteractionGaussian() {
    add_parameters({
        make_autoparameter(&CoreInteraction::eps, "eps"),
        make_autoparameter(&CoreInteraction::sig, "sig"),
        make_autoparameter(&CoreInteraction::cut, "cutoff"),
    });
  }

private:
  void make_new_instance(VariantMap const &params) override {
    m_handle = make_shared_from_args<CoreInteraction, double, double, double>(
        params, "eps", "sig", "cutoff");
  }
};
#endif // ESPRESSO_GAUSSIAN

#ifdef ESPRESSO_BMHTF_NACL
class InteractionBMHTF
    : public InteractionPotentialInterface<::BMHTF_Parameters> {
protected:
  CoreInteractionIAParamsPtr get_ptr_offset() const override {
    return &::IA_parameters::bmhtf;
  }

public:
  InteractionBMHTF() {
    add_parameters({
        make_autoparameter(&CoreInteraction::A, "a"),
        make_autoparameter(&CoreInteraction::B, "b"),
        make_autoparameter(&CoreInteraction::C, "c"),
        make_autoparameter(&CoreInteraction::D, "d"),
        make_autoparameter(&CoreInteraction::sig, "sig"),
        make_autoparameter(&CoreInteraction::cut, "cutoff"),
    });
  }

private:
  void make_new_instance(VariantMap const &params) override {
    m_handle = make_shared_from_args<CoreInteraction, double, double, double,
                                     double, double, double>(
        params, "a", "b", "c", "d", "sig", "cutoff");
  }
};
#endif // ESPRESSO_BMHTF_NACL

#ifdef ESPRESSO_MORSE
class InteractionMorse
    : public InteractionPotentialInterface<::Morse_Parameters> {
protected:
  CoreInteractionIAParamsPtr get_ptr_offset() const override {
    return &::IA_parameters::morse;
  }

public:
  InteractionMorse() {
    add_parameters({
        make_autoparameter(&CoreInteraction::eps, "eps"),
        make_autoparameter(&CoreInteraction::alpha, "alpha"),
        make_autoparameter(&CoreInteraction::rmin, "rmin"),
        make_autoparameter(&CoreInteraction::cut, "cutoff"),
    });
  }

private:
  void make_new_instance(VariantMap const &params) override {
    m_handle =
        make_shared_from_args<CoreInteraction, double, double, double, double>(
            params, "eps", "alpha", "rmin", "cutoff");
  }
};
#endif // ESPRESSO_MORSE

#ifdef ESPRESSO_BUCKINGHAM
class InteractionBuckingham
    : public InteractionPotentialInterface<::Buckingham_Parameters> {
protected:
  CoreInteractionIAParamsPtr get_ptr_offset() const override {
    return &::IA_parameters::buckingham;
  }

public:
  InteractionBuckingham() {
    add_parameters({
        make_autoparameter(&CoreInteraction::A, "a"),
        make_autoparameter(&CoreInteraction::B, "b"),
        make_autoparameter(&CoreInteraction::C, "c"),
        make_autoparameter(&CoreInteraction::D, "d"),
        make_autoparameter(&CoreInteraction::cut, "cutoff"),
        make_autoparameter(&CoreInteraction::discont, "discont"),
        make_autoparameter(&CoreInteraction::shift, "shift"),
    });
  }

private:
  void make_new_instance(VariantMap const &params) override {
    m_handle = make_shared_from_args<CoreInteraction, double, double, double,
                                     double, double, double, double>(
        params, "a", "b", "c", "d", "cutoff", "discont", "shift");
  }
};
#endif // ESPRESSO_BUCKINGHAM

#ifdef ESPRESSO_SOFT_SPHERE
class InteractionSoftSphere
    : public InteractionPotentialInterface<::SoftSphere_Parameters> {
protected:
  CoreInteractionIAParamsPtr get_ptr_offset() const override {
    return &::IA_parameters::soft_sphere;
  }

public:
  InteractionSoftSphere() {
    add_parameters({
        make_autoparameter(&CoreInteraction::a, "a"),
        make_autoparameter(&CoreInteraction::n, "n"),
        make_autoparameter(&CoreInteraction::cut, "cutoff"),
        make_autoparameter(&CoreInteraction::offset, "offset"),
    });
  }

private:
  void make_new_instance(VariantMap const &params) override {
    m_handle =
        make_shared_from_args<CoreInteraction, double, double, double, double>(
            params, "a", "n", "cutoff", "offset");
  }
};
#endif // ESPRESSO_SOFT_SPHERE

#ifdef ESPRESSO_HAT
class InteractionHat : public InteractionPotentialInterface<::Hat_Parameters> {
protected:
  CoreInteractionIAParamsPtr get_ptr_offset() const override {
    return &::IA_parameters::hat;
  }

public:
  InteractionHat() {
    add_parameters({
        make_autoparameter(&CoreInteraction::Fmax, "F_max"),
        make_autoparameter(&CoreInteraction::r, "cutoff"),
    });
  }

private:
  void make_new_instance(VariantMap const &params) override {
    m_handle = make_shared_from_args<CoreInteraction, double, double>(
        params, "F_max", "cutoff");
  }
};
#endif // ESPRESSO_HAT

#ifdef ESPRESSO_GAY_BERNE
class InteractionGayBerne
    : public InteractionPotentialInterface<::GayBerne_Parameters> {
protected:
  CoreInteractionIAParamsPtr get_ptr_offset() const override {
    return &::IA_parameters::gay_berne;
  }

public:
  InteractionGayBerne() {
    add_parameters({
        make_autoparameter(&CoreInteraction::eps, "eps"),
        make_autoparameter(&CoreInteraction::sig, "sig"),
        make_autoparameter(&CoreInteraction::cut, "cut"),
        make_autoparameter(&CoreInteraction::k1, "k1"),
        make_autoparameter(&CoreInteraction::k2, "k2"),
        make_autoparameter(&CoreInteraction::mu, "mu"),
        make_autoparameter(&CoreInteraction::nu, "nu"),
    });
  }

private:
  std::string inactive_parameter() const override { return "cut"; }

  void make_new_instance(VariantMap const &params) override {
    m_handle = make_shared_from_args<CoreInteraction, double, double, double,
                                     double, double, double, double>(
        params, "eps", "sig", "cut", "k1", "k2", "mu", "nu");
  }
};
#endif // ESPRESSO_GAY_BERNE

#ifdef ESPRESSO_TABULATED
class InteractionTabulated
    : public InteractionPotentialInterface<::TabulatedPotential> {
protected:
  CoreInteractionIAParamsPtr get_ptr_offset() const override {
    return &::IA_parameters::tab;
  }

public:
  InteractionTabulated() {
    add_parameters({
        make_autoparameter(&CoreInteraction::minval, "min"),
        make_autoparameter(&CoreInteraction::maxval, "max"),
        make_autoparameter(&CoreInteraction::force_tab, "force"),
        make_autoparameter(&CoreInteraction::energy_tab, "energy"),
    });
  }

private:
  std::string inactive_parameter() const override { return "max"; }

  void make_new_instance(VariantMap const &params) override {
    m_handle = make_shared_from_args<CoreInteraction, double, double,
                                     std::vector<double>, std::vector<double>>(
        params, "min", "max", "force", "energy");
  }

public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "get_cutoff") {
      return m_handle.get()->cutoff();
    }
    return InteractionPotentialInterface<CoreInteraction>::do_call_method(
        name, params);
  }
};
#endif // ESPRESSO_TABULATED

#ifdef ESPRESSO_DPD
class InteractionDPD : public InteractionPotentialInterface<::DPD_Parameters> {
protected:
  CoreInteractionIAParamsPtr get_ptr_offset() const override {
    return &::IA_parameters::dpd;
  }

public:
  InteractionDPD() {
    add_parameters({
        {"weight_function", AutoParameter::read_only,
         [this]() { return m_handle.get()->radial.wf; }},
        {"gamma", AutoParameter::read_only,
         [this]() { return m_handle.get()->radial.gamma; }},
        {"k", AutoParameter::read_only,
         [this]() { return m_handle.get()->radial.k; }},
        {"r_cut", AutoParameter::read_only,
         [this]() { return m_handle.get()->radial.cutoff; }},
        {"trans_weight_function", AutoParameter::read_only,
         [this]() { return m_handle.get()->trans.wf; }},
        {"trans_gamma", AutoParameter::read_only,
         [this]() { return m_handle.get()->trans.gamma; }},
        {"trans_r_cut", AutoParameter::read_only,
         [this]() { return m_handle.get()->trans.cutoff; }},
    });
    std::ignore = get_ptr_offset(); // for code coverage
  }

private:
  std::string inactive_parameter() const override { return "r_cut"; }

  void make_new_instance(VariantMap const &params) override {
    m_handle = make_shared_from_args<CoreInteraction, double, double, double,
                                     double, double, double, double>(
        params, "gamma", "k", "r_cut", "weight_function", "trans_gamma",
        "trans_r_cut", "trans_weight_function");
    if (m_handle->radial.wf != 0 and m_handle->radial.wf != 1) {
      throw std::domain_error(
          "DPDInteraction parameter 'weight_function' must be 0 or 1");
    }
    if (m_handle->trans.wf != 0 and m_handle->trans.wf != 1) {
      throw std::domain_error(
          "DPDInteraction parameter 'trans_weight_function' must be 0 or 1");
    }
  }
};
#endif // ESPRESSO_DPD

#ifdef ESPRESSO_THOLE
class InteractionThole
    : public InteractionPotentialInterface<::Thole_Parameters> {
protected:
  CoreInteractionIAParamsPtr get_ptr_offset() const override {
    return &::IA_parameters::thole;
  }

public:
  InteractionThole() {
    add_parameters({
        make_autoparameter(&CoreInteraction::scaling_coeff, "scaling_coeff"),
        make_autoparameter(&CoreInteraction::q1q2, "q1q2"),
    });
  }

private:
  std::string inactive_parameter() const override { return "scaling_coeff"; }
  double get_inactive_cutoff() const override { return 0.; }

  void make_new_instance(VariantMap const &params) override {
    m_handle = make_shared_from_args<CoreInteraction, double, double>(
        params, "scaling_coeff", "q1q2");
  }
};
#endif // ESPRESSO_THOLE

#ifdef ESPRESSO_SMOOTH_STEP
class InteractionSmoothStep
    : public InteractionPotentialInterface<::SmoothStep_Parameters> {
protected:
  CoreInteractionIAParamsPtr get_ptr_offset() const override {
    return &::IA_parameters::smooth_step;
  }

public:
  InteractionSmoothStep() {
    add_parameters({
        make_autoparameter(&CoreInteraction::eps, "eps"),
        make_autoparameter(&CoreInteraction::sig, "sig"),
        make_autoparameter(&CoreInteraction::cut, "cutoff"),
        make_autoparameter(&CoreInteraction::d, "d"),
        make_autoparameter(&CoreInteraction::n, "n"),
        make_autoparameter(&CoreInteraction::k0, "k0"),
    });
  }

private:
  void make_new_instance(VariantMap const &params) override {
    m_handle = make_shared_from_args<CoreInteraction, double, double, double,
                                     double, int, double>(
        params, "eps", "sig", "cutoff", "d", "n", "k0");
  }
};
#endif // ESPRESSO_SMOOTH_STEP

class NonBondedInteractionHandle
    : public AutoParameters<NonBondedInteractionHandle> {
  std::shared_ptr<::IA_parameters> m_handle;
  std::shared_ptr<NonBondedInteractionHandle *> m_self;
  std::weak_ptr<std::function<void()>> m_notify_cutoff_change;
#ifdef ESPRESSO_WCA
  std::shared_ptr<InteractionWCA> m_wca;
#endif
#ifdef ESPRESSO_LENNARD_JONES
  std::shared_ptr<InteractionLJ> m_lj;
#endif
#ifdef ESPRESSO_LENNARD_JONES_GENERIC
  std::shared_ptr<InteractionLJGen> m_ljgen;
#endif
#ifdef ESPRESSO_LJCOS
  std::shared_ptr<InteractionLJcos> m_ljcos;
#endif
#ifdef ESPRESSO_LJCOS2
  std::shared_ptr<InteractionLJcos2> m_ljcos2;
#endif
#ifdef ESPRESSO_HERTZIAN
  std::shared_ptr<InteractionHertzian> m_hertzian;
#endif
#ifdef ESPRESSO_GAUSSIAN
  std::shared_ptr<InteractionGaussian> m_gaussian;
#endif
#ifdef ESPRESSO_BMHTF_NACL
  std::shared_ptr<InteractionBMHTF> m_bmhtf;
#endif
#ifdef ESPRESSO_MORSE
  std::shared_ptr<InteractionMorse> m_morse;
#endif
#ifdef ESPRESSO_BUCKINGHAM
  std::shared_ptr<InteractionBuckingham> m_buckingham;
#endif
#ifdef ESPRESSO_SOFT_SPHERE
  std::shared_ptr<InteractionSoftSphere> m_soft_sphere;
#endif
#ifdef ESPRESSO_HAT
  std::shared_ptr<InteractionHat> m_hat;
#endif
#ifdef ESPRESSO_GAY_BERNE
  std::shared_ptr<InteractionGayBerne> m_gay_berne;
#endif
#ifdef ESPRESSO_TABULATED
  std::shared_ptr<InteractionTabulated> m_tabulated;
#endif
#ifdef ESPRESSO_DPD
  std::shared_ptr<InteractionDPD> m_dpd;
#endif
#ifdef ESPRESSO_THOLE
  std::shared_ptr<InteractionThole> m_thole;
#endif
#ifdef ESPRESSO_SMOOTH_STEP
  std::shared_ptr<InteractionSmoothStep> m_smooth_step;
#endif

public:
  NonBondedInteractionHandle() {
    m_self = std::make_shared<NonBondedInteractionHandle *>(this);
    std::vector<AutoParameter> params;
    apply([this, &params]<typename T>(std::shared_ptr<T> &member,
                                      std::string const &name,
                                      std::string const &) {
      auto const setter = [this, &member](Variant const &v) {
        member = get_value<std::shared_ptr<T>>(v);
        member->attach(m_self, m_handle);
        // copying the new value to the core may queue a runtime error message,
        // but we can't detect it to roll back to the last valid state
        member->update_core();
      };
      params.emplace_back(name.c_str(), setter, [&member]() { return member; });
    });
    add_parameters(std::move(params));
  }

public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "reset") {
      if (not context()->is_head_node()) {
        return {};
      }
      auto const new_params = VariantMap{{"notify", false}};
      apply([&new_params](auto &so, std::string const &, std::string const &) {
        so->call_method("deactivate", new_params);
      });
      if (get_value_or<bool>(params, "notify", true)) {
        call_method("on_non_bonded_ia_change", {});
      }
      return {};
    }
    if (name == "on_non_bonded_ia_change") {
      on_non_bonded_ia_change();
      return {};
    }
    return {};
  }

  void do_construct(VariantMap const &params) override {
    m_handle = std::make_shared<::IA_parameters>();
    if (not context()->is_head_node()) {
      return;
    }
    apply([this, &params]<typename T>(std::shared_ptr<T> const &,
                                      std::string const &name,
                                      std::string const &so_name) {
      auto so = params.contains(name)
                    ? get_value<std::shared_ptr<T>>(params.at(name))
                    : std::dynamic_pointer_cast<T>(
                          context()->make_shared(so_name, {}));
      set_parameter(name, so);
    });
  }

  void attach(
      std::function<void(std::shared_ptr<::IA_parameters> const &)> cb_register,
      std::weak_ptr<std::function<void()>> cb_notify_cutoff_change) {
    cb_register(m_handle);
    m_notify_cutoff_change = cb_notify_cutoff_change;
  }

private:
  void apply(auto const &&fun) {
#ifdef ESPRESSO_WCA
    fun(m_wca, "wca", "Interactions::InteractionWCA");
#endif
#ifdef ESPRESSO_LENNARD_JONES
    fun(m_lj, "lennard_jones", "Interactions::InteractionLJ");
#endif
#ifdef ESPRESSO_LENNARD_JONES_GENERIC
    fun(m_ljgen, "generic_lennard_jones", "Interactions::InteractionLJGen");
#endif
#ifdef ESPRESSO_LJCOS
    fun(m_ljcos, "lennard_jones_cos", "Interactions::InteractionLJcos");
#endif
#ifdef ESPRESSO_LJCOS2
    fun(m_ljcos2, "lennard_jones_cos2", "Interactions::InteractionLJcos2");
#endif
#ifdef ESPRESSO_HERTZIAN
    fun(m_hertzian, "hertzian", "Interactions::InteractionHertzian");
#endif
#ifdef ESPRESSO_GAUSSIAN
    fun(m_gaussian, "gaussian", "Interactions::InteractionGaussian");
#endif
#ifdef ESPRESSO_BMHTF_NACL
    fun(m_bmhtf, "bmhtf", "Interactions::InteractionBMHTF");
#endif
#ifdef ESPRESSO_MORSE
    fun(m_morse, "morse", "Interactions::InteractionMorse");
#endif
#ifdef ESPRESSO_BUCKINGHAM
    fun(m_buckingham, "buckingham", "Interactions::InteractionBuckingham");
#endif
#ifdef ESPRESSO_SOFT_SPHERE
    fun(m_soft_sphere, "soft_sphere", "Interactions::InteractionSoftSphere");
#endif
#ifdef ESPRESSO_HAT
    fun(m_hat, "hat", "Interactions::InteractionHat");
#endif
#ifdef ESPRESSO_GAY_BERNE
    fun(m_gay_berne, "gay_berne", "Interactions::InteractionGayBerne");
#endif
#ifdef ESPRESSO_TABULATED
    fun(m_tabulated, "tabulated", "Interactions::InteractionTabulated");
#endif
#ifdef ESPRESSO_DPD
    fun(m_dpd, "dpd", "Interactions::InteractionDPD");
#endif
#ifdef ESPRESSO_THOLE
    fun(m_thole, "thole", "Interactions::InteractionThole");
#endif
#ifdef ESPRESSO_SMOOTH_STEP
    fun(m_smooth_step, "smooth_step", "Interactions::InteractionSmoothStep");
#endif
  }

public:
  void on_non_bonded_ia_change() {
    if (auto callback = m_notify_cutoff_change.lock()) {
      (*callback)();
    }
  }
};

template <class CoreIA>
void InteractionPotentialInterface<CoreIA>::update_core(bool notify) {
  assert(m_handle);
  if (auto core_struct = m_core_struct.lock()) {
    core_struct.get()->*get_ptr_offset() = *m_handle;
    if (notify) {
      if (auto si_struct = m_si_struct.lock()) {
        (**si_struct).on_non_bonded_ia_change();
      }
    }
  }
}

} // namespace Interactions
} // namespace ScriptInterface
