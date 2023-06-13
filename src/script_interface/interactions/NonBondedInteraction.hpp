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

/** @file
 *  The ScriptInterface counterparts of the non-bonded interactions parameters
 *  structs from the core are defined here.
 */

#ifndef SCRIPT_INTERFACE_INTERACTIONS_NONBONDED_INTERACTION_HPP
#define SCRIPT_INTERFACE_INTERACTIONS_NONBONDED_INTERACTION_HPP

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/get_value.hpp"

#include "core/event.hpp"
#include "core/nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace ScriptInterface {
namespace Interactions {

template <class CoreIA>
class InteractionPotentialInterface
    : public AutoParameters<InteractionPotentialInterface<CoreIA>> {
  /** @brief Particle type pair. */
  std::array<int, 2> m_types = {-1, -1};

public:
  using CoreInteraction = CoreIA;

protected:
  using BaseClass = AutoParameters<InteractionPotentialInterface<CoreIA>>;
  using BaseClass::context;
  using BaseClass::get_valid_parameters;

  /** @brief Managed object. */
  std::shared_ptr<CoreInteraction> m_ia_si;
  /** @brief Pointer to the corresponding member in a handle. */
  virtual CoreInteraction IA_parameters::*get_ptr_offset() const = 0;
  /** @brief Create a new instance using the constructor with range checks. */
  virtual void make_new_instance(VariantMap const &params) = 0;
  /** @brief Which parameter indicates whether the potential is inactive. */
  virtual std::string inactive_parameter() const { return "cutoff"; }
  /** @brief Which magic value indicates the potential is inactive. */
  virtual double inactive_cutoff() const { return INACTIVE_CUTOFF; }

  template <typename T>
  auto make_autoparameter(T CoreInteraction::*ptr, char const *name) {
    return AutoParameter{name, AutoParameter::read_only,
                         [this, ptr]() { return m_ia_si.get()->*ptr; }};
  }

private:
  void check_valid_parameters(VariantMap const &params) const {
    auto const keys = get_valid_parameters();
    for (auto const &key : keys) {
      if (params.count(std::string(key)) == 0) {
        throw std::runtime_error("Parameter '" + key + "' is missing");
      }
    }
    for (auto const &kv : params) {
      if (std::find(keys.begin(), keys.end(), kv.first) == keys.end()) {
        throw std::runtime_error("Parameter '" + kv.first +
                                 "' is not recognized");
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
      if (m_types[0] != -1) {
        copy_si_to_core();
        on_non_bonded_ia_change();
      }
      return {};
    }
    if (name == "deactivate") {
      m_ia_si = std::make_shared<CoreInteraction>();
      if (m_types[0] != -1) {
        copy_si_to_core();
        on_non_bonded_ia_change();
      }
      return {};
    }
    if (name == "is_registered") {
      return m_types[0] != -1;
    }
    if (name == "bind_types") {
      auto types = get_value<std::vector<int>>(params, "_types");
      if (types[0] > types[1]) {
        std::swap(types[0], types[1]);
      }
      if (m_types[0] == -1 or
          (m_types[0] == types[0] and m_types[1] == types[1])) {
        m_types[0] = types[0];
        m_types[1] = types[1];
      } else {
        context()->parallel_try_catch([this]() {
          throw std::runtime_error(
              "Non-bonded interaction is already bound to interaction pair [" +
              std::to_string(m_types[0]) + ", " + std::to_string(m_types[1]) +
              "]");
        });
      }
      return {};
    }
    return {};
  }

  void do_construct(VariantMap const &params) final {
    if (params.count("_types") != 0) {
      do_call_method("bind_types", params);
      m_ia_si = std::make_shared<CoreInteraction>();
      copy_core_to_si();
    } else {
      if (std::abs(get_value<double>(params, inactive_parameter()) -
                   inactive_cutoff()) < 1e-9) {
        m_ia_si = std::make_shared<CoreInteraction>();
      } else {
        context()->parallel_try_catch([this, &params]() {
          check_valid_parameters(params);
          make_new_instance(params);
        });
      }
    }
  }

  void copy_si_to_core() {
    assert(m_ia_si != nullptr);
    auto const key = get_ia_param_key(m_types[0], m_types[1]);
    assert(key < ::nonbonded_ia_params.size());
    ::nonbonded_ia_params[key].get()->*get_ptr_offset() = *m_ia_si;
  }

  void copy_core_to_si() {
    assert(m_ia_si != nullptr);
    auto const key = get_ia_param_key(m_types[0], m_types[1]);
    assert(key < ::nonbonded_ia_params.size());
    *m_ia_si = ::nonbonded_ia_params[key].get()->*get_ptr_offset();
  }
};

#ifdef WCA
class InteractionWCA : public InteractionPotentialInterface<::WCA_Parameters> {
protected:
  CoreInteraction IA_parameters::*get_ptr_offset() const override {
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
  double inactive_cutoff() const override { return 0.; }

  void make_new_instance(VariantMap const &params) override {
    m_ia_si = make_shared_from_args<CoreInteraction, double, double>(
        params, "epsilon", "sigma");
  }

public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "get_cutoff") {
      return m_ia_si.get()->cut;
    }
    return InteractionPotentialInterface<CoreInteraction>::do_call_method(
        name, params);
  }
};
#endif // WCA

#ifdef LENNARD_JONES
class InteractionLJ : public InteractionPotentialInterface<::LJ_Parameters> {
protected:
  CoreInteraction IA_parameters::*get_ptr_offset() const override {
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
    auto const *shift_string = boost::get<std::string>(&params.at("shift"));
    if (shift_string != nullptr) {
      if (*shift_string != "auto") {
        throw std::invalid_argument(
            "LJ parameter 'shift' has to be 'auto' or a float");
      }
      new_params["shift"] = 0.;
    }
    m_ia_si = make_shared_from_args<CoreInteraction, double, double, double,
                                    double, double, double>(
        new_params, "epsilon", "sigma", "cutoff", "offset", "min", "shift");
    if (shift_string != nullptr) {
      m_ia_si->shift = m_ia_si->get_auto_shift();
    }
  }
};
#endif // LENNARD_JONES

#ifdef LENNARD_JONES_GENERIC
class InteractionLJGen
    : public InteractionPotentialInterface<::LJGen_Parameters> {
protected:
  CoreInteraction IA_parameters::*get_ptr_offset() const override {
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
#ifdef LJGEN_SOFTCORE
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
    auto const *shift_string = boost::get<std::string>(&params.at("shift"));
    if (shift_string != nullptr) {
      if (*shift_string != "auto") {
        throw std::invalid_argument(
            "Generic LJ parameter 'shift' has to be 'auto' or a float");
      }
      new_params["shift"] = 0.;
    }
    m_ia_si = make_shared_from_args<CoreInteraction, double, double, double,
                                    double, double,
#ifdef LJGEN_SOFTCORE
                                    double, double,
#endif
                                    double, double, double, double>(
        new_params, "epsilon", "sigma", "cutoff", "shift", "offset",
#ifdef LJGEN_SOFTCORE
        "lam", "delta",
#endif
        "e1", "e2", "b1", "b2");
    if (shift_string != nullptr) {
      m_ia_si->shift = m_ia_si->get_auto_shift();
    }
  }
};
#endif // LENNARD_JONES_GENERIC

#ifdef LJCOS
class InteractionLJcos
    : public InteractionPotentialInterface<::LJcos_Parameters> {
protected:
  CoreInteraction IA_parameters::*get_ptr_offset() const override {
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
    m_ia_si =
        make_shared_from_args<CoreInteraction, double, double, double, double>(
            params, "epsilon", "sigma", "cutoff", "offset");
  }
};
#endif // LJCOS

#ifdef LJCOS2
class InteractionLJcos2
    : public InteractionPotentialInterface<::LJcos2_Parameters> {
protected:
  CoreInteraction IA_parameters::*get_ptr_offset() const override {
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
  double inactive_cutoff() const override { return 0.; }

  void make_new_instance(VariantMap const &params) override {
    m_ia_si =
        make_shared_from_args<CoreInteraction, double, double, double, double>(
            params, "epsilon", "sigma", "offset", "width");
  }

public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "get_cutoff") {
      return m_ia_si.get()->cut;
    }
    return InteractionPotentialInterface<CoreInteraction>::do_call_method(
        name, params);
  }
};
#endif // LJCOS2

#ifdef HERTZIAN
class InteractionHertzian
    : public InteractionPotentialInterface<::Hertzian_Parameters> {
protected:
  CoreInteraction IA_parameters::*get_ptr_offset() const override {
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
    m_ia_si = make_shared_from_args<CoreInteraction, double, double>(
        params, "eps", "sig");
  }
};
#endif // HERTZIAN

#ifdef GAUSSIAN
class InteractionGaussian
    : public InteractionPotentialInterface<::Gaussian_Parameters> {
protected:
  CoreInteraction IA_parameters::*get_ptr_offset() const override {
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
    m_ia_si = make_shared_from_args<CoreInteraction, double, double, double>(
        params, "eps", "sig", "cutoff");
  }
};
#endif // GAUSSIAN

#ifdef BMHTF_NACL
class InteractionBMHTF
    : public InteractionPotentialInterface<::BMHTF_Parameters> {
protected:
  CoreInteraction IA_parameters::*get_ptr_offset() const override {
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
    m_ia_si = make_shared_from_args<CoreInteraction, double, double, double,
                                    double, double, double>(
        params, "a", "b", "c", "d", "sig", "cutoff");
  }
};
#endif // BMHTF_NACL

#ifdef MORSE
class InteractionMorse
    : public InteractionPotentialInterface<::Morse_Parameters> {
protected:
  CoreInteraction IA_parameters::*get_ptr_offset() const override {
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
    m_ia_si =
        make_shared_from_args<CoreInteraction, double, double, double, double>(
            params, "eps", "alpha", "rmin", "cutoff");
  }
};
#endif // MORSE

#ifdef BUCKINGHAM
class InteractionBuckingham
    : public InteractionPotentialInterface<::Buckingham_Parameters> {
protected:
  CoreInteraction IA_parameters::*get_ptr_offset() const override {
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
    m_ia_si = make_shared_from_args<CoreInteraction, double, double, double,
                                    double, double, double, double>(
        params, "a", "b", "c", "d", "cutoff", "discont", "shift");
  }
};
#endif // BUCKINGHAM

#ifdef SOFT_SPHERE
class InteractionSoftSphere
    : public InteractionPotentialInterface<::SoftSphere_Parameters> {
protected:
  CoreInteraction IA_parameters::*get_ptr_offset() const override {
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
    m_ia_si =
        make_shared_from_args<CoreInteraction, double, double, double, double>(
            params, "a", "n", "cutoff", "offset");
  }
};
#endif // SOFT_SPHERE

#ifdef HAT
class InteractionHat : public InteractionPotentialInterface<::Hat_Parameters> {
protected:
  CoreInteraction IA_parameters::*get_ptr_offset() const override {
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
    m_ia_si = make_shared_from_args<CoreInteraction, double, double>(
        params, "F_max", "cutoff");
  }
};
#endif // HAT

#ifdef GAY_BERNE
class InteractionGayBerne
    : public InteractionPotentialInterface<::GayBerne_Parameters> {
protected:
  CoreInteraction IA_parameters::*get_ptr_offset() const override {
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
    m_ia_si = make_shared_from_args<CoreInteraction, double, double, double,
                                    double, double, double, double>(
        params, "eps", "sig", "cut", "k1", "k2", "mu", "nu");
  }
};
#endif // GAY_BERNE

#ifdef TABULATED
class InteractionTabulated
    : public InteractionPotentialInterface<::TabulatedPotential> {
protected:
  CoreInteraction IA_parameters::*get_ptr_offset() const override {
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
    m_ia_si = make_shared_from_args<CoreInteraction, double, double,
                                    std::vector<double>, std::vector<double>>(
        params, "min", "max", "force", "energy");
  }

public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "get_cutoff") {
      return m_ia_si.get()->cutoff();
    }
    return InteractionPotentialInterface<CoreInteraction>::do_call_method(
        name, params);
  }
};
#endif // TABULATED

#ifdef DPD
class InteractionDPD : public InteractionPotentialInterface<::DPD_Parameters> {
protected:
  CoreInteraction IA_parameters::*get_ptr_offset() const override {
    return &::IA_parameters::dpd;
  }

public:
  InteractionDPD() {
    add_parameters({
        {"weight_function", AutoParameter::read_only,
         [this]() { return m_ia_si.get()->radial.wf; }},
        {"gamma", AutoParameter::read_only,
         [this]() { return m_ia_si.get()->radial.gamma; }},
        {"k", AutoParameter::read_only,
         [this]() { return m_ia_si.get()->radial.k; }},
        {"r_cut", AutoParameter::read_only,
         [this]() { return m_ia_si.get()->radial.cutoff; }},
        {"trans_weight_function", AutoParameter::read_only,
         [this]() { return m_ia_si.get()->trans.wf; }},
        {"trans_gamma", AutoParameter::read_only,
         [this]() { return m_ia_si.get()->trans.gamma; }},
        {"trans_r_cut", AutoParameter::read_only,
         [this]() { return m_ia_si.get()->trans.cutoff; }},
    });
    std::ignore = get_ptr_offset(); // for code coverage
  }

private:
  std::string inactive_parameter() const override { return "r_cut"; }

  void make_new_instance(VariantMap const &params) override {
    m_ia_si = make_shared_from_args<CoreInteraction, double, double, double,
                                    double, double, double, double>(
        params, "gamma", "k", "r_cut", "weight_function", "trans_gamma",
        "trans_r_cut", "trans_weight_function");
    if (m_ia_si->radial.wf != 0 and m_ia_si->radial.wf != 1) {
      throw std::domain_error(
          "DPDInteraction parameter 'weight_function' must be 0 or 1");
    }
    if (m_ia_si->trans.wf != 0 and m_ia_si->trans.wf != 1) {
      throw std::domain_error(
          "DPDInteraction parameter 'trans_weight_function' must be 0 or 1");
    }
  }
};
#endif // DPD

#ifdef THOLE
class InteractionThole
    : public InteractionPotentialInterface<::Thole_Parameters> {
protected:
  CoreInteraction IA_parameters::*get_ptr_offset() const override {
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
  double inactive_cutoff() const override { return 0.; }

  void make_new_instance(VariantMap const &params) override {
    m_ia_si = make_shared_from_args<CoreInteraction, double, double>(
        params, "scaling_coeff", "q1q2");
  }
};
#endif // THOLE

#ifdef SMOOTH_STEP
class InteractionSmoothStep
    : public InteractionPotentialInterface<::SmoothStep_Parameters> {
protected:
  CoreInteraction IA_parameters::*get_ptr_offset() const override {
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
    m_ia_si = make_shared_from_args<CoreInteraction, double, double, double,
                                    double, int, double>(
        params, "eps", "sig", "cutoff", "d", "n", "k0");
  }
};
#endif // SMOOTH_STEP

class NonBondedInteractionHandle
    : public AutoParameters<NonBondedInteractionHandle> {
  std::array<int, 2> m_types = {-1, -1};
  std::shared_ptr<::IA_parameters> m_interaction;
#ifdef WCA
  std::shared_ptr<InteractionWCA> m_wca;
#endif
#ifdef LENNARD_JONES
  std::shared_ptr<InteractionLJ> m_lj;
#endif
#ifdef LENNARD_JONES_GENERIC
  std::shared_ptr<InteractionLJGen> m_ljgen;
#endif
#ifdef LJCOS
  std::shared_ptr<InteractionLJcos> m_ljcos;
#endif
#ifdef LJCOS2
  std::shared_ptr<InteractionLJcos2> m_ljcos2;
#endif
#ifdef HERTZIAN
  std::shared_ptr<InteractionHertzian> m_hertzian;
#endif
#ifdef GAUSSIAN
  std::shared_ptr<InteractionGaussian> m_gaussian;
#endif
#ifdef BMHTF_NACL
  std::shared_ptr<InteractionBMHTF> m_bmhtf;
#endif
#ifdef MORSE
  std::shared_ptr<InteractionMorse> m_morse;
#endif
#ifdef BUCKINGHAM
  std::shared_ptr<InteractionBuckingham> m_buckingham;
#endif
#ifdef SOFT_SPHERE
  std::shared_ptr<InteractionSoftSphere> m_soft_sphere;
#endif
#ifdef HAT
  std::shared_ptr<InteractionHat> m_hat;
#endif
#ifdef GAY_BERNE
  std::shared_ptr<InteractionGayBerne> m_gay_berne;
#endif
#ifdef TABULATED
  std::shared_ptr<InteractionTabulated> m_tabulated;
#endif
#ifdef DPD
  std::shared_ptr<InteractionDPD> m_dpd;
#endif
#ifdef THOLE
  std::shared_ptr<InteractionThole> m_thole;
#endif
#ifdef SMOOTH_STEP
  std::shared_ptr<InteractionSmoothStep> m_smooth_step;
#endif

  template <class T>
  auto make_autoparameter(std::shared_ptr<T> &member, const char *key) const {
    auto const setter = [this, &member](Variant const &v) {
      member = get_value<std::shared_ptr<T>>(v);
      if (m_types[0] != -1) {
        auto const types = Variant{std::vector<int>{{m_types[0], m_types[1]}}};
        member->do_call_method("bind_types", VariantMap{{"_types", types}});
        member->copy_si_to_core();
        on_non_bonded_ia_change();
      }
    };
    return AutoParameter{key, setter, [&member]() { return member; }};
  }

public:
  NonBondedInteractionHandle() {
    add_parameters({
#ifdef WCA
        make_autoparameter(m_wca, "wca"),
#endif
#ifdef LENNARD_JONES
        make_autoparameter(m_lj, "lennard_jones"),
#endif
#ifdef LENNARD_JONES_GENERIC
        make_autoparameter(m_ljgen, "generic_lennard_jones"),
#endif
#ifdef LJCOS
        make_autoparameter(m_ljcos, "lennard_jones_cos"),
#endif
#ifdef LJCOS2
        make_autoparameter(m_ljcos2, "lennard_jones_cos2"),
#endif
#ifdef HERTZIAN
        make_autoparameter(m_hertzian, "hertzian"),
#endif
#ifdef GAUSSIAN
        make_autoparameter(m_gaussian, "gaussian"),
#endif
#ifdef BMHTF_NACL
        make_autoparameter(m_bmhtf, "bmhtf"),
#endif
#ifdef MORSE
        make_autoparameter(m_morse, "morse"),
#endif
#ifdef BUCKINGHAM
        make_autoparameter(m_buckingham, "buckingham"),
#endif
#ifdef SOFT_SPHERE
        make_autoparameter(m_soft_sphere, "soft_sphere"),
#endif
#ifdef HAT
        make_autoparameter(m_hat, "hat"),
#endif
#ifdef GAY_BERNE
        make_autoparameter(m_gay_berne, "gay_berne"),
#endif
#ifdef TABULATED
        make_autoparameter(m_tabulated, "tabulated"),
#endif
#ifdef DPD
        make_autoparameter(m_dpd, "dpd"),
#endif
#ifdef THOLE
        make_autoparameter(m_thole, "thole"),
#endif
#ifdef SMOOTH_STEP
        make_autoparameter(m_smooth_step, "smooth_step"),
#endif
    });
  }

private:
  template <class T>
  void set_member(std::shared_ptr<T> &member, std::string key,
                  std::string so_name, VariantMap const &params) {
    auto const ia_types = VariantMap{{"_types", params.at("_types")}};
    if (params.count(key) != 0) {
      member = get_value<std::shared_ptr<T>>(params.at(key));
      member->do_call_method("bind_types", ia_types);
      member->copy_si_to_core();
    } else {
      auto so_object = context()->make_shared_local(so_name, ia_types);
      member = std::dynamic_pointer_cast<T>(so_object);
      member->copy_core_to_si();
    }
  }

public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "get_types") {
      return std::vector<int>{{m_types[0], m_types[1]}};
    }
    return {};
  }

  void do_construct(VariantMap const &params) override {
    assert(params.count("_types") != 0);
    auto const types = get_value<std::vector<int>>(params.at("_types"));
    m_types[0] = std::min(types[0], types[1]);
    m_types[1] = std::max(types[0], types[1]);
    make_particle_type_exist(m_types[1]);
    // create interface objects
    auto const key = get_ia_param_key(m_types[0], m_types[1]);
    m_interaction = ::nonbonded_ia_params[key];
#ifdef WCA
    set_member(m_wca, "wca", "Interactions::InteractionWCA", params);
#endif
#ifdef LENNARD_JONES
    set_member(m_lj, "lennard_jones", "Interactions::InteractionLJ", params);
#endif
#ifdef LENNARD_JONES_GENERIC
    set_member(m_ljgen, "generic_lennard_jones",
               "Interactions::InteractionLJGen", params);
#endif
#ifdef LJCOS
    set_member(m_ljcos, "lennard_jones_cos", "Interactions::InteractionLJcos",
               params);
#endif
#ifdef LJCOS2
    set_member(m_ljcos2, "lennard_jones_cos2",
               "Interactions::InteractionLJcos2", params);
#endif
#ifdef HERTZIAN
    set_member(m_hertzian, "hertzian", "Interactions::InteractionHertzian",
               params);
#endif
#ifdef GAUSSIAN
    set_member(m_gaussian, "gaussian", "Interactions::InteractionGaussian",
               params);
#endif
#ifdef BMHTF_NACL
    set_member(m_bmhtf, "bmhtf", "Interactions::InteractionBMHTF", params);
#endif
#ifdef MORSE
    set_member(m_morse, "morse", "Interactions::InteractionMorse", params);
#endif
#ifdef BUCKINGHAM
    set_member(m_buckingham, "buckingham",
               "Interactions::InteractionBuckingham", params);
#endif
#ifdef SOFT_SPHERE
    set_member(m_soft_sphere, "soft_sphere",
               "Interactions::InteractionSoftSphere", params);
#endif
#ifdef HAT
    set_member(m_hat, "hat", "Interactions::InteractionHat", params);
#endif
#ifdef GAY_BERNE
    set_member(m_gay_berne, "gay_berne", "Interactions::InteractionGayBerne",
               params);
#endif
#ifdef TABULATED
    set_member(m_tabulated, "tabulated", "Interactions::InteractionTabulated",
               params);
#endif
#ifdef DPD
    set_member(m_dpd, "dpd", "Interactions::InteractionDPD", params);
#endif
#ifdef THOLE
    set_member(m_thole, "thole", "Interactions::InteractionThole", params);
#endif
#ifdef SMOOTH_STEP
    set_member(m_smooth_step, "smooth_step",
               "Interactions::InteractionSmoothStep", params);
#endif
  }

  auto get_ia() const { return m_interaction; }
};

} // namespace Interactions
} // namespace ScriptInterface

#endif
