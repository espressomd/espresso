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

#include "core/PropagationMode.hpp"
#include "core/bonded_interactions/bonded_interaction_data.hpp"
#include "core/thermostat.hpp"

#include <script_interface/Context.hpp>
#include <script_interface/auto_parameters/AutoParameters.hpp>
#include <script_interface/system/Leaf.hpp>
#ifdef WALBERLA
#include <script_interface/walberla/LBFluid.hpp>
#endif

#include <cassert>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

namespace ScriptInterface {
namespace Thermostat {

template <typename CoreClass>
class Interface : public AutoParameters<Interface<CoreClass>, System::Leaf> {
  using BaseClass = AutoParameters<Interface<CoreClass>, System::Leaf>;

public:
  using CoreThermostat = CoreClass;
  using BaseClass::do_set_parameter;
  using BaseClass::get_parameter;
  using System::Leaf::bind_system;

protected:
  using BaseClass::add_parameters;
  using BaseClass::context;
  using BaseClass::get_parameter_insertion_order;
  using System::Leaf::m_system;

  bool is_active = false;
  std::shared_ptr<CoreThermostat> m_handle;
  /** @brief Basic lock mechanism that follows RAII. */
  std::weak_ptr<bool> m_edit_lock;

  void check_lock() {
    if (m_edit_lock.expired()) {
      throw AutoParameter::WriteError{};
    }
  }

  void on_bind_system(::System::System &system) override {
    get_member_handle(*system.thermostat) = m_handle;
    system.on_thermostat_param_change();
    is_active = true;
  }

  void on_detach_system(::System::System &system) override {
    get_member_handle(*system.thermostat).reset();
    system.on_thermostat_param_change();
    is_active = false;
  }

  void sanity_checks_positive(double value, std::string const &name) const {
    if (value < 0.) {
      throw std::domain_error("Parameter '" + name + "' cannot be negative");
    }
  }

  void sanity_checks_positive(Utils::Vector3d const &value,
                              std::string const &name) const {
    if (not(value >= Utils::Vector3d::broadcast(0.))) {
      throw std::domain_error("Parameter '" + name + "' cannot be negative");
    }
  }

  virtual bool invalid_rng_state(VariantMap const &params) const {
    return (not params.count("seed") or is_none(params.at("seed"))) and
           is_seed_required();
  }

private:
  virtual std::shared_ptr<CoreThermostat> &
  get_member_handle(::Thermostat::Thermostat &thermostat) = 0;

  void set_new_parameters(VariantMap const &params) {
    if (params.count("__check_rng_state") and invalid_rng_state(params)) {
      context()->parallel_try_catch([]() {
        throw std::invalid_argument("Parameter 'seed' is needed on first "
                                    "activation of the thermostat");
      });
    }
    for (auto const &key : get_parameter_insertion_order()) {
      if (params.count(key)) {
        auto const &v = params.at(key);
        if (key == "is_active") {
          is_active = get_value<bool>(v);
        } else {
          do_set_parameter(key.c_str(), v);
        }
      }
    }
  }

protected:
  template <typename T>
  auto make_autoparameter(T CoreThermostat::*member, char const *name) {
    return AutoParameter{
        name,
        [this, member, name = std::string(name)](Variant const &v) {
          check_lock();
          auto const value = get_value<T>(v);
          context()->parallel_try_catch(
              [&]() { sanity_checks_positive(value, name); });
          m_handle.get()->*member = std::move(value);
        },
        [this, member]() { return m_handle.get()->*member; }};
  }

  template <typename T>
  auto make_autogamma(T CoreThermostat::*member, char const *name) {
    return AutoParameter{
        name,
        [this, member, name = std::string(name)](Variant const &v) {
          check_lock();
          if (is_none(v)) {
            return;
          }
#ifdef PARTICLE_ANISOTROPY
          static_assert(std::is_same_v<T, Utils::Vector3d>);
          T gamma{};
          try {
            gamma = T::broadcast(get_value<double>(v));
          } catch (...) {
            gamma = get_value<T>(v);
          }
#else
          auto const gamma = get_value<T>(v);
#endif // PARTICLE_ANISOTROPY
          context()->parallel_try_catch(
              [&]() { sanity_checks_positive(gamma, name); });
          m_handle.get()->*member = gamma;
        },
        [this, member]() {
          auto constexpr gamma_null = ::Thermostat::gamma_null;
          auto const gamma = m_handle.get()->*member;
          return (gamma >= gamma_null) ? Variant{gamma} : Variant{None{}};
        }};
  }

  virtual VariantMap extend_parameters(VariantMap const &parameters) const {
    auto params = parameters;
    if (not is_seed_required()) {
      for (auto key : {std::string("seed"), std::string("philox_counter")}) {
        if (params.count(key) == 0ul) {
          params[key] = get_parameter(key);
        }
      }
    }
    return params;
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "override_philox_counter") {
      // only call this method if you know what you are doing
      set_rng_counter(params.at("counter"));
      return {};
    }
    return {};
  }

public:
  Interface() {
    add_parameters({
        {"seed",
         [this](Variant const &v) {
           check_lock();
           context()->parallel_try_catch([&]() {
             if (not is_none(v))
               set_rng_seed(v);
           });
         },
         [this]() {
           return m_handle->is_seed_required() ? Variant{None{}}
                                               : Variant{get_rng_seed()};
         }},
        {"philox_counter",
         [this](Variant const &v) {
           check_lock();
           context()->parallel_try_catch([&]() {
             if (not is_none(v))
               set_rng_counter(v);
           });
         },
         [this]() { return get_rng_counter(); }},
        {"is_active", AutoParameter::read_only, [this]() { return is_active; }},
    });
  }

  virtual std::optional<double> extract_kT(VariantMap const &params) const {
    if (params.count("kT")) {
      auto const value = get_value<double>(params, "kT");
      sanity_checks_positive(value, "kT");
      return value;
    }
    return {std::nullopt};
  }

  auto release_lock() {
    auto lock = std::make_shared<bool>(false);
    m_edit_lock = lock;
    return lock;
  }

  auto is_activated() const { return is_active; }

  virtual bool is_seed_required() const { return m_handle->is_seed_required(); }

  auto get_rng_seed() const {
    auto const seed = m_handle->rng_seed();
    assert(seed <= static_cast<uint32_t>(std::numeric_limits<int>::max()));
    return static_cast<int>(seed);
  }

  auto get_rng_counter() const {
    auto const counter = m_handle->rng_counter();
    assert(counter <= static_cast<uint64_t>(std::numeric_limits<int>::max()));
    return static_cast<int>(counter);
  }

  void set_rng_seed(Variant const &value) {
    auto const seed = get_value<int>(value);
    if (seed < 0) {
      throw std::domain_error("Parameter 'seed' must be a positive integer");
    }
    assert(static_cast<uint64_t>(seed) <=
           static_cast<uint64_t>(std::numeric_limits<uint32_t>::max()));
    m_handle->rng_initialize(static_cast<uint32_t>(seed));
  }

  void set_rng_counter(Variant const &value) {
    auto const counter = get_value<int>(value);
    assert(counter >= 0);
    assert(static_cast<uint64_t>(counter) <=
           std::numeric_limits<uint64_t>::max());
    m_handle->set_rng_counter(static_cast<uint64_t>(counter));
  }

  void do_construct(VariantMap const &params) override {
    m_handle = std::make_shared<CoreThermostat>();
    if (not params.empty()) {
      auto const read_write_lock = release_lock();
      set_new_parameters(params);
    }
  }

  void update_and_bind(VariantMap const &params, bool was_active,
                       std::shared_ptr<::System::System> system) {
    auto const old_handle = m_handle;
    auto new_params = extend_parameters(params);
    new_params["__global_kT"] = system->thermostat->kT;
    new_params["__check_rng_state"] = true;
    try {
      m_handle = std::make_shared<CoreThermostat>();
      set_new_parameters(new_params);
      bind_system(system);
    } catch (...) {
      assert(not is_active);
      m_handle = old_handle;
      if (was_active) {
        bind_system(system);
      }
      throw;
    }
  }

  virtual ::ThermostatFlags get_thermo_flag() const = 0;
};

class Langevin : public Interface<::LangevinThermostat> {
  std::shared_ptr<CoreThermostat> &
  get_member_handle(::Thermostat::Thermostat &thermostat) override {
    return thermostat.langevin;
  }

public:
  Langevin() {
    add_parameters({
        make_autogamma(&CoreThermostat::gamma, "gamma"),
#ifdef ROTATION
        make_autogamma(&CoreThermostat::gamma_rotation, "gamma_rotation"),
#endif
    });
  }

  ::ThermostatFlags get_thermo_flag() const final { return THERMO_LANGEVIN; }

protected:
  VariantMap extend_parameters(VariantMap const &parameters) const override {
    auto params =
        Interface<::LangevinThermostat>::extend_parameters(parameters);
#ifdef ROTATION
    // If gamma_rotation is not set explicitly, use the translational one.
    if (params.count("gamma_rotation") == 0ul and params.count("gamma")) {
      params["gamma_rotation"] = params.at("gamma");
    }
#endif // ROTATION
    return params;
  }
};

class Brownian : public Interface<::BrownianThermostat> {
  std::shared_ptr<CoreThermostat> &
  get_member_handle(::Thermostat::Thermostat &thermostat) override {
    return thermostat.brownian;
  }

public:
  Brownian() {
    add_parameters({
        make_autogamma(&CoreThermostat::gamma, "gamma"),
#ifdef ROTATION
        make_autogamma(&CoreThermostat::gamma_rotation, "gamma_rotation"),
#endif
    });
  }

  ::ThermostatFlags get_thermo_flag() const final { return THERMO_BROWNIAN; }

protected:
  VariantMap extend_parameters(VariantMap const &parameters) const override {
    auto params =
        Interface<::BrownianThermostat>::extend_parameters(parameters);
#ifdef ROTATION
    // If gamma_rotation is not set explicitly, use the translational one.
    if (params.count("gamma_rotation") == 0ul and params.count("gamma")) {
      params["gamma_rotation"] = params.at("gamma");
    }
#endif // ROTATION
    return params;
  }
};

#ifdef NPT
class IsotropicNpt : public Interface<::IsotropicNptThermostat> {
  std::shared_ptr<CoreThermostat> &
  get_member_handle(::Thermostat::Thermostat &thermostat) override {
    return thermostat.npt_iso;
  }

public:
  IsotropicNpt() {
    add_parameters({
        make_autoparameter(&CoreThermostat::gamma0, "gamma0"),
        make_autoparameter(&CoreThermostat::gammav, "gammav"),
    });
  }

  ::ThermostatFlags get_thermo_flag() const final { return THERMO_NPT_ISO; }
};
#endif // NPT

#ifdef WALBERLA
class LBThermostat : public Interface<::LBThermostat> {
  std::shared_ptr<CoreThermostat> &
  get_member_handle(::Thermostat::Thermostat &thermostat) override {
    return thermostat.lb;
  }

public:
  LBThermostat() {
    add_parameters({
        {"gamma",
         [this](Variant const &v) {
           check_lock();
           if (is_none(v)) {
             return;
           }
           auto const gamma = get_value<double>(v);
           context()->parallel_try_catch(
               [&]() { sanity_checks_positive(gamma, "gamma"); });
           m_handle->gamma = gamma;
         },
         [this]() {
           auto const gamma = m_handle->gamma;
           return (gamma >= 0.) ? Variant{gamma} : Variant{None{}};
         }},
    });
  }

  ::ThermostatFlags get_thermo_flag() const final { return THERMO_LB; }

  std::optional<double> extract_kT(VariantMap const &params) const override {
    auto const obj =
        get_value<std::shared_ptr<walberla::LBFluid>>(params, "LB_fluid");
    auto const value = get_value<double>(obj->get_parameter("kT"));
    sanity_checks_positive(value, "kT");
    return value;
  }

protected:
  bool invalid_rng_state(VariantMap const &params) const override {
    return (not params.count("seed") or is_none(params.at("seed"))) and
           params.count("__global_kT") and is_seed_required() and
           get_value<double>(params, "__global_kT") != 0.;
  }
};
#endif // WALBERLA

#ifdef DPD
class DPDThermostat : public Interface<::DPDThermostat> {
  std::shared_ptr<CoreThermostat> &
  get_member_handle(::Thermostat::Thermostat &thermostat) override {
    return thermostat.dpd;
  }

public:
  ::ThermostatFlags get_thermo_flag() const final { return THERMO_DPD; }
};
#endif

#ifdef STOKESIAN_DYNAMICS
class Stokesian : public Interface<::StokesianThermostat> {
  std::shared_ptr<CoreThermostat> &
  get_member_handle(::Thermostat::Thermostat &thermostat) override {
    return thermostat.stokesian;
  }

public:
  ::ThermostatFlags get_thermo_flag() const final { return THERMO_SD; }
};
#endif

class ThermalizedBond : public Interface<::ThermalizedBondThermostat> {
  std::shared_ptr<CoreThermostat> &
  get_member_handle(::Thermostat::Thermostat &thermostat) override {
    return thermostat.thermalized_bond;
  }

public:
  ::ThermostatFlags get_thermo_flag() const final { return THERMO_BOND; }
};

class Thermostat : public AutoParameters<Thermostat, System::Leaf> {
  std::shared_ptr<Langevin> langevin;
  std::shared_ptr<Brownian> brownian;
#ifdef NPT
  std::shared_ptr<IsotropicNpt> npt_iso;
#endif
#ifdef WALBERLA
  std::shared_ptr<LBThermostat> lb;
#endif
#ifdef DPD
  std::shared_ptr<DPDThermostat> dpd;
#endif
#ifdef STOKESIAN_DYNAMICS
  std::shared_ptr<Stokesian> stokesian;
#endif
  std::shared_ptr<ThermalizedBond> thermalized_bond;
  std::shared_ptr<::Thermostat::Thermostat> m_handle;
  std::unique_ptr<VariantMap> m_params;

  template <typename Fun> void apply(Fun fun) {
    fun(*langevin);
    fun(*brownian);
#ifdef NPT
    fun(*npt_iso);
#endif
#ifdef WALBERLA
    fun(*lb);
#endif
#ifdef DPD
    fun(*dpd);
#endif
#ifdef STOKESIAN_DYNAMICS
    fun(*stokesian);
#endif
    fun(*thermalized_bond);
  }

protected:
  template <typename T>
  auto make_autoparameter(T Thermostat::*member, char const *name) {
    return AutoParameter{
        name,
        [this, member, name = std::string(name)](Variant const &v) {
          auto &thermostat = this->*member;
          if (thermostat) {
            throw WriteError{name};
          }
          thermostat = get_value<T>(v);
        },
        [this, member]() { return this->*member; }};
  }

  template <typename T>
  void setup_thermostat(std::shared_ptr<T> &thermostat,
                        VariantMap const &params) {
    auto const original_kT = m_handle->kT;
    std::optional<double> new_kT;
    context()->parallel_try_catch(
        [&]() { new_kT = thermostat->extract_kT(params); });
    auto const thermo_flag = thermostat->get_thermo_flag();
    if (new_kT) {
      context()->parallel_try_catch(
          [&]() { update_global_kT(original_kT, *new_kT, thermo_flag); });
    }
    auto const was_active = thermostat->is_activated();
    turn_thermostat_off(*thermostat);
    auto read_write_lock = thermostat->release_lock();
    context()->parallel_try_catch([&]() {
      try {
        if (new_kT) {
          m_handle->kT = *new_kT;
        }
        thermostat->update_and_bind(params, was_active, m_system.lock());
        m_handle->thermo_switch |= thermo_flag;
      } catch (...) {
        auto success = false;
        try {
          m_handle->kT = original_kT;
          if (was_active) {
            m_handle->thermo_switch |= thermo_flag;
            thermostat->bind_system(m_system.lock());
          }
          success = true;
          throw success;
        } catch (...) {
          assert(success &&
                 "An exception occurred when setting up the thermostat. "
                 "An exception also occurred when attempting to restore the "
                 "original thermostat. The system is now in an invalid state.");
        }
        throw;
      }
    });
  }

  template <typename T> void turn_thermostat_off(T &thermostat) {
    auto const thermo_flag = thermostat.get_thermo_flag();
    if (m_handle->thermo_switch & thermo_flag) {
      thermostat.detach_system();
      m_handle->thermo_switch &= ~thermo_flag;
      if (m_handle->thermo_switch == 0) {
        m_handle->kT = -1.;
      }
    }
  }

  void update_global_kT(double old_kT, double new_kT, int thermo_flag) {
    if (new_kT >= 0.) {
      auto const same_kT = ::Thermostat::are_kT_equal(old_kT, new_kT);
      auto const thermo_switch = m_handle->thermo_switch;
      if (thermo_switch != THERMO_OFF and thermo_switch != thermo_flag and
          thermo_switch != THERMO_BOND and not same_kT and old_kT >= 0.) {
        throw std::runtime_error(
            "Cannot set parameter 'kT' to " + std::to_string(new_kT) +
            ": there are currently active thermostats with kT=" +
            std::to_string(old_kT));
      }
      get_system().check_kT(new_kT);
      if (not same_kT) {
        m_handle->kT = new_kT;
        get_system().on_temperature_change();
      }
    }
  }

public:
  Thermostat() {
    add_parameters({
        {"kT", AutoParameter::read_only,
         [this]() {
           return (m_handle->kT >= 0.) ? Variant{m_handle->kT}
                                       : Variant{None{}};
         }},
        make_autoparameter(&Thermostat::langevin, "langevin"),
        make_autoparameter(&Thermostat::brownian, "brownian"),
#ifdef NPT
        make_autoparameter(&Thermostat::npt_iso, "npt_iso"),
#endif
#ifdef WALBERLA
        make_autoparameter(&Thermostat::lb, "lb"),
#endif
#ifdef DPD
        make_autoparameter(&Thermostat::dpd, "dpd"),
#endif
#ifdef STOKESIAN_DYNAMICS
        make_autoparameter(&Thermostat::stokesian, "stokesian"),
#endif
        make_autoparameter(&Thermostat::thermalized_bond, "thermalized_bond"),
    });
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "set_langevin") {
      setup_thermostat(langevin, params);
      return {};
    }
    if (name == "set_brownian") {
      setup_thermostat(brownian, params);
      return {};
    }
#ifdef NPT
    if (name == "set_npt") {
      setup_thermostat(npt_iso, params);
      return {};
    }
#endif // NPT
#ifdef WALBERLA
    if (name == "set_lb") {
      setup_thermostat(lb, params);
      return {};
    }
#endif // WALBERLA
#ifdef DPD
    if (name == "set_dpd") {
      setup_thermostat(dpd, params);
      return {};
    }
#endif // DPD
#ifdef STOKESIAN_DYNAMICS
    if (name == "set_stokesian") {
      setup_thermostat(stokesian, params);
      return {};
    }
#endif // STOKESIAN_DYNAMICS
    if (name == "set_thermalized_bond") {
      setup_thermostat(thermalized_bond, params);
      return {};
    }
    if (name == "turn_off") {
      apply([this](auto &thermostat) { turn_thermostat_off(thermostat); });
      assert(m_handle->thermo_switch == THERMO_OFF);
      get_system().on_temperature_change();
      return {};
    }
    return {};
  }

  void do_construct(VariantMap const &params) override {
    m_params = std::make_unique<VariantMap>(params);
  }

  void on_bind_system(::System::System &system) override {
    assert(m_params != nullptr);
    m_handle = system.thermostat;
    auto const &params = *m_params;
    if (not params.empty()) {
      reload_checkpointed_thermostats(params);
      m_params.reset();
      return;
    }
    m_params.reset();
    if (not context()->is_head_node()) {
      return;
    }
    make_default_constructed_thermostats();
  }

private:
  /**
   * @brief Reload thermostats from checkpointed data.
   */
  void reload_checkpointed_thermostats(VariantMap const &params) {
    for (auto const &key : get_parameter_insertion_order()) {
      if (key != "kT") {
        auto const &v = params.at(key);
        do_set_parameter(key.c_str(), v);
      }
    }
    if (not is_none(params.at("kT"))) {
      m_handle->kT = get_value<double>(params, "kT");
    }
    apply([this](auto &thermostat) {
      if (get_value<bool>(thermostat.get_parameter("is_active"))) {
        thermostat.bind_system(m_system.lock());
        m_handle->thermo_switch |= thermostat.get_thermo_flag();
      }
    });
    get_system().on_thermostat_param_change();
  }

  /**
   * @brief Instantiate default-contructed thermostats.
   * Can only be run on the head node!
   */
  void make_default_constructed_thermostats() {
    assert(context()->is_head_node());
    auto const make_thermostat = [this](char const *name, char const *so_name) {
      set_parameter(name, Variant{context()->make_shared(so_name, {})});
    };
    make_thermostat("langevin", "Thermostat::Langevin");
    make_thermostat("brownian", "Thermostat::Brownian");
#ifdef NPT
    make_thermostat("npt_iso", "Thermostat::IsotropicNpt");
#endif
#ifdef WALBERLA
    make_thermostat("lb", "Thermostat::LB");
#endif
#ifdef DPD
    make_thermostat("dpd", "Thermostat::DPD");
#endif
#ifdef STOKESIAN_DYNAMICS
    make_thermostat("stokesian", "Thermostat::Stokesian");
#endif
    make_thermostat("thermalized_bond", "Thermostat::ThermalizedBond");
  }
};

} // namespace Thermostat
} // namespace ScriptInterface
