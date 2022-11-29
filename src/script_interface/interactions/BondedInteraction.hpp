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

/** @file
 *  The ScriptInterface counterparts of the bonded interactions parameters
 *  structs from the core are defined here.
 *
 */

#ifndef SCRIPT_INTERFACE_INTERACTIONS_BONDED_INTERACTION_HPP
#define SCRIPT_INTERFACE_INTERACTIONS_BONDED_INTERACTION_HPP

#include "core/bonded_interactions/bonded_interaction_data.hpp"
#include "core/immersed_boundaries.hpp"
#include "core/thermostat.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/get_value.hpp"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/variant.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Interactions {

class BondedInteraction : public AutoParameters<BondedInteraction> {
protected:
  std::shared_ptr<::Bonded_IA_Parameters> m_bonded_ia;

public:
  std::shared_ptr<::Bonded_IA_Parameters> bonded_ia() { return m_bonded_ia; }
  std::shared_ptr<const ::Bonded_IA_Parameters> bonded_ia() const {
    return m_bonded_ia;
  }

protected:
  using AutoParameters<BondedInteraction>::context;
  using AutoParameters<BondedInteraction>::valid_parameters;

  virtual std::set<std::string> get_valid_parameters() const {
    auto const vec = valid_parameters();
    auto valid_keys = std::set<std::string>();
    std::transform(vec.begin(), vec.end(),
                   std::inserter(valid_keys, valid_keys.begin()),
                   [](auto const &key) { return std::string{key}; });
    return valid_keys;
  }

private:
  void check_valid_parameters(VariantMap const &params) const {
    auto const valid_keys = get_valid_parameters();
    for (auto const &key : valid_keys) {
      if (params.count(std::string(key)) == 0) {
        throw std::runtime_error("Parameter '" + key + "' is missing");
      }
    }
    for (auto const &kv : params) {
      if (valid_keys.count(kv.first) == 0) {
        throw std::runtime_error("Parameter '" + kv.first +
                                 "' is not recognized");
      }
    }
  }

  void do_construct(VariantMap const &params) override {
    // Check if initialization "by id" or "by parameters"
    if (params.find("bond_id") != params.end()) {
      auto const bond_id = get_value<int>(params, "bond_id");
      context()->parallel_try_catch([&]() {
        if (not ::bonded_ia_params.contains(bond_id)) {
          throw std::runtime_error("No bond with id " +
                                   std::to_string(bond_id) +
                                   " exists in the ESPResSo core");
        }
      });
      m_bonded_ia = ::bonded_ia_params.at(bond_id);
    } else {
      context()->parallel_try_catch([&]() {
        check_valid_parameters(params);
        construct_bond(params);
      });
    }
  }

  virtual void construct_bond(VariantMap const &params) = 0;

public:
  bool operator==(BondedInteraction const &other) {
    return m_bonded_ia == other.m_bonded_ia;
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    // this feature is needed to compare bonds
    if (name == "is_same_bond") {
      auto const bond_so =
          get_value<std::shared_ptr<BondedInteraction>>(params, "bond");
      return *this == *bond_so;
    }
    if (name == "get_num_partners") {
      return number_of_partners(*bonded_ia());
    }
    if (name == "get_zero_based_type") {
      auto const bond_id = get_value<int>(params, "bond_id");
      return ::bonded_ia_params.get_zero_based_type(bond_id);
    }

    return {};
  }
};

template <class CoreIA> class BondedInteractionImpl : public BondedInteraction {

public:
  using CoreBondedInteraction = CoreIA;
  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }
};

class FeneBond : public BondedInteractionImpl<::FeneBond> {
public:
  FeneBond() {
    add_parameters({
        {"k", AutoParameter::read_only, [this]() { return get_struct().k; }},
        {"d_r_max", AutoParameter::read_only,
         [this]() { return get_struct().drmax; }},
        {"r_0", AutoParameter::read_only, [this]() { return get_struct().r0; }},
    });
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia = std::make_shared<::Bonded_IA_Parameters>(
        CoreBondedInteraction(get_value<double>(params, "k"),
                              get_value<double>(params, "d_r_max"),
                              get_value<double>(params, "r_0")));
  }
};

class HarmonicBond : public BondedInteractionImpl<::HarmonicBond> {
public:
  HarmonicBond() {
    add_parameters({
        {"k", AutoParameter::read_only, [this]() { return get_struct().k; }},
        {"r_0", AutoParameter::read_only, [this]() { return get_struct().r; }},
        {"r_cut", AutoParameter::read_only,
         [this]() { return get_struct().r_cut; }},
    });
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction(
            get_value<double>(params, "k"), get_value<double>(params, "r_0"),
            get_value<double>(params, "r_cut")));
  }
};

class QuarticBond : public BondedInteractionImpl<::QuarticBond> {
public:
  QuarticBond() {
    add_parameters({
        {"k0", AutoParameter::read_only, [this]() { return get_struct().k0; }},
        {"k1", AutoParameter::read_only, [this]() { return get_struct().k1; }},
        {"r", AutoParameter::read_only, [this]() { return get_struct().r; }},
        {"r_cut", AutoParameter::read_only,
         [this]() { return get_struct().r_cut; }},
    });
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction(
            get_value<double>(params, "k0"), get_value<double>(params, "k1"),
            get_value<double>(params, "r"),
            get_value<double>(params, "r_cut")));
  }
};

class BondedCoulomb : public BondedInteractionImpl<::BondedCoulomb> {
public:
  BondedCoulomb() {
    add_parameters({
        {"prefactor", AutoParameter::read_only,
         [this]() { return get_struct().prefactor; }},
    });
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia = std::make_shared<::Bonded_IA_Parameters>(
        CoreBondedInteraction(get_value<double>(params, "prefactor")));
  }
};

class BondedCoulombSR : public BondedInteractionImpl<::BondedCoulombSR> {
public:
  BondedCoulombSR() {
    add_parameters({
        {"q1q2", AutoParameter::read_only,
         [this]() { return get_struct().q1q2; }},
    });
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia = std::make_shared<::Bonded_IA_Parameters>(
        CoreBondedInteraction(get_value<double>(params, "q1q2")));
  }
};

class AngleHarmonicBond : public BondedInteractionImpl<::AngleHarmonicBond> {
public:
  AngleHarmonicBond() {
    add_parameters({
        {"bend", AutoParameter::read_only,
         [this]() { return get_struct().bend; }},
        {"phi0", AutoParameter::read_only,
         [this]() { return get_struct().phi0; }},
    });
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia = std::make_shared<::Bonded_IA_Parameters>(
        CoreBondedInteraction(get_value<double>(params, "bend"),
                              get_value<double>(params, "phi0")));
  }
};

class AngleCosineBond : public BondedInteractionImpl<::AngleCosineBond> {
public:
  AngleCosineBond() {
    add_parameters({
        {"bend", AutoParameter::read_only,
         [this]() { return get_struct().bend; }},
        {"phi0", AutoParameter::read_only,
         [this]() { return get_struct().phi0; }},
    });
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia = std::make_shared<::Bonded_IA_Parameters>(
        CoreBondedInteraction(get_value<double>(params, "bend"),
                              get_value<double>(params, "phi0")));
  }
};

class AngleCossquareBond : public BondedInteractionImpl<::AngleCossquareBond> {
public:
  AngleCossquareBond() {
    add_parameters({
        {"bend", AutoParameter::read_only,
         [this]() { return get_struct().bend; }},
        {"phi0", AutoParameter::read_only,
         [this]() { return get_struct().phi0; }},
    });
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia = std::make_shared<::Bonded_IA_Parameters>(
        CoreBondedInteraction(get_value<double>(params, "bend"),
                              get_value<double>(params, "phi0")));
  }
};

class DihedralBond : public BondedInteractionImpl<::DihedralBond> {
public:
  DihedralBond() {
    add_parameters({
        {"mult", AutoParameter::read_only,
         [this]() { return get_struct().mult; }},
        {"bend", AutoParameter::read_only,
         [this]() { return get_struct().bend; }},
        {"phase", AutoParameter::read_only,
         [this]() { return get_struct().phase; }},
    });
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction(
            get_value<int>(params, "mult"), get_value<double>(params, "bend"),
            get_value<double>(params, "phase")));
  }
};

class TabulatedDistanceBond
    : public BondedInteractionImpl<::TabulatedDistanceBond> {
public:
  TabulatedDistanceBond() {
    add_parameters({
        {"min", AutoParameter::read_only,
         [this]() { return get_struct().pot->minval; }},
        {"max", AutoParameter::read_only,
         [this]() { return get_struct().pot->maxval; }},
        {"energy", AutoParameter::read_only,
         [this]() { return get_struct().pot->energy_tab; }},
        {"force", AutoParameter::read_only,
         [this]() { return get_struct().pot->force_tab; }},
    });
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction(
            get_value<double>(params, "min"), get_value<double>(params, "max"),
            get_value<std::vector<double>>(params, "energy"),
            get_value<std::vector<double>>(params, "force")));
  }
};

class TabulatedAngleBond : public BondedInteractionImpl<::TabulatedAngleBond> {
public:
  TabulatedAngleBond() {
    add_parameters({
        {"min", AutoParameter::read_only,
         [this]() { return get_struct().pot->minval; }},
        {"max", AutoParameter::read_only,
         [this]() { return get_struct().pot->maxval; }},
        {"energy", AutoParameter::read_only,
         [this]() { return get_struct().pot->energy_tab; }},
        {"force", AutoParameter::read_only,
         [this]() { return get_struct().pot->force_tab; }},
    });
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction(
            get_value<double>(params, "min"), get_value<double>(params, "max"),
            get_value<std::vector<double>>(params, "energy"),
            get_value<std::vector<double>>(params, "force")));
  }
};

class TabulatedDihedralBond
    : public BondedInteractionImpl<::TabulatedDihedralBond> {
public:
  TabulatedDihedralBond() {
    add_parameters({
        {"min", AutoParameter::read_only,
         [this]() { return get_struct().pot->minval; }},
        {"max", AutoParameter::read_only,
         [this]() { return get_struct().pot->maxval; }},
        {"energy", AutoParameter::read_only,
         [this]() { return get_struct().pot->energy_tab; }},
        {"force", AutoParameter::read_only,
         [this]() { return get_struct().pot->force_tab; }},
    });
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction(
            get_value<double>(params, "min"), get_value<double>(params, "max"),
            get_value<std::vector<double>>(params, "energy"),
            get_value<std::vector<double>>(params, "force")));
  }
};

class ThermalizedBond : public BondedInteractionImpl<::ThermalizedBond> {
public:
  ThermalizedBond() {
    add_parameters({
        {"temp_com", AutoParameter::read_only,
         [this]() { return get_struct().temp_com; }},
        {"gamma_com", AutoParameter::read_only,
         [this]() { return get_struct().gamma_com; }},
        {"temp_distance", AutoParameter::read_only,
         [this]() { return get_struct().temp_distance; }},
        {"gamma_distance", AutoParameter::read_only,
         [this]() { return get_struct().gamma_distance; }},
        {"r_cut", AutoParameter::read_only,
         [this]() { return get_struct().r_cut; }},
        {"seed", AutoParameter::read_only,
         []() { return static_cast<size_t>(thermalized_bond.rng_seed()); }},
    });
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "get_rng_state") {
      auto const state = ::thermalized_bond.rng_counter();
      // check it's safe to forward the current state as an integer
      assert(state < static_cast<uint64_t>(std::numeric_limits<int>::max()));
      return static_cast<int>(state);
    }
    return BondedInteraction::do_call_method(name, params);
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia = std::make_shared<::Bonded_IA_Parameters>(
        CoreBondedInteraction(get_value<double>(params, "temp_com"),
                              get_value<double>(params, "gamma_com"),
                              get_value<double>(params, "temp_distance"),
                              get_value<double>(params, "gamma_distance"),
                              get_value<double>(params, "r_cut")));

    if (is_none(params.at("seed"))) {
      if (::thermalized_bond.is_seed_required()) {
        throw std::invalid_argument("A parameter 'seed' has to be given on "
                                    "first activation of a thermalized bond");
      }
    } else {
      auto const seed = get_value<int>(params, "seed");
      if (seed < 0) {
        throw std::domain_error("Parameter 'seed' must be >= 0");
      }
      ::thermalized_bond.rng_initialize(static_cast<uint32_t>(seed));
    }

    // handle checkpointing
    if (params.count("rng_state") and not is_none(params.at("rng_state"))) {
      auto const state = get_value<int>(params, "rng_state");
      assert(state >= 0);
      ::thermalized_bond.set_rng_counter(static_cast<uint64_t>(state));
    }
  }

  std::set<std::string> get_valid_parameters() const override {
    auto names =
        BondedInteractionImpl<CoreBondedInteraction>::get_valid_parameters();
    names.insert("rng_state");
    return names;
  }
};

class RigidBond : public BondedInteractionImpl<::RigidBond> {
public:
  RigidBond() {
    add_parameters({
        {"r", AutoParameter::read_only,
         [this]() { return std::sqrt(get_struct().d2); }},
        {"ptol", AutoParameter::read_only,
         [this]() { return 0.5 * get_struct().p_tol; }},
        {"vtol", AutoParameter::read_only,
         [this]() { return get_struct().v_tol; }},
    });
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction(
            get_value<double>(params, "r"), get_value<double>(params, "ptol"),
            get_value<double>(params, "vtol")));
  }
};

class IBMTriel : public BondedInteractionImpl<::IBMTriel> {
public:
  IBMTriel() {
    add_parameters({
        {"k1", AutoParameter::read_only, [this]() { return get_struct().k1; }},
        {"k2", AutoParameter::read_only, [this]() { return get_struct().k2; }},
        {"maxDist", AutoParameter::read_only,
         [this]() { return get_struct().maxDist; }},
        {"elasticLaw", AutoParameter::read_only,
         [this]() {
           if (get_struct().elasticLaw == tElasticLaw::NeoHookean) {
             return std::string("NeoHookean");
           }
           return std::string("Skalak");
         }},
    });
  }

private:
  void construct_bond(VariantMap const &params) override {
    auto const law_name = get_value<std::string>(params, "elasticLaw");
    tElasticLaw elastic_law;
    if (law_name == "NeoHookean") {
      elastic_law = tElasticLaw::NeoHookean;
    } else if (law_name == "Skalak") {
      elastic_law = tElasticLaw::Skalak;
    } else {
      throw std::invalid_argument(
          "Invalid value for parameter 'elasticLaw': '" + law_name + "'");
    }
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction(
            get_value<int>(params, "ind1"), get_value<int>(params, "ind2"),
            get_value<int>(params, "ind3"),
            get_value<double>(params, "maxDist"), elastic_law,
            get_value<double>(params, "k1"), get_value<double>(params, "k2")));
  }

  std::set<std::string> get_valid_parameters() const override {
    auto names =
        BondedInteractionImpl<CoreBondedInteraction>::get_valid_parameters();
    names.insert("ind1");
    names.insert("ind2");
    names.insert("ind3");
    return names;
  }
};

class IBMVolCons : public BondedInteractionImpl<::IBMVolCons> {
public:
  IBMVolCons() {
    add_parameters({
        {"softID", AutoParameter::read_only,
         [this]() { return get_struct().softID; }},
        {"kappaV", AutoParameter::read_only,
         [this]() { return get_struct().kappaV; }},
    });
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "current_volume") {
      return ::immersed_boundaries.get_current_volume(get_struct().softID);
    }
    return BondedInteraction::do_call_method(name, params);
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia = std::make_shared<::Bonded_IA_Parameters>(
        CoreBondedInteraction(get_value<int>(params, "softID"),
                              get_value<double>(params, "kappaV")));
  }
};

class IBMTribend : public BondedInteractionImpl<::IBMTribend> {
public:
  IBMTribend() {
    add_parameters({
        {"kb", AutoParameter::read_only, [this]() { return get_struct().kb; }},
        {"refShape", AutoParameter::read_only,
         [this]() {
           return (m_flat) ? std::string("Flat") : std::string("Initial");
         }},
        {"theta0", AutoParameter::read_only,
         [this]() { return get_struct().theta0; }},
    });
  }

private:
  bool m_flat;
  void construct_bond(VariantMap const &params) override {
    auto const shape_name = get_value<std::string>(params, "refShape");
    if (shape_name == "Flat") {
      m_flat = true;
    } else if (shape_name == "Initial") {
      m_flat = false;
    } else {
      throw std::invalid_argument("Invalid value for parameter 'refShape': '" +
                                  shape_name + "'");
    }
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction(
            get_value<int>(params, "ind1"), get_value<int>(params, "ind2"),
            get_value<int>(params, "ind3"), get_value<int>(params, "ind4"),
            get_value<double>(params, "kb"), m_flat));
  }

  std::set<std::string> get_valid_parameters() const override {
    auto names =
        BondedInteractionImpl<CoreBondedInteraction>::get_valid_parameters();
    names.erase("theta0");
    names.insert("ind1");
    names.insert("ind2");
    names.insert("ind3");
    names.insert("ind4");
    return names;
  }
};

class OifGlobalForcesBond
    : public BondedInteractionImpl<::OifGlobalForcesBond> {
public:
  OifGlobalForcesBond() {
    add_parameters({
        {"A0_g", AutoParameter::read_only,
         [this]() { return get_struct().A0_g; }},
        {"ka_g", AutoParameter::read_only,
         [this]() { return get_struct().ka_g; }},
        {"V0", AutoParameter::read_only, [this]() { return get_struct().V0; }},
        {"kv", AutoParameter::read_only, [this]() { return get_struct().kv; }},
    });
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction(
            get_value<double>(params, "A0_g"),
            get_value<double>(params, "ka_g"), get_value<double>(params, "V0"),
            get_value<double>(params, "kv")));
  }
};

class OifLocalForcesBond : public BondedInteractionImpl<::OifLocalForcesBond> {
public:
  OifLocalForcesBond() {
    add_parameters({
        {"r0", AutoParameter::read_only, [this]() { return get_struct().r0; }},
        {"ks", AutoParameter::read_only, [this]() { return get_struct().ks; }},
        {"kslin", AutoParameter::read_only,
         [this]() { return get_struct().kslin; }},
        {"phi0", AutoParameter::read_only,
         [this]() { return get_struct().phi0; }},
        {"kb", AutoParameter::read_only, [this]() { return get_struct().kb; }},
        {"A01", AutoParameter::read_only,
         [this]() { return get_struct().A01; }},
        {"A02", AutoParameter::read_only,
         [this]() { return get_struct().A02; }},
        {"kal", AutoParameter::read_only,
         [this]() { return get_struct().kal; }},
        {"kvisc", AutoParameter::read_only,
         [this]() { return get_struct().kvisc; }},
    });
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction(
            get_value<double>(params, "r0"), get_value<double>(params, "ks"),
            get_value<double>(params, "kslin"),
            get_value<double>(params, "phi0"), get_value<double>(params, "kb"),
            get_value<double>(params, "A01"), get_value<double>(params, "A02"),
            get_value<double>(params, "kal"),
            get_value<double>(params, "kvisc")));
  }
};

class VirtualBond : public BondedInteractionImpl<::VirtualBond> {
public:
  VirtualBond() {
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction());
  }

private:
  void construct_bond(VariantMap const &) override {}
};

} // namespace Interactions
} // namespace ScriptInterface

#endif
