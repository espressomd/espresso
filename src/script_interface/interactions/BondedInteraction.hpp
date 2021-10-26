/*
 * Copyright (C) 2021 The ESPResSo project
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
#include "core/thermostat.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/get_value.hpp"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/variant.hpp>

#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Interactions {

template <class T> class BondedInteractionInterface {
protected:
  std::shared_ptr<::Bonded_IA_Parameters> m_bonded_ia;

public:
  std::shared_ptr<::Bonded_IA_Parameters> bonded_ia() { return m_bonded_ia; }
  std::shared_ptr<const ::Bonded_IA_Parameters> bonded_ia() const {
    return m_bonded_ia;
  }
};

class BondedInteraction : public AutoParameters<BondedInteraction>,
                          public BondedInteractionInterface<BondedInteraction> {
  void do_construct(VariantMap const &params) override {
    // Check if initialization "by id" or "by parameters"
    if (params.find("bond_id") != params.end()) {
      m_bonded_ia = ::bonded_ia_params.at(get_value<int>(params, "bond_id"));
    } else {
      construct_bond(params);
    }
  }

private:
  virtual void construct_bond(VariantMap const &params) {}

public:
  template <typename T>
  bool operator==(const BondedInteractionInterface<T> &other) {
    return m_bonded_ia == other.m_bonded_ia;
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    // this feature is needed to compare bonds
    if (name == "get_address") {
      return reinterpret_cast<std::size_t>(bonded_ia().get());
    }
    if (name == "get_num_partners") {
      return number_of_partners(*bonded_ia());
    }

    return {};
  }
};

class FeneBond : public BondedInteraction {
  using CoreBondedInteraction = ::FeneBond;

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

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }
};

class HarmonicBond : public BondedInteraction {
  using CoreBondedInteraction = ::HarmonicBond;

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

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }
};

class QuarticBond : public BondedInteraction {
  using CoreBondedInteraction = ::QuarticBond;

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

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }
};

class BondedCoulomb : public BondedInteraction {
  using CoreBondedInteraction = ::BondedCoulomb;

public:
  BondedCoulomb() {
    add_parameters({
        {"prefactor", AutoParameter::read_only,
         [this]() { return get_struct().prefactor; }},
    });
  }

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia = std::make_shared<::Bonded_IA_Parameters>(
        CoreBondedInteraction(get_value<double>(params, "prefactor")));
  }
};

class BondedCoulombSR : public BondedInteraction {
  using CoreBondedInteraction = ::BondedCoulombSR;

public:
  BondedCoulombSR() {
    add_parameters({
        {"q1q2", AutoParameter::read_only,
         [this]() { return get_struct().q1q2; }},
    });
  }

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia = std::make_shared<::Bonded_IA_Parameters>(
        CoreBondedInteraction(get_value<double>(params, "q1q2")));
  }
};

class AngleHarmonicBond : public BondedInteraction {
  using CoreBondedInteraction = ::AngleHarmonicBond;

public:
  AngleHarmonicBond() {
    add_parameters({
        {"bend", AutoParameter::read_only,
         [this]() { return get_struct().bend; }},
        {"phi0", AutoParameter::read_only,
         [this]() { return get_struct().phi0; }},
    });
  }

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia = std::make_shared<::Bonded_IA_Parameters>(
        CoreBondedInteraction(get_value<double>(params, "bend"),
                              get_value<double>(params, "phi0")));
  }
};

class AngleCosineBond : public BondedInteraction {
  using CoreBondedInteraction = ::AngleCosineBond;

public:
  AngleCosineBond() {
    add_parameters({
        {"bend", AutoParameter::read_only,
         [this]() { return get_struct().bend; }},
        {"phi0", AutoParameter::read_only,
         [this]() { return get_struct().phi0; }},
    });
  }

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia = std::make_shared<::Bonded_IA_Parameters>(
        CoreBondedInteraction(get_value<double>(params, "bend"),
                              get_value<double>(params, "phi0")));
  }
};

class AngleCossquareBond : public BondedInteraction {
  using CoreBondedInteraction = ::AngleCossquareBond;

public:
  AngleCossquareBond() {
    add_parameters({
        {"bend", AutoParameter::read_only,
         [this]() { return get_struct().bend; }},
        {"phi0", AutoParameter::read_only,
         [this]() { return get_struct().phi0; }},
    });
  }

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia = std::make_shared<::Bonded_IA_Parameters>(
        CoreBondedInteraction(get_value<double>(params, "bend"),
                              get_value<double>(params, "phi0")));
  }
};

class DihedralBond : public BondedInteraction {
  using CoreBondedInteraction = ::DihedralBond;

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

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction(
            get_value<int>(params, "mult"), get_value<double>(params, "bend"),
            get_value<double>(params, "phase")));
  }
};

class TabulatedDistanceBond : public BondedInteraction {
  using CoreBondedInteraction = ::TabulatedDistanceBond;

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

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
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

class TabulatedAngleBond : public BondedInteraction {
  using CoreBondedInteraction = ::TabulatedAngleBond;

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

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
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

class TabulatedDihedralBond : public BondedInteraction {
  using CoreBondedInteraction = ::TabulatedDihedralBond;

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

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
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

class ThermalizedBond : public BondedInteraction {
  using CoreBondedInteraction = ::ThermalizedBond;

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

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia = std::make_shared<::Bonded_IA_Parameters>(
        CoreBondedInteraction(get_value<double>(params, "temp_com"),
                              get_value<double>(params, "gamma_com"),
                              get_value<double>(params, "temp_distance"),
                              get_value<double>(params, "gamma_distance"),
                              get_value<double>(params, "r_cut")));
    thermalized_bond.rng_initialize(
        static_cast<uint32_t>(get_value<int>(params, "seed")));
  }
};

class RigidBond : public BondedInteraction {
  using CoreBondedInteraction = ::RigidBond;

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

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction(
            get_value<double>(params, "r"), get_value<double>(params, "ptol"),
            get_value<double>(params, "vtol")));
  }
};

class IBMTriel : public BondedInteraction {
  using CoreBondedInteraction = ::IBMTriel;

private:
  tElasticLaw str2elastic_law(std::string el) {
    if (boost::iequals(el, "NeoHookean")) {
      return tElasticLaw::NeoHookean;
    }
    return tElasticLaw::Skalak;
  }

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

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction(
            get_value<int>(params, "ind1"), get_value<int>(params, "ind2"),
            get_value<int>(params, "ind3"),
            get_value<double>(params, "maxDist"),
            str2elastic_law(get_value<std::string>(params, "elasticLaw")),
            get_value<double>(params, "k1"), get_value<double>(params, "k2")));
  }
};

class IBMVolCons : public BondedInteraction {
  using CoreBondedInteraction = ::IBMVolCons;

public:
  IBMVolCons() {
    add_parameters({
        {"softID", AutoParameter::read_only,
         [this]() { return get_struct().softID; }},
        {"kappaV", AutoParameter::read_only,
         [this]() { return get_struct().kappaV; }},
    });
  }

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }

private:
  void construct_bond(VariantMap const &params) override {
    m_bonded_ia = std::make_shared<::Bonded_IA_Parameters>(
        CoreBondedInteraction(get_value<int>(params, "softID"),
                              get_value<double>(params, "kappaV")));
  }
};

class IBMTribend : public BondedInteraction {
  using CoreBondedInteraction = ::IBMTribend;

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

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }

private:
  bool m_flat;
  void construct_bond(VariantMap const &params) override {
    auto const &refShape = get_value<std::string>(params, "refShape");
    m_flat = boost::iequals(refShape, "Flat");
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction(
            get_value<int>(params, "ind1"), get_value<int>(params, "ind2"),
            get_value<int>(params, "ind3"), get_value<int>(params, "ind4"),
            get_value<double>(params, "kb"), m_flat));
  }
};

class OifGlobalForcesBond : public BondedInteraction {
  using CoreBondedInteraction = ::OifGlobalForcesBond;

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

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
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

class OifLocalForcesBond : public BondedInteraction {
  using CoreBondedInteraction = ::OifLocalForcesBond;

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

  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
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

class VirtualBond : public BondedInteraction {
  using CoreBondedInteraction = ::VirtualBond;

public:
  VirtualBond() {
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction());
  }

public:
  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }
};

} // namespace Interactions
} // namespace ScriptInterface

#endif
