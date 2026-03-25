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

#include <config/config.hpp>

#include "ParticleHandle.hpp"

#include "script_interface/Variant.hpp"
#include "script_interface/cell_system/CellSystem.hpp"
#include "script_interface/get_value.hpp"
#include "script_interface/interactions/BondedInteractions.hpp"
#include "script_interface/system/System.hpp"

#include "core/BoxGeometry.hpp"
#include "core/PropagationMode.hpp"
#include "core/bonded_interactions/bonded_interaction_data.hpp"
#include "core/bonds.hpp"
#include "core/cell_system/CellStructure.hpp"
#include "core/exclusions.hpp"
#include "core/nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "core/particle_node.hpp"
#include "core/propagation.hpp"
#include "core/rotation.hpp"
#include "core/system/System.hpp"
#include "core/virtual_sites.hpp"
#include "core/virtual_sites/com.hpp"

#include <utils/Vector.hpp>
#include <utils/mpi/reduce_optional.hpp>

#include <boost/format.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/communicator.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <vector>

namespace ScriptInterface {
namespace Particles {

#ifdef ESPRESSO_ROTATION
static auto constexpr contradicting_arguments_quat = std::to_array<
    std::array<std::string_view, 3>>({
    {"dip", "dipm",
     "Setting 'dip' is sufficient as it defines the scalar dipole moment."},
    {"quat", "director",
     "Setting 'quat' is sufficient as it defines the director."},
    {"dip", "quat",
     "Setting 'dip' would overwrite 'quat'. Set 'quat' and 'dipm' instead."},
    {"dip", "director",
     "Setting 'dip' would overwrite 'director'. Set 'director' and "
     "'dipm' instead."},
});

static void sanity_checks_rotation(VariantMap const &params) {
  // if we are not constructing a particle from a checkpoint file,
  // check the quaternion is not accidentally set twice by the user
  if (not params.contains("__cpt_sentinel")) {
    auto formatter =
        boost::format("Contradicting particle attributes: '%s' and '%s'. %s");
    for (auto const &[prop1, prop2, reason] : contradicting_arguments_quat) {
      if (params.contains(std::string{prop1}) and
          params.contains(std::string{prop2})) {
        auto const err_msg = boost::str(formatter % prop1 % prop2 % reason);
        throw std::invalid_argument(err_msg);
      }
    }
  }
}
#endif // ESPRESSO_ROTATION

#if defined(ESPRESSO_ROTATION) or defined(ESPRESSO_EXTERNAL_FORCES)
static uint8_t bitfield_from_flag(Utils::Vector3i const &flag) {
  auto bitfield = static_cast<uint8_t>(0u);
  if (flag[0])
    bitfield |= static_cast<uint8_t>(1u);
  if (flag[1])
    bitfield |= static_cast<uint8_t>(2u);
  if (flag[2])
    bitfield |= static_cast<uint8_t>(4u);
  return bitfield;
}
#endif

#ifdef ESPRESSO_ROTATION
static auto quat2vector(Utils::Quaternion<double> const &q) {
  return Utils::Vector4d{{q[0], q[1], q[2], q[3]}};
}

static auto get_quaternion_safe(std::string const &name, Variant const &value) {
  auto const q = get_value<Utils::Vector4d>(value);
  if (q.norm2() == 0.) {
    throw std::domain_error(error_msg(name, "must be non-zero"));
  }
  return Utils::Quaternion<double>{{q[0], q[1], q[2], q[3]}};
}
#endif // ESPRESSO_ROTATION

#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
static auto get_gamma_safe(Variant const &value) {
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  try {
    return Utils::Vector3d::broadcast(get_value<double>(value));
  } catch (...) {
    return get_value<Utils::Vector3d>(value);
  }
#else  // ESPRESSO_PARTICLE_ANISOTROPY
  return get_value<double>(value);
#endif // ESPRESSO_PARTICLE_ANISOTROPY
}
#endif // ESPRESSO_THERMOSTAT_PER_PARTICLE

template <typename T, class F>
T ParticleHandle::get_particle_property(F const &fun) const {
  auto &cell_structure = get_cell_structure()->get_cell_structure();
  auto const &comm = context()->get_comm();
  auto const ptr = const_cast<Particle const *>(
      get_real_particle(comm, m_pid, cell_structure));
  std::optional<T> ret;
  if (ptr == nullptr) {
    ret = {};
  } else {
    ret = {fun(*ptr)};
  }
  return Utils::Mpi::reduce_optional(comm, ret);
}

template <typename T>
T ParticleHandle::get_particle_property(T const &(Particle::*getter)()
                                            const) const {
  return get_particle_property<T>(
      [getter](Particle const &p) { return (p.*getter)(); });
}

template <class F>
void ParticleHandle::set_particle_property(F const &fun) const {
  auto &cell_structure = get_cell_structure()->get_cell_structure();
  auto const &comm = context()->get_comm();
  auto const ptr = get_real_particle(comm, m_pid, cell_structure);
  if (ptr != nullptr) {
    fun(*ptr);
  }
  get_system()->on_particle_change();
}

template <typename T>
void ParticleHandle::set_particle_property(T &(Particle::*setter)(),
                                           Variant const &value) const {
  set_particle_property(
      [&value, setter](Particle &p) { (p.*setter)() = get_value<T>(value); });
}

ParticleHandle::ParticleHandle() {
  /* Warning: the order of particle property setters matters! Some properties
   * override each other, e.g. quat/director/dip or dip/dipm.
   * This is relevant during checkpointing: all particle properties are set at
   * once, in the order specified in the following call to `add_parameters()`.
   */
  add_parameters({
      {"id", AutoParameter::read_only, [this]() { return m_pid; }},
      {"type",
       [this](Variant const &value) {
         auto const old_type = get_particle_property(&Particle::type);
         auto const new_type = get_value<int>(value);
         if (new_type < 0) {
           throw std::domain_error(
               error_msg("type", "must be an integer >= 0"));
         }
         get_system()->nonbonded_ias->make_particle_type_exist(new_type);
         on_particle_type_change(m_pid, old_type, new_type);
         set_particle_property(&Particle::type, value);
       },
       [this]() { return get_particle_data(m_pid).type(); }},
      {"pos",
       [this](Variant const &value) {
         auto const pos = get_value<Utils::Vector3d>(value);
         particle_checks(m_pid, pos);
         set_particle_pos(m_pid, pos);
       },
       [this]() {
         auto const p = get_particle_data(m_pid);
         auto const pos = p.pos();
         auto const image_box = p.image_box();
         return get_system()->box_geo->unfolded_position(pos, image_box);
       }},
      {"v",
       [this](Variant const &value) {
         set_particle_property(&Particle::v, value);
       },
       [this]() { return get_particle_data(m_pid).v(); }},
      {"f",
       [this](Variant const &value) {
         set_particle_property(&Particle::force, value);
       },
       [this]() { return get_particle_data(m_pid).force(); }},
      {"mass",
#ifdef ESPRESSO_MASS
       [this](Variant const &value) {
         if (get_value<double>(value) <= 0.) {
           throw std::domain_error(error_msg("mass", "must be a float > 0"));
         }
         set_particle_property(&Particle::mass, value);
       },
#else  // ESPRESSO_MASS
       [](Variant const &value) {
         auto const default_mass = Particle().mass();
         if (std::abs(get_value<double>(value) - default_mass) > 1e-10) {
           throw std::runtime_error("Feature MASS not compiled in");
         }
       },
#endif // ESPRESSO_MASS
       [this]() { return get_particle_data(m_pid).mass(); }},
      {"q",
#ifdef ESPRESSO_ELECTROSTATICS
       [this](Variant const &value) {
         set_particle_property(&Particle::q, value);
       },
#else  // ESPRESSO_ELECTROSTATICS
       [](Variant const &value) {
         if (get_value<double>(value) != 0.) {
           throw std::runtime_error("Feature ELECTROSTATICS not compiled in");
         }
       },
#endif // ESPRESSO_ELECTROSTATICS
       [this]() { return get_particle_data(m_pid).q(); }},
#ifdef ESPRESSO_DIPOLES
      {"dip",
       [this](Variant const &value) {
         set_particle_property([&value](Particle &p) {
           auto const dip = get_value<Utils::Vector3d>(value);
           std::tie(p.quat(), p.dipm()) = convert_dip_to_quat(dip);
         });
       },
       [this]() { return get_particle_data(m_pid).calc_dip(); }},
      {"dipm",
       [this](Variant const &value) {
         set_particle_property(&Particle::dipm, value);
       },
       [this]() { return get_particle_data(m_pid).dipm(); }},
#endif // ESPRESSO_DIPOLES
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
      {"dip_fld",
       [this](Variant const &value) {
         set_particle_property(&Particle::dip_fld, value);
       },
       [this]() { return get_particle_data(m_pid).dip_fld(); }},
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
      {"magnetodynamics",
       [this](Variant const &value) {
         set_particle_property([&value](Particle &p) {
           auto const dict = get_value<VariantMap>(value);
           if (dict.contains("is_enabled"))
             p.stoner_wohlfarth_is_enabled() =
                 get_value<bool>(dict.at("is_enabled"));
           if (dict.contains("sw_phi_0"))
             p.stoner_wohlfarth_phi_0() =
                 get_value<double>(dict.at("sw_phi_0"));
           if (dict.contains("sat_mag"))
             p.saturation_magnetization() =
                 get_value<double>(dict.at("sat_mag"));
           if (dict.contains("anisotropy_field_inv"))
             p.magnetic_anisotropy_field_inv() =
                 get_value<double>(dict.at("anisotropy_field_inv"));
           if (dict.contains("anisotropy_energy"))
             p.magnetic_anisotropy_energy() =
                 get_value<double>(dict.at("anisotropy_energy"));
           if (dict.contains("sw_tau0_inv"))
             p.stoner_wohlfarth_tau0_inv() =
                 get_value<double>(dict.at("sw_tau0_inv"));
           if (dict.contains("sw_dt_incr"))
             p.stoner_wohlfarth_dt_incr() =
                 get_value<double>(dict.at("sw_dt_incr"));
         });
       },
       [this]() {
         auto const &p = get_particle_data(m_pid);
         return VariantMap{
             {"is_enabled", p.stoner_wohlfarth_is_enabled()},
             {"sw_phi_0", p.stoner_wohlfarth_phi_0()},
             {"sat_mag", p.saturation_magnetization()},
             {"anisotropy_field_inv", p.magnetic_anisotropy_field_inv()},
             {"anisotropy_energy", p.magnetic_anisotropy_energy()},
             {"sw_tau0_inv", p.stoner_wohlfarth_tau0_inv()},
             {"sw_dt_incr", p.stoner_wohlfarth_dt_incr()},
         };
       }},
#endif // ESPRESSO_THERMAL_STONER_WOHLFARTH
#ifdef ESPRESSO_ROTATION
      {"director",
       [this](Variant const &value) {
         set_particle_property([&value](Particle &p) {
           auto const director = get_value<Utils::Vector3d>(value).normalized();
           p.quat() = Utils::convert_director_to_quaternion(director);
         });
       },
       [this]() {
         auto const quat = get_particle_data(m_pid).quat();
         return Utils::convert_quaternion_to_director(quat);
       }},
      {"quat",
       [this](Variant const &value) {
         auto const quat = get_quaternion_safe("quat", value);
         set_particle_property([&quat](Particle &p) { p.quat() = quat; });
       },
       [this]() { return quat2vector(get_particle_data(m_pid).quat()); }},
      {"omega_body",
       [this](Variant const &value) {
         set_particle_property(&Particle::omega, value);
       },
       [this]() { return get_particle_data(m_pid).omega(); }},
      {"rotation",
       [this](Variant const &value) {
         set_particle_property([&value](Particle &p) {
           auto const rotation_flag =
               Utils::Vector3i{get_value<Utils::Vector3b>(value)};
           p.rotation() = bitfield_from_flag(rotation_flag);
         });
       },
       [this]() {
         auto const rotation_bits = get_particle_data(m_pid).rotation();
         return Utils::Vector3b{{::detail::get_nth_bit(rotation_bits, 0),
                                 ::detail::get_nth_bit(rotation_bits, 1),
                                 ::detail::get_nth_bit(rotation_bits, 2)}};
       }},
      {"omega_lab",
       [this](Variant const &value) {
         set_particle_property([&value](Particle &p) {
           auto const omega = get_value<Utils::Vector3d>(value);
           p.omega() = convert_vector_space_to_body(p, omega);
         });
       },
       [this]() {
         auto &p = get_particle_data(m_pid);
         return convert_vector_body_to_space(p, p.omega());
       }},
      {"torque_lab",
       [this](Variant const &value) {
         set_particle_property([&value](Particle &p) {
           auto const torque = get_value<Utils::Vector3d>(value);
           p.torque() = convert_vector_space_to_body(p, torque);
         });
       },
       [this]() {
         auto &p = get_particle_data(m_pid);
         return convert_vector_body_to_space(p, p.torque());
       }},
#endif // ESPRESSO_ROTATION
#ifdef ESPRESSO_ROTATIONAL_INERTIA
      {"rinertia",
       [this](Variant const &value) {
         set_particle_property(&Particle::rinertia, value);
       },
       [this]() { return get_particle_data(m_pid).rinertia(); }},
#endif // ESPRESSO_ROTATIONAL_INERTIA
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
      {"mu_E",
       [this](Variant const &value) {
         set_particle_property(&Particle::mu_E, value);
       },
       [this]() { return get_particle_data(m_pid).mu_E(); }},
#endif // ESPRESSO_LB_ELECTROHYDRODYNAMICS
#ifdef ESPRESSO_EXTERNAL_FORCES
      {"fix",
       [this](Variant const &value) {
         set_particle_property([&value](Particle &p) {
           auto const fix_flag =
               Utils::Vector3i{(get_value<Utils::Vector3b>(value))};
           p.fixed() = bitfield_from_flag(fix_flag);
         });
       },
       [this]() {
         auto const fixed = get_particle_data(m_pid).fixed();
         return Utils::Vector3b{{::detail::get_nth_bit(fixed, 0),
                                 ::detail::get_nth_bit(fixed, 1),
                                 ::detail::get_nth_bit(fixed, 2)}};
       }},
      {"ext_force",
       [this](Variant const &value) {
         set_particle_property(&Particle::ext_force, value);
       },
       [this]() { return get_particle_data(m_pid).ext_force(); }},
#ifdef ESPRESSO_ROTATION
      {"ext_torque",
       [this](Variant const &value) {
         set_particle_property(&Particle::ext_torque, value);
       },
       [this]() { return get_particle_data(m_pid).ext_torque(); }},
#endif // ESPRESSO_ROTATION
#endif // ESPRESSO_EXTERNAL_FORCES
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
      {"gamma",
       [this](Variant const &value) {
         set_particle_property(&Particle::gamma,
                               Variant{get_gamma_safe(value)});
       },
       [this]() { return get_particle_data(m_pid).gamma(); }},
#ifdef ESPRESSO_ROTATION
      {"gamma_rot",
       [this](Variant const &value) {
         set_particle_property(&Particle::gamma_rot,
                               Variant{get_gamma_safe(value)});
       },
       [this]() { return get_particle_data(m_pid).gamma_rot(); }},
#endif // ESPRESSO_ROTATION
#endif // ESPRESSO_THERMOSTAT_PER_PARTICLE
      {"pos_folded", AutoParameter::read_only,
       [this]() {
         auto const &box_geo = *get_system()->box_geo;
         return box_geo.folded_position(get_particle_data(m_pid).pos());
       }},
      {"lees_edwards_offset",
       [this](Variant const &value) {
         set_particle_property(&Particle::lees_edwards_offset, value);
       },
       [this]() { return get_particle_data(m_pid).lees_edwards_offset(); }},
      {"lees_edwards_flag", AutoParameter::read_only,
       [this]() { return get_particle_data(m_pid).lees_edwards_flag(); }},
      {"image_box", AutoParameter::read_only,
       [this]() {
         auto const &box_geo = *get_system()->box_geo;
         auto const p = get_particle_data(m_pid);
         return box_geo.folded_image_box(p.pos(), p.image_box());
       }},
      {"node", AutoParameter::read_only,
       [this]() {
         return (context()->is_head_node()) ? get_particle_node(m_pid) : -1;
       }},
      {"mol_id",
       [this](Variant const &value) {
         auto const mol_id = get_value<int>(value);
         if (mol_id < 0) {
           throw std::domain_error(
               error_msg("mol_id", "must be an integer >= 0"));
         }
         set_particle_property(&Particle::mol_id, Variant{mol_id});
       },
       [this]() { return get_particle_data(m_pid).mol_id(); }},
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
      {"vs_quat",
       [this](Variant const &value) {
         auto const quat = get_quaternion_safe("vs_quat", value);
         set_particle_property(
             [&quat](Particle &p) { p.vs_relative().quat = quat; });
       },
       [this]() {
         return quat2vector(get_particle_data(m_pid).vs_relative().quat);
       }},
      {"vs_relative",
       [this](Variant const &value) {
         ParticleProperties::VirtualSitesRelativeParameters vs_relative{};
         try {
           auto const array = get_value<std::vector<Variant>>(value);
           if (array.size() != 3) {
             throw 0;
           }
           vs_relative.distance = get_value<double>(array[1]);
           vs_relative.to_particle_id = get_value<int>(array[0]);
           vs_relative.rel_orientation =
               get_quaternion_safe("vs_relative", array[2]);
         } catch (...) {
           throw std::invalid_argument(error_msg(
               "vs_relative", "must take the form [id, distance, quaternion]"));
         }
         set_particle_property(
             [&vs_relative](Particle &p) { p.vs_relative() = vs_relative; });
       },
       [this]() {
         auto const vs_rel = get_particle_data(m_pid).vs_relative();
         return std::vector<Variant>{{vs_rel.to_particle_id, vs_rel.distance,
                                      quat2vector(vs_rel.rel_orientation)}};
       }},
#endif // ESPRESSO_VIRTUAL_SITES_RELATIVE
      {"propagation",
       [this](Variant const &value) {
         auto const propagation = get_value<int>(value);
         if (!is_valid_propagation_combination(propagation)) {
           throw std::domain_error(error_msg(
               "propagation", "propagation combination not accepted: " +
                                  propagation_bitmask_to_string(propagation)));
         }
         set_particle_property(&Particle::propagation, value);
       },
       [this]() { return get_particle_data(m_pid).propagation(); }},
#ifdef ESPRESSO_ENGINE
      {"swimming",
       [this](Variant const &value) {
         set_particle_property([&value](Particle &p) {
           ParticleParametersSwimming swim{};
           swim.swimming = true;
           auto const dict = get_value<VariantMap>(value);
           if (dict.contains("f_swim")) {
             swim.f_swim = get_value<double>(dict.at("f_swim"));
           }
           if (dict.contains("is_engine_force_on_fluid")) {
             auto const is_engine_force_on_fluid =
                 get_value<bool>(dict.at("is_engine_force_on_fluid"));
             swim.is_engine_force_on_fluid = is_engine_force_on_fluid;
           }
           p.swimming() = swim;
         });
       },
       [this]() {
         auto const swim = get_particle_data(m_pid).swimming();
         return VariantMap{
             {"f_swim", swim.f_swim},
             {"is_engine_force_on_fluid", swim.is_engine_force_on_fluid},
         };
       }},
#endif // ESPRESSO_ENGINE
  });
}

Variant ParticleHandle::do_call_method(std::string const &name,
                                       VariantMap const &params) {
  if (name == "set_param_parallel") {
    auto const param_name = get_value<std::string>(params, "name");
    if (not params.contains("value")) {
      throw Exception("Parameter '" + param_name + "' is missing.");
    }
    auto const &value = params.at("value");
    context()->parallel_try_catch(
        [&]() { do_set_parameter(param_name, value); });
    return {};
  }
  if (name == "update_params") {
    // Set new properties
    context()->parallel_try_catch([&]() {
#ifdef ESPRESSO_ROTATION
      sanity_checks_rotation(params);
#endif
      for (auto const &name : get_parameter_insertion_order()) {
        if (params.contains(name) and name != "bonds") {
          do_set_parameter(name, params.at(name));
        }
      }
    });

    // Set bonds
    if (params.contains("bonds_ids")) {
      // Remove old bonds
      set_particle_property([&](Particle &p) { p.bonds().clear(); });
      // Add new bonds
      auto const bonds_ids = get_value<std::vector<int>>(params, "bonds_ids");
      auto const bonds_partner_ids =
          get_value<std::vector<std::vector<int>>>(params, "bonds_parts");
      for (std::size_t i = 0; i < bonds_ids.size(); i += 1) {
        std::vector<int> particle_ids = {m_pid};
        std::ranges::copy(bonds_partner_ids[i],
                          std::back_inserter(particle_ids));
        ::add_bond(*get_system(), bonds_ids[i], particle_ids);
        get_system()->on_particle_change();
      }
    }
#ifdef ESPRESSO_EXCLUSIONS
    // set exclusions
    if (params.contains("exclusions")) {
      std::vector<int> exclusion_list;
      if (is_type<int>(params.at("exclusions"))) {
        exclusion_list.emplace_back(get_value<int>(params, "exclusions"));
      } else {
        exclusion_list = get_value<std::vector<int>>(params, "exclusions");
      }
      context()->parallel_try_catch([&]() {
        auto &cell_structure = get_cell_structure()->get_cell_structure();
        for (auto const pid : exclusion_list) {
          particle_exclusion_sanity_checks(m_pid, pid, cell_structure,
                                           context()->get_comm());
        }
      });
      set_particle_property([this, &exclusion_list](Particle &p) {
        auto &cell_structure = get_cell_structure()->get_cell_structure();
        for (auto const pid : p.exclusions()) {
          local_remove_exclusion(m_pid, pid, cell_structure);
        }
        for (auto const pid : exclusion_list) {
          if (!p.has_exclusion(pid)) {
            local_add_exclusion(m_pid, pid, cell_structure);
          }
        }
      });
    }
#endif // ESPRESSO_EXCLUSIONS
  }
  if (name == "get_bond_by_id") {
    if (not context()->is_head_node()) {
      return {};
    }
    return get_bonded_ias()->call_method("get_bond", params);
  }
  if (name == "get_bonds_view") {
    if (not context()->is_head_node()) {
      return {};
    }
    auto const bond_list = get_particle_data(m_pid).bonds();
    std::vector<std::vector<Variant>> bonds_flat;
    for (auto const &&bond_view : bond_list) {
      std::vector<Variant> bond_flat;
      bond_flat.emplace_back(bond_view.bond_id());
      for (auto const pid : bond_view.partner_ids()) {
        bond_flat.emplace_back(pid);
      }
      bonds_flat.emplace_back(std::move(bond_flat));
    }
    return make_vector_of_variants(bonds_flat);
  }
  if (name == "add_bond") {
    auto const bond_id = get_value<int>(params, "bond_id");
    auto const partner_ids = get_value<std::vector<int>>(params, "part_id");
    std::vector<int> particle_ids = {m_pid};
    std::ranges::copy(partner_ids, std::back_inserter(particle_ids));
    ::add_bond(*get_system(), bond_id, particle_ids);
    get_system()->on_particle_change();
  } else if (name == "del_bond") {
    set_particle_property([&params](Particle &p) {
      auto const bond_id = get_value<int>(params, "bond_id");
      auto const part_id = get_value<std::vector<int>>(params, "part_id");
      auto const bond_view =
          BondView(bond_id, {part_id.data(), part_id.size()});
      auto &bond_list = p.bonds();
      auto it = std::find(bond_list.begin(), bond_list.end(), bond_view);
      if (it != bond_list.end()) {
        bond_list.erase(it);
      }
    });
  } else if (name == "delete_all_bonds") {
    set_particle_property([&](Particle &p) { p.bonds().clear(); });
  } else if (name == "is_valid_bond_id") {
    auto const bond_id = get_value<int>(params, "bond_id");
    return get_system()->bonded_ias->get_zero_based_type(bond_id) != 0;
  }
  if (name == "remove_particle") {
    context()->parallel_try_catch([&]() {
      auto &cell_structure = get_cell_structure()->get_cell_structure();
      std::ignore =
          get_real_particle(context()->get_comm(), m_pid, cell_structure);
      remove_particle(m_pid);
    });
  } else if (name == "is_virtual") {
    if (not context()->is_head_node()) {
      return {};
    }
    return get_particle_data(m_pid).is_virtual();
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  } else if (name == "vs_auto_relate_to") {
    if (not context()->is_head_node()) {
      return {};
    }
    auto const other_pid = get_value<int>(params, "pid");
    auto const override_cutoff_check =
        get_value<bool>(params, "override_cutoff_check");
    if (m_pid == other_pid) {
      throw std::invalid_argument("A virtual site cannot relate to itself");
    }
    if (other_pid < 0) {
      throw std::domain_error("Invalid particle id: " +
                              std::to_string(other_pid));
    }
    auto const system = get_system();
    /* note: this code can be rewritten as parallel code, but only with a call
     * to `cells_update_ghosts(DATA_PART_POSITION | DATA_PART_PROPERTIES)`, as
     * there is no guarantee the virtual site has visibility of the relative
     * particle through the ghost layer during particle creation. However,
     * ghost updates can scramble the particle ordering in the local cells,
     * which is an issue for checkpointing: the H5MD writer will use the
     * scrambled ordering before writing to a checkpoint file and the
     * non-scrambled ordering after reloading from a checkpoint file.
     */
    auto const &p_current = get_particle_data(m_pid);
    auto const &p_relate_to = get_particle_data(other_pid);
    auto const [quat, dist] = calculate_vs_relate_to_params(
        p_current, p_relate_to, *system->box_geo, system->get_min_global_cut(),
        override_cutoff_check);
    set_parameter("vs_relative", Variant{std::vector<Variant>{
                                     {other_pid, dist, quat2vector(quat)}}});
    set_parameter("propagation",
                  Variant{static_cast<int>(PropagationMode::TRANS_VS_RELATIVE |
                                           PropagationMode::ROT_VS_RELATIVE)});
#endif // ESPRESSO_VIRTUAL_SITES_RELATIVE
#ifdef ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS
  } else if (name == "vs_com_relate_to") {
    auto &cell_structure = get_cell_structure()->get_cell_structure();
    auto const molid = get_value<int>(params, "molid");
    auto const maybe_exists_vs = get_pid_for_vs_com(cell_structure, molid);
    if (not context()->is_head_node()) {
      return {};
    }
    if (molid < 0) {
      throw std::domain_error("Invalid molecule id: " + std::to_string(molid));
    }
    if (maybe_exists_vs) {
      throw std::runtime_error(
          "Molecule id: " + std::to_string(molid) +
          " is already tracked by virtual site with particle id: " +
          std::to_string(*maybe_exists_vs));
    }
    set_parameter("mol_id", params.at("molid"));
    set_parameter(
        "propagation",
        Variant{static_cast<int>(PropagationMode::TRANS_VS_CENTER_OF_MASS)});
#endif // ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS
#ifdef ESPRESSO_EXCLUSIONS
  } else if (name == "has_exclusion") {
    auto const other_pid = get_value<int>(params, "pid");
    auto &cell_structure = get_cell_structure()->get_cell_structure();
    auto const p =
        get_real_particle(context()->get_comm(), m_pid, cell_structure);
    if (p != nullptr) {
      return p->has_exclusion(other_pid);
    }
  }
  if (name == "add_exclusion") {
    auto const other_pid = get_value<int>(params, "pid");
    auto &cell_structure = get_cell_structure()->get_cell_structure();
    context()->parallel_try_catch([&]() {
      particle_exclusion_sanity_checks(m_pid, other_pid, cell_structure,
                                       context()->get_comm());
    });
    local_add_exclusion(m_pid, other_pid, cell_structure);
    get_system()->on_particle_change();
  } else if (name == "del_exclusion") {
    auto const other_pid = get_value<int>(params, "pid");
    auto &cell_structure = get_cell_structure()->get_cell_structure();
    context()->parallel_try_catch([&]() {
      particle_exclusion_sanity_checks(m_pid, other_pid, cell_structure,
                                       context()->get_comm());
    });
    local_remove_exclusion(m_pid, other_pid, cell_structure);
    get_system()->on_particle_change();
  } else if (name == "set_exclusions") {
    auto &cell_structure = get_cell_structure()->get_cell_structure();
    std::vector<int> exclusion_list;
    try {
      auto const pid = get_value<int>(params, "p_ids");
      exclusion_list.push_back(pid);
    } catch (...) {
      exclusion_list = get_value<std::vector<int>>(params, "p_ids");
    }
    context()->parallel_try_catch([&]() {
      for (auto const pid : exclusion_list) {
        particle_exclusion_sanity_checks(m_pid, pid, cell_structure,
                                         context()->get_comm());
      }
    });
    set_particle_property([this, &exclusion_list](Particle &p) {
      auto &cell_structure = get_cell_structure()->get_cell_structure();
      for (auto const pid : p.exclusions()) {
        local_remove_exclusion(m_pid, pid, cell_structure);
      }
      for (auto const pid : exclusion_list) {
        if (!p.has_exclusion(pid)) {
          local_add_exclusion(m_pid, pid, cell_structure);
        }
      }
    });
  } else if (name == "get_exclusions") {
    if (not context()->is_head_node()) {
      return {};
    }
    auto const excl_list = get_particle_data(m_pid).exclusions();
    return Variant{std::vector<int>{excl_list.begin(), excl_list.end()}};
#endif // ESPRESSO_EXCLUSIONS
#ifdef ESPRESSO_ROTATION
  }
  if (name == "rotate_particle") {
    set_particle_property([&params](Particle &p) {
      auto const axis = get_value<Utils::Vector3d>(params, "axis");
      auto const angle = get_value<double>(params, "angle");
      local_rotate_particle(p, axis, angle);
    });
  }
  if (name == "convert_vector_body_to_space") {
    return get_particle_property<std::vector<double>>(
        [&params](Particle const &p) {
          auto const vec = get_value<Utils::Vector3d>(params, "vec");
          return convert_vector_body_to_space(p, vec).as_vector();
        });
  }
  if (name == "convert_vector_space_to_body") {
    return get_particle_property<std::vector<double>>(
        [&params](Particle const &p) {
          auto const vec = get_value<Utils::Vector3d>(params, "vec");
          return convert_vector_space_to_body(p, vec).as_vector();
        });
#endif // ESPRESSO_ROTATION
  }
  return {};
}

std::size_t ParticleHandle::setup_hidden_args(VariantMap const &params) {
  auto n_extra_args = params.size() - params.count("id");
  if (params.contains("__cell_structure")) {
    auto so = get_value<std::shared_ptr<CellSystem::CellSystem>>(
        params, "__cell_structure");
    so->configure(*this);
    m_cell_structure = so;
    --n_extra_args;
  }
  if (params.contains("__bonded_ias")) {
    m_bonded_ias = get_value<std::shared_ptr<Interactions::BondedInteractions>>(
        params, "__bonded_ias");
    --n_extra_args;
  }
  return n_extra_args;
}

void ParticleHandle::do_construct(VariantMap const &params) {
  auto const n_extra_args = setup_hidden_args(params);
  m_pid = (params.contains("id")) ? get_value<int>(params, "id")
                                  : get_maximal_particle_id() + 1;

#ifndef NDEBUG
  if (not params.contains("id")) {
    auto head_node_reference = m_pid;
    boost::mpi::broadcast(context()->get_comm(), head_node_reference, 0);
    assert(m_pid == head_node_reference && "global max_seen_pid has diverged");
  }
#endif

  // create a new particle if extra arguments were passed
  if (n_extra_args == 0) {
    return;
  }

  auto const pos = get_value<Utils::Vector3d>(params, "pos");
  context()->parallel_try_catch([&]() {
    particle_checks(m_pid, pos);
    auto &cell_structure = get_cell_structure()->get_cell_structure();
    auto ptr = cell_structure.get_local_particle(m_pid);
    if (ptr != nullptr) {
      throw std::invalid_argument("Particle " + std::to_string(m_pid) +
                                  " already exists");
    }
  });

#ifdef ESPRESSO_ROTATION
  context()->parallel_try_catch([&]() { sanity_checks_rotation(params); });
#endif // ESPRESSO_ROTATION

  // create a default-constructed particle
  make_new_particle(m_pid, pos);

  try {
    context()->parallel_try_catch([&]() {
      /* clang-format off */
      // set particle properties (filter out read-only and deferred properties)
      std::set<std::string_view> const skip = {
          "pos_folded", "pos", "id", "exclusions", "node", "image_box", "bonds",
          "lees_edwards_flag", "__cpt_sentinel",
      };
      /* clang-format on */
      for (auto const &name : get_parameter_insertion_order()) {
        if (params.contains(name) and not skip.contains(name)) {
          do_set_parameter(name, params.at(name));
        }
      }
      if (not params.contains("type")) {
        do_set_parameter("type", 0);
      }
    });
  } catch (...) {
    remove_particle(m_pid);
    throw;
  }
}

} // namespace Particles
} // namespace ScriptInterface
