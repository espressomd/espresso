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

#include "config/config.hpp"

#include "ParticleHandle.hpp"

#include "script_interface/Variant.hpp"
#include "script_interface/communication.hpp"
#include "script_interface/get_value.hpp"

#include "core/bonded_interactions/bonded_interaction_data.hpp"
#include "core/cells.hpp"
#include "core/event.hpp"
#include "core/exclusions.hpp"
#include "core/grid.hpp"
#include "core/nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "core/particle_node.hpp"
#include "core/rotation.hpp"
#include "core/virtual_sites.hpp"

#include <utils/Vector.hpp>

#include <boost/format.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/communicator.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

namespace ScriptInterface {
namespace Particles {

static void particle_checks(int p_id, Utils::Vector3d const &pos) {
  if (p_id < 0) {
    throw std::domain_error("Invalid particle id: " + std::to_string(p_id));
  }
#ifndef __FAST_MATH__
  if (std::isnan(pos[0]) or std::isnan(pos[1]) or std::isnan(pos[2]) or
      std::isinf(pos[0]) or std::isinf(pos[1]) or std::isinf(pos[2])) {
    throw std::domain_error("Particle position must be finite");
  }
#endif // __FAST_MATH__
}

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

static auto error_msg(std::string const &name, std::string const &reason) {
  std::stringstream msg;
  msg << "attribute '" << name << "' of 'ParticleHandle' " << reason;
  return msg.str();
}

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

#ifdef THERMOSTAT_PER_PARTICLE
static auto get_gamma_safe(Variant const &value) {
#ifdef PARTICLE_ANISOTROPY
  try {
    return Utils::Vector3d::broadcast(get_value<double>(value));
  } catch (...) {
    return get_value<Utils::Vector3d>(value);
  }
#else  // PARTICLE_ANISOTROPY
  return get_value<double>(value);
#endif // PARTICLE_ANISOTROPY
}
#endif // THERMOSTAT_PER_PARTICLE

static auto get_real_particle(boost::mpi::communicator const &comm, int p_id) {
  if (p_id < 0) {
    throw std::domain_error("Invalid particle id: " + std::to_string(p_id));
  }
  auto ptr = ::cell_structure.get_local_particle(p_id);
  if (ptr != nullptr and ptr->is_ghost()) {
    ptr = nullptr;
  }
  auto const n_found = boost::mpi::all_reduce(
      comm, static_cast<int>(ptr != nullptr), std::plus<>());
  if (n_found == 0) {
    throw std::runtime_error("Particle with id " + std::to_string(p_id) +
                             " not found");
  }
  return ptr;
}

template <typename T, class F>
T ParticleHandle::get_particle_property(F const &fun) const {
  auto const &comm = context()->get_comm();
  auto const ptr = const_cast<Particle const *>(get_real_particle(comm, m_pid));
  boost::optional<T> ret;
  if (ptr == nullptr) {
    ret = {};
  } else {
    ret = {fun(*ptr)};
  }
  return mpi_reduce_optional(comm, ret);
}

template <typename T>
T ParticleHandle::get_particle_property(T const &(Particle::*getter)()
                                            const) const {
  return get_particle_property<T>(
      [getter](Particle const &p) { return (p.*getter)(); });
}

template <class F>
void ParticleHandle::set_particle_property(F const &fun) const {
  auto const &comm = context()->get_comm();
  auto const ptr = get_real_particle(comm, m_pid);
  if (ptr != nullptr) {
    fun(*ptr);
  }
  on_particle_change();
}

template <typename T>
void ParticleHandle::set_particle_property(T &(Particle::*setter)(),
                                           Variant const &value) const {
  set_particle_property(
      [&value, setter](Particle &p) { (p.*setter)() = get_value<T>(value); });
}

ParticleHandle::ParticleHandle() {
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
         make_particle_type_exist(new_type);
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
         return unfolded_position(pos, image_box, ::box_geo.length());
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
#ifdef MASS
       [this](Variant const &value) {
         if (get_value<double>(value) <= 0.) {
           throw std::domain_error(error_msg("mass", "must be a float > 0"));
         }
         set_particle_property(&Particle::mass, value);
       },
#else  // MASS
       [](Variant const &value) {
         auto const default_mass = Particle().mass();
         if (std::abs(get_value<double>(value) - default_mass) > 1e-10) {
           throw std::runtime_error("Feature MASS not compiled in");
         }
       },
#endif // MASS
       [this]() { return get_particle_data(m_pid).mass(); }},
      {"q",
#ifdef ELECTROSTATICS
       [this](Variant const &value) {
         set_particle_property(&Particle::q, value);
       },
#else  // ELECTROSTATICS
       [](Variant const &value) {
         if (get_value<double>(value) != 0.) {
           throw std::runtime_error("Feature ELECTROSTATICS not compiled in");
         }
       },
#endif // ELECTROSTATICS
       [this]() { return get_particle_data(m_pid).q(); }},
      {"virtual",
#ifdef VIRTUAL_SITES
       [this](Variant const &value) {
         set_particle_property(&Particle::virtual_flag, value);
       },
#else  // VIRTUAL_SITES
       [](Variant const &value) {
         if (get_value<bool>(value)) {
           throw std::runtime_error("Feature VIRTUAL_SITES not compiled in");
         }
       },
#endif // VIRTUAL_SITES
       [this]() { return get_particle_data(m_pid).is_virtual(); }},
#ifdef ROTATION
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
               Utils::Vector3i{(get_value<Utils::Vector3b>(value))};
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
#endif // ROTATION
#ifdef DIPOLES
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
#endif // DIPOLES
#ifdef DIPOLE_FIELD_TRACKING
      {"dip_fld",
       [this](Variant const &value) {
         set_particle_property(&Particle::dip_fld, value);
       },
       [this]() { return get_particle_data(m_pid).dip_fld(); }},
#endif
#ifdef ROTATIONAL_INERTIA
      {"rinertia",
       [this](Variant const &value) {
         set_particle_property(&Particle::rinertia, value);
       },
       [this]() { return get_particle_data(m_pid).rinertia(); }},
#endif // ROTATIONAL_INERTIA
#ifdef LB_ELECTROHYDRODYNAMICS
      {"mu_E",
       [this](Variant const &value) {
         set_particle_property(&Particle::mu_E, value);
       },
       [this]() { return get_particle_data(m_pid).mu_E(); }},
#endif // LB_ELECTROHYDRODYNAMICS
#ifdef EXTERNAL_FORCES
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
#ifdef ROTATION
      {"ext_torque",
       [this](Variant const &value) {
         set_particle_property(&Particle::ext_torque, value);
       },
       [this]() { return get_particle_data(m_pid).ext_torque(); }},
#endif // ROTATION
#endif // EXTERNAL_FORCES
#ifdef THERMOSTAT_PER_PARTICLE
      {"gamma",
       [this](Variant const &value) {
         set_particle_property(&Particle::gamma,
                               Variant{get_gamma_safe(value)});
       },
       [this]() { return get_particle_data(m_pid).gamma(); }},
#ifdef ROTATION
      {"gamma_rot",
       [this](Variant const &value) {
         set_particle_property(&Particle::gamma_rot,
                               Variant{get_gamma_safe(value)});
       },
       [this]() { return get_particle_data(m_pid).gamma_rot(); }},
#endif // ROTATION
#endif // THERMOSTAT_PER_PARTICLE
      {"pos_folded", AutoParameter::read_only,
       [this]() {
         return folded_position(get_particle_data(m_pid).pos(), ::box_geo);
       }},

      {"lees_edwards_offset",
       [this](Variant const &value) {
         set_particle_property(&Particle::lees_edwards_offset, value);
       },
       [this]() { return get_particle_data(m_pid).lees_edwards_offset(); }},
      {"lees_edwards_flag", AutoParameter::read_only,
       [this]() { return get_particle_data(m_pid).lees_edwards_flag(); }},
      {"image_box", AutoParameter::read_only,
       [this]() { return get_particle_data(m_pid).image_box(); }},
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
#ifdef VIRTUAL_SITES_RELATIVE
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
#endif // VIRTUAL_SITES_RELATIVE
#ifdef ENGINE
      {"swimming",
       [this](Variant const &value) {
         set_particle_property([&value](Particle &p) {
           ParticleParametersSwimming swim{};
           swim.swimming = true;
           auto const dict = get_value<VariantMap>(value);
           if (dict.count("f_swim") != 0) {
             swim.f_swim = get_value<double>(dict.at("f_swim"));
           }
           if (dict.count("is_engine_force_on_fluid") != 0) {
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
#endif // ENGINE
  });
}

#ifdef EXCLUSIONS
/**
 * @brief Locally add an exclusion to a particle.
 * @param pid1 the identity of the first exclusion partner
 * @param pid2 the identity of the second exclusion partner
 */
static void local_add_exclusion(int pid1, int pid2) {
  if (auto p1 = ::cell_structure.get_local_particle(pid1)) {
    add_exclusion(*p1, pid2);
  }
  if (auto p2 = ::cell_structure.get_local_particle(pid2)) {
    add_exclusion(*p2, pid1);
  }
}

/**
 * @brief Locally remove an exclusion to a particle.
 * @param pid1 the identity of the first exclusion partner
 * @param pid2 the identity of the second exclusion partner
 */
static void local_remove_exclusion(int pid1, int pid2) {
  if (auto p1 = ::cell_structure.get_local_particle(pid1)) {
    delete_exclusion(*p1, pid2);
  }
  if (auto p2 = ::cell_structure.get_local_particle(pid2)) {
    delete_exclusion(*p2, pid1);
  }
}

void ParticleHandle::particle_exclusion_sanity_checks(int pid1,
                                                      int pid2) const {
  if (pid1 == pid2) {
    throw std::runtime_error("Particles cannot exclude themselves (id " +
                             std::to_string(pid1) + ")");
  }
  std::ignore = get_real_particle(context()->get_comm(), pid1);
  std::ignore = get_real_particle(context()->get_comm(), pid2);
}
#endif // EXCLUSIONS

Variant ParticleHandle::do_call_method(std::string const &name,
                                       VariantMap const &params) {
  if (name == "set_param_parallel") {
    auto const param_name = get_value<std::string>(params, "name");
    if (params.count("value") == 0) {
      throw Exception("Parameter '" + param_name + "' is missing.");
    }
    auto const &value = params.at("value");
    context()->parallel_try_catch(
        [&]() { do_set_parameter(param_name, value); });
    return {};
  }
  if (name == "get_bonds_view") {
    if (not context()->is_head_node()) {
      return {};
    }
    auto const bond_list = get_particle_data(m_pid).bonds();
    std::vector<std::vector<int>> bonds_flat;
    for (auto const &&bond_view : bond_list) {
      std::vector<int> bond_flat;
      bond_flat.emplace_back(bond_view.bond_id());
      for (auto const pid : bond_view.partner_ids()) {
        bond_flat.emplace_back(pid);
      }
      bonds_flat.emplace_back(std::move(bond_flat));
    }
    return make_vector_of_variants(bonds_flat);
  }
  if (name == "add_bond") {
    set_particle_property([&params](Particle &p) {
      auto const bond_id = get_value<int>(params, "bond_id");
      auto const part_id = get_value<std::vector<int>>(params, "part_id");
      auto const bond_view =
          BondView(bond_id, {part_id.data(), part_id.size()});
      p.bonds().insert(bond_view);
    });
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
    return ::bonded_ia_params.get_zero_based_type(bond_id) != 0;
  }
  if (name == "remove_particle") {
    context()->parallel_try_catch([&]() {
      std::ignore = get_real_particle(context()->get_comm(), m_pid);
      remove_particle(m_pid);
    });
#ifdef VIRTUAL_SITES_RELATIVE
  } else if (name == "vs_relate_to") {
    if (not context()->is_head_node()) {
      return {};
    }
    auto const other_pid = get_value<int>(params, "pid");
    if (m_pid == other_pid) {
      throw std::invalid_argument("A virtual site cannot relate to itself");
    }
    if (other_pid < 0) {
      throw std::domain_error("Invalid particle id: " +
                              std::to_string(other_pid));
    }
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
    auto const [quat, dist] =
        calculate_vs_relate_to_params(p_current, p_relate_to);
    set_parameter("vs_relative", Variant{std::vector<Variant>{
                                     {other_pid, dist, quat2vector(quat)}}});
    set_parameter("virtual", true);
#endif // VIRTUAL_SITES_RELATIVE
#ifdef EXCLUSIONS
  } else if (name == "has_exclusion") {
    auto const other_pid = get_value<int>(params, "pid");
    auto const p = get_real_particle(context()->get_comm(), m_pid);
    if (p != nullptr) {
      return p->has_exclusion(other_pid);
    }
  }
  if (name == "add_exclusion") {
    auto const other_pid = get_value<int>(params, "pid");
    context()->parallel_try_catch(
        [&]() { particle_exclusion_sanity_checks(m_pid, other_pid); });
    local_add_exclusion(m_pid, other_pid);
    on_particle_change();
  } else if (name == "del_exclusion") {
    auto const other_pid = get_value<int>(params, "pid");
    context()->parallel_try_catch(
        [&]() { particle_exclusion_sanity_checks(m_pid, other_pid); });
    local_remove_exclusion(m_pid, other_pid);
    on_particle_change();
  } else if (name == "set_exclusions") {
    std::vector<int> exclusion_list;
    try {
      auto const pid = get_value<int>(params, "p_ids");
      exclusion_list.push_back(pid);
    } catch (...) {
      exclusion_list = get_value<std::vector<int>>(params, "p_ids");
    }
    context()->parallel_try_catch([&]() {
      for (auto const pid : exclusion_list) {
        particle_exclusion_sanity_checks(m_pid, pid);
      }
    });
    set_particle_property([this, &exclusion_list](Particle &p) {
      for (auto const pid : p.exclusions()) {
        local_remove_exclusion(m_pid, pid);
      }
      for (auto const pid : exclusion_list) {
        if (!p.has_exclusion(pid)) {
          local_add_exclusion(m_pid, pid);
        }
      }
    });
  } else if (name == "get_exclusions") {
    if (not context()->is_head_node()) {
      return {};
    }
    auto const excl_list = get_particle_data(m_pid).exclusions();
    return Variant{std::vector<int>{excl_list.begin(), excl_list.end()}};
#endif // EXCLUSIONS
#ifdef ROTATION
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
#endif // ROTATION
  }
  return {};
}

#ifdef ROTATION
static auto const contradicting_arguments_quat = std::vector<
    std::array<std::string, 3>>{{
    {{"dip", "dipm",
      "Setting 'dip' is sufficient as it defines the scalar dipole moment."}},
    {{"quat", "director",
      "Setting 'quat' is sufficient as it defines the director."}},
    {{"dip", "quat",
      "Setting 'dip' would overwrite 'quat'. Set 'quat' and 'dipm' instead."}},
    {{"dip", "director",
      "Setting 'dip' would overwrite 'director'. Set 'director' and "
      "'dipm' instead."}},
}};
#endif // ROTATION

void ParticleHandle::do_construct(VariantMap const &params) {
  auto const n_extra_args = params.size() - params.count("id");
  auto const has_param = [&params](std::string const key) {
    return params.count(key) == 1;
  };
  m_pid = (has_param("id")) ? get_value<int>(params, "id")
                            : get_maximal_particle_id() + 1;

#ifndef NDEBUG
  if (!has_param("id")) {
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
    auto ptr = ::cell_structure.get_local_particle(m_pid);
    if (ptr != nullptr) {
      throw std::invalid_argument("Particle " + std::to_string(m_pid) +
                                  " already exists");
    }
  });

#ifdef ROTATION
  context()->parallel_try_catch([&]() {
    // if we are not constructing a particle from a checkpoint file,
    // check the quaternion is not accidentally set twice by the user
    if (not has_param("__cpt_sentinel")) {
      auto formatter =
          boost::format("Contradicting particle attributes: '%s' and '%s'. %s");
      for (auto const &[prop1, prop2, reason] : contradicting_arguments_quat) {
        if (has_param(prop1) and has_param(prop2)) {
          auto const err_msg = boost::str(formatter % prop1 % prop2 % reason);
          throw std::invalid_argument(err_msg);
        }
      }
    }
  });
#endif // ROTATION

  // create a default-constructed particle
  make_new_particle(m_pid, pos);

  context()->parallel_try_catch([&]() {
    // set particle properties (filter out read-only and deferred properties)
    std::set<std::string> const skip = {
        "pos_folded", "pos", "quat", "director",  "id",    "lees_edwards_flag",
        "exclusions", "dip", "node", "image_box", "bonds", "__cpt_sentinel",
    };
#ifdef ROTATION
    // multiple parameters can potentially set the quaternion, but only one
    // can be allowed to; these conditionals are required to handle a reload
    // from a checkpoint file, where all properties exist (avoids accidentally
    // overwriting the quaternion by the default-constructed dipole moment)
    for (std::string name : {"quat", "director", "dip"}) {
      if (has_param(name)) {
        do_set_parameter(name, params.at(name));
        break;
      }
    }
#endif // ROTATION
    std::vector<std::string> sorted_param_names = {};
    std::for_each(params.begin(), params.end(), [&](auto const &kv) {
      if (skip.count(kv.first) == 0) {
        sorted_param_names.push_back(kv.first);
      }
    });
    std::sort(sorted_param_names.begin(), sorted_param_names.end());
    for (auto const &name : sorted_param_names) {
      do_set_parameter(name, params.at(name));
    }
    if (not has_param("type")) {
      do_set_parameter("type", 0);
    }
  });
}

} // namespace Particles
} // namespace ScriptInterface
