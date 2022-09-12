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

#include "config.hpp"

#include "ParticleHandle.hpp"

#include "script_interface/get_value.hpp"

#include "core/bonded_interactions/bonded_interaction_data.hpp"
#include "core/grid.hpp"
#include "core/particle_data.hpp"
#include "core/particle_node.hpp"
#include "core/rotation.hpp"
#include "core/virtual_sites.hpp"

#include <utils/Vector.hpp>

#include <algorithm>
#include <cmath>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Particles {

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

static auto get_bond_vector(VariantMap const &params) {
  auto const bond_id = get_value<int>(params, "bond_id");
  auto const part_id = get_value<std::vector<int>>(params, "part_id");
  std::vector<int> bond_view;
  bond_view.emplace_back(bond_id);
  for (auto const pid : part_id) {
    bond_view.emplace_back(pid);
  }
  return bond_view;
}

ParticleHandle::ParticleHandle() {
  add_parameters({
      {"id", AutoParameter::read_only, [this]() { return m_pid; }},
      {"type",
       [this](Variant const &value) {
         auto const p_type = get_value<int>(value);
         if (p_type < 0) {
           throw std::domain_error(
               error_msg("type", "must be an integer >= 0"));
         }
         set_particle_type(m_pid, p_type);
       },
       [this]() { return particle().type(); }},
      {"pos",
       [this](Variant const &value) {
         place_particle(m_pid, get_value<Utils::Vector3d>(value));
       },
       [this]() {
         auto const &p = particle();
         return unfolded_position(p.pos(), p.image_box(), ::box_geo.length());
       }},
      {"v",
       [this](Variant const &value) {
         set_particle_v(m_pid, get_value<Utils::Vector3d>(value));
       },
       [this]() { return particle().v(); }},
      {"f",
       [this](Variant const &value) {
         set_particle_f(m_pid, get_value<Utils::Vector3d>(value));
       },
       [this]() { return particle().force(); }},
      {"mass",
#ifdef MASS
       [this](Variant const &value) {
         set_particle_mass(m_pid, get_value<double>(value));
       },
#else  // MASS
       [](Variant const &value) {
         if (std::abs(get_value<double>(value) - 1.) > 1e-10) {
           throw std::runtime_error("Feature MASS not compiled in");
         }
       },
#endif // MASS
       [this]() { return particle().mass(); }},
      {"q",
#ifdef ELECTROSTATICS
       [this](Variant const &value) {
         set_particle_q(m_pid, get_value<double>(value));
       },
#else  // ELECTROSTATICS
       [](Variant const &value) {
         if (get_value<double>(value) != 0.) {
           throw std::runtime_error("Feature ELECTROSTATICS not compiled in");
         }
       },
#endif // ELECTROSTATICS
       [this]() { return particle().q(); }},
      {"virtual",
#ifdef VIRTUAL_SITES
       [this](Variant const &value) {
         set_particle_virtual(m_pid, get_value<bool>(value));
       },
#else  // VIRTUAL_SITES
       [](Variant const &value) {
         if (get_value<bool>(value)) {
           throw std::runtime_error("Feature VIRTUAL_SITES not compiled in");
         }
       },
#endif // VIRTUAL_SITES
       [this]() { return particle().is_virtual(); }},
#ifdef ROTATION
      {"director",
       [this](Variant const &value) {
         set_particle_director(m_pid, get_value<Utils::Vector3d>(value));
       },
       [this]() { return particle().calc_director(); }},
      {"quat",
       [this](Variant const &value) {
         set_particle_quat(m_pid, get_quaternion_safe("quat", value));
       },
       [this]() { return quat2vector(particle().quat()); }},
      {"omega_body",
       [this](Variant const &value) {
         set_particle_omega_body(m_pid, get_value<Utils::Vector3d>(value));
       },
       [this]() { return particle().omega(); }},
      {"rotation",
       [this](Variant const &value) {
         set_particle_rotation(
             m_pid, Utils::Vector3i{get_value<Utils::Vector3b>(value)});
       },
       [this]() {
         auto const &p = particle();
         return Utils::Vector3b{{p.can_rotate_around(0), p.can_rotate_around(1),
                                 p.can_rotate_around(2)}};
       }},
      {"omega_lab",
       [this](Variant const &value) {
         set_particle_omega_lab(m_pid, get_value<Utils::Vector3d>(value));
       },
       [this]() {
         auto const &p = particle();
         return convert_vector_body_to_space(p, p.omega());
       }},
      {"torque_lab",
       [this](Variant const &value) {
         set_particle_torque_lab(m_pid, get_value<Utils::Vector3d>(value));
       },
       [this]() {
         auto const &p = particle();
         return convert_vector_body_to_space(p, p.torque());
       }},
#endif // ROTATION
#ifdef DIPOLES
      {"dip",
       [this](Variant const &value) {
         set_particle_dip(m_pid, get_value<Utils::Vector3d>(value));
       },
       [this]() { return particle().calc_dip(); }},
      {"dipm",
       [this](Variant const &value) {
         set_particle_dipm(m_pid, get_value<double>(value));
       },
       [this]() { return particle().dipm(); }},
#endif // DIPOLES
#ifdef ROTATIONAL_INERTIA
      {"rinertia",
       [this](Variant const &value) {
         set_particle_rotational_inertia(m_pid,
                                         get_value<Utils::Vector3d>(value));
       },
       [this]() { return particle().rinertia(); }},
#endif // ROTATIONAL_INERTIA
#ifdef LB_ELECTROHYDRODYNAMICS
      {"mu_E",
       [this](Variant const &value) {
         set_particle_mu_E(m_pid, get_value<Utils::Vector3d>(value));
       },
       [this]() { return particle().mu_E(); }},
#endif // LB_ELECTROHYDRODYNAMICS
#ifdef EXTERNAL_FORCES
      {"fix",
       [this](Variant const &value) {
         set_particle_fix(m_pid,
                          Utils::Vector3i{get_value<Utils::Vector3b>(value)});
       },
       [this]() {
         auto const &p = particle();
         return Utils::Vector3b{
             {p.is_fixed_along(0), p.is_fixed_along(1), p.is_fixed_along(2)}};
       }},
      {"ext_force",
       [this](Variant const &value) {
         set_particle_ext_force(m_pid, get_value<Utils::Vector3d>(value));
       },
       [this]() { return particle().ext_force(); }},
#ifdef ROTATION
      {"ext_torque",
       [this](Variant const &value) {
         set_particle_ext_torque(m_pid, get_value<Utils::Vector3d>(value));
       },
       [this]() { return particle().ext_torque(); }},
#endif // ROTATION
#endif // EXTERNAL_FORCES
#ifdef THERMOSTAT_PER_PARTICLE
      {"gamma",
       [this](Variant const &value) {
         set_particle_gamma(m_pid, get_gamma_safe(value));
       },
       [this]() { return particle().gamma(); }},
#ifdef ROTATION
      {"gamma_rot",
       [this](Variant const &value) {
         set_particle_gamma_rot(m_pid, get_gamma_safe(value));
       },
       [this]() { return particle().gamma_rot(); }},
#endif // ROTATION
#endif // THERMOSTAT_PER_PARTICLE
      {"pos_folded", AutoParameter::read_only,
       [this]() { return folded_position(particle().pos(), ::box_geo); }},

      {"lees_edwards_offset",
       [this](Variant const &value) {
         set_particle_lees_edwards_offset(m_pid, get_value<double>(value));
       },
       [this]() { return particle().lees_edwards_offset(); }},
      {"lees_edwards_flag", AutoParameter::read_only,
       [this]() { return particle().lees_edwards_flag(); }},
      {"image_box", AutoParameter::read_only,
       [this]() { return particle().image_box(); }},
      {"node", AutoParameter::read_only,
       [this]() { return get_particle_node(m_pid); }},
      {"mol_id",
       [this](Variant const &value) {
         auto const mol_id = get_value<int>(value);
         if (mol_id < 0) {
           throw std::domain_error(
               error_msg("mol_id", "must be an integer >= 0"));
         }
         set_particle_mol_id(m_pid, mol_id);
       },
       [this]() { return particle().mol_id(); }},
#ifdef VIRTUAL_SITES_RELATIVE
      {"vs_quat",
       [this](Variant const &value) {
         set_particle_vs_quat(m_pid, get_quaternion_safe("vs_quat", value));
       },
       [this]() { return quat2vector(particle().vs_relative().quat); }},
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
         set_particle_vs_relative(m_pid, vs_relative.to_particle_id,
                                  vs_relative.distance,
                                  vs_relative.rel_orientation);
       },
       [this]() {
         auto const &p = particle();
         return std::vector<Variant>{
             {p.vs_relative().to_particle_id, p.vs_relative().distance,
              quat2vector(p.vs_relative().rel_orientation)}};
       }},
#endif // VIRTUAL_SITES_RELATIVE
#ifdef ENGINE
      {"swimming",
       [this](Variant const &value) {
         ParticleParametersSwimming swim{};
         swim.swimming = true;
         auto const dict = get_value<VariantMap>(value);
         if (dict.count("f_swim") != 0) {
           swim.f_swim = get_value<double>(dict.at("f_swim"));
         }
         if (dict.count("v_swim") != 0) {
           swim.v_swim = get_value<double>(dict.at("v_swim"));
         }
         if (swim.f_swim != 0. and swim.v_swim != 0.) {
           throw std::invalid_argument(error_msg(
               "swimming",
               "cannot be set with 'v_swim' and 'f_swim' at the same time"));
         }
         if (dict.count("mode") != 0) {
           auto const mode = get_value<std::string>(dict.at("mode"));
           if (mode == "pusher") {
             swim.push_pull = -1;
           } else if (mode == "puller") {
             swim.push_pull = +1;
           } else if (mode == "N/A") {
             swim.push_pull = 0;
           } else {
             throw std::invalid_argument(
                 error_msg("swimming.mode",
                           "has to be either 'pusher', 'puller' or 'N/A'"));
           }
         }
         if (dict.count("dipole_length") != 0) {
           swim.dipole_length = get_value<double>(dict.at("dipole_length"));
         }
         set_particle_swimming(m_pid, swim);
       },
       [this]() {
         auto const &p = particle();
         auto const &swim = p.swimming();
         std::string mode;
         if (swim.push_pull == -1) {
           mode = "pusher";
         } else if (swim.push_pull == 1) {
           mode = "puller";
         } else {
           mode = "N/A";
         }
         return VariantMap{{{"mode", mode},
                            {"v_swim", swim.v_swim},
                            {"f_swim", swim.f_swim},
                            {"dipole_length", swim.dipole_length}}};
       }},
#endif // ENGINE
  });
}

Variant ParticleHandle::do_call_method(std::string const &name,
                                       VariantMap const &params) {
  if (name == "get_bonds_view") {
    std::vector<std::vector<int>> bonds_flat;
    for (auto const &bond_view : get_particle_bonds(m_pid)) {
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
    add_particle_bond(m_pid, get_bond_vector(params));
  } else if (name == "del_bond") {
    delete_particle_bond(m_pid, get_bond_vector(params));
  } else if (name == "delete_all_bonds") {
    delete_particle_bonds(m_pid);
  } else if (name == "is_valid_bond_id") {
    auto const bond_id = get_value<int>(params, "bond_id");
    return ::bonded_ia_params.get_zero_based_type(bond_id) != 0;
  }
  if (name == "remove_particle") {
    remove_particle(m_pid);
#ifdef VIRTUAL_SITES_RELATIVE
  } else if (name == "vs_relate_to") {
    vs_relate_to(m_pid, get_value<int>(params, "pid"));
#endif // VIRTUAL_SITES_RELATIVE
#ifdef EXCLUSIONS
  } else if (name == "has_exclusion") {
    auto const &p = get_particle_data(m_pid);
    return p.has_exclusion(get_value<int>(params, "pid"));
  }
  if (name == "add_exclusion") {
    add_particle_exclusion(m_pid, get_value<int>(params, "pid"));
  } else if (name == "del_exclusion") {
    remove_particle_exclusion(m_pid, get_value<int>(params, "pid"));
  } else if (name == "set_exclusions") {
    auto const &p = particle();
    for (auto const pid : p.exclusions_as_vector()) {
      remove_particle_exclusion(m_pid, pid);
    }
    std::vector<int> exclusion_list;
    try {
      auto const pid = get_value<int>(params, "p_ids");
      exclusion_list.push_back(pid);
    } catch (...) {
      exclusion_list = get_value<std::vector<int>>(params, "p_ids");
    }
    for (auto const pid : exclusion_list) {
      if (!p.has_exclusion(pid)) {
        add_particle_exclusion(m_pid, pid);
      }
    }
  } else if (name == "get_exclusions") {
    return particle().exclusions_as_vector();
#endif // EXCLUSIONS
#ifdef ROTATION
  }
  if (name == "rotate_particle") {
    rotate_particle(m_pid, get_value<Utils::Vector3d>(params, "axis"),
                    get_value<double>(params, "angle"));
  }
  if (name == "convert_vector_body_to_space") {
    auto const &p = get_particle_data(m_pid);
    return convert_vector_body_to_space(
               p, get_value<Utils::Vector3d>(params, "vec"))
        .as_vector();
  }
  if (name == "convert_vector_space_to_body") {
    auto const &p = get_particle_data(m_pid);
    return convert_vector_space_to_body(
               p, get_value<Utils::Vector3d>(params, "vec"))
        .as_vector();
#endif // ROTATION
  }
  return {};
}

void ParticleHandle::do_construct(VariantMap const &params) {
  m_pid = (params.count("id")) ? get_value<int>(params, "id")
                               : get_maximal_particle_id() + 1;

  // create a new particle if extra arguments were passed
  if ((params.size() - params.count("id")) > 0) {
    if (particle_exists(m_pid)) {
      throw std::invalid_argument("Particle " + std::to_string(m_pid) +
                                  " already exists");
    }

    // create a default-constructed particle
    auto const pos = get_value<Utils::Vector3d>(params, "pos");
    place_particle(m_pid, pos);

    // set particle properties (filter out read-only and deferred properties)
    std::vector<std::string> skip = {
        "id",        "pos",  "dipm",  "pos_folded", "lees_edwards_flag",
        "image_box", "node", "bonds", "exclusions",
    };
    for (auto const &kv : params) {
      if (std::find(skip.begin(), skip.end(), kv.first) == skip.end()) {
        do_set_parameter(kv.first, kv.second);
      }
    }
  }
}

} // namespace Particles
} // namespace ScriptInterface
