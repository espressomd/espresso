/*
 * Copyright (C) 2013-2022 The ESPResSo project
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

#include "System.hpp"

#include "config/config.hpp"

#include "core/BoxGeometry.hpp"
#include "core/Particle.hpp"
#include "core/ParticleRange.hpp"
#include "core/cell_system/CellStructure.hpp"
#include "core/cell_system/CellStructureType.hpp"
#include "core/cells.hpp"
#include "core/communication.hpp"
#include "core/nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "core/object-in-fluid/oif_global_forces.hpp"
#include "core/particle_node.hpp"
#include "core/propagation.hpp"
#include "core/rotation.hpp"
#include "core/system/System.hpp"
#include "core/system/System.impl.hpp"

#include "script_interface/analysis/Analysis.hpp"
#include "script_interface/bond_breakage/BreakageSpecs.hpp"
#include "script_interface/cell_system/CellSystem.hpp"
#include "script_interface/electrostatics/Container.hpp"
#include "script_interface/galilei/ComFixed.hpp"
#include "script_interface/galilei/Galilei.hpp"
#include "script_interface/integrators/IntegratorHandle.hpp"
#include "script_interface/lees_edwards/LeesEdwards.hpp"
#include "script_interface/magnetostatics/Container.hpp"
#include "script_interface/thermostat/thermostat.hpp"

#include <utils/Vector.hpp>
#include <utils/math/vec_rotate.hpp>

#include <boost/mpi/collectives.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace ScriptInterface {
namespace System {

static bool system_created = false;

/** @brief Container for leaves of the system class. */
struct System::Leaves {
  Leaves() = default;
  std::shared_ptr<CellSystem::CellSystem> cell_system;
  std::shared_ptr<Integrators::IntegratorHandle> integrator;
  std::shared_ptr<Thermostat::Thermostat> thermostat;
  std::shared_ptr<Analysis::Analysis> analysis;
  std::shared_ptr<Galilei::ComFixed> comfixed;
  std::shared_ptr<Galilei::Galilei> galilei;
  std::shared_ptr<BondBreakage::BreakageSpecs> bond_breakage;
  std::shared_ptr<LeesEdwards::LeesEdwards> lees_edwards;
#ifdef ELECTROSTATICS
  std::shared_ptr<Coulomb::Container> electrostatics;
#endif
#ifdef DIPOLES
  std::shared_ptr<Dipoles::Container> magnetostatics;
#endif
};

System::System() : m_instance{}, m_leaves{std::make_shared<Leaves>()} {
  auto const add_parameter =
      [this, ptr = m_leaves.get()](std::string key, auto(Leaves::*member)) {
        add_parameters({AutoParameter(
            key.c_str(),
            [this, ptr, member, key](Variant const &val) {
              auto &dst = ptr->*member;
              if (dst != nullptr) {
                throw WriteError(key);
              }
              dst = get_value<std::remove_reference_t<decltype(dst)>>(val);
              dst->bind_system(m_instance);
            },
            [ptr, member]() { return ptr->*member; })});
      };

  add_parameters({
      {"box_l",
       [this](Variant const &v) {
         context()->parallel_try_catch([&]() {
           auto const new_value = get_value<Utils::Vector3d>(v);
           if (not(new_value > Utils::Vector3d::broadcast(0.))) {
             throw std::domain_error("Attribute 'box_l' must be > 0");
           }
           m_instance->veto_boxl_change();
           m_instance->box_geo->set_length(new_value);
           m_instance->on_boxl_change();
         });
       },
       [this]() { return m_instance->box_geo->length(); }},
      {"periodicity",
       [this](Variant const &v) {
         auto const periodicity = get_value<Utils::Vector3b>(v);
         for (unsigned int i = 0u; i < 3u; ++i) {
           m_instance->box_geo->set_periodic(i, periodicity[i]);
         }
         context()->parallel_try_catch(
             [&]() { m_instance->on_periodicity_change(); });
       },
       [this]() {
         return Utils::Vector3b{m_instance->box_geo->periodic(0),
                                m_instance->box_geo->periodic(1),
                                m_instance->box_geo->periodic(2)};
       }},
      {"min_global_cut",
       [this](Variant const &v) {
         context()->parallel_try_catch([&]() {
           auto const new_value = get_value<double>(v);
           if (new_value < 0. and new_value != INACTIVE_CUTOFF) {
             throw std::domain_error("Attribute 'min_global_cut' must be >= 0");
           }
           m_instance->set_min_global_cut(new_value);
         });
       },
       [this]() { return m_instance->get_min_global_cut(); }},
      {"max_oif_objects", ::max_oif_objects},
  });
  add_parameter("cell_system", &Leaves::cell_system);
  add_parameter("integrator", &Leaves::integrator);
  add_parameter("thermostat", &Leaves::thermostat);
  add_parameter("analysis", &Leaves::analysis);
  add_parameter("comfixed", &Leaves::comfixed);
  add_parameter("galilei", &Leaves::galilei);
  add_parameter("bond_breakage", &Leaves::bond_breakage);
  add_parameter("lees_edwards", &Leaves::lees_edwards);
#ifdef ELECTROSTATICS
  add_parameter("electrostatics", &Leaves::electrostatics);
#endif
#ifdef DIPOLES
  add_parameter("magnetostatics", &Leaves::magnetostatics);
#endif
}

void System::do_construct(VariantMap const &params) {
  /* When reloading the system state from a checkpoint file,
   * the order of global variables instantiation matters.
   * The @c box_l must be set before any other global variable.
   * All these globals re-initialize the cell system, and we
   * cannot use the default-constructed @c box_geo when e.g.
   * long-range interactions exist in the system, otherwise
   * runtime errors about the local geometry being smaller
   * than the interaction range would be raised.
   */
  m_instance = ::System::System::create();
  ::System::set_system(m_instance);

  // domain decomposition can only be set after box_l is set
  m_instance->set_cell_structure_topology(CellStructureType::NSQUARE);
  do_set_parameter("box_l", params.at("box_l"));
  m_instance->set_cell_structure_topology(CellStructureType::REGULAR);

  m_instance->lb.bind_system(m_instance);
  m_instance->ek.bind_system(m_instance);

  for (auto const &key : get_parameter_insertion_order()) {
    if (key != "box_l" and params.count(key) != 0ul) {
      do_set_parameter(key, params.at(key));
    }
  }
}

static void rotate_system(CellStructure &cell_structure, double phi,
                          double theta, double alpha) {
  auto const particles = cell_structure.local_particles();

  // Calculate center of mass
  Utils::Vector3d local_com{};
  double local_mass = 0.0;

  for (auto const &p : particles) {
    if (not p.is_virtual()) {
      local_com += p.mass() * p.pos();
      local_mass += p.mass();
    }
  }

  auto const total_mass =
      boost::mpi::all_reduce(comm_cart, local_mass, std::plus<>());
  auto const com =
      boost::mpi::all_reduce(comm_cart, local_com, std::plus<>()) / total_mass;

  // Rotation axis in Cartesian coordinates
  Utils::Vector3d axis;
  axis[0] = std::sin(theta) * std::cos(phi);
  axis[1] = std::sin(theta) * std::sin(phi);
  axis[2] = std::cos(theta);

  // Rotate particle coordinates
  for (auto &p : particles) {
    // Move the center of mass of the system to the origin
    p.pos() = com + Utils::vec_rotate(axis, alpha, p.pos() - com);
#ifdef ROTATION
    local_rotate_particle(p, axis, alpha);
#endif
  }

  cell_structure.set_resort_particles(Cells::RESORT_GLOBAL);
  cell_structure.update_ghosts_and_resort_particle(Cells::DATA_PART_POSITION |
                                                   Cells::DATA_PART_PROPERTIES);
}

/** Rescale all particle positions in direction @p dir by a factor @p scale.
 *  @param cell_structure cell structure
 *  @param dir   direction to scale (0/1/2 = x/y/z, 3 = x+y+z isotropically)
 *  @param scale factor by which to rescale (>1: stretch, <1: contract)
 */
static void rescale_particles(CellStructure &cell_structure, int dir,
                              double scale) {
  for (auto &p : cell_structure.local_particles()) {
    if (dir < 3)
      p.pos()[dir] *= scale;
    else {
      p.pos() *= scale;
    }
  }
}

Variant System::do_call_method(std::string const &name,
                               VariantMap const &parameters) {
  if (name == "is_system_created") {
    return system_created;
  }
  if (name == "lock_system_creation") {
    system_created = true;
    return {};
  }
  if (name == "rescale_boxl") {
    auto &box_geo = *m_instance->box_geo;
    auto const coord = get_value<int>(parameters, "coord");
    auto const length = get_value<double>(parameters, "length");
    assert(coord >= 0);
    assert(coord != 3 or ((box_geo.length()[0] == box_geo.length()[1]) and
                          (box_geo.length()[1] == box_geo.length()[2])));
    auto const scale = (coord == 3) ? length * box_geo.length_inv()[0]
                                    : length * box_geo.length_inv()[coord];
    context()->parallel_try_catch([&]() {
      if (length <= 0.) {
        throw std::domain_error("Parameter 'd_new' must be > 0");
      }
      m_instance->veto_boxl_change(true);
    });
    auto new_value = Utils::Vector3d{};
    if (coord == 3) {
      new_value = Utils::Vector3d::broadcast(length);
    } else {
      new_value = box_geo.length();
      new_value[static_cast<unsigned>(coord)] = length;
    }
    // when shrinking, rescale the particles first
    if (scale <= 1.) {
      rescale_particles(*m_instance->cell_structure, coord, scale);
      m_instance->on_particle_change();
    }
    m_instance->box_geo->set_length(new_value);
    m_instance->on_boxl_change();
    if (scale > 1.) {
      rescale_particles(*m_instance->cell_structure, coord, scale);
      m_instance->on_particle_change();
    }
    return {};
  }
  if (name == "setup_type_map") {
    auto const types = get_value<std::vector<int>>(parameters, "type_list");
    for (auto const type : types) {
      ::init_type_map(type);
    }
    return {};
  }
  if (name == "number_of_particles") {
    auto const type = get_value<int>(parameters, "type");
    return ::number_of_particles_with_type(type);
  }
  if (name == "velocity_difference") {
    auto const pos1 = get_value<Utils::Vector3d>(parameters, "pos1");
    auto const pos2 = get_value<Utils::Vector3d>(parameters, "pos2");
    auto const v1 = get_value<Utils::Vector3d>(parameters, "v1");
    auto const v2 = get_value<Utils::Vector3d>(parameters, "v2");
    return m_instance->box_geo->velocity_difference(pos2, pos1, v2, v1);
  }
  if (name == "distance_vec") {
    auto const pos1 = get_value<Utils::Vector3d>(parameters, "pos1");
    auto const pos2 = get_value<Utils::Vector3d>(parameters, "pos2");
    return m_instance->box_geo->get_mi_vector(pos2, pos1);
  }

  if (name == "cutoff_by_types") {
    auto types = get_value<std::vector<int>>(parameters, "types");
    return m_instance->maximal_cutoff(types);
  }

  if (name == "rotate_system") {
    rotate_system(*m_instance->cell_structure,
                  get_value<double>(parameters, "phi"),
                  get_value<double>(parameters, "theta"),
                  get_value<double>(parameters, "alpha"));
    m_instance->on_particle_change();
    m_instance->update_dependent_particles();
    return {};
  }
  if (name == "get_propagation_modes_enum") {
    return make_unordered_map_of_variants(propagation_flags_map());
  }
  if (name == "session_shutdown") {
    if (m_instance) {
      if (&::System::get_system() == m_instance.get()) {
        ::System::reset_system();
      }
      assert(m_instance.use_count() == 1l);
      m_instance.reset();
    }
    return {};
  }
  return {};
}

} // namespace System
} // namespace ScriptInterface
