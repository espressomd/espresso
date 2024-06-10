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

#include "config/config.hpp"

#include "ObservableStat.hpp"

#include "core/bonded_interactions/bonded_interaction_data.hpp"
#include "core/nonbonded_interactions/nonbonded_interaction_data.hpp"

#include "core/Observable_stat.hpp"

#include <utils/Span.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Analysis {

/**
 * @brief Generate an observable summary.
 * @param[in] system   The system to analyze.
 * @param[in] obs      The observable handle.
 * @param[in] calc_sp  Whether to compute a scalar pressure.
 */
static auto get_summary(::System::System const &system,
                        Observable_stat const &obs, bool const calc_sp) {
  auto const obs_dim = obs.get_chunk_size();

  auto const get_obs_contribs = [obs_dim,
                                 calc_sp](Utils::Span<double> const views) {
    if (obs_dim == 1) {
      return std::vector<Variant>(views.begin(), views.end());
    }
    assert(obs_dim == 9ul);
    assert(views.size() % 9ul == 0ul);
    std::vector<Variant> out;
    for (std::size_t i = 0ul; i < views.size() / 9ul; ++i) {
      auto const view = Utils::Span<double>{views.data() + i * 9ul, 9ul};
      if (calc_sp) {
        auto const trace = view[0] + view[4] + view[8];
        out.emplace_back(trace / 3.);
      } else {
        auto const flat_matrix = std::vector<double>(view.begin(), view.end());
        out.emplace_back(flat_matrix);
      }
    }
    return out;
  };

  auto const get_obs_contrib =
      [&get_obs_contribs](Utils::Span<double> const views) -> Variant {
    return get_obs_contribs(views)[0];
  };

  std::unordered_map<std::string, Variant> dict;
  dict["kinetic"] = get_obs_contrib(obs.kinetic);
  dict["external_fields"] = get_obs_contrib(obs.external_fields);

  {
    auto values = std::vector<double>(obs_dim);
    for (std::size_t i = 0ul; i < obs_dim; ++i) {
      values[i] = obs.accumulate(0., i);
    }
    dict["total"] = get_obs_contrib({values.data(), obs_dim});
  }

  auto const n_bonds = static_cast<int>(::bonded_ia_params.get_next_key());
  for (int bond_id = 0; bond_id < n_bonds; ++bond_id) {
    if (::bonded_ia_params.get_zero_based_type(bond_id) != 0) {
      dict["bonded," + std::to_string(bond_id)] =
          get_obs_contrib(obs.bonded_contribution(bond_id));
    }
  }

  auto const n_nonbonded =
      system.nonbonded_ias->get_max_seen_particle_type() + 1;
  for (int i = 0; i < n_nonbonded; ++i) {
    for (int j = i; j < n_nonbonded; ++j) {
      auto const indices = std::to_string(i) + "," + std::to_string(j);
      dict["non_bonded_intra," + indices] =
          get_obs_contrib(obs.non_bonded_intra_contribution(i, j));
      dict["non_bonded_inter," + indices] =
          get_obs_contrib(obs.non_bonded_inter_contribution(i, j));
    }
  }

#ifdef ELECTROSTATICS
  {
    auto const values = get_obs_contribs(obs.coulomb);
    for (std::size_t i = 0ul; i < values.size(); ++i) {
      dict["coulomb," + std::to_string(i)] = values[i];
    }
  }
#endif // ELECTROSTATICS

#ifdef DIPOLES
  {
    auto const values = get_obs_contribs(obs.dipolar);
    for (std::size_t i = 0ul; i < values.size(); ++i) {
      dict["dipolar," + std::to_string(i)] = values[i];
    }
  }
#endif // DIPOLES

#ifdef VIRTUAL_SITES
  {
    auto const values = get_obs_contribs(obs.virtual_sites);
    for (std::size_t i = 0ul; i < values.size(); ++i) {
      dict["virtual_sites," + std::to_string(i)] = values[i];
    }
  }
#endif // VIRTUAL_SITES

  return dict;
}

Variant ObservableStat::do_call_method(std::string const &name,
                                       VariantMap const &parameters) {
  auto &system = get_system();
  if (name == "calculate_energy") {
    auto const obs = system.calculate_energy();
    return get_summary(system, *obs, false);
  }
  if (name == "calculate_scalar_pressure") {
    auto const obs = system.calculate_pressure();
    return get_summary(system, *obs, true);
  }
  if (name == "calculate_pressure_tensor") {
    auto const obs = system.calculate_pressure();
    return get_summary(system, *obs, false);
  }
  return {};
}

} // namespace Analysis
} // namespace ScriptInterface
