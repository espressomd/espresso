/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
/** \file
 *  Statistical tools to analyze simulations.
 *
 *  The corresponding header file is statistics.hpp.
 */

#include "statistics.hpp"

#include "Particle.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "communication.hpp"
#include "energy.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "integrate.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "npt.hpp"
#include "partCfg_global.hpp"
#include "pressure.hpp"
#include "short_range_loop.hpp"

#include <utils/NoOp.hpp>
#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/contains.hpp>

#include <cstdlib>
#include <limits>

/** Previous particle configurations (needed for offline analysis and
 *  correlation analysis)
 */
std::vector<std::vector<Utils::Vector3d>> configs;

int get_n_configs() { return static_cast<int>(configs.size()); }

int get_n_part_conf() {
  return (configs.size()) ? static_cast<int>(configs[0].size()) : 0;
}

/****************************************************************************************
 *                                 basic observables calculation
 ****************************************************************************************/

double mindist(PartCfg &partCfg, const std::vector<int> &set1,
               const std::vector<int> &set2) {
  using Utils::contains;

  auto mindist2 = std::numeric_limits<double>::infinity();

  for (auto jt = partCfg.begin(); jt != partCfg.end(); ++jt) {
    /* check which sets particle j belongs to (bit 0: set1, bit1: set2) */
    auto in_set = 0u;
    if (set1.empty() || contains(set1, jt->p.type))
      in_set = 1u;
    if (set2.empty() || contains(set2, jt->p.type))
      in_set |= 2u;
    if (in_set == 0)
      continue;

    for (auto it = std::next(jt); it != partCfg.end(); ++it)
      /* accept a pair if particle j is in set1 and particle i in set2 or vice
       * versa. */
      if (((in_set & 1u) && (set2.empty() || contains(set2, it->p.type))) ||
          ((in_set & 2u) && (set1.empty() || contains(set1, it->p.type))))
        mindist2 = std::min(mindist2,
                            get_mi_vector(jt->r.p, it->r.p, box_geo).norm2());
  }

  return std::sqrt(mindist2);
}

Utils::Vector3d local_particle_momentum() {
  auto const particles = cell_structure.local_particles();
  auto const momentum =
      std::accumulate(particles.begin(), particles.end(), Utils::Vector3d{},
                      [](Utils::Vector3d &m, Particle const &p) {
                        return m + p.p.mass * p.m.v;
                      });

  return momentum;
}

REGISTER_CALLBACK_REDUCTION(local_particle_momentum,
                            std::plus<Utils::Vector3d>())

Utils::Vector3d calc_linear_momentum(int include_particles,
                                     int include_lbfluid) {
  Utils::Vector3d linear_momentum{};
  if (include_particles) {
    linear_momentum +=
        mpi_call(::Communication::Result::reduction,
                 std::plus<Utils::Vector3d>(), local_particle_momentum);
  }
  if (include_lbfluid) {
    linear_momentum += lb_lbfluid_calc_fluid_momentum();
  }
  return linear_momentum;
}

Utils::Vector3d centerofmass(PartCfg &partCfg, int type) {
  Utils::Vector3d com{};
  double mass = 0.0;

  for (auto const &p : partCfg) {
    if ((p.p.type == type) || (type == -1))
      if (not p.p.is_virtual) {
        com += p.r.p * p.p.mass;
        mass += p.p.mass;
      }
  }
  com /= mass;
  return com;
}

Utils::Vector3d angularmomentum(PartCfg &partCfg, int type) {
  Utils::Vector3d am{};

  for (auto const &p : partCfg) {
    if ((p.p.type == type) || (type == -1))
      if (not p.p.is_virtual) {
        am += p.p.mass * vector_product(p.r.p, p.m.v);
      }
  }
  return am;
}

void momentofinertiamatrix(PartCfg &partCfg, int type, double *MofImatrix) {
  int i, count;
  double massi;
  Utils::Vector3d p1{};
  count = 0;

  for (i = 0; i < 9; i++)
    MofImatrix[i] = 0.;

  auto const com = centerofmass(partCfg, type);
  for (auto const &p : partCfg) {
    if (type == p.p.type and (not p.p.is_virtual)) {
      count++;
      p1 = p.r.p - com;
      massi = p.p.mass;
      MofImatrix[0] += massi * (p1[1] * p1[1] + p1[2] * p1[2]);
      MofImatrix[4] += massi * (p1[0] * p1[0] + p1[2] * p1[2]);
      MofImatrix[8] += massi * (p1[0] * p1[0] + p1[1] * p1[1]);
      MofImatrix[1] -= massi * (p1[0] * p1[1]);
      MofImatrix[2] -= massi * (p1[0] * p1[2]);
      MofImatrix[5] -= massi * (p1[1] * p1[2]);
    }
  }
  /* use symmetry */
  MofImatrix[3] = MofImatrix[1];
  MofImatrix[6] = MofImatrix[2];
  MofImatrix[7] = MofImatrix[5];
}

std::vector<int> nbhood(PartCfg &partCfg, const Utils::Vector3d &pos,
                        double r_catch, const Utils::Vector3i &planedims) {
  std::vector<int> ids;

  auto const r2 = r_catch * r_catch;
  auto const pt = Utils::Vector3d{pos[0], pos[1], pos[2]};

  Utils::Vector3d d;

  for (auto const &p : partCfg) {
    if ((planedims[0] + planedims[1] + planedims[2]) == 3) {
      d = get_mi_vector(pt, p.r.p, box_geo);
    } else {
      /* Calculate the in plane distance */
      for (int j = 0; j < 3; j++) {
        d[j] = planedims[j] * (p.r.p[j] - pt[j]);
      }
    }

    if (d.norm2() < r2) {
      ids.push_back(p.p.identity);
    }
  }

  return ids;
}

double distto(PartCfg &partCfg, const Utils::Vector3d &pos, int pid) {
  auto mindist = std::numeric_limits<double>::infinity();

  for (auto const &part : partCfg) {
    if (pid != part.p.identity) {
      auto const d = get_mi_vector({pos[0], pos[1], pos[2]}, part.r.p, box_geo);
      mindist = std::min(mindist, d.norm2());
    }
  }
  return std::sqrt(mindist);
}

void calc_part_distribution(PartCfg &partCfg, std::vector<int> const &p1_types,
                            std::vector<int> const &p2_types, double r_min,
                            double r_max, int r_bins, bool log_flag,
                            double *low, double *dist) {
  int ind, cnt = 0;
  double inv_bin_width = 0.0;
  double min_dist, min_dist2 = 0.0, start_dist2;

  start_dist2 = Utils::sqr(box_geo.length()[0] + box_geo.length()[1] +
                           box_geo.length()[2]);
  /* bin preparation */
  *low = 0.0;
  for (int i = 0; i < r_bins; i++)
    dist[i] = 0.0;
  if (log_flag)
    inv_bin_width = (double)r_bins / (log(r_max) - log(r_min));
  else
    inv_bin_width = (double)r_bins / (r_max - r_min);

  /* particle loop: p1_types */
  for (auto const &p1 : partCfg) {
    for (int t1 : p1_types) {
      if (p1.p.type == t1) {
        min_dist2 = start_dist2;
        /* particle loop: p2_types */
        for (auto const &p2 : partCfg) {
          if (p1 != p2) {
            for (int t2 : p2_types) {
              if (p2.p.type == t2) {
                auto const act_dist2 =
                    get_mi_vector(p1.r.p, p2.r.p, box_geo).norm2();
                if (act_dist2 < min_dist2) {
                  min_dist2 = act_dist2;
                }
              }
            }
          }
        }
        min_dist = sqrt(min_dist2);
        if (min_dist <= r_max) {
          if (min_dist >= r_min) {
            /* calculate bin index */
            if (log_flag)
              ind = (int)((log(min_dist) - log(r_min)) * inv_bin_width);
            else
              ind = (int)((min_dist - r_min) * inv_bin_width);
            if (ind >= 0 && ind < r_bins) {
              dist[ind] += 1.0;
            }
          } else {
            *low += 1.0;
          }
        }
        cnt++;
      }
    }
  }
  if (cnt == 0)
    return;

  /* normalization */
  *low /= (double)cnt;
  for (int i = 0; i < r_bins; i++)
    dist[i] /= (double)cnt;
}

void calc_rdf(PartCfg &partCfg, std::vector<int> const &p1_types,
              std::vector<int> const &p2_types, double r_min, double r_max,
              int r_bins, std::vector<double> &rdf) {
  calc_rdf(partCfg, &p1_types[0], p1_types.size(), &p2_types[0],
           p2_types.size(), r_min, r_max, r_bins, &rdf[0]);
}

void calc_rdf(PartCfg &partCfg, int const *p1_types, int n_p1,
              int const *p2_types, int n_p2, double r_min, double r_max,
              int r_bins, double *rdf) {
  long int cnt = 0;
  int ind;
  bool mixed_flag = false;
  if (n_p1 == n_p2) {
    for (int i = 0; i < n_p1; i++)
      if (p1_types[i] != p2_types[i])
        mixed_flag = true;
  } else {
    mixed_flag = true;
  }

  auto const bin_width = (r_max - r_min) / (double)r_bins;
  auto const inv_bin_width = 1.0 / bin_width;
  for (int i = 0; i < r_bins; i++)
    rdf[i] = 0.0;
  /* particle loop: p1_types */
  for (auto it = partCfg.begin(); it != partCfg.end(); ++it) {
    for (int t1 = 0; t1 < n_p1; t1++) {
      if (it->p.type == p1_types[t1]) {
        /* distinguish mixed and identical rdf's */
        auto jt = mixed_flag ? partCfg.begin() : std::next(it);

        /* particle loop: p2_types */
        for (; jt != partCfg.end(); ++jt) {
          for (int t2 = 0; t2 < n_p2; t2++) {
            if (jt->p.type == p2_types[t2]) {
              auto const dist = get_mi_vector(it->r.p, jt->r.p, box_geo).norm();
              if (dist > r_min && dist < r_max) {
                ind = (int)((dist - r_min) * inv_bin_width);
                rdf[ind]++;
              }
              cnt++;
            }
          }
        }
      }
    }
  }
  if (cnt == 0)
    return;

  /* normalization */
  auto const volume = box_geo.volume();
  for (int i = 0; i < r_bins; i++) {
    auto const r_in = i * bin_width + r_min;
    auto const r_out = r_in + bin_width;
    auto const bin_volume = (4.0 / 3.0) * Utils::pi() *
                            ((r_out * r_out * r_out) - (r_in * r_in * r_in));
    rdf[i] *= volume / (bin_volume * static_cast<double>(cnt));
  }
}

void calc_rdf_av(PartCfg &partCfg, std::vector<int> const &p1_types,
                 std::vector<int> const &p2_types, double r_min, double r_max,
                 int r_bins, std::vector<double> &rdf, int n_conf) {
  calc_rdf_av(partCfg, &p1_types[0], p1_types.size(), &p2_types[0],
              p2_types.size(), r_min, r_max, r_bins, &rdf[0], n_conf);
}

void calc_rdf_av(PartCfg &partCfg, int const *p1_types, int n_p1,
                 int const *p2_types, int n_p2, double r_min, double r_max,
                 int r_bins, double *rdf, int n_conf) {
  long int cnt = 0;
  int cnt_conf = 1;
  bool mixed_flag = false;
  std::vector<double> rdf_tmp(r_bins);

  if (n_p1 == n_p2) {
    for (int i = 0; i < n_p1; i++)
      if (p1_types[i] != p2_types[i])
        mixed_flag = true;
  } else
    mixed_flag = true;

  auto const bin_width = (r_max - r_min) / (double)r_bins;
  auto const inv_bin_width = 1.0 / bin_width;
  auto const volume = box_geo.volume();
  for (int l = 0; l < r_bins; l++)
    rdf_tmp[l] = rdf[l] = 0.0;

  while (cnt_conf <= n_conf) {
    for (int l = 0; l < r_bins; l++)
      rdf_tmp[l] = 0.0;
    cnt = 0;
    auto const k = configs.size() - cnt_conf;
    int i = 0;
    for (auto it = partCfg.begin(); it != partCfg.end(); ++it) {
      for (int t1 = 0; t1 < n_p1; t1++) {
        if (it->p.type == p1_types[t1]) {
          /* distinguish mixed and identical rdf's */
          auto jt = mixed_flag ? partCfg.begin() : std::next(it);
          int j = mixed_flag ? 0 : i + 1;

          // particle loop: p2_types
          for (; jt != partCfg.end(); ++jt) {
            for (int t2 = 0; t2 < n_p2; t2++) {
              if (jt->p.type == p2_types[t2]) {
                auto const dist =
                    get_mi_vector(configs[k][i], configs[k][j], box_geo).norm();
                if (dist > r_min && dist < r_max) {
                  auto const ind =
                      static_cast<int>((dist - r_min) * inv_bin_width);
                  rdf_tmp[ind]++;
                }
                cnt++;
              }
            }
            j++;
          }
        }
      }
      i++;
    }
    // normalization

    for (int i = 0; i < r_bins; i++) {
      auto const r_in = i * bin_width + r_min;
      auto const r_out = r_in + bin_width;
      auto const bin_volume = (4.0 / 3.0) * Utils::pi() *
                              ((r_out * r_out * r_out) - (r_in * r_in * r_in));
      rdf[i] += rdf_tmp[i] * volume / (bin_volume * static_cast<double>(cnt));
    }

    cnt_conf++;
  } // cnt_conf loop
  for (int i = 0; i < r_bins; i++) {
    rdf[i] /= (cnt_conf - 1);
  }
}

std::vector<double> calc_structurefactor(PartCfg &partCfg,
                                         std::vector<int> const &p_types,
                                         int order) {
  auto const order2 = order * order;
  std::vector<double> ff;
  ff.resize(2 * order2);
  ff[2 * order2] = 0;
  auto const twoPI_L = 2 * Utils::pi() / box_geo.length()[0];

  if (order < 1) {
    fprintf(stderr,
            "WARNING: parameter \"order\" has to be a whole positive number");
    fflush(nullptr);
    errexit();
  } else {
    for (int qi = 0; qi < 2 * order2; qi++) {
      ff[qi] = 0.0;
    }
    for (int i = 0; i <= order; i++) {
      for (int j = -order; j <= order; j++) {
        for (int k = -order; k <= order; k++) {
          auto const n = i * i + j * j + k * k;
          if ((n <= order2) && (n >= 1)) {
            double C_sum = 0.0, S_sum = 0.0;
            for (auto const &p : partCfg) {
              for (int t : p_types) {
                if (p.p.type == t) {
                  auto const qr =
                      twoPI_L * (Utils::Vector3i{{i, j, k}} * p.r.p);
                  C_sum += cos(qr);
                  S_sum += sin(qr);
                }
              }
            }
            ff[2 * n - 2] += C_sum * C_sum + S_sum * S_sum;
            ff[2 * n - 1]++;
          }
        }
      }
    }
    int n = 0;
    for (auto const &p : partCfg) {
      for (int t : p_types) {
        if (p.p.type == t)
          n++;
      }
    }
    for (int qi = 0; qi < order2; qi++)
      if (ff[2 * qi + 1] != 0)
        ff[2 * qi] /= n * ff[2 * qi + 1];
  }
  return ff;
}

std::vector<std::vector<double>> modify_stucturefactor(int order,
                                                       double const *sf) {
  int length = 0;

  for (int i = 0; i < order * order; i++) {
    if (sf[2 * i + 1] > 0) {
      length++;
    }
  }

  auto const qfak = 2.0 * Utils::pi() / box_geo.length()[0];
  std::vector<double> intern;
  intern.assign(2, 0.0);
  std::vector<std::vector<double>> structure_factor;
  structure_factor.assign(length, intern);

  int cnt = 0;
  for (int i = 0; i < order * order; i++) {
    if (sf[2 * i + 1] > 0) {
      structure_factor[cnt][0] = qfak * sqrt(i + 1);
      structure_factor[cnt][1] = sf[2 * i];
      cnt++;
    }
  }

  return structure_factor;
}

/****************************************************************************************
 *                                 config storage functions
 ****************************************************************************************/

void analyze_append(PartCfg &partCfg) {
  std::vector<Utils::Vector3d> config;
  for (auto const &p : partCfg) {
    config.emplace_back(p.r.p);
  }
  configs.emplace_back(config);
}

/****************************************************************************************
 *                                 Observables handling
 ****************************************************************************************/

void obsstat_realloc_and_clear(Observable_stat *stat, int n_pre, int n_bonded,
                               int n_non_bonded, int n_coulomb, int n_dipolar,
                               int n_vs, int c_size) {

  // Number of doubles to store pressure in
  const int total =
      c_size *
      (n_pre + static_cast<int>(bonded_ia_params.size()) + n_non_bonded +
       n_coulomb + n_dipolar + n_vs + Observable_stat::n_external_field);

  // Allocate mem for the double list
  stat->data.resize(total);

  // Number of doubles per interaction (pressure=1, stress tensor=9,...)
  stat->chunk_size = c_size;

  // Number of chunks for different interaction types
  stat->n_coulomb = n_coulomb;
  stat->n_dipolar = n_dipolar;
  stat->n_virtual_sites = n_vs;
  // Pointers to the start of different contributions
  stat->bonded = stat->data.data() + c_size * n_pre;
  stat->non_bonded = stat->bonded + c_size * bonded_ia_params.size();
  stat->coulomb = stat->non_bonded + c_size * n_non_bonded;
  stat->dipolar = stat->coulomb + c_size * n_coulomb;
  stat->virtual_sites = stat->dipolar + c_size * n_dipolar;
  stat->external_fields = stat->virtual_sites + c_size * n_vs;

  // Set all observables to zero
  for (int i = 0; i < total; i++)
    stat->data[i] = 0.0;
}

void obsstat_realloc_and_clear_non_bonded(Observable_stat_non_bonded *stat_nb,
                                          int n_nonbonded, int c_size) {
  auto const total = c_size * (n_nonbonded + n_nonbonded);

  stat_nb->data_nb.resize(total);
  stat_nb->chunk_size_nb = c_size;
  stat_nb->non_bonded_intra = stat_nb->data_nb.data();
  stat_nb->non_bonded_inter = stat_nb->non_bonded_intra + c_size * n_nonbonded;

  for (int i = 0; i < total; i++)
    stat_nb->data_nb[i] = 0.0;
}

void invalidate_obs() {
  total_energy.init_status = 0;
  total_pressure.init_status = 0;
  total_p_tensor.init_status = 0;
}

void update_pressure(int v_comp) {
  Utils::Vector3d p_vel;
  /* if desired (v_comp==1) replace ideal component with instantaneous one */
  if (total_pressure.init_status != 1 + v_comp) {
    init_virials(&total_pressure);
    init_p_tensor(&total_p_tensor);

    init_virials_non_bonded(&total_pressure_non_bonded);
    init_p_tensor_non_bonded(&total_p_tensor_non_bonded);

    if (v_comp && (integ_switch == INTEG_METHOD_NPT_ISO) &&
        !(nptiso.invalidate_p_vel)) {
      if (total_pressure.init_status == 0)
        master_pressure_calc(0);
      total_pressure.data[0] = 0.0;
      MPI_Reduce(nptiso.p_vel.data(), p_vel.data(), 3, MPI_DOUBLE, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      for (int i = 0; i < 3; i++)
        if (nptiso.geometry & nptiso.nptgeom_dir[i])
          total_pressure.data[0] += p_vel[i];
      total_pressure.data[0] /= (nptiso.dimension * nptiso.volume);
      total_pressure.init_status = 1 + v_comp;
    } else
      master_pressure_calc(v_comp);
  }
}
