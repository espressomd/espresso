/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/** \file
 *  Statistical tools to analyze simulations.
 *
 *  The corresponding header file is statistics.hpp.
 */

#include "statistics.hpp"

#include "communication.hpp"
#include "energy.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "npt.hpp"
#include "partCfg_global.hpp"
#include "particle_data.hpp"
#include "pressure.hpp"
#include "short_range_loop.hpp"
#include "statistics_chain.hpp"
#include "virtual_sites.hpp"

#include <utils/NoOp.hpp>
#include <utils/constants.hpp>
#include <utils/contains.hpp>

#include <cstdlib>
#include <cstring>
#include <limits>

/** Previous particle configurations (needed for offline analysis and
    correlation analysis) */
double **configs = nullptr;
int n_configs = 0;
int n_part_conf = 0;

/****************************************************************************************
 *                                 helper functions
 ****************************************************************************************/
/****************************************************************************************
 *                                 basic observables calculation
 ****************************************************************************************/

double mindist(PartCfg &partCfg, IntList const &set1, IntList const &set2) {
  int in_set;

  auto mindist2 = std::numeric_limits<double>::infinity();

  for (auto jt = partCfg.begin(); jt != (--partCfg.end()); ++jt) {
    /* check which sets particle j belongs to
       bit 0: set1, bit1: set2
    */
    in_set = 0;
    if (set1.empty() || contains(set1, jt->p.type))
      in_set = 1;
    if (set2.empty() || contains(set2, jt->p.type))
      in_set |= 2;
    if (in_set == 0)
      continue;

    for (auto it = std::next(jt); it != partCfg.end(); ++it)
      /* accept a pair if particle j is in set1 and particle i in set2 or vice
       * versa. */
      if (((in_set & 1) && (set2.empty() || contains(set2, it->p.type))) ||
          ((in_set & 2) && (set1.empty() || contains(set1, it->p.type))))
        mindist2 = std::min(mindist2, min_distance2(jt->r.p, it->r.p));
  }

  return std::sqrt(mindist2);
}

void predict_momentum_particles(double *result) {
  double momentum[3] = {0.0, 0.0, 0.0};

  for (auto const &p : local_cells.particles()) {
    auto const mass = p.p.mass;

    momentum[0] += mass * (p.m.v[0] + p.f.f[0] * 0.5 * time_step / p.p.mass);
    momentum[1] += mass * (p.m.v[1] + p.f.f[1] * 0.5 * time_step / p.p.mass);
    momentum[2] += mass * (p.m.v[2] + p.f.f[2] * 0.5 * time_step / p.p.mass);
  }

  MPI_Reduce(momentum, result, 3, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
}

Utils::Vector3d calc_linear_momentum(int include_particles,
                                     int include_lbfluid) {
  Utils::Vector3d linear_momentum{};
  if (include_particles) {
    Utils::Vector3d momentum_particles{};
    mpi_gather_stats(4, momentum_particles.data(), nullptr, nullptr, nullptr);
    linear_momentum += momentum_particles;
  }
  if (include_lbfluid) {
#if defined(LB) or defined(LB_GPU)
    linear_momentum += lb_lbfluid_calc_fluid_momentum();
#endif
  }
  return linear_momentum;
}

Utils::Vector3d centerofmass(PartCfg &partCfg, int type) {
  Utils::Vector3d com{};
  double mass = 0.0;

  for (auto const &p : partCfg) {
    if ((p.p.type == type) || (type == -1)) {
      for (int j = 0; j < 3; j++) {
        com[j] += p.r.p[j] * (p).p.mass;
      }
      mass += (p).p.mass;
    }
  }
  for (int j = 0; j < 3; j++)
    com[j] /= mass;
  return com;
}

void angularmomentum(PartCfg &partCfg, int type, double *com) {
  com[0] = com[1] = com[2] = 0.;

  for (auto const &p : partCfg) {
    if (type == p.p.type) {
      auto const tmp = vector_product(p.r.p, p.m.v);
      for (int i = 0; i < 3; i++) {
        com[i] += tmp[i] * p.p.mass;
      }
    }
  }
}

void momentofinertiamatrix(PartCfg &partCfg, int type, double *MofImatrix) {
  int i, count;
  double p1[3], massi;
  count = 0;

  for (i = 0; i < 9; i++)
    MofImatrix[i] = 0.;

  auto const com = centerofmass(partCfg, type);
  for (auto const &p : partCfg) {
    if (type == p.p.type) {
      count++;
      for (i = 0; i < 3; i++) {
        p1[i] = p.r.p[i] - com[i];
      }
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

IntList nbhood(PartCfg &partCfg, double pt[3], double r,
               int const planedims[3]) {
  IntList ids;
  Utils::Vector3d d;

  auto const r2 = r * r;

  for (auto const &p : partCfg) {
    if ((planedims[0] + planedims[1] + planedims[2]) == 3) {
      d = get_mi_vector(pt, p.r.p);
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

double distto(PartCfg &partCfg, double p[3], int pid) {
  auto mindist = std::numeric_limits<double>::infinity();

  for (auto const &part : partCfg) {
    if (pid != part.p.identity) {
      auto const d = get_mi_vector(p, part.r.p);
      mindist = std::min(mindist, d.norm2());
    }
  }
  return std::sqrt(mindist);
}

void calc_part_distribution(PartCfg &partCfg, int const *p1_types, int n_p1,
                            int const *p2_types, int n_p2, double r_min,
                            double r_max, int r_bins, int log_flag, double *low,
                            double *dist) {
  int t1, t2, ind, cnt = 0;
  double inv_bin_width = 0.0;
  double min_dist, min_dist2 = 0.0, start_dist2;

  start_dist2 = Utils::sqr(box_l[0] + box_l[1] + box_l[2]);
  /* bin preparation */
  *low = 0.0;
  for (int i = 0; i < r_bins; i++)
    dist[i] = 0.0;
  if (log_flag == 1)
    inv_bin_width = (double)r_bins / (log(r_max) - log(r_min));
  else
    inv_bin_width = (double)r_bins / (r_max - r_min);

  /* particle loop: p1_types*/
  for (auto const &p1 : partCfg) {
    for (t1 = 0; t1 < n_p1; t1++) {
      if (p1.p.type == p1_types[t1]) {
        min_dist2 = start_dist2;
        /* particle loop: p2_types*/
        for (auto const &p2 : partCfg) {
          if (p1 != p2) {
            for (t2 = 0; t2 < n_p2; t2++) {
              if (p2.p.type == p2_types[t2]) {
                auto const act_dist2 = get_mi_vector(p1.r.p, p2.r.p).norm2();
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
            if (log_flag == 1)
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
  int i, t1, t2, ind;
  int mixed_flag = 0;
  double inv_bin_width = 0.0, bin_width = 0.0;
  double volume, bin_volume, r_in, r_out;

  if (n_p1 == n_p2) {
    for (i = 0; i < n_p1; i++)
      if (p1_types[i] != p2_types[i])
        mixed_flag = 1;
  } else
    mixed_flag = 1;

  bin_width = (r_max - r_min) / (double)r_bins;
  inv_bin_width = 1.0 / bin_width;
  for (i = 0; i < r_bins; i++)
    rdf[i] = 0.0;
  /* particle loop: p1_types*/
  for (auto it = partCfg.begin(); it != partCfg.end(); ++it) {
    for (t1 = 0; t1 < n_p1; t1++) {
      if (it->p.type == p1_types[t1]) {
        /* distinguish mixed and identical rdf's */
        auto jt = (mixed_flag == 1) ? partCfg.begin() : std::next(it);

        /* particle loop: p2_types*/
        for (; jt != partCfg.end(); ++jt) {
          for (t2 = 0; t2 < n_p2; t2++) {
            if (jt->p.type == p2_types[t2]) {
              auto const dist = get_mi_vector(it->r.p, jt->r.p).norm();
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

  /* normalization */
  volume = box_l[0] * box_l[1] * box_l[2];
  for (i = 0; i < r_bins; i++) {
    r_in = i * bin_width + r_min;
    r_out = r_in + bin_width;
    bin_volume = (4.0 / 3.0) * Utils::pi() *
                 ((r_out * r_out * r_out) - (r_in * r_in * r_in));
    rdf[i] *= volume / (bin_volume * cnt);
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
  int mixed_flag = 0;
  double inv_bin_width = 0.0, bin_width = 0.0;
  double volume, bin_volume, r_in, r_out;
  double *rdf_tmp;

  rdf_tmp = (double *)Utils::malloc(r_bins * sizeof(double));

  if (n_p1 == n_p2) {
    for (int i = 0; i < n_p1; i++)
      if (p1_types[i] != p2_types[i])
        mixed_flag = 1;
  } else
    mixed_flag = 1;

  bin_width = (r_max - r_min) / (double)r_bins;
  inv_bin_width = 1.0 / bin_width;
  volume = box_l[0] * box_l[1] * box_l[2];
  for (int l = 0; l < r_bins; l++)
    rdf_tmp[l] = rdf[l] = 0.0;

  while (cnt_conf <= n_conf) {
    for (int l = 0; l < r_bins; l++)
      rdf_tmp[l] = 0.0;
    cnt = 0;
    auto const k = n_configs - cnt_conf;
    int i = 0;
    for (auto it = partCfg.begin(); it != partCfg.end(); ++it) {
      for (int t1 = 0; t1 < n_p1; t1++) {
        if (it->p.type == p1_types[t1]) {
          /* distinguish mixed and identical rdf's */
          auto jt = (mixed_flag == 1) ? partCfg.begin() : std::next(it);
          int j = (mixed_flag == 1) ? 0 : i + 1;

          // particle loop: p2_types
          for (; jt != partCfg.end(); ++jt) {
            for (int t2 = 0; t2 < n_p2; t2++) {
              if (jt->p.type == p2_types[t2]) {
                using Utils::make_const_span;

                auto const dist =
                    get_mi_vector(make_const_span(configs[k] + 3 * i, 3),
                                  make_const_span(configs[k] + 3 * j, 3))
                        .norm();
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
      r_in = i * bin_width + r_min;
      r_out = r_in + bin_width;
      bin_volume = (4.0 / 3.0) * Utils::pi() *
                   ((r_out * r_out * r_out) - (r_in * r_in * r_in));
      rdf[i] += rdf_tmp[i] * volume / (bin_volume * cnt);
    }

    cnt_conf++;
  } // cnt_conf loop
  for (int i = 0; i < r_bins; i++) {
    rdf[i] /= (cnt_conf - 1);
  }
  free(rdf_tmp);
}

void calc_structurefactor(PartCfg &partCfg, int const *p_types, int n_types,
                          int order, double **_ff) {
  int i, j, k, n, qi, t, order2;
  double qr, twoPI_L, C_sum, S_sum, *ff = nullptr;

  order2 = order * order;
  *_ff = ff = Utils::realloc(ff, 2 * order2 * sizeof(double));
  ff[2 * order2] = 0;
  twoPI_L = 2 * Utils::pi() / box_l[0];

  if ((n_types < 0) || (n_types > max_seen_particle_type)) {
    fprintf(stderr, "WARNING: Wrong number of particle types!");
    fflush(nullptr);
    errexit();
  } else if (order < 1) {
    fprintf(stderr,
            "WARNING: parameter \"order\" has to be a whole positive number");
    fflush(nullptr);
    errexit();
  } else {
    for (qi = 0; qi < 2 * order2; qi++) {
      ff[qi] = 0.0;
    }
    for (i = 0; i <= order; i++) {
      for (j = -order; j <= order; j++) {
        for (k = -order; k <= order; k++) {
          n = i * i + j * j + k * k;
          if ((n <= order2) && (n >= 1)) {
            C_sum = S_sum = 0.0;
            for (auto const &p : partCfg) {
              for (t = 0; t < n_types; t++) {
                if (p.p.type == p_types[t]) {
                  qr = twoPI_L * (i * p.r.p[0] + j * p.r.p[1] + k * p.r.p[2]);
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
    n = 0;
    for (auto const &p : partCfg) {
      for (t = 0; t < n_types; t++) {
        if (p.p.type == p_types[t])
          n++;
      }
    }
    for (qi = 0; qi < order2; qi++)
      if (ff[2 * qi + 1] != 0)
        ff[2 * qi] /= n * ff[2 * qi + 1];
  }
}

std::vector<std::vector<double>> modify_stucturefactor(int order,
                                                       double const *sf) {
  int length = 0;

  for (int i = 0; i < order * order; i++) {
    if (sf[2 * i + 1] > 0) {
      length++;
    }
  }

  double qfak = 2.0 * Utils::pi() / box_l[0];
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

int calc_cylindrical_average(
    PartCfg &partCfg, std::vector<double> const &center_,
    std::vector<double> const &direction_, double length, double radius,
    int bins_axial, int bins_radial, std::vector<int> types,
    std::map<std::string, std::vector<std::vector<std::vector<double>>>>
        &distribution) {
  int index_axial;
  int index_radial;
  double binwd_axial = length / bins_axial;
  double binwd_radial = radius / bins_radial;

  auto center = Utils::Vector3d{center_};
  auto direction = Utils::Vector3d{direction_};

  // Select all particle types if the only entry in types is -1
  bool all_types = false;
  if (types.size() == 1 && types[0] == -1)
    all_types = true;

  distribution.insert(
      std::pair<std::string, std::vector<std::vector<std::vector<double>>>>(
          "density",
          std::vector<std::vector<std::vector<double>>>(types.size())));
  distribution.insert(
      std::pair<std::string, std::vector<std::vector<std::vector<double>>>>(
          "v_r", std::vector<std::vector<std::vector<double>>>(types.size())));
  distribution.insert(
      std::pair<std::string, std::vector<std::vector<std::vector<double>>>>(
          "v_t", std::vector<std::vector<std::vector<double>>>(types.size())));

  for (unsigned int type = 0; type < types.size(); type++) {
    distribution["density"][type].resize(bins_radial);
    distribution["v_r"][type].resize(bins_radial);
    distribution["v_t"][type].resize(bins_radial);
    for (int index_radial = 0; index_radial < bins_radial; index_radial++) {
      distribution["density"][type][index_radial].assign(bins_axial, 0.0);
      distribution["v_r"][type][index_radial].assign(bins_axial, 0.0);
      distribution["v_t"][type][index_radial].assign(bins_axial, 0.0);
    }
  }

  auto const norm_direction = direction.norm();

  for (auto const &p : partCfg) {
    for (unsigned int type_id = 0; type_id < types.size(); type_id++) {
      if (types[type_id] == p.p.type || all_types) {
        auto const pos = folded_position(p);

        Utils::Vector3d vel{p.m.v};

        auto const diff = pos - center;

        // Find the height of the particle above the axis (height) and
        // the distance from the center point (dist)
        auto const hat = vector_product(direction, diff);
        auto const height = hat.norm();
        auto const dist = direction * diff / norm_direction;

        // Determine the components of the velocity parallel and
        // perpendicular to the direction vector
        double v_radial;
        if (height == 0)
          v_radial = vector_product(vel, direction).norm() / norm_direction;
        else
          v_radial = vel * hat / height;

        auto const v_axial = vel * direction / norm_direction;

        // Work out relevant indices for x and y
        index_radial = static_cast<int>(floor(height / binwd_radial));
        index_axial =
            static_cast<int>(floor((dist + 0.5 * length) / binwd_axial));

        if ((index_radial < bins_radial && index_radial >= 0) &&
            (index_axial < bins_axial && index_axial >= 0)) {
          distribution["density"][type_id][index_radial][index_axial] += 1;
          distribution["v_r"][type_id][index_radial][index_axial] += v_radial;
          distribution["v_t"][type_id][index_radial][index_axial] += v_axial;
        }
      }
    }
  }

  // Now we turn the counts into densities by dividing by one radial
  // bin (binvolume).  We also divide the velocities by the counts.
  double binvolume;
  for (unsigned int type_id = 0; type_id < types.size(); type_id++) {
    for (int index_radial = 0; index_radial < bins_radial; index_radial++) {
      // All bins are cylindrical shells of thickness binwd_radial.
      // The volume is thus: binvolume = pi*(r_outer - r_inner)^2 * length
      if (index_radial == 0)
        binvolume = M_PI * binwd_radial * binwd_radial * length;
      else
        binvolume = M_PI * (index_radial * index_radial + 2 * index_radial) *
                    binwd_radial * binwd_radial * length;
      for (int index_axial = 0; index_axial < bins_axial; index_axial++) {
        if (distribution["density"][type_id][index_radial][index_axial] != 0) {
          distribution["v_r"][type_id][index_radial][index_axial] /=
              distribution["density"][type_id][index_radial][index_axial];
          distribution["v_t"][type_id][index_radial][index_axial] /=
              distribution["density"][type_id][index_radial][index_axial];
          distribution["density"][type_id][index_radial][index_axial] /=
              binvolume;
        }
      }
    }
  }

  return ES_OK;
}

/****************************************************************************************
 *                                 config storage functions
 ****************************************************************************************/

void analyze_append(PartCfg &partCfg) {
  n_part_conf = partCfg.size();
  configs = Utils::realloc(configs, (n_configs + 1) * sizeof(double *));
  configs[n_configs] =
      (double *)Utils::malloc(3 * n_part_conf * sizeof(double));
  int i = 0;
  for (auto const &p : partCfg) {
    configs[n_configs][3 * i + 0] = p.r.p[0];
    configs[n_configs][3 * i + 1] = p.r.p[1];
    configs[n_configs][3 * i + 2] = p.r.p[2];
    i++;
  }
  n_configs++;
}

void analyze_configs(double const *tmp_config, int count) {
  int i;
  n_part_conf = count;
  configs = Utils::realloc(configs, (n_configs + 1) * sizeof(double *));
  configs[n_configs] =
      (double *)Utils::malloc(3 * n_part_conf * sizeof(double));
  for (i = 0; i < n_part_conf; i++) {
    configs[n_configs][3 * i] = tmp_config[3 * i];
    configs[n_configs][3 * i + 1] = tmp_config[3 * i + 1];
    configs[n_configs][3 * i + 2] = tmp_config[3 * i + 2];
  }
  n_configs++;
}

/****************************************************************************************
 *                                 Observables handling
 ****************************************************************************************/

void obsstat_realloc_and_clear(Observable_stat *stat, int n_pre, int n_bonded,
                               int n_non_bonded, int n_coulomb, int n_dipolar,
                               int n_vs, int c_size) {

  // Number of doubles to store pressure in
  const int total =
      c_size * (n_pre + bonded_ia_params.size() + n_non_bonded + n_coulomb +
                n_dipolar + n_vs + Observable_stat::n_external_field);

  // Allocate mem for the double list
  stat->data.resize(total);

  // Number of doubles per interaction (pressure=1, stress tensor=9,...)
  stat->chunk_size = c_size;

  // Number of chunks for different interaction types
  stat->n_coulomb = n_coulomb;
  stat->n_dipolar = n_dipolar;
  stat->n_non_bonded = n_non_bonded;
  stat->n_virtual_sites = n_vs;
  // Pointers to the start of different contributions
  stat->bonded = stat->data.e + c_size * n_pre;
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
  int i, total = c_size * (n_nonbonded + n_nonbonded);

  stat_nb->data_nb.resize(total);
  stat_nb->chunk_size_nb = c_size;
  stat_nb->n_nonbonded = n_nonbonded;
  stat_nb->non_bonded_intra = stat_nb->data_nb.e;
  stat_nb->non_bonded_inter = stat_nb->non_bonded_intra + c_size * n_nonbonded;

  for (i = 0; i < total; i++)
    stat_nb->data_nb[i] = 0.0;
}

void invalidate_obs() {
  total_energy.init_status = 0;
  total_pressure.init_status = 0;
  total_p_tensor.init_status = 0;
}

void update_pressure(int v_comp) {
  int i;
  double p_vel[3];
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
      total_pressure.data.e[0] = 0.0;
      MPI_Reduce(nptiso.p_vel, p_vel, 3, MPI_DOUBLE, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      for (i = 0; i < 3; i++)
        if (nptiso.geometry & nptiso.nptgeom_dir[i])
          total_pressure.data.e[0] += p_vel[i];
      total_pressure.data.e[0] /= (nptiso.dimension * nptiso.volume);
      total_pressure.init_status = 1 + v_comp;
    } else
      master_pressure_calc(v_comp);
  }
}
