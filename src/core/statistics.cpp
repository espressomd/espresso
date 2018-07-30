/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file statistics.cpp
    This is the place for analysis (so far...).
    Implementation of statistics.hpp
*/
#include "statistics.hpp"
#include "communication.hpp"
#include "energy.hpp"
#include "grid.hpp"
#include "initialize.hpp"
#include "interaction_data.hpp"
#include "lb.hpp"
#include "npt.hpp"
#include "partCfg_global.hpp"
#include "particle_data.hpp"
#include "pressure.hpp"
#include "short_range_loop.hpp"
#include "statistics_chain.hpp"
#include "statistics_cluster.hpp"
#include "statistics_fluid.hpp"
#include "utils.hpp"
#include "utils/NoOp.hpp"
#include "utils/list_contains.hpp"
#include "virtual_sites.hpp"

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
  double pt[3];
  int in_set;

  auto mindist2 = std::numeric_limits<double>::infinity();

  for (auto jt = partCfg.begin(); jt != (--partCfg.end()); ++jt) {
    pt[0] = jt->r.p[0];
    pt[1] = jt->r.p[1];
    pt[2] = jt->r.p[2];
    /* check which sets particle j belongs to
       bit 0: set1, bit1: set2
    */
    in_set = 0;
    if (set1.empty() || list_contains(set1, jt->p.type))
      in_set = 1;
    if (set2.empty() || list_contains(set2, jt->p.type))
      in_set |= 2;
    if (in_set == 0)
      continue;

    for (auto it = std::next(jt); it != partCfg.end(); ++it)
      /* accept a pair if particle j is in set1 and particle i in set2 or vice
       * versa. */
      if (((in_set & 1) && (set2.empty() || list_contains(set2, it->p.type))) ||
          ((in_set & 2) && (set1.empty() || list_contains(set1, it->p.type))))
        mindist2 = std::min(mindist2, min_distance2(pt, it->r.p));
  }

  return std::sqrt(mindist2);
}

void merge_aggregate_lists(int *head_list, int *agg_id_list, int p1molid,
                           int p2molid, int *link_list) {
  int target1, target2, head_p1;
  /* merge list containing p2molid into list containing p1molid*/
  target1 = head_list[agg_id_list[p2molid]];
  head_list[agg_id_list[p2molid]] = -2;
  head_p1 = head_list[agg_id_list[p1molid]];
  head_list[agg_id_list[p1molid]] = target1;
  agg_id_list[target1] = agg_id_list[p1molid];
  target2 = link_list[target1];
  while (target2 != -1) {
    target1 = target2;
    target2 = link_list[target1];
    agg_id_list[target1] = agg_id_list[p1molid];
  }
  agg_id_list[target1] = agg_id_list[p1molid];
  link_list[target1] = head_p1;
}

int aggregation(double dist_criteria2, int min_contact, int s_mol_id,
                int f_mol_id, int *head_list, int *link_list, int *agg_id_list,
                int *agg_num, int *agg_size, int *agg_max, int *agg_min,
                int *agg_avg, int *agg_std, int charge) {
  int target1;
  int *contact_num, ind;

  if (min_contact > 1) {
    contact_num =
        (int *)Utils::malloc(topology.size() * topology.size() * sizeof(int));
    for (int i = 0; i < topology.size() * topology.size(); i++)
      contact_num[i] = 0;
  } else {
    contact_num = (int *)0; /* Just to keep the compiler happy */
  }

  on_observable_calc();

  for (int i = s_mol_id; i <= f_mol_id; i++) {
    head_list[i] = i;
    link_list[i] = -1;
    agg_id_list[i] = i;
    agg_size[i] = 0;
  }

  short_range_loop(Utils::NoOp{},
                   [&](Particle &p1, Particle &p2, Distance &d) {
                     auto p1molid = p1.p.mol_id;
                     auto p2molid = p2.p.mol_id;
                     if (((p1molid <= f_mol_id) && (p1molid >= s_mol_id)) &&
                         ((p2molid <= f_mol_id) && (p2molid >= s_mol_id))) {
                       if (agg_id_list[p1molid] != agg_id_list[p2molid]) {
#ifdef ELECTROSTATICS
                         if (charge && (p1.p.q * p2.p.q >= 0)) {
                           return;
                         }
#endif
                         if (d.dist2 < dist_criteria2) {
                           if (p1molid > p2molid) {
                             ind = p1molid * topology.size() + p2molid;
                           } else {
                             ind = p2molid * topology.size() + p1molid;
                           }
                           if (min_contact > 1) {
                             contact_num[ind]++;
                             if (contact_num[ind] >= min_contact) {
                               merge_aggregate_lists(head_list, agg_id_list,
                                                     p1molid, p2molid,
                                                     link_list);
                             }
                           } else {
                             merge_aggregate_lists(head_list, agg_id_list,
                                                   p1molid, p2molid, link_list);
                           }
                         }
                       }
                     }
                   });

  /* count number of aggregates
     find aggregate size
     find max and find min size, and std */
  for (int i = s_mol_id; i <= f_mol_id; i++) {
    if (head_list[i] != -2) {
      (*agg_num)++;
      agg_size[*agg_num - 1]++;
      target1 = head_list[i];
      while (link_list[target1] != -1) {
        target1 = link_list[target1];
        agg_size[*agg_num - 1]++;
      }
    }
  }

  for (int i = 0; i < *agg_num; i++) {
    *agg_avg += agg_size[i];
    *agg_std += agg_size[i] * agg_size[i];
    if (*agg_min > agg_size[i]) {
      *agg_min = agg_size[i];
    }
    if (*agg_max < agg_size[i]) {
      *agg_max = agg_size[i];
    }
  }

  return 0;
}

/** Calculate momentum of all particles in the local domain
 * @param result Result for this processor (Output)
 */
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

/** Calculate total momentum of the system (particles & LB fluid)
 * inputs are bools to include particles and fluid in the linear momentum
 * calculation
 * @param momentum Result for this processor (Output)
 */
std::vector<double> calc_linear_momentum(int include_particles,
                                         int include_lbfluid) {
  double momentum_particles[3] = {0., 0., 0.};
  std::vector<double> linear_momentum(3, 0.0);
  if (include_particles) {
    mpi_gather_stats(4, momentum_particles, nullptr, nullptr, nullptr);
    linear_momentum[0] += momentum_particles[0];
    linear_momentum[1] += momentum_particles[1];
    linear_momentum[2] += momentum_particles[2];
  }
  if (include_lbfluid) {
    double momentum_fluid[3] = {0., 0., 0.};
#ifdef LB
    if (lattice_switch & LATTICE_LB) {
      mpi_gather_stats(6, momentum_fluid, nullptr, nullptr, nullptr);
    }
#endif
#ifdef LB_GPU
    if (lattice_switch & LATTICE_LB_GPU) {
      lb_calc_fluid_momentum_GPU(momentum_fluid);
    }
#endif
    linear_momentum[0] += momentum_fluid[0];
    linear_momentum[1] += momentum_fluid[1];
    linear_momentum[2] += momentum_fluid[2];
  }
  return linear_momentum;
}

std::vector<double> centerofmass(PartCfg &partCfg, int type) {
  std::vector<double> com(3);
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

std::vector<double> centerofmass_vel(PartCfg &partCfg, int type) {
  /*center of mass velocity scaled with time_step*/
  std::vector<double> com_vel(3);
  int count = 0;

  for (auto const &p : partCfg) {
    if (type == p.p.type) {
      for (int i = 0; i < 3; i++) {
        com_vel[i] += p.m.v[i];
      }
      count++;
    }
  }

  for (int i = 0; i < 3; i++) {
    com_vel[i] /= count;
  }
  return com_vel;
}

void angularmomentum(PartCfg &partCfg, int type, double *com) {
  double tmp[3];
  com[0] = com[1] = com[2] = 0.;

  for (auto const &p : partCfg) {
    if (type == p.p.type) {
      vector_product(p.r.p, p.m.v, tmp);
      for (int i = 0; i < 3; i++) {
        com[i] += tmp[i] * p.p.mass;
      }
    }
  }
  return;
}

void momentofinertiamatrix(PartCfg &partCfg, int type, double *MofImatrix) {
  int i, count;
  double p1[3], massi;
  std::vector<double> com(3);
  count = 0;

  for (i = 0; i < 9; i++)
    MofImatrix[i] = 0.;
  com = centerofmass(partCfg, type);
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

IntList nbhood(PartCfg &partCfg, double pt[3], double r, int planedims[3]) {
  IntList ids;
  Vector3d d;

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

void calc_part_distribution(PartCfg &partCfg, int *p1_types, int n_p1,
                            int *p2_types, int n_p2, double r_min, double r_max,
                            int r_bins, int log_flag, double *low,
                            double *dist) {
  int t1, t2, ind, cnt = 0;
  double inv_bin_width = 0.0;
  double min_dist, min_dist2 = 0.0, start_dist2, act_dist2;

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
                act_dist2 = min_distance2(p1.r.p, p2.r.p);
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

void calc_rdf(PartCfg &partCfg, std::vector<int> &p1_types,
              std::vector<int> &p2_types, double r_min, double r_max,
              int r_bins, std::vector<double> &rdf) {
  calc_rdf(partCfg, &p1_types[0], p1_types.size(), &p2_types[0],
           p2_types.size(), r_min, r_max, r_bins, &rdf[0]);
}

void calc_rdf(PartCfg &partCfg, int *p1_types, int n_p1, int *p2_types,
              int n_p2, double r_min, double r_max, int r_bins, double *rdf) {
  long int cnt = 0;
  int i, t1, t2, ind;
  int mixed_flag = 0;
  double inv_bin_width = 0.0, bin_width = 0.0, dist;
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
              dist = min_distance(it->r.p, jt->r.p);
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
    bin_volume =
        (4.0 / 3.0) * PI * ((r_out * r_out * r_out) - (r_in * r_in * r_in));
    rdf[i] *= volume / (bin_volume * cnt);
  }
}

void calc_rdf_av(PartCfg &partCfg, std::vector<int> &p1_types,
                 std::vector<int> &p2_types, double r_min, double r_max,
                 int r_bins, std::vector<double> &rdf, int n_conf) {
  calc_rdf_av(partCfg, &p1_types[0], p1_types.size(), &p2_types[0],
              p2_types.size(), r_min, r_max, r_bins, &rdf[0], n_conf);
}

void calc_rdf_av(PartCfg &partCfg, int *p1_types, int n_p1, int *p2_types,
                 int n_p2, double r_min, double r_max, int r_bins, double *rdf,
                 int n_conf) {
  long int cnt = 0;
  int cnt_conf = 1;
  int mixed_flag = 0;
  double inv_bin_width = 0.0, bin_width = 0.0;
  double volume, bin_volume, r_in, r_out;
  double *rdf_tmp, p1[3], p2[3];

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
                p1[0] = configs[k][3 * i + 0];
                p1[1] = configs[k][3 * i + 1];
                p1[2] = configs[k][3 * i + 2];
                p2[0] = configs[k][3 * j + 0];
                p2[1] = configs[k][3 * j + 1];
                p2[2] = configs[k][3 * j + 2];
                auto const dist = min_distance(p1, p2);
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
      bin_volume =
          (4.0 / 3.0) * PI * ((r_out * r_out * r_out) - (r_in * r_in * r_in));
      rdf[i] += rdf_tmp[i] * volume / (bin_volume * cnt);
    }

    cnt_conf++;
  } // cnt_conf loop
  for (int i = 0; i < r_bins; i++) {
    rdf[i] /= (cnt_conf - 1);
  }
  free(rdf_tmp);
}

void calc_structurefactor(PartCfg &partCfg, int *p_types, int n_types,
                          int order, double **_ff) {
  int i, j, k, n, qi, t, order2;
  double qr, twoPI_L, C_sum, S_sum, *ff = nullptr;

  order2 = order * order;
  *_ff = ff = Utils::realloc(ff, 2 * order2 * sizeof(double));
  ff[2 * order2] = 0;
  twoPI_L = 2 * PI / box_l[0];

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

std::vector<std::vector<double>> modify_stucturefactor(int order, double *sf) {
  int length = 0;

  for (int i = 0; i < order * order; i++) {
    if (sf[2 * i + 1] > 0) {
      length++;
    }
  }

  double qfak = 2.0 * PI / box_l[0];
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

// calculates average density profile in dir direction over last n_conf
// configurations
void density_profile_av(PartCfg &partCfg, int n_conf, int n_bin, double density,
                        int dir, double *rho_ave, int type) {
  int i, j, k, m, n;
  double r;
  double r_bin;
  double pos[3];
  int image_box[3];

  // calculation over last n_conf configurations

  // bin width
  r_bin = box_l[dir] / (double)(n_bin);

  for (i = 0; i < n_bin; i++)
    rho_ave[i] = 0;

  k = n_configs - n_conf;

  while (k < n_configs) {
    r = 0;
    j = 0;
    while (r < box_l[dir]) {
      n = 0;
      for (auto const &p : partCfg) {
        // com particles
        if (p.p.type == type) {
          for (m = 0; m < 3; m++) {
            pos[m] = configs[k][3 * i + m];
            image_box[m] = 0;
          }
          fold_coordinate(pos, image_box, dir);
          if (pos[dir] <= r + r_bin && pos[dir] > r)
            n++;
        }
      }

      rho_ave[j] += (double)(n) / (box_l[1] * box_l[2] * r_bin) / density;
      j++;
      r += r_bin;
    }
    k++;
  } // k loop

  // normalization
  for (i = 0; i < n_bin; i++)
    rho_ave[i] /= n_conf;
}

int calc_cylindrical_average(
    PartCfg &partCfg, std::vector<double> center_,
    std::vector<double> direction_, double length, double radius,
    int bins_axial, int bins_radial, std::vector<int> types,
    std::map<std::string, std::vector<std::vector<std::vector<double>>>>
        &distribution) {
  int index_axial;
  int index_radial;
  double binwd_axial = length / bins_axial;
  double binwd_radial = radius / bins_radial;

  auto center = Vector3d{std::move(center_)};
  auto direction = Vector3d{std::move(direction_)};

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

        Vector3d vel{p.m.v};

        auto const diff = pos - center;

        // Find the height of the particle above the axis (height) and
        // the distance from the center point (dist)
        auto const hat = direction.cross(diff);
        auto const height = hat.norm();
        auto const dist = direction.dot(diff) / norm_direction;

        // Determine the components of the velocity parallel and
        // perpendicular to the direction vector
        double v_radial;
        if (height == 0)
          v_radial = vel.cross(direction).norm() / norm_direction;
        else
          v_radial = vel.dot(hat) / height;

        auto const v_axial = vel.dot(direction) / norm_direction;

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
  // bin (binvolume).  We also divide the velocites by the counts.
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

int calc_radial_density_map(PartCfg &partCfg, int xbins, int ybins,
                            int thetabins, double xrange, double yrange,
                            double axis[3], double center[3], IntList *beadids,
                            DoubleList *density_map,
                            DoubleList *density_profile) {
  int i, j, t;
  int bi;
  int nbeadtypes;
  int beadcount;
  double vectprod[3];
  double pvector[3];
  double xdist, ydist, rdist, xav, yav, theta;
  double xbinwidth, ybinwidth, binvolume;
  double thetabinwidth;
  double *thetaradii;
  int *thetacounts;
  int xindex, yindex, tindex;
  xbinwidth = xrange / (double)(xbins);
  ybinwidth = yrange / (double)(ybins);

  nbeadtypes = beadids->n;

  beadcount = 0;
  xav = 0.0;
  yav = 0.0;

  for (auto const &pi : partCfg) {
    for (bi = 0; bi < nbeadtypes; bi++) {
      if (beadids->e[bi] == pi.p.type) {
        /* Find the vector from the point to the center */
        vecsub(center, folded_position(pi), pvector);

        /* Work out x and y coordinates with respect to rotation axis */

        /* Find the minimum distance of the point from the axis */
        vector_product(axis, pvector, vectprod);
        xdist = sqrt(sqrlen(vectprod) / sqrlen(axis));

        /* Find the projection of the vector from the point to the center
           onto the axis vector */
        ydist = scalar(axis, pvector) / sqrt(sqrlen(axis));

        /* Work out relevant indices for x and y */
        xindex = (int)(floor(xdist / xbinwidth));
        yindex = (int)(floor((ydist + yrange * 0.5) / ybinwidth));

        /* Check array bounds */
        if ((xindex < xbins && xindex > 0) && (yindex < ybins && yindex > 0)) {
          density_map[bi].e[ybins * xindex + yindex] += 1;
          xav += xdist;
          yav += ydist;
          beadcount += 1;
        }
      }
    }
  }

  /* Now turn counts into densities for the density map */
  for (bi = 0; bi < nbeadtypes; bi++) {
    for (i = 0; i < xbins; i++) {
      /* All bins are cylinders and therefore constant in yindex */
      binvolume = PI * (2 * i * xbinwidth + xbinwidth * xbinwidth) * yrange;
      for (j = 0; j < ybins; j++) {
        density_map[bi].e[ybins * i + j] /= binvolume;
      }
    }
  }

  /* if required calculate the theta density profile */
  if (thetabins > 0) {
    /* Convert the center to an output of the density center */
    xav = xav / (double)(beadcount);
    yav = yav / (double)(beadcount);
    thetabinwidth = 2 * PI / (double)(thetabins);
    thetaradii =
        (double *)Utils::malloc(thetabins * nbeadtypes * sizeof(double));
    thetacounts = (int *)Utils::malloc(thetabins * nbeadtypes * sizeof(int));
    for (bi = 0; bi < nbeadtypes; bi++) {
      for (t = 0; t < thetabins; t++) {
        thetaradii[bi * thetabins + t] = 0.0;
        thetacounts[bi * thetabins + t] = 0.0;
      }
    }
    /* Maybe there is a nicer way to do this but now I will just repeat the loop
     * over all particles */
    for (auto const &pi : partCfg) {
      for (bi = 0; bi < nbeadtypes; bi++) {
        if (beadids->e[bi] == pi.p.type) {
          vecsub(center, folded_position(pi), pvector);
          vector_product(axis, pvector, vectprod);
          xdist = sqrt(sqrlen(vectprod) / sqrlen(axis));
          ydist = scalar(axis, pvector) / sqrt(sqrlen(axis));
          /* Center the coordinates */

          xdist = xdist - xav;
          ydist = ydist - yav;
          rdist = sqrt(xdist * xdist + ydist * ydist);
          if (ydist >= 0) {
            theta = acos(xdist / rdist);
          } else {
            theta = 2 * PI - acos(xdist / rdist);
          }
          tindex = (int)(floor(theta / thetabinwidth));
          thetaradii[bi * thetabins + tindex] += xdist + xav;
          thetacounts[bi * thetabins + tindex] += 1;
          if (tindex >= thetabins) {
            fprintf(stderr, "ERROR: outside density_profile array bounds in "
                            "calc_radial_density_map");
            fflush(nullptr);
            errexit();
          } else {
            density_profile[bi].e[tindex] += 1;
          }
        }
      }
    }

    /* normalize the theta densities*/
    for (bi = 0; bi < nbeadtypes; bi++) {
      for (t = 0; t < thetabins; t++) {
        rdist = thetaradii[bi * thetabins + t] /
                (double)(thetacounts[bi * thetabins + t]);
        density_profile[bi].e[t] /= rdist * rdist;
      }
    }

    free(thetaradii);
    free(thetacounts);
  }

  return ES_OK;
}

int calc_vanhove(PartCfg &partCfg, int ptype, double rmin, double rmax,
                 int rbins, int tmax, double *msd, double **vanhove) {
  int c1, c3, c3_max, ind;
  double p1[3], p2[3], dist;
  double bin_width, inv_bin_width;
  std::vector<int> ids;

  for (auto const &p : partCfg) {
    if (p.p.type == ptype) {
      ids.push_back(p.p.identity);
    }
  }

  if (ids.empty()) {
    return 0;
  }

  /* preparation */
  bin_width = (rmax - rmin) / (double)rbins;
  inv_bin_width = 1.0 / bin_width;

  /* calculate msd and store distribution in vanhove */
  for (c1 = 0; c1 < n_configs; c1++) {
    c3_max = (c1 + tmax + 1) > n_configs ? n_configs : c1 + tmax + 1;
    for (c3 = (c1 + 1); c3 < c3_max; c3++) {
      for (auto const &id : ids) {
        p1[0] = configs[c1][3 * id];
        p1[1] = configs[c1][3 * id + 1];
        p1[2] = configs[c1][3 * id + 2];
        p2[0] = configs[c3][3 * id];
        p2[1] = configs[c3][3 * id + 1];
        p2[2] = configs[c3][3 * id + 2];
        dist = distance(p1, p2);
        if (dist > rmin && dist < rmax) {
          ind = (int)((dist - rmin) * inv_bin_width);
          vanhove[(c3 - c1 - 1)][ind]++;
        }
        msd[(c3 - c1 - 1)] += dist * dist;
      }
    }
  }

  /* normalize */
  for (c1 = 0; c1 < (tmax); c1++) {
    for (int i = 0; i < rbins; i++) {
      vanhove[c1][i] /= (double)(n_configs - c1 - 1) * ids.size();
    }
    msd[c1] /= (double)(n_configs - c1 - 1) * ids.size();
  }

  return ids.size();
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

void analyze_push(PartCfg &partCfg) {
  n_part_conf = partCfg.size();
  free(configs[0]);
  for (int i = 0; i < n_configs - 1; i++) {
    configs[i] = configs[i + 1];
  }
  configs[n_configs - 1] =
      (double *)Utils::malloc(3 * n_part_conf * sizeof(double));

  int i = 0;
  for (auto const &p : partCfg) {
    configs[n_configs - 1][3 * i + 0] = p.r.p[0];
    configs[n_configs - 1][3 * i + 1] = p.r.p[1];
    configs[n_configs - 1][3 * i + 2] = p.r.p[2];

    i++;
  }
}

void analyze_replace(PartCfg &partCfg, int ind) {
  n_part_conf = partCfg.size();

  int i = 0;
  for (auto const &p : partCfg) {
    configs[ind][3 * i + 0] = p.r.p[0];
    configs[ind][3 * i + 1] = p.r.p[1];
    configs[ind][3 * i + 2] = p.r.p[2];

    i++;
  }
}

void analyze_remove(int ind) {
  int i;
  free(configs[ind]);
  for (i = ind; i < n_configs - 1; i++) {
    configs[i] = configs[i + 1];
  }
  n_configs--;
  configs = Utils::realloc(configs, n_configs * sizeof(double *));
  if (n_configs == 0)
    n_part_conf = 0;
}

void analyze_configs(double *tmp_config, int count) {
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

void analyze_activate(PartCfg &partCfg, int ind) {
  int i;
  double pos[3];
  n_part_conf = partCfg.size();

  for (i = 0; i < n_part_conf; i++) {
    pos[0] = configs[ind][3 * i];
    pos[1] = configs[ind][3 * i + 1];
    pos[2] = configs[ind][3 * i + 2];
    if (place_particle(i, pos) == ES_ERROR) {
      runtimeErrorMsg() << "failed upon replacing particle " << i
                        << "  in Espresso";
    }
  }
}

/****************************************************************************************
 *                                 Observables handling
 ****************************************************************************************/

void obsstat_realloc_and_clear(Observable_stat *stat, int n_pre, int n_bonded,
                               int n_non_bonded, int n_coulomb, int n_dipolar,
                               int n_vs, int c_size) {

  int i;
  // Number of doubles to store pressure in
  int total = c_size * (n_pre + bonded_ia_params.size() + n_non_bonded + n_coulomb +
                        n_dipolar + n_vs);

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

  // Set all obseravables to zero
  for (i = 0; i < total; i++)
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
