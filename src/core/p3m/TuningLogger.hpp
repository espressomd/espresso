/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#ifndef ESPRESSO_SRC_CORE_P3M_TUNING_LOGGER_HPP
#define ESPRESSO_SRC_CORE_P3M_TUNING_LOGGER_HPP

#include "config/config.hpp"

#if defined(P3M) || defined(DP3M)

#include <utils/Vector.hpp>

#include <cstdio>
#include <string>

class TuningLogger {
public:
  enum class Mode { Coulomb, Dipolar };
  TuningLogger(bool verbose, std::string name, Mode mode)
      : m_verbose{verbose}, m_name{std::move(name)}, m_mode{mode} {}

  void log_tuning_start() const {
    if (m_verbose) {
      std::printf("mesh cao r_cut_iL    alpha_L     err       "
                  "rs_err    ks_err    time [ms]\n");
    }
  }

  template <typename... Types>
  void log_success(double time, Types... parameter_set) const {
    if (m_verbose) {
      row(parameter_set...);
      std::printf(" %-8.2f\n", time);
    }
  }

  template <typename... Types>
  void log_skip(std::string reason, Types... parameter_set) const {
    if (m_verbose) {
      row(parameter_set...);
      std::printf(" %s\n", reason.c_str());
    }
  }

  void log_cao_too_large(int mesh, int cao) const {
    if (m_verbose) {
      std::printf("%-4d %-3d cao too large for this mesh\n", mesh, cao);
    }
  }

  void tuning_goals(double accuracy, double prefactor, double box_l,
                    int n_particles, double sum_prop) const {
    if (m_verbose) {
      std::string particle_trait;
      std::string particle_property;
      switch (m_mode) {
      case Mode::Coulomb:
        particle_trait = "charged";
        particle_property = "Sum[q_i^2]";
        break;
      case Mode::Dipolar:
        particle_trait = "magnetic";
        particle_property = "Sum[mu_i^2]";
        break;
      }
      std::printf("%s tune parameters: Accuracy goal = %.5e prefactor = %.5e\n"
                  "System: box_l = %.5e # %s part = %d %s = %.5e\n",
                  m_name.c_str(), accuracy, prefactor, box_l,
                  particle_trait.c_str(), n_particles,
                  particle_property.c_str(), sum_prop);
    }
  }

  void tuning_results(Utils::Vector3i const &mesh, int cao, double r_cut_iL,
                      double alpha_L, double accuracy, double time) const {
    if (m_verbose) {
      std::printf(
          "\nresulting parameters: mesh: (%d, %d, %d), cao: %d, r_cut_iL: %.4e,"
          "\n                      alpha_L: %.4e, accuracy: %.4e, time: %.2f\n",
          mesh[0], mesh[1], mesh[2], cao, r_cut_iL, alpha_L, accuracy, time);
    }
  }

  void report_fixed_cao(int cao) const {
    if (m_verbose) {
      std::printf("fixed cao %d\n", cao);
    }
  }

  void report_fixed_r_cut_iL(double r_cut_iL) const {
    if (m_verbose) {
      std::printf("fixed r_cut_iL %f\n", r_cut_iL);
    }
  }

  void report_fixed_mesh(Utils::Vector3i const &mesh) const {
    if (m_verbose) {
      std::printf("fixed mesh (%d, %d, %d)\n", mesh[0], mesh[1], mesh[2]);
    }
  }

  auto get_name() const { return m_name; }

private:
  bool m_verbose;
  std::string m_name;
  Mode m_mode;

  void row(int mesh, int cao, double r_cut_iL, double alpha_L, double accuracy,
           double rs_err, double ks_err) const {
    std::printf("%-4d %-3d %.5e %.5e %.3e %.3e %.3e", mesh, cao, r_cut_iL,
                alpha_L, accuracy, rs_err, ks_err);
  }
};

#endif // P3M or DP3M

#endif
