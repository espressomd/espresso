/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include <fcs.h>
#include <list>
#include <mpi.h>
#include <string>
#include <vector>

namespace Scafacos {

/** \brief Abstraction of a method from the scafacos library */

struct Scafacos {
  Scafacos(std::string method, MPI_Comm comm, const std::string &parameters);
  ~Scafacos();
  /** Parse parameter string */
  void parse_parameters(const std::string &s);
  /** Get the parameters from the library */
  std::string get_parameters();
  /** Get active method name */
  std::string get_method();

  /** Set parameters common to all methods */
  void set_common_parameters(const double *box_l, const int *periodicity,
                             int total_particles);
  /** Calculate short range pair force if supported by the method */
  inline double pair_force(double dist) const {
    if (has_near) {
      fcs_float field;
      fcs_compute_near_field(handle, dist, &field);
      return field;
    }

    return 0.0;
  }

  /** Get pair energy for near field contrib */
  inline double pair_energy(double dist) const {
    if (has_near) {
      fcs_float potential;
      fcs_compute_near_potential(handle, dist, &potential);
      return potential;
    }

    return 0.0;
  }

  /** Calculate the fields and potentials for charges */
  void run(std::vector<double> &charges, std::vector<double> &positions,
           std::vector<double> &fields, std::vector<double> &potentials);

/** Calculate fields and potentials for dipoles */
#ifdef FCS_ENABLE_DIPOLES
  void run_dipolar(std::vector<double> &dipoles, std::vector<double> &positions,
                   std::vector<double> &fields,
                   std::vector<double> &potentials);

#endif
  /** Tune parameters */
  void tune(std::vector<double> &charges, std::vector<double> &positions);
  /** Get shortrange cutoff (0.0 if not supported) */
  double r_cut() const;
  /** Set cutoff */
  void set_r_cut(double r_cut);

  /** Get a list of methods supported by the library */
  static std::list<std::string> available_methods();

  /** Scafacos used for dipolar ia */
  bool dipolar() { return m_dipolar; }

  /** Switch scafacos to dipolar ia */
  void set_dipolar(bool d);

  /** Handle from the library */
  FCS handle;
  /** Whether the method supports near field delegation */
  bool has_near;
  /** The scafacos method name of this instance */
  const std::string method;
  /** The last parameters set */
  std::string m_last_parameters;

  /** Scafacos used for dipolar ia */
  bool m_dipolar;
};

} // namespace Scafacos
