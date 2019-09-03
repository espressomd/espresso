/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
#ifndef REPART_HPP_INCLUDED
#define REPART_HPP_INCLUDED

#include <functional>
#include <string>
#include <vector>

namespace generic_dd {
namespace repart {

/** Represents a linear combination of possible weighting functions
 * for cells and cell pairs.
 * For cells, the weights are supposed to indicate some sort of load associated
 * with the cell. For cell pairs, the weights should indicate some sort of cost
 * that results from placing each of the two cells on a different process.
 */
struct Metric {
  using metric_func = std::function<void(std::vector<double> &)>;
  using cc_metric_func = std::function<double(int, int)>;

  Metric() = default;
  Metric(const std::string &desc) { set_metric(desc); }

  /** Set "desc" as metric. Might throw a std::invalid_argument exception if
   * desc is not understood. Metric description strings are linear combinations
   * of single metrics.
   * E.g. "2.0*ncells +1.7*nghostpart"
   * The space after the metric name ("ncell") is mandatory.
   * Factor, multiplication and addition sign are mandatory. Negative constants
   * are only allowed for the first factor. For further use subtraction instead
   * of addition, e.g. "-1.0*ncells -1.7*npart".
   * Single metric names are also acceptable and interpreted as "1.0<name>".
   * Valid metrics are: ncells, npart, ndistpairs, bondedia, and rand.
   * @param desc string to describe the metric
   */
  void set_metric(const std::string &desc);

  /** Returns cell weights.
   * @return vector of weights. One double for each local cell.
   */
  std::vector<double> operator()() const;

  /** Returns true is this metric is set to calculate cell-cell weights
   * (e.g. transfer counts).
   */
  bool has_cell_cell_metric() const { return !ccmdesc.empty(); }

  /** Returns cell-cell metric weights if has_cell_cell_metric() == true. Else
   * returns 1.0 always.
   * @param i valid cell index
   * @param j valid cell index
   */
  double cell_cell_weight(int i, int j) const;

  /** Returns the load accumulated over all cells for this process.
   */
  double curload() const;

  /** Returns the average of curload() over all processes. Needs to be called
   * by all processes.
   */
  double paverage() const;

  /** Returns the maximum of curload() over all processes. Needs to be called
   * by all processes.
   */
  double pmax() const;

  /** Returns the imbalance of curload() over all processes. Needs to be called
   * by all processes.
   */
  double pimbalance() const;

private:
  void parse_metric_desc(const std::string &desc);
  void parse_cell_metric_desc(const std::string &desc);
  void parse_cell_cell_metric_desc(const std::string &desc);

  std::vector<std::pair<double, metric_func>> mdesc;
  std::vector<std::pair<double, cc_metric_func>> ccmdesc;
};

} // namespace repart
} // namespace generic_dd

#endif
