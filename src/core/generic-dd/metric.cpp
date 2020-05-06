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
#include <algorithm>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/operations.hpp>
#include <chrono>
#include <functional>
#include <random>
#include <regex>
#include <stdexcept>

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cells.hpp"
#include "communication.hpp" // comm_cart
#include "metric.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "short_range_loop.hpp"

/** Fills weights with a constant.
 */
static void metric_ncells(std::vector<double> &weights) {
  std::fill(weights.begin(), weights.end(), 1.0);
}

/** Fills weights with the number of particles per cell.
 */
static void metric_npart(std::vector<double> &weights) {
  std::transform(cell_structure.m_local_cells.begin(),
                 cell_structure.m_local_cells.end(), weights.begin(),
                 [](const Cell *c) { return c->particles().size(); });
}

/** Returns the number of distance pairs calculated for cell "c".
 */
static int cell_ndistpairs(Cell *c) {
  int nnp = std::accumulate(c->m_neighbors.red().begin(),
                            c->m_neighbors.red().begin(), 0,
                            [](int acc, const Cell *neigh) {
                              return acc + neigh->particles().size();
                            });
  return c->particles().size() * nnp;
}

/** Fills weights with the number of distance pairs per cell.
 */
static void metric_ndistpairs(std::vector<double> &weights) {
  std::transform(cell_structure.m_local_cells.begin(),
                 cell_structure.m_local_cells.end(), weights.begin(),
                 cell_ndistpairs);
}

/** Returns the number of bond partners for all particles of cell "c".
 */
int cell_nbondedia(Cell *cell) {
  int nbondedia = 0;
  for (const Particle &p : cell->particles()) {
    for (int i = 0; i < p.bl.n;) {
      int type_num = p.bl.e[i++];
      Bonded_ia_parameters *iaparams = &bonded_ia_params[type_num];
      // int type = iaparams->type;

      // This could be incremented conditionally if "type" has a specific value
      // to only count bonded_ia of a certain type.
      nbondedia++;
      i += iaparams->num; // Skip the all partner particle numbers
    }
  }
  return nbondedia;
}

/** Fills weights with the number of bond partners for all particles in each
 * cell.
 */
static void metric_nbondedia(std::vector<double> &weights) {
  std::transform(cell_structure.m_local_cells.begin(),
                 cell_structure.m_local_cells.end(), weights.begin(),
                 cell_nbondedia);
}

namespace {
/** Generator for random integers
 */
struct Randintgen {
  Randintgen()
      : mt(std::chrono::high_resolution_clock::now()
               .time_since_epoch()
               .count()),
        d(1, 1000) {}
  Randintgen(Randintgen &&other)
      : mt(std::move(other.mt)), d(std::move(other.d)) {}
  int operator()() { return d(mt); }

private:
  std::mt19937 mt;
  std::uniform_int_distribution<int> d;
};
} // namespace

/** Fills weights with random integers between 1 and 1000.
 */
static void metric_rand(std::vector<double> &weights) {
  std::generate(weights.begin(), weights.end(), Randintgen{});
}

/** Returns 1.0
 */
static double cc_metric_uniform(int i, int j) {
  (void)i;
  (void)j;
  return 1.0;
}

/** Returns the product of the particle numbers of both cells.
 */
static double cc_metric_multiply_npart(int i, int j) {
  return cells[i].particles().size() * cells[j].particles().size();
}

/** Returns the sum of the particle number of both cells.
 */
static double cc_metric_add_npart(int i, int j) {
  return cells[i].particles().size() + cells[j].particles().size();
}

namespace generic_dd {
namespace repart {

/** Get the appropriate metric function described in string "desc".
 * Throws a std::invalid_argument exception if no such metric is available
 */
static Metric::metric_func get_single_metric_func(const std::string &desc) {
  using desc_func_pair = std::pair<std::string, Metric::metric_func>;
  static const std::vector<desc_func_pair> mets = {
      {"ncells", metric_ncells},
      {"npart", metric_npart},
      {"ndistpairs", metric_ndistpairs},
      {"nbondedia", metric_nbondedia},
      {"rand", metric_rand}};

  for (const auto &t : mets) {
    if (desc == std::get<0>(t)) {
      return std::get<1>(t);
    }
  }

  throw std::invalid_argument(std::string("No such metric available: ") + desc);
}

/** Get the appropriate cell pair metric function described in string "desc".
 * Throws a std::invalid_argument exception if no such metric is available
 */
static Metric::cc_metric_func
get_single_cc_metric_func(const std::string &desc) {
  using desc_func_pair = std::pair<std::string, Metric::cc_metric_func>;
  static const std::vector<desc_func_pair> mets = {
      {"uniform", cc_metric_uniform},
      {"add_npart", cc_metric_add_npart},
      {"multiply_npart", cc_metric_multiply_npart}};
  for (const auto &t : mets) {
    if (desc == std::get<0>(t)) {
      return std::get<1>(t);
    }
  }

  throw std::invalid_argument(std::string("No such metric available: ") + desc);
}

void Metric::set_metric(const std::string &desc) {
  mdesc.clear();
  ccmdesc.clear();
  parse_metric_desc(desc);
}

static std::string remove_whitespace(const std::string &s) {
  static const std::regex ws("\\s");
  return std::regex_replace(s, ws, "");
}

void Metric::parse_metric_desc(const std::string &desc) {
  std::string::size_type pos;
  // Cell cell metric denoted in square brackets at the end.
  if ((pos = desc.find('[')) != std::string::npos) {
    auto cdesc = desc.substr(0, pos - 1);
    auto ccdesc = desc.substr(pos + 1, desc.find(']') - pos - 1);
    // std::cout << "Parsing '" << cdesc << "' as metric and '" << ccdesc << "'
    // as cell cell metric" << std::endl;
    parse_cell_metric_desc(cdesc);
    parse_cell_cell_metric_desc(ccdesc);
  } else {
    parse_cell_metric_desc(desc);
  }
}

/** Parses a string representing a linear combination of words interpreted via
 * "parse_metric". The results are added to "mvec".
 */
template <typename MetricVec, typename StringParseFunc>
static void parse_linear_combination(const std::string &desc, MetricVec &mvec,
                                     StringParseFunc parse_metric) {
  static const std::regex term_re(
      "\\s*([\\+-]?\\s*\\d*(\\.\\d*)?)?\\s*\\*?\\s*(\\w+)");
  static const auto match_to_term = [parse_metric](const std::smatch &match) {
    // There is whitespace in the numbers (first group in term_re), so we have
    // to remove it before calling atof.
    std::string num = remove_whitespace(match.str(1));
    // Number might be omitted, i.e. "foo+bar". In this case 'num' only holds a
    // sign or nothing at all (for the very first term).
    if (num == "+" || num == "-" || num.empty())
      num += "1.0";
    try {
      return std::make_pair(stod(num), parse_metric(match.str(3)));
    } catch (...) {
      std::cout << "Error in metric parsing at: " << match.str(0) << std::endl;
      throw;
    }
  };

  auto t_begin =
      std::sregex_iterator(std::begin(desc), std::end(desc), term_re);
  auto t_end = decltype(t_begin){};

  std::transform(t_begin, t_end, std::back_inserter(mvec), match_to_term);
}

void Metric::parse_cell_cell_metric_desc(const std::string &desc) {
  parse_linear_combination(desc, ccmdesc, get_single_cc_metric_func);
}

void Metric::parse_cell_metric_desc(const std::string &desc) {
  parse_linear_combination(desc, mdesc, get_single_metric_func);
}

/** Evaluates the represented linear combination of metric functions and returns
 * the resulting weights.
 */
std::vector<double> Metric::operator()() const {
  std::vector<double> w(cell_structure.m_local_cells.size(), 0.0),
      tmp(cell_structure.m_local_cells.size());

  for (const auto &t : mdesc) {
    double factor = std::get<0>(t);
    metric_func func = std::get<1>(t);

    func(tmp);
    // Multiply add to w.
    std::transform(
        w.begin(), w.end(), tmp.begin(), w.begin(),
        [factor](double acc, double val) { return acc + factor * val; });
  }
  return w;
}

double Metric::curload() const {
  const auto ws = (*this)();
  return std::accumulate(std::begin(ws), std::end(ws), 0.0);
}

double Metric::paverage() const {
  const double l = curload();
  return boost::mpi::all_reduce(comm_cart, l, std::plus<void>{}) /
         comm_cart.size();
}

double Metric::pmax() const {
  const double l = curload();
  return boost::mpi::all_reduce(comm_cart, l, boost::mpi::maximum<double>{});
}

double Metric::pimbalance() const { return pmax() / paverage(); }

double Metric::cell_cell_weight(int i, int j) const {
  if (!has_cell_cell_metric()) {
    return 1.0;
  }

  double res = 0.0;
  for (const auto &t : ccmdesc) {
    double factor = std::get<0>(t);
    cc_metric_func func = std::get<1>(t);
    res += factor * func(i, j);
  }

  return res;
}

} // namespace repart
} // namespace generic_dd
