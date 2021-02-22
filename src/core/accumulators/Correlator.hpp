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
/** @file
 *
 * This module computes correlations (and other two time averages) on
 * the fly and from files.
 *
 * The basic idea is that the user can write arbitrary functions A and B that
 * can depend on e.g. particle coordinates or whatever state of the MD box.
 * These functions can be vector-valued, always indicated by dim_A and dim_B.
 *
 * The way they can be correlated is can be formulated as a (vector-valued)
 * operation on A and B. One example would be to calculate the product
 * component by component, and if A and B are both the particle velocities,
 * then one would obtain
 * @f[ <v_1^x v_1^x(t)>, <v_1^y v_1^y(t)>, <v_1^z v_1^z(t)>,
 *     <v_2^x v_2^x(t)>, <v_2^y v_2^y(t)>, <v_2^z v_2^z(t)>, \ldots @f]
 *
 * The idea of the algorithm is to not keep all As and Bs in memory, but
 * feed the algorithm (or the container object) successively by new As and Bs
 * and new correlation estimates are calculated on the fly. Those As and Bs
 * which have exceeded a certain age are first compressed (see below) and then
 * discarded.
 *
 * To save memory, increase statistics and make the calculation possible for
 * many orders of magnitude in time, the order-<i>n</i> blocking algorithm in
 * @cite frenkel02b (algorithm 9, section 4.4.2 p. 95)
 * is used. Thus not all As and Bs of the whole "past" are stored but
 * some of them are blocked. In this implementation, a blocking based on 2 is
 * always applied: all As and Bs not older than a certain @c tau_lin are stored
 * as they were, those which are older are stored as average values. All As
 * and Bs older than <tt>2*tau_lin</tt> are compressed in blocks of four, etc.
 *
 * This leads to a hierarchical "history": on level 0 the last @c tau_lin values
 * are stored. This is done in a cyclic array: the newest is always appended
 * at the end of the array, and if the array length is reached, values
 * are appended at the beginning, overwriting older values. We therefore
 * have to carry the index <tt>newest[0]</tt> which tells, what is the index
 * of the newest value of A and B.
 *
 * As soon as the first row of As and Bs is full, two of them are
 * compressed by calculating the arithmetic mean, and stored in the
 * first element on level 1. Then we can overwrite these two values
 * on level 0, because we have saved them. If necessary,
 * we produce space on level 0 by compressing to level 1. On
 * level 1 an array with @c tau_lin entries is available to store
 * the level-1-compressed values. It is successively filled
 * and also cyclic. When it is filled, the values are stored
 * on level 2 and so on.
 *
 * This allows to have a "history" over many orders of magnitude
 * in time, without the full memory effort.
 *
 * Correlations are only calculated on each level. For
 * <tt>tau=1,2,..,tau_lin</tt> the values are taken from level 1.
 * For <tt>tau=tau_lin, tau_lin+2, .., 2*tau_lin</tt> we take the values
 * from level 2. On level 2 we also have values for <tt>0,..tau_lin-2</tt>,
 * but these are discarded as we have "better" estimates on level 1.
 *
 * The functions A and B can get extra arguments. This can e.g. be an
 * integer describing the "type" of the particles to be considered,
 * or in case of files, it is a file_data_source pointer, which tells
 * the function from where to read. These arguments have to be a single
 * pointer to an object that carries all relevant information
 * to obtain A and B (from the MD Box or somewhere else).
 *
 * The correlation has to be initialized with all necessary information, i.e.
 * all function pointers, the dimensions of A and B and their dimensions, etc.
 *
 *
 * TODO: There is a lot of stuff to do:
 * - Expand the file_data_source so that one can specify which
 *  columns of the file are to be processed
 * - calculate an estimate of average values. This might be
 *  even necessary to calculate
 *  @f$\left<\left(A(\tau)-\left<A\right>\right)\left(B(\tau)-\left<B\right>\right)\right>@f$,
 *  which is often probably what people want
 * - Use the A_args to calculate As and Bs only for particular
 *  particle types (especially and example, so that other people can follow)
 * - Use the A_args to calculate molecular stuff in combination with
 *  the topology concept
 * - Write a destructor
 */
#ifndef _STATISTICS_CORRELATION_H
#define _STATISTICS_CORRELATION_H

#include "AccumulatorBase.hpp"
#include "integrate.hpp"
#include "observables/Observable.hpp"

#include <utils/Vector.hpp>

#include <boost/multi_array.hpp>
#include <boost/serialization/access.hpp>

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace Accumulators {

/** The main correlator class
 *
 *  Data organization:
 *  We use a ring-like way to manage the data: at the beginning we have a
 *  linear array, which we fill from index 0 to @c tau_lin. The index
 *  <tt>newest[i]</tt> always indicates the latest entry of the hierarchic
 *  "past" For every new entry in is incremented and if @c tau_lin is reached,
 *  it starts again from the beginning.
 */
class Correlator : public AccumulatorBase {
  using obs_ptr = std::shared_ptr<Observables::Observable>;

public:
  /** The initialization procedure for the correlation object. All important
   *  parameters have to be specified at the same time. They cannot be changed
   *  later, so every instance of the correlation class has to be fed with
   *  correct data from the very beginning.
   *
   *  @param delta_N The number of time steps between subsequent updates
   *  @param tau_lin The linear part of the correlation function.
   *  @param tau_max maximal time delay tau to sample
   *  @param obs1 First observable to correlate
   *  @param obs2 Second observable to correlate
   *  @param corr_operation how to correlate the two observables A and B
   *      (this has no default)
   *  @param compress1_ how the A values should be compressed (usually
   *      the linear compression method)
   *  @param compress2_ how the B values should be compressed (usually
   *      the linear compression method)
   *  @param correlation_args_ optional arguments for the correlation function
   *      (currently only used when @p corr_operation is "fcs_acf")
   *
   */
  Correlator(int tau_lin, double tau_max, int delta_N, std::string compress1_,
             std::string compress2_, std::string corr_operation, obs_ptr obs1,
             obs_ptr obs2, Utils::Vector3d correlation_args_ = {})
      : AccumulatorBase(delta_N), finalized(false), t(0),
        m_correlation_args(correlation_args_), m_tau_lin(tau_lin),
        m_dt(delta_N * time_step), m_tau_max(tau_max),
        compressA_name(std::move(compress1_)),
        compressB_name(std::move(compress2_)),
        corr_operation_name(std::move(corr_operation)), A_obs(std::move(obs1)),
        B_obs(std::move(obs2)) {
    initialize();
  }

private:
  void initialize();

public:
  /** The function to process a new datapoint of A and B
   *
   *  First the function finds out if it is necessary to make some space for
   *  the new entries of A and B. Then, if necessary, it compresses old values
   *  of A and B to make room for the new value. Finally, the new values of A
   *  and B are stored in <tt>A[newest[0]]</tt> and <tt>B[newest[0]]</tt>,
   *  where the <tt>newest</tt> indices have been increased before. Finally,
   *  the correlation estimate is updated.
   *  TODO: Not all correlation estimates have to be updated.
   */
  void update() override;

  /** At the end of data collection, go through the whole hierarchy and
   *  correlate data left there.
   */
  int finalize();

  /** Return correlation result */
  std::vector<double> get_correlation();
  size_t n_values() const {
    return m_tau_lin + 1 + (m_tau_lin + 1) / 2 * (m_hierarchy_depth - 1);
  }
  std::vector<size_t> shape() const override {
    std::vector<size_t> shape = m_shape;
    shape.insert(shape.begin(), n_values());
    return shape;
  }
  std::vector<int> get_samples_sizes() const {
    return std::vector<int>(n_sweeps.begin(), n_sweeps.end());
  }
  std::vector<double> get_lag_times() const;

  int tau_lin() const { return m_tau_lin; }
  double tau_max() const { return m_tau_max; }
  double dt() const { return m_dt; }

  Utils::Vector3d const &correlation_args() const { return m_correlation_args; }
  void set_correlation_args(Utils::Vector3d const &args) {
    m_correlation_args = args;
  }

  std::string const &compress1() const { return compressA_name; }
  std::string const &compress2() const { return compressB_name; }
  std::string const &correlation_operation() const {
    return corr_operation_name;
  }

  /** Partial serialization of state that is not accessible via the interface.
   */
  std::string get_internal_state() const;
  void set_internal_state(std::string const &);

private:
  bool finalized; ///< whether the correlation is finalized
  unsigned int t; ///< global time in number of frames

  Utils::Vector3d m_correlation_args; ///< additional arguments, which the
                                      ///< correlation may need (currently
                                      ///< only used by fcs_acf)

  int m_hierarchy_depth; ///< maximum level of data compression
  int m_tau_lin;         ///< number of frames in the linear correlation
  size_t m_dim_corr;     ///< number of columns for the correlation
  double m_dt;           ///< time interval at which samples arrive
  double m_tau_max;      ///< maximum time for which the correlation should be
                         ///< calculated

  std::string compressA_name;
  std::string compressB_name;
  std::string corr_operation_name; ///< Name of the correlation operator

  std::shared_ptr<Observables::Observable> A_obs;
  std::shared_ptr<Observables::Observable> B_obs;

  std::vector<int> tau; ///< time differences
  boost::multi_array<std::vector<double>, 2> A;
  boost::multi_array<std::vector<double>, 2> B;

  boost::multi_array<double, 2> result; ///< output quantity

  /// number of correlation sweeps at a particular value of tau
  std::vector<size_t> n_sweeps;
  /// number of data values already present at a particular value of tau
  std::vector<unsigned> n_vals;
  /// index of the newest entry in each hierarchy level
  std::vector<size_t> newest;

  std::vector<double> A_accumulated_average; ///< all A values are added up here
  std::vector<double> B_accumulated_average; ///< all B values are added up here
  size_t n_data; ///< a counter for calculated averages and variances

  size_t dim_A;                ///< dimensionality of A
  size_t dim_B;                ///< dimensionality of B
  std::vector<size_t> m_shape; ///< dimensionality of the correlation

  using correlation_operation_type = std::vector<double> (*)(
      std::vector<double> const &, std::vector<double> const &,
      Utils::Vector3d const &);

  correlation_operation_type corr_operation;

  using compression_function = std::vector<double> (*)(
      std::vector<double> const &A1, std::vector<double> const &A2);

  // compressing functions
  compression_function compressA;
  compression_function compressB;
};

} // namespace Accumulators
#endif
