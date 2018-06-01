/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project

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
/* * @File Header file for the correlation class
 *
 * This module allow to compute correlations (and other two time averages) on
 * the fly an from files.
 *
 * The basic idea is that the user can write arbitrary function A and B that
 * can depend on e.g. particle coordinates or whatever state of the MD box.
 * These function can be vector valued, always indicated by dim_A and dim_B.
 *
 * The way they can be correlated is can be formulated as a (vector valued)
 * operation on A and B. One example would be to calculate the product
 * component by component, and if A and B are both the particle velocities,
 * then one would obtain
 * { <v_1x v_1x(t)>  <v_1y v_1y(t)>  <v_1z v_1z(t)>  <v_2x v_2x(t)>  <v_2y
 * v_2y(t)>  <v_2z v_2z(t)> ... }
 *
 * The idea of the algorithm is not to keep all As and Bs in memory, but
 * feed the algorithm (or the container object) successively by new As and Bs
 * and new correlation estimates are calculated on the fly. Those As and Bs
 * which have exceeded a certain age are (first compressed, see below) and then
 * discarded.
 *
 * To save memory, increase statistics and make the calculation possible for
 * many orders of magnitude in time, the blocking algorithm in Frenkels book
 * (p.91)
 * is applied. Thus not all As and Bs of the whole "past" are stored but
 * some of them are blocked. In this implementation always a blocking based
 * on 2 is applied: All As and Bs not older than a certain tau_lin are stored
 * as they were, those which are older are not entire stored, but only
 * their compressed average value. All As and Bs older than 2*tau_lin
 * are compressed in blocks of four etc.
 *
 * This leads to a hiearchical "history": on level 0 the last tau_lin values
 * are stored. This is done in a cyclic array: The newest is always appended
 * at the end of the array, and if the array length is reached, values
 * are appended at the beginning, overwriting older values. We therefore
 * have to carry the index newest[0] which tells, what is the index of the
 * newest value of A and B.
 *
 * As soon as the first row of As and Bs is full, two of them are
 * compressed by calculating the arithmetic mean, and stored in the
 * first element on level 1. Then we can overwrite these two values
 * on level 0, because we have saved them. Always if necessary
 * we produce space on level 0 by compressing to level 1. On
 * level 1 also an array with tau_lin entries is available to store
 * the level-1-compresssed values. It is sucessivly filled
 * and also cyclic. When it is filled, the values are stored
 * on level 2 and so on.
 *
 * This allows to have a "history" over many orders of magnitude
 * in time, without the full memory effort.
 *
 * Correlations are only calculated on each level: For
 * tau=1,2,..,tau_lin the values are taken from level 1.
 * For tau=tau_lin, tau_lin+2, .., 2*tau_lin we take the values
 * from level 2. On level 2 we halso have values for 0,..tau_lin-2,
 * but these are discared as we have "better" estimates on level 1.
 *
 * The functions A and B can get extra arguments. This can e.g. be an
 * integer describing to the "type" of the particles to be considered,
 * or in case of files, it is a file_data_source pointer, which tells
 * the function from where to read. These arguments have to be a single
 * pointer to an object that carries all relevant information
 * to obtain A and B (from the MD Box or somewhere else).
 *
 * The correlation has to be initialized with all necessary information,
 * i.e. all function pointers, the dimensions of A and B and their dimensions
 * etc.
 * When new A and B values should be processed, a single call of
 * double_correlation_get_data(..) with the correlation object as an
 * argument is enough to update A and B and the correlation estimate.
 *
 * Eventually the correlation can be printed.
 *
 *
 * There is a lot of stuff to do:
 * ==============================
 * -> Write input functions that take TCL arrays for A and B to
 *  make the method available for TCL-coded variables
 * -> Expand the file_data_source so that one can specify which
 *  columns of the file are to be processed
 * -> Write more correlation operations (scalar product)
 * -> Write more observable
 * -> calculate an estimate of average values. This might be
 *  even necessary to calculate <(A-<A>)(B(tau)-<B>), which
 *  is often probably what people want
 * -> Use the A_args to calculate As and Bs only for particular
 *  particle types (especially and example, so that other people can follow)
 * -> Use the A_args to calculate molecular stuff in combination with
 *  the topology concept
 * -> Tidy up the parsers (just a bit)
 * -> Write some nice samples
 * -> Write a destructor
 * -> Finally: write the users guide
 *
 *
 */
#ifndef _STATISTICS_CORRELATION_H
#define _STATISTICS_CORRELATION_H

#include "Vector.hpp"
#include "observables/Observable.hpp"
#include "utils.hpp"

#include <boost/multi_array.hpp>
#include <boost/serialization/access.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>

namespace Correlators {

// IDs of different correlations

void autoupdate_correlations();

/** The struct that is used to calculate correlations *
 *
 * Data organization:
 * We use a ring-like way to manage the data: at the beginning we have a linear
 * array
 * which we fill from index 0 to tau_lin. The index newest[i] always indicates
 * the latest
 * entry of the hierarchic "past" For every new entry in is incremented and if
 * tau_lin is reached,
 * it starts again from the beginning.
 */
class Correlator {
  using obs_ptr = std::shared_ptr<Observables::Observable>;

public:
  /**
   * The initialization procedure for the correlation object. All important
   * parameters have to be speciefied
   * at the same time. They can not be change later, so every instance of the
   * correlation class
   * has to be fed with correct data from the very beginning.
   *
   * @param _dt: The time interval between subsequent updates
   * @param _tau_lin: The linear part of the correlation function.
   * @param _tau_max: maximal time delay tau to sample
   * @param _window_distance: The distance in time domain between update of the
   * correlation estimate
   * @param _dim_A: The dimension of the A vector
   * @param _dim_B: The dimension of the B vector
   * @param _A: First observable to correlate
   * @param _B: Second observable to correlate
   * @param _dim_corr: The dimension of the correlation to be calculated
   * @param _corr_operation_name: how to correlate the two observables A and B
   *     (this has no default)
   * @param _compressA_name: how the A values should be compressed (usually
   *     the linear compression method)
   * @param _compressB_name: how the B values should be compressed (usually
   *     the linear compression method)
   * @param _args: the parameters of the observables
   *
   */
  Correlator(int tau_lin, double tau_max, double dt,
             std::string const &compress1_, std::string const &compress2_,
             std::string const &corr_operation, const obs_ptr &obs1, const obs_ptr &obs2)
      : autoupdate(0), finalized(0), t(0), m_tau_lin(tau_lin),
        m_dt(dt), m_tau_max(tau_max), compressA_name(compress1_),
        compressB_name(compress2_), corr_operation_name(corr_operation),
        A_obs(obs1), B_obs(obs2) {
    initialize();
  }

private:
  void initialize();

public:
  /** The function to process a new datapoint of A and B
   *
   * First the function finds out if it necessary to make some space for the new
   * entries of A and B.
   * Then, if necessary, it compresses old Values of A and B to make for the new
   * value. Finally
   * The new values of A and B are stored in A[newest[0]] and B[newest[0]],
   * where the newest indices
   * have been increased before. Finally the correlation estimate is updated.
   * TODO: Not all
   * the correlation estimates have to be updated.
   *
   */
  int get_data();

  /** At the end of data collection, go through the whole hierarchy and
   * correlate data left there
   *
   * This works pretty much the same as get_data, but does not feed on new data,
   * just uses what
   * is already available.
   *
   */
  int finalize();

  /** Return an estimate of the integrated correlation time
   *
   *  We calculate the correlation time for each dim_corr by normalizing the
   * correlation,
   * integrating it and finding out where C(tau)=tau;
   *
   */

  int get_correlation_time(double *correlation_time);

  /** The function to process a new datapoint of A and B
   *
   * First the function finds out if it necessary to make some space for the new
   * entries of A and B.
   * Then, if necessary, it compresses old Values of A and B to make for the new
   * value. Finally
   * The new values of A and B are stored in A[newest[0]] and B[newest[0]],
   * where the newest indices
   * have been increased before. Finally the correlation estimate is updated.
   * TODO: Not all
   * the correlation estimates have to be updated.
   *
   */

  /** Return correlation result */
  std::vector<double> get_correlation();

  int tau_lin() const { return m_tau_lin; }
  double tau_max() const { return m_tau_max; }
  double last_update() const { return m_last_update; }
  double dt() const { return m_dt; }
  int dim_corr() const { return m_dim_corr; }
  int n_result() const { return m_n_result; }

  Vector3d const &correlation_args() const { return m_correlation_args; }
  void set_correlation_args(Vector3d const &args) { m_correlation_args = args; }

  std::string const &compress1() const { return compressA_name; }
  std::string const &compress2() const { return compressB_name; }
  std::string const &correlation_operation() const {
    return corr_operation_name;
  }

  void start_auto_update();
  void stop_auto_update();
  int autoupdate;

  /* Partial serialization of state that is not accessible
     via the interface. */
  std::string get_internal_state() const;
  void set_internal_state(std::string const &);

private:
  unsigned int finalized; // non-zero of correlation is finialized
  unsigned int t;         // global time in number of frames

  Vector3d
      m_correlation_args; // additional arguments, which the correlation may
                          // need (currently only used by fcs_acf)

  int hierarchy_depth; // maximum level of data compression
  int m_tau_lin;       // number of frames in the linear correlation
  int m_dim_corr;
  double m_dt; // time interval at which samples arrive
  double
      m_tau_max; // maximum time, for which the correlation should be calculated
  int update_frequency; // time distance between updates in MD timesteps

  std::string compressA_name;
  std::string compressB_name;
  std::string corr_operation_name;

  int m_n_result; // the total number of result values

  std::shared_ptr<Observables::Observable> A_obs;
  std::shared_ptr<Observables::Observable> B_obs;

  std::vector<int> tau; // time differences
  boost::multi_array<std::vector<double>, 2> A;
  boost::multi_array<std::vector<double>, 2> B;

  boost::multi_array<double, 2> result; // output quantity

  // The actual allocated storage space
  std::vector<unsigned int>
      n_sweeps; // number of correlation sweeps at a particular value of tau
  std::vector<unsigned int>
      n_vals; // number of data values already present at a particular
              // value of tau
  std::vector<unsigned int>
      newest; // index of the newest entry in each hierarchy level

  std::vector<double> A_accumulated_average; // all A values are added up here
  std::vector<double> B_accumulated_average; // all B values are added up here
  unsigned int n_data; // a counter to calculated averages and variances

  double m_last_update;

  unsigned int dim_A; // dimensionality of A
  unsigned int dim_B;

  using correlation_operation_type = std::vector<double> (*)(
      std::vector<double> const &, std::vector<double> const &, Vector3d);

  correlation_operation_type corr_operation;

  using compression_function = std::vector<double> (*)(
      std::vector<double> const &A1, std::vector<double> const &A2);

  // compressing functions
  compression_function compressA;
  compression_function compressB;
};

extern int correlations_autoupdate;

} // Namespace correlators
#endif
