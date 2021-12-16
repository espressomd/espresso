//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \\file LeesEdwards_Collision.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "LeesEdwards_Collision.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#define FUNC_PREFIX

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning push
#pragma warning(disable : 1599)
#endif

using namespace std;

namespace walberla {
namespace pystencils {

namespace internal_leesedwards_collision_leesedwards_collision {
static FUNC_PREFIX void leesedwards_collision_leesedwards_collision(
    double *RESTRICT const _data_force, double *RESTRICT _data_pdfs,
    double *RESTRICT const _data_velocity, int64_t const _size_force_0,
    int64_t const _size_force_1, int64_t const _size_force_2,
    int64_t const _stride_force_0, int64_t const _stride_force_1,
    int64_t const _stride_force_2, int64_t const _stride_force_3,
    int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1,
    int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3,
    int64_t const _stride_velocity_0, int64_t const _stride_velocity_1,
    int64_t const _stride_velocity_2, int64_t const _stride_velocity_3,
    double omega) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    double *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    double *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    double *RESTRICT _data_velocity_20_30 =
        _data_velocity + _stride_velocity_2 * ctr_2;
    double *RESTRICT _data_velocity_20_31 =
        _data_velocity + _stride_velocity_2 * ctr_2 + _stride_velocity_3;
    double *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    double *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    double *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    double *RESTRICT _data_velocity_20_32 =
        _data_velocity + _stride_velocity_2 * ctr_2 + 2 * _stride_velocity_3;
    double *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      double *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      double *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      double *RESTRICT _data_velocity_20_30_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_30;
      double *RESTRICT _data_velocity_20_31_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_31;
      double *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      double *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      double *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      double *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      double *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      double *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      double *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      double *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      double *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      double *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      double *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      double *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      double *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      double *RESTRICT _data_velocity_20_32_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_32;
      double *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      double *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      double *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      double *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      double *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      double *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      double *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      for (int64_t ctr_0 = 0; ctr_0 < _size_force_0; ctr_0 += 1) {
        const double xi_1 = _data_force_20_30_10[_stride_force_0 * ctr_0];
        const double xi_2 = _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0];
        const double xi_3 = _data_velocity_20_30_10[_stride_velocity_0 * ctr_0];
        const double xi_4 = _data_velocity_20_31_10[_stride_velocity_0 * ctr_0];
        const double xi_5 = _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0];
        const double xi_6 = _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0];
        const double xi_7 = _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0];
        const double xi_8 = _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0];
        const double xi_9 = _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const double xi_10 = _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0];
        const double xi_11 = _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0];
        const double xi_12 = _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0];
        const double xi_13 = _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0];
        const double xi_14 = _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0];
        const double xi_15 = _data_force_20_31_10[_stride_force_0 * ctr_0];
        const double xi_16 = _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0];
        const double xi_17 = _data_force_20_32_10[_stride_force_0 * ctr_0];
        const double xi_18 =
            _data_velocity_20_32_10[_stride_velocity_0 * ctr_0];
        const double xi_19 = _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0];
        const double xi_20 = _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0];
        const double xi_21 = _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0];
        const double xi_22 = _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0];
        const double xi_23 = _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0];
        const double xi_24 = _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0];
        const double xi_25 = _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0];
        const double rho = xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_16 +
                           xi_19 + xi_2 + xi_20 + xi_21 + xi_22 + xi_23 +
                           xi_24 + xi_25 + xi_5 + xi_6 + xi_7 + xi_8 + xi_9;
        const double u1Pu2 = xi_18 + xi_4;
        const double u1Mu2 = -xi_18 + xi_4;
        _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0] =
            omega * (rho * (xi_18 * xi_18) * -0.333333333333333 +
                     rho * (xi_3 * xi_3) * -0.333333333333333 +
                     rho * (xi_4 * xi_4) * -0.333333333333333 +
                     rho * 0.333333333333333 - xi_9) +
            xi_9 +
            (omega * -0.5 + 1.0) *
                (-xi_1 * xi_3 - xi_15 * xi_4 - xi_17 * xi_18);
        _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0] =
            omega * (rho * (xi_18 * xi_18) * -0.166666666666667 +
                     rho * (xi_3 * xi_3) * -0.166666666666667 +
                     rho * (xi_4 * xi_4) * -0.166666666666667 +
                     rho * ((xi_4 * xi_4) * 0.333333333333333 +
                            xi_4 * 0.166666666666667 - 0.111111111111111) +
                     rho * 0.166666666666667 - xi_20) +
            xi_20 +
            (omega * -0.5 + 1.0) *
                (xi_1 * xi_3 * -0.166666666666667 +
                 xi_15 * (xi_4 * 2.0 + 1.0) * 0.166666666666667 +
                 xi_17 * xi_18 * -0.166666666666667);
        _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0] =
            omega * (rho * (xi_18 * xi_18) * -0.166666666666667 +
                     rho * (xi_3 * xi_3) * -0.166666666666667 +
                     rho * (xi_4 * xi_4) * -0.166666666666667 +
                     rho * ((xi_4 * xi_4) * 0.333333333333333 +
                            xi_4 * -0.166666666666667 - 0.111111111111111) +
                     rho * 0.166666666666667 - xi_25) +
            xi_25 +
            (omega * -0.5 + 1.0) *
                (xi_1 * xi_3 * -0.166666666666667 +
                 xi_15 * (xi_4 * 2.0 - 1.0) * 0.166666666666667 +
                 xi_17 * xi_18 * -0.166666666666667);
        _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0] =
            omega * (rho * (xi_18 * xi_18) * -0.166666666666667 +
                     rho * (xi_3 * xi_3) * -0.166666666666667 +
                     rho * (xi_4 * xi_4) * -0.166666666666667 +
                     rho * ((xi_3 * xi_3) * 0.333333333333333 +
                            xi_3 * -0.166666666666667 - 0.111111111111111) +
                     rho * 0.166666666666667 - xi_22) +
            xi_22 +
            (omega * -0.5 + 1.0) *
                (xi_1 * (xi_3 * 2.0 - 1.0) * 0.166666666666667 +
                 xi_15 * xi_4 * -0.166666666666667 +
                 xi_17 * xi_18 * -0.166666666666667);
        _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0] =
            omega * (rho * (xi_18 * xi_18) * -0.166666666666667 +
                     rho * (xi_3 * xi_3) * -0.166666666666667 +
                     rho * (xi_4 * xi_4) * -0.166666666666667 +
                     rho * ((xi_3 * xi_3) * 0.333333333333333 +
                            xi_3 * 0.166666666666667 - 0.111111111111111) +
                     rho * 0.166666666666667 - xi_24) +
            xi_24 +
            (omega * -0.5 + 1.0) *
                (xi_1 * (xi_3 * 2.0 + 1.0) * 0.166666666666667 +
                 xi_15 * xi_4 * -0.166666666666667 +
                 xi_17 * xi_18 * -0.166666666666667);
        _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0] =
            omega * (rho * (xi_18 * xi_18) * -0.166666666666667 +
                     rho * (xi_4 * xi_4) * -0.166666666666667 +
                     rho *
                         ((xi_3 + ((ctr_1 <= 0) ? (1.0) : (0.0)) *
                                      0.10000000000000001) *
                          (xi_3 + ((ctr_1 <= 0) ? (1.0) : (0.0)) *
                                      0.10000000000000001)) *
                         -0.166666666666667 +
                     rho * ((xi_18 * xi_18) * 0.333333333333333 +
                            xi_18 * 0.166666666666667 - 0.111111111111111) +
                     rho * 0.166666666666667 - xi_21) +
            xi_21 +
            (omega * -0.5 + 1.0) *
                (xi_1 *
                     (xi_3 +
                      ((ctr_1 <= 0) ? (1.0) : (0.0)) * 0.10000000000000001) *
                     -0.166666666666667 +
                 xi_15 * xi_4 * -0.166666666666667 +
                 xi_17 * (xi_18 * 2.0 + 1.0) * 0.166666666666667);
        _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0] =
            omega * (rho * (xi_18 * xi_18) * -0.166666666666667 +
                     rho * (xi_4 * xi_4) * -0.166666666666667 +
                     rho *
                         ((xi_3 + ((ctr_1 >= 63) ? (-1.0) : (0.0)) *
                                      0.10000000000000001) *
                          (xi_3 + ((ctr_1 >= 63) ? (-1.0) : (0.0)) *
                                      0.10000000000000001)) *
                         -0.166666666666667 +
                     rho * ((xi_18 * xi_18) * 0.333333333333333 +
                            xi_18 * -0.166666666666667 - 0.111111111111111) +
                     rho * 0.166666666666667 - xi_8) +
            xi_8 +
            (omega * -0.5 + 1.0) *
                (xi_1 *
                     (xi_3 +
                      ((ctr_1 >= 63) ? (-1.0) : (0.0)) * 0.10000000000000001) *
                     -0.166666666666667 +
                 xi_15 * xi_4 * -0.166666666666667 +
                 xi_17 * (xi_18 * 2.0 - 1.0) * 0.166666666666667);
        _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0] =
            omega *
                (rho * (xi_18 * xi_18) * -0.0416666666666667 +
                 rho * (xi_3 * xi_3) * -0.0416666666666667 +
                 rho * (xi_4 * xi_4) * -0.0416666666666667 +
                 rho * ((xi_18 * xi_18) * 0.0416666666666667 +
                        xi_3 * -0.0833333333333333 + xi_4 * 0.0833333333333333 +
                        ((xi_3 - xi_4) * (xi_3 - xi_4)) * 0.125 -
                        0.0138888888888889) +
                 rho * 0.0416666666666667 - xi_16) +
            xi_16 +
            (omega * -0.5 + 1.0) *
                (xi_1 * (xi_3 * -2.0 + xi_4 * 3.0 + 1.0) * -0.0833333333333333 +
                 xi_15 * (xi_3 * -3.0 + xi_4 * 2.0 + 1.0) * 0.0833333333333333 +
                 xi_17 * xi_18 * -0.0833333333333333);
        _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0] =
            omega *
                (rho * (xi_18 * xi_18) * -0.0416666666666667 +
                 rho * (xi_3 * xi_3) * -0.0416666666666667 +
                 rho * (xi_4 * xi_4) * -0.0416666666666667 +
                 rho * ((xi_18 * xi_18) * 0.0416666666666667 +
                        xi_3 * 0.0833333333333333 + xi_4 * 0.0833333333333333 +
                        ((xi_3 + xi_4) * (xi_3 + xi_4)) * 0.125 -
                        0.0138888888888889) +
                 rho * 0.0416666666666667 - xi_11) +
            xi_11 +
            (omega * -0.5 + 1.0) *
                (xi_1 * (xi_3 * 2.0 + xi_4 * 3.0 + 1.0) * 0.0833333333333333 +
                 xi_15 * (xi_3 * 3.0 + xi_4 * 2.0 + 1.0) * 0.0833333333333333 +
                 xi_17 * xi_18 * -0.0833333333333333);
        _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0] =
            omega * (rho * (xi_18 * xi_18) * -0.0416666666666667 +
                     rho * (xi_3 * xi_3) * -0.0416666666666667 +
                     rho * (xi_4 * xi_4) * -0.0416666666666667 +
                     rho * ((xi_18 * xi_18) * 0.0416666666666667 +
                            xi_3 * -0.0833333333333333 +
                            xi_4 * -0.0833333333333333 +
                            ((xi_3 + xi_4) * (xi_3 + xi_4)) * 0.125 -
                            0.0138888888888889) +
                     rho * 0.0416666666666667 - xi_14) +
            xi_14 +
            (omega * -0.5 + 1.0) *
                (xi_1 * (xi_3 * 2.0 + xi_4 * 3.0 - 1.0) * 0.0833333333333333 +
                 xi_15 * (xi_3 * 3.0 + xi_4 * 2.0 - 1.0) * 0.0833333333333333 +
                 xi_17 * xi_18 * -0.0833333333333333);
        _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0] =
            omega *
                (rho * (xi_18 * xi_18) * -0.0416666666666667 +
                 rho * (xi_3 * xi_3) * -0.0416666666666667 +
                 rho * (xi_4 * xi_4) * -0.0416666666666667 +
                 rho * ((xi_18 * xi_18) * 0.0416666666666667 +
                        xi_3 * 0.0833333333333333 + xi_4 * -0.0833333333333333 +
                        ((xi_3 - xi_4) * (xi_3 - xi_4)) * 0.125 -
                        0.0138888888888889) +
                 rho * 0.0416666666666667 - xi_13) +
            xi_13 +
            (omega * -0.5 + 1.0) *
                (xi_1 * (xi_3 * 2.0 + xi_4 * -3.0 + 1.0) * 0.0833333333333333 +
                 xi_15 * (xi_3 * 3.0 + xi_4 * -2.0 + 1.0) *
                     -0.0833333333333333 +
                 xi_17 * xi_18 * -0.0833333333333333);
        _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0] =
            omega *
                (rho * (xi_18 * xi_18) * -0.0416666666666667 +
                 rho * (xi_4 * xi_4) * -0.0416666666666667 +
                 rho *
                     ((xi_3 +
                       ((ctr_1 <= 0) ? (1.0) : (0.0)) * 0.10000000000000001) *
                      (xi_3 +
                       ((ctr_1 <= 0) ? (1.0) : (0.0)) * 0.10000000000000001)) *
                     -0.0416666666666667 +
                 rho * ((u1Pu2 * u1Pu2) * 0.125 + u1Pu2 * 0.0833333333333333 +
                        ((xi_3 + ((ctr_1 <= 0) ? (1.0) : (0.0)) *
                                     0.10000000000000001) *
                         (xi_3 + ((ctr_1 <= 0) ? (1.0) : (0.0)) *
                                     0.10000000000000001)) *
                            0.0416666666666667 -
                        0.0138888888888889) +
                 rho * 0.0416666666666667 - xi_2) +
            xi_2 +
            (omega * -0.5 + 1.0) *
                (xi_1 *
                     (xi_3 +
                      ((ctr_1 <= 0) ? (1.0) : (0.0)) * 0.10000000000000001) *
                     -0.0833333333333333 +
                 xi_15 * (xi_18 * 3.0 + xi_4 * 2.0 + 1.0) * 0.0833333333333333 +
                 xi_17 * (xi_18 * 2.0 + xi_4 * 3.0 + 1.0) * 0.0833333333333333);
        _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0] =
            omega *
                (rho * (xi_18 * xi_18) * -0.0416666666666667 +
                 rho * (xi_4 * xi_4) * -0.0416666666666667 +
                 rho *
                     ((xi_3 +
                       ((ctr_1 <= 0) ? (1.0) : (0.0)) * 0.10000000000000001) *
                      (xi_3 +
                       ((ctr_1 <= 0) ? (1.0) : (0.0)) * 0.10000000000000001)) *
                     -0.0416666666666667 +
                 rho * ((u1Mu2 * u1Mu2) * 0.125 + u1Mu2 * -0.0833333333333333 +
                        ((xi_3 + ((ctr_1 <= 0) ? (1.0) : (0.0)) *
                                     0.10000000000000001) *
                         (xi_3 + ((ctr_1 <= 0) ? (1.0) : (0.0)) *
                                     0.10000000000000001)) *
                            0.0416666666666667 -
                        0.0138888888888889) +
                 rho * 0.0416666666666667 - xi_6) +
            xi_6 +
            (omega * -0.5 + 1.0) * (xi_1 *
                                        (xi_3 + ((ctr_1 <= 0) ? (1.0) : (0.0)) *
                                                    0.10000000000000001) *
                                        -0.0833333333333333 +
                                    xi_15 * (xi_18 * 3.0 + xi_4 * -2.0 + 1.0) *
                                        -0.0833333333333333 +
                                    xi_17 * (xi_18 * 2.0 + xi_4 * -3.0 + 1.0) *
                                        0.0833333333333333);
        _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0] =
            omega *
                (rho * (xi_18 * xi_18) * -0.0416666666666667 +
                 rho * (xi_4 * xi_4) * -0.0416666666666667 +
                 rho *
                     ((xi_3 +
                       ((ctr_1 <= 0) ? (1.0) : (0.0)) * 0.10000000000000001) *
                      (xi_3 +
                       ((ctr_1 <= 0) ? (1.0) : (0.0)) * 0.10000000000000001)) *
                     -0.0416666666666667 +
                 rho *
                     (xi_18 * 0.0833333333333333 + xi_3 * -0.0833333333333333 +
                      (xi_4 * xi_4) * 0.0416666666666667 +
                      ((-xi_18 + xi_3 +
                        ((ctr_1 <= 0) ? (1.0) : (0.0)) * 0.10000000000000001) *
                       (-xi_18 + xi_3 +
                        ((ctr_1 <= 0) ? (1.0) : (0.0)) * 0.10000000000000001)) *
                          0.125 +
                      ((ctr_1 <= 0) ? (1.0) : (0.0)) * -0.0083333333333333332 -
                      0.0138888888888889) +
                 rho * 0.0416666666666667 - xi_10) +
            xi_10 +
            (omega * -0.5 + 1.0) *
                (xi_1 *
                     (xi_18 * 3.0 + xi_3 * -2.0 +
                      ((ctr_1 <= 0) ? (1.0) : (0.0)) * -0.20000000000000001 +
                      1.0) *
                     -0.0833333333333333 +
                 xi_15 * xi_4 * -0.0833333333333333 +
                 xi_17 *
                     (xi_18 * 2.0 + xi_3 * -3.0 +
                      ((ctr_1 <= 0) ? (1.0) : (0.0)) * -0.30000000000000004 +
                      1.0) *
                     0.0833333333333333);
        _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0] =
            omega *
                (rho * (xi_18 * xi_18) * -0.0416666666666667 +
                 rho * (xi_4 * xi_4) * -0.0416666666666667 +
                 rho *
                     ((xi_3 +
                       ((ctr_1 <= 0) ? (1.0) : (0.0)) * 0.10000000000000001) *
                      (xi_3 +
                       ((ctr_1 <= 0) ? (1.0) : (0.0)) * 0.10000000000000001)) *
                     -0.0416666666666667 +
                 rho *
                     (xi_18 * 0.0833333333333333 + xi_3 * 0.0833333333333333 +
                      (xi_4 * xi_4) * 0.0416666666666667 +
                      ((xi_18 + xi_3 +
                        ((ctr_1 <= 0) ? (1.0) : (0.0)) * 0.10000000000000001) *
                       (xi_18 + xi_3 +
                        ((ctr_1 <= 0) ? (1.0) : (0.0)) * 0.10000000000000001)) *
                          0.125 +
                      ((ctr_1 <= 0) ? (1.0) : (0.0)) * 0.0083333333333333332 -
                      0.0138888888888889) +
                 rho * 0.0416666666666667 - xi_7) +
            xi_7 +
            (omega * -0.5 + 1.0) *
                (xi_1 *
                     (xi_18 * 3.0 + xi_3 * 2.0 +
                      ((ctr_1 <= 0) ? (1.0) : (0.0)) * 0.20000000000000001 +
                      1.0) *
                     0.0833333333333333 +
                 xi_15 * xi_4 * -0.0833333333333333 +
                 xi_17 *
                     (xi_18 * 2.0 + xi_3 * 3.0 +
                      ((ctr_1 <= 0) ? (1.0) : (0.0)) * 0.30000000000000004 +
                      1.0) *
                     0.0833333333333333);
        _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0] =
            omega *
                (rho * (xi_18 * xi_18) * -0.0416666666666667 +
                 rho * (xi_4 * xi_4) * -0.0416666666666667 +
                 rho *
                     ((xi_3 +
                       ((ctr_1 >= 63) ? (-1.0) : (0.0)) * 0.10000000000000001) *
                      (xi_3 + ((ctr_1 >= 63) ? (-1.0) : (0.0)) *
                                  0.10000000000000001)) *
                     -0.0416666666666667 +
                 rho * ((u1Mu2 * u1Mu2) * 0.125 + u1Mu2 * 0.0833333333333333 +
                        ((xi_3 + ((ctr_1 >= 63) ? (-1.0) : (0.0)) *
                                     0.10000000000000001) *
                         (xi_3 + ((ctr_1 >= 63) ? (-1.0) : (0.0)) *
                                     0.10000000000000001)) *
                            0.0416666666666667 -
                        0.0138888888888889) +
                 rho * 0.0416666666666667 - xi_5) +
            xi_5 +
            (omega * -0.5 + 1.0) *
                (xi_1 *
                     (xi_3 +
                      ((ctr_1 >= 63) ? (-1.0) : (0.0)) * 0.10000000000000001) *
                     -0.0833333333333333 +
                 xi_15 * (xi_18 * -3.0 + xi_4 * 2.0 + 1.0) *
                     0.0833333333333333 +
                 xi_17 * (xi_18 * -2.0 + xi_4 * 3.0 + 1.0) *
                     -0.0833333333333333);
        _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0] =
            omega *
                (rho * (xi_18 * xi_18) * -0.0416666666666667 +
                 rho * (xi_4 * xi_4) * -0.0416666666666667 +
                 rho *
                     ((xi_3 +
                       ((ctr_1 >= 63) ? (-1.0) : (0.0)) * 0.10000000000000001) *
                      (xi_3 + ((ctr_1 >= 63) ? (-1.0) : (0.0)) *
                                  0.10000000000000001)) *
                     -0.0416666666666667 +
                 rho * ((u1Pu2 * u1Pu2) * 0.125 + u1Pu2 * -0.0833333333333333 +
                        ((xi_3 + ((ctr_1 >= 63) ? (-1.0) : (0.0)) *
                                     0.10000000000000001) *
                         (xi_3 + ((ctr_1 >= 63) ? (-1.0) : (0.0)) *
                                     0.10000000000000001)) *
                            0.0416666666666667 -
                        0.0138888888888889) +
                 rho * 0.0416666666666667 - xi_12) +
            xi_12 +
            (omega * -0.5 + 1.0) *
                (xi_1 *
                     (xi_3 +
                      ((ctr_1 >= 63) ? (-1.0) : (0.0)) * 0.10000000000000001) *
                     -0.0833333333333333 +
                 xi_15 * (xi_18 * 3.0 + xi_4 * 2.0 - 1.0) * 0.0833333333333333 +
                 xi_17 * (xi_18 * 2.0 + xi_4 * 3.0 - 1.0) * 0.0833333333333333);
        _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0] =
            omega * (rho * (xi_18 * xi_18) * -0.0416666666666667 +
                     rho * (xi_4 * xi_4) * -0.0416666666666667 +
                     rho *
                         ((xi_3 + ((ctr_1 >= 63) ? (-1.0) : (0.0)) *
                                      0.10000000000000001) *
                          (xi_3 + ((ctr_1 >= 63) ? (-1.0) : (0.0)) *
                                      0.10000000000000001)) *
                         -0.0416666666666667 +
                     rho * (xi_18 * -0.0833333333333333 +
                            xi_3 * -0.0833333333333333 +
                            (xi_4 * xi_4) * 0.0416666666666667 +
                            ((xi_18 + xi_3 +
                              ((ctr_1 >= 63) ? (-1.0) : (0.0)) *
                                  0.10000000000000001) *
                             (xi_18 + xi_3 +
                              ((ctr_1 >= 63) ? (-1.0) : (0.0)) *
                                  0.10000000000000001)) *
                                0.125 +
                            ((ctr_1 >= 63) ? (-1.0) : (0.0)) *
                                -0.0083333333333333332 -
                            0.0138888888888889) +
                     rho * 0.0416666666666667 - xi_23) +
            xi_23 +
            (omega * -0.5 + 1.0) *
                (xi_1 *
                     (xi_18 * 3.0 + xi_3 * 2.0 +
                      ((ctr_1 >= 63) ? (-1.0) : (0.0)) * 0.20000000000000001 -
                      1.0) *
                     0.0833333333333333 +
                 xi_15 * xi_4 * -0.0833333333333333 +
                 xi_17 *
                     (xi_18 * 2.0 + xi_3 * 3.0 +
                      ((ctr_1 >= 63) ? (-1.0) : (0.0)) * 0.30000000000000004 -
                      1.0) *
                     0.0833333333333333);
        _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0] =
            omega *
                (rho * (xi_18 * xi_18) * -0.0416666666666667 +
                 rho * (xi_4 * xi_4) * -0.0416666666666667 +
                 rho *
                     ((xi_3 +
                       ((ctr_1 >= 63) ? (-1.0) : (0.0)) * 0.10000000000000001) *
                      (xi_3 + ((ctr_1 >= 63) ? (-1.0) : (0.0)) *
                                  0.10000000000000001)) *
                     -0.0416666666666667 +
                 rho *
                     (xi_18 * -0.0833333333333333 + xi_3 * 0.0833333333333333 +
                      (xi_4 * xi_4) * 0.0416666666666667 +
                      ((-xi_18 + xi_3 +
                        ((ctr_1 >= 63) ? (-1.0) : (0.0)) *
                            0.10000000000000001) *
                       (-xi_18 + xi_3 +
                        ((ctr_1 >= 63) ? (-1.0) : (0.0)) *
                            0.10000000000000001)) *
                          0.125 +
                      ((ctr_1 >= 63) ? (-1.0) : (0.0)) * 0.0083333333333333332 -
                      0.0138888888888889) +
                 rho * 0.0416666666666667 - xi_19) +
            xi_19 +
            (omega * -0.5 + 1.0) *
                (xi_1 *
                     (xi_18 * -3.0 + xi_3 * 2.0 +
                      ((ctr_1 >= 63) ? (-1.0) : (0.0)) * 0.20000000000000001 +
                      1.0) *
                     0.0833333333333333 +
                 xi_15 * xi_4 * -0.0833333333333333 +
                 xi_17 *
                     (xi_18 * -2.0 + xi_3 * 3.0 +
                      ((ctr_1 >= 63) ? (-1.0) : (0.0)) * 0.30000000000000004 +
                      1.0) *
                     -0.0833333333333333);
      }
    }
  }
}
} // namespace internal_leesedwards_collision_leesedwards_collision

void LeesEdwards_Collision::run(IBlock *block) {
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);
  auto velocity = block->getData<field::GhostLayerField<double, 3>>(velocityID);
  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);

  auto &omega = this->omega_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  double *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  double *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(velocity->nrOfGhostLayers()));
  double *RESTRICT const _data_velocity = velocity->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->xSize()) + 0));
  const int64_t _size_force_0 = int64_t(cell_idx_c(force->xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->ySize()) + 0));
  const int64_t _size_force_1 = int64_t(cell_idx_c(force->ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->zSize()) + 0));
  const int64_t _size_force_2 = int64_t(cell_idx_c(force->zSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  const int64_t _stride_force_0 = int64_t(force->xStride());
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
  const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
  const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
  const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
  internal_leesedwards_collision_leesedwards_collision::
      leesedwards_collision_leesedwards_collision(
          _data_force, _data_pdfs, _data_velocity, _size_force_0, _size_force_1,
          _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2,
          _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2,
          _stride_pdfs_3, _stride_velocity_0, _stride_velocity_1,
          _stride_velocity_2, _stride_velocity_3, omega);
}

void LeesEdwards_Collision::runOnCellInterval(
    const shared_ptr<StructuredBlockStorage> &blocks,
    const CellInterval &globalCellInterval, cell_idx_t ghostLayers,
    IBlock *block) {
  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);
  auto velocity = block->getData<field::GhostLayerField<double, 3>>(velocityID);
  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);

  auto &omega = this->omega_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force->nrOfGhostLayers()));
  double *RESTRICT const _data_force =
      force->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
  double *RESTRICT _data_pdfs =
      pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(velocity->nrOfGhostLayers()));
  double *RESTRICT const _data_velocity =
      velocity->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.xSize()) + 0));
  const int64_t _size_force_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.ySize()) + 0));
  const int64_t _size_force_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.zSize()) + 0));
  const int64_t _size_force_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  const int64_t _stride_force_0 = int64_t(force->xStride());
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
  const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
  const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
  const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
  internal_leesedwards_collision_leesedwards_collision::
      leesedwards_collision_leesedwards_collision(
          _data_force, _data_pdfs, _data_velocity, _size_force_0, _size_force_1,
          _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2,
          _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2,
          _stride_pdfs_3, _stride_velocity_0, _stride_velocity_1,
          _stride_velocity_2, _stride_velocity_3, omega);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif