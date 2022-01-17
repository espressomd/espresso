// kernel generated with pystencils v0.4.4, lbmpy v0.4.3+6.g13df23c,
// lbmpy_walberla/pystencils_walberla from commit ref: refs/heads/le_ghost_vel

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
//! \\file CollideSweepDoublePrecisionLeesEdwards.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "CollideSweepDoublePrecisionLeesEdwards.h"
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

namespace internal_607d8a5c7ac58c25acf09ad94bb82cf4 {
static FUNC_PREFIX void
collidesweepdoubleprecisionleesedwards_collidesweepdoubleprecisionleesedwards(
    double *RESTRICT const _data_force, double *RESTRICT _data_pdfs,
    int64_t const _size_force_0, int64_t const _size_force_1,
    int64_t const _size_force_2, int64_t const _stride_force_0,
    int64_t const _stride_force_1, int64_t const _stride_force_2,
    int64_t const _stride_force_3, int64_t const _stride_pdfs_0,
    int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
    int64_t const _stride_pdfs_3, double grid_size, double omega_shear,
    double v_s) {
  const double xi_0 = 1 / (omega_shear * -0.25000000000000000 + 2.0);
  const double rr_0 = xi_0 * (omega_shear * -2.0 + 4.0);
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    double *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    double *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    double *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    double *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      double *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      double *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      double *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      double *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      double *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      double *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      double *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      double *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      double *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      double *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      double *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      double *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      double *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      double *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      double *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      double *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      double *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      double *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      double *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      double *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      double *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      for (int64_t ctr_0 = 0; ctr_0 < _size_force_0; ctr_0 += 1) {
        const double xi_25 = _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0];
        const double xi_26 = _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0];
        const double xi_27 = _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0];
        const double xi_28 = _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0];
        const double xi_29 = _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0];
        const double xi_30 = _data_force_20_31_10[_stride_force_0 * ctr_0];
        const double xi_31 = _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0];
        const double xi_32 = _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0];
        const double xi_33 = _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0];
        const double xi_34 = _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0];
        const double xi_35 = _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const double xi_36 = _data_force_20_30_10[_stride_force_0 * ctr_0];
        const double xi_37 = _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0];
        const double xi_38 = _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0];
        const double xi_39 = _data_force_20_32_10[_stride_force_0 * ctr_0];
        const double xi_40 = _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0];
        const double xi_41 = _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0];
        const double xi_42 = _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0];
        const double xi_43 = _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0];
        const double xi_44 = _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0];
        const double xi_45 = _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0];
        const double xi_46 = _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0];
        const double xi_3 = xi_26;
        const double xi_4 = xi_45;
        const double xi_5 = xi_27;
        const double xi_6 = xi_28;
        const double xi_7 = xi_29;
        const double xi_8 = xi_30;
        const double xi_9 = xi_31;
        const double xi_10 = xi_32;
        const double xi_11 = xi_33;
        const double xi_12 = xi_34;
        const double xi_13 = xi_35;
        const double xi_14 = xi_36;
        const double xi_15 = xi_38;
        const double xi_16 = xi_37;
        const double xi_17 = xi_39;
        const double xi_18 = xi_40;
        const double xi_19 = xi_41;
        const double xi_20 = xi_42;
        const double xi_21 = xi_43;
        const double xi_22 = xi_44;
        const double xi_23 = xi_25;
        const double xi_24 = xi_46;
        const double vel0Term = xi_16 + xi_18 + xi_19 + xi_22 + xi_6;
        const double vel1Term = xi_10 + xi_11 + xi_15 + xi_5;
        const double vel2Term = xi_20 + xi_24 + xi_3;
        const double rho = vel0Term + vel1Term + vel2Term + xi_12 + xi_13 +
                           xi_21 + xi_23 + xi_4 + xi_7 + xi_9;
        const double xi_1 = 1 / (rho);
        const double u_0 =
            xi_1 * xi_14 * 0.50000000000000000 +
            xi_1 * (vel0Term - xi_11 - xi_12 - xi_23 - xi_3 - xi_4);
        const double u_1 =
            xi_1 * xi_8 * 0.50000000000000000 +
            xi_1 * (vel1Term + xi_16 - xi_19 - xi_21 - xi_23 - xi_24 - xi_9);
        const double u_2 = xi_1 * xi_17 * 0.50000000000000000 +
                           xi_1 * (vel2Term - xi_10 + xi_15 - xi_21 - xi_22 -
                                   xi_4 + xi_6 - xi_7);
        const double forceTerm_0 =
            omega_shear * u_0 * xi_14 * 0.50000000000000000 +
            omega_shear * u_1 * xi_8 * 0.50000000000000000 +
            omega_shear * u_2 * xi_17 * 0.50000000000000000 - u_0 * xi_14 -
            u_1 * xi_8 - u_2 * xi_17;
        const double forceTerm_1 =
            omega_shear * u_0 * xi_14 * 0.083333333333333333 +
            omega_shear * u_1 * xi_8 * -0.16666666666666667 +
            omega_shear * u_2 * xi_17 * 0.083333333333333333 +
            rr_0 * xi_8 * -0.083333333333333333 +
            u_0 * xi_14 * -0.16666666666666667 +
            u_1 * xi_8 * 0.33333333333333333 +
            u_2 * xi_17 * -0.16666666666666667 + xi_8 * 0.16666666666666667;
        const double forceTerm_2 =
            omega_shear * u_0 * xi_14 * 0.083333333333333333 +
            omega_shear * u_1 * xi_8 * -0.16666666666666667 +
            omega_shear * u_2 * xi_17 * 0.083333333333333333 +
            rr_0 * xi_8 * 0.083333333333333333 +
            u_0 * xi_14 * -0.16666666666666667 +
            u_1 * xi_8 * 0.33333333333333333 +
            u_2 * xi_17 * -0.16666666666666667 + xi_8 * -0.16666666666666667;
        const double forceTerm_3 =
            omega_shear * u_0 * xi_14 * -0.16666666666666667 +
            omega_shear * u_1 * xi_8 * 0.083333333333333333 +
            omega_shear * u_2 * xi_17 * 0.083333333333333333 +
            rr_0 * xi_14 * 0.083333333333333333 +
            u_0 * xi_14 * 0.33333333333333333 +
            u_1 * xi_8 * -0.16666666666666667 +
            u_2 * xi_17 * -0.16666666666666667 + xi_14 * -0.16666666666666667;
        const double forceTerm_4 =
            omega_shear * u_0 * xi_14 * -0.16666666666666667 +
            omega_shear * u_1 * xi_8 * 0.083333333333333333 +
            omega_shear * u_2 * xi_17 * 0.083333333333333333 +
            rr_0 * xi_14 * -0.083333333333333333 +
            u_0 * xi_14 * 0.33333333333333333 +
            u_1 * xi_8 * -0.16666666666666667 +
            u_2 * xi_17 * -0.16666666666666667 + xi_14 * 0.16666666666666667;
        const double forceTerm_5 =
            omega_shear * u_0 * xi_14 * 0.083333333333333333 +
            omega_shear * u_1 * xi_8 * 0.083333333333333333 +
            omega_shear * u_2 * xi_17 * -0.16666666666666667 +
            rr_0 * xi_17 * -0.083333333333333333 +
            u_0 * xi_14 * -0.16666666666666667 +
            u_1 * xi_8 * -0.16666666666666667 +
            u_2 * xi_17 * 0.33333333333333333 + xi_17 * 0.16666666666666667;
        const double forceTerm_6 =
            omega_shear * u_0 * xi_14 * 0.083333333333333333 +
            omega_shear * u_1 * xi_8 * 0.083333333333333333 +
            omega_shear * u_2 * xi_17 * -0.16666666666666667 +
            rr_0 * xi_17 * 0.083333333333333333 +
            u_0 * xi_14 * -0.16666666666666667 +
            u_1 * xi_8 * -0.16666666666666667 +
            u_2 * xi_17 * 0.33333333333333333 + xi_17 * -0.16666666666666667;
        const double forceTerm_7 =
            omega_shear * u_0 * xi_14 * -0.083333333333333333 +
            omega_shear * u_0 * xi_8 * 0.12500000000000000 +
            omega_shear * u_1 * xi_14 * 0.12500000000000000 +
            omega_shear * u_1 * xi_8 * -0.083333333333333333 +
            omega_shear * u_2 * xi_17 * 0.041666666666666667 +
            rr_0 * xi_14 * 0.041666666666666667 +
            rr_0 * xi_8 * -0.041666666666666667 +
            u_0 * xi_14 * 0.16666666666666667 +
            u_0 * xi_8 * -0.25000000000000000 +
            u_1 * xi_14 * -0.25000000000000000 +
            u_1 * xi_8 * 0.16666666666666667 +
            u_2 * xi_17 * -0.083333333333333333 +
            xi_14 * -0.083333333333333333 + xi_8 * 0.083333333333333333;
        const double forceTerm_8 =
            omega_shear * u_0 * xi_14 * -0.083333333333333333 +
            omega_shear * u_0 * xi_8 * -0.12500000000000000 +
            omega_shear * u_1 * xi_14 * -0.12500000000000000 +
            omega_shear * u_1 * xi_8 * -0.083333333333333333 +
            omega_shear * u_2 * xi_17 * 0.041666666666666667 +
            rr_0 * xi_14 * -0.041666666666666667 +
            rr_0 * xi_8 * -0.041666666666666667 +
            u_0 * xi_14 * 0.16666666666666667 +
            u_0 * xi_8 * 0.25000000000000000 +
            u_1 * xi_14 * 0.25000000000000000 +
            u_1 * xi_8 * 0.16666666666666667 +
            u_2 * xi_17 * -0.083333333333333333 + xi_14 * 0.083333333333333333 +
            xi_8 * 0.083333333333333333;
        const double forceTerm_9 =
            omega_shear * u_0 * xi_14 * -0.083333333333333333 +
            omega_shear * u_0 * xi_8 * -0.12500000000000000 +
            omega_shear * u_1 * xi_14 * -0.12500000000000000 +
            omega_shear * u_1 * xi_8 * -0.083333333333333333 +
            omega_shear * u_2 * xi_17 * 0.041666666666666667 +
            rr_0 * xi_14 * 0.041666666666666667 +
            rr_0 * xi_8 * 0.041666666666666667 +
            u_0 * xi_14 * 0.16666666666666667 +
            u_0 * xi_8 * 0.25000000000000000 +
            u_1 * xi_14 * 0.25000000000000000 +
            u_1 * xi_8 * 0.16666666666666667 +
            u_2 * xi_17 * -0.083333333333333333 +
            xi_14 * -0.083333333333333333 + xi_8 * -0.083333333333333333;
        const double forceTerm_10 =
            omega_shear * u_0 * xi_14 * -0.083333333333333333 +
            omega_shear * u_0 * xi_8 * 0.12500000000000000 +
            omega_shear * u_1 * xi_14 * 0.12500000000000000 +
            omega_shear * u_1 * xi_8 * -0.083333333333333333 +
            omega_shear * u_2 * xi_17 * 0.041666666666666667 +
            rr_0 * xi_14 * -0.041666666666666667 +
            rr_0 * xi_8 * 0.041666666666666667 +
            u_0 * xi_14 * 0.16666666666666667 +
            u_0 * xi_8 * -0.25000000000000000 +
            u_1 * xi_14 * -0.25000000000000000 +
            u_1 * xi_8 * 0.16666666666666667 +
            u_2 * xi_17 * -0.083333333333333333 + xi_14 * 0.083333333333333333 +
            xi_8 * -0.083333333333333333;
        const double forceTerm_11 =
            omega_shear * u_0 * xi_14 * 0.041666666666666667 +
            omega_shear * u_1 * xi_17 * -0.12500000000000000 +
            omega_shear * u_1 * xi_8 * -0.083333333333333333 +
            omega_shear * u_2 * xi_17 * -0.083333333333333333 +
            omega_shear * u_2 * xi_8 * -0.12500000000000000 +
            rr_0 * xi_17 * -0.041666666666666667 +
            rr_0 * xi_8 * -0.041666666666666667 +
            u_0 * xi_14 * -0.083333333333333333 +
            u_1 * xi_17 * 0.25000000000000000 +
            u_1 * xi_8 * 0.16666666666666667 +
            u_2 * xi_17 * 0.16666666666666667 +
            u_2 * xi_8 * 0.25000000000000000 + xi_17 * 0.083333333333333333 +
            xi_8 * 0.083333333333333333;
        const double forceTerm_12 =
            omega_shear * u_0 * xi_14 * 0.041666666666666667 +
            omega_shear * u_1 * xi_17 * 0.12500000000000000 +
            omega_shear * u_1 * xi_8 * -0.083333333333333333 +
            omega_shear * u_2 * xi_17 * -0.083333333333333333 +
            omega_shear * u_2 * xi_8 * 0.12500000000000000 +
            rr_0 * xi_17 * -0.041666666666666667 +
            rr_0 * xi_8 * 0.041666666666666667 +
            u_0 * xi_14 * -0.083333333333333333 +
            u_1 * xi_17 * -0.25000000000000000 +
            u_1 * xi_8 * 0.16666666666666667 +
            u_2 * xi_17 * 0.16666666666666667 +
            u_2 * xi_8 * -0.25000000000000000 + xi_17 * 0.083333333333333333 +
            xi_8 * -0.083333333333333333;
        const double forceTerm_13 =
            omega_shear * u_0 * xi_14 * -0.083333333333333333 +
            omega_shear * u_0 * xi_17 * 0.12500000000000000 +
            omega_shear * u_1 * xi_8 * 0.041666666666666667 +
            omega_shear * u_2 * xi_14 * 0.12500000000000000 +
            omega_shear * u_2 * xi_17 * -0.083333333333333333 +
            rr_0 * xi_14 * 0.041666666666666667 +
            rr_0 * xi_17 * -0.041666666666666667 +
            u_0 * xi_14 * 0.16666666666666667 +
            u_0 * xi_17 * -0.25000000000000000 +
            u_1 * xi_8 * -0.083333333333333333 +
            u_2 * xi_14 * -0.25000000000000000 +
            u_2 * xi_17 * 0.16666666666666667 + xi_14 * -0.083333333333333333 +
            xi_17 * 0.083333333333333333;
        const double forceTerm_14 =
            omega_shear * u_0 * xi_14 * -0.083333333333333333 +
            omega_shear * u_0 * xi_17 * -0.12500000000000000 +
            omega_shear * u_1 * xi_8 * 0.041666666666666667 +
            omega_shear * u_2 * xi_14 * -0.12500000000000000 +
            omega_shear * u_2 * xi_17 * -0.083333333333333333 +
            rr_0 * xi_14 * -0.041666666666666667 +
            rr_0 * xi_17 * -0.041666666666666667 +
            u_0 * xi_14 * 0.16666666666666667 +
            u_0 * xi_17 * 0.25000000000000000 +
            u_1 * xi_8 * -0.083333333333333333 +
            u_2 * xi_14 * 0.25000000000000000 +
            u_2 * xi_17 * 0.16666666666666667 + xi_14 * 0.083333333333333333 +
            xi_17 * 0.083333333333333333;
        const double forceTerm_15 =
            omega_shear * u_0 * xi_14 * 0.041666666666666667 +
            omega_shear * u_1 * xi_17 * 0.12500000000000000 +
            omega_shear * u_1 * xi_8 * -0.083333333333333333 +
            omega_shear * u_2 * xi_17 * -0.083333333333333333 +
            omega_shear * u_2 * xi_8 * 0.12500000000000000 +
            rr_0 * xi_17 * 0.041666666666666667 +
            rr_0 * xi_8 * -0.041666666666666667 +
            u_0 * xi_14 * -0.083333333333333333 +
            u_1 * xi_17 * -0.25000000000000000 +
            u_1 * xi_8 * 0.16666666666666667 +
            u_2 * xi_17 * 0.16666666666666667 +
            u_2 * xi_8 * -0.25000000000000000 + xi_17 * -0.083333333333333333 +
            xi_8 * 0.083333333333333333;
        const double forceTerm_16 =
            omega_shear * u_0 * xi_14 * 0.041666666666666667 +
            omega_shear * u_1 * xi_17 * -0.12500000000000000 +
            omega_shear * u_1 * xi_8 * -0.083333333333333333 +
            omega_shear * u_2 * xi_17 * -0.083333333333333333 +
            omega_shear * u_2 * xi_8 * -0.12500000000000000 +
            rr_0 * xi_17 * 0.041666666666666667 +
            rr_0 * xi_8 * 0.041666666666666667 +
            u_0 * xi_14 * -0.083333333333333333 +
            u_1 * xi_17 * 0.25000000000000000 +
            u_1 * xi_8 * 0.16666666666666667 +
            u_2 * xi_17 * 0.16666666666666667 +
            u_2 * xi_8 * 0.25000000000000000 + xi_17 * -0.083333333333333333 +
            xi_8 * -0.083333333333333333;
        const double forceTerm_17 =
            omega_shear * u_0 * xi_14 * -0.083333333333333333 +
            omega_shear * u_0 * xi_17 * -0.12500000000000000 +
            omega_shear * u_1 * xi_8 * 0.041666666666666667 +
            omega_shear * u_2 * xi_14 * -0.12500000000000000 +
            omega_shear * u_2 * xi_17 * -0.083333333333333333 +
            rr_0 * xi_14 * 0.041666666666666667 +
            rr_0 * xi_17 * 0.041666666666666667 +
            u_0 * xi_14 * 0.16666666666666667 +
            u_0 * xi_17 * 0.25000000000000000 +
            u_1 * xi_8 * -0.083333333333333333 +
            u_2 * xi_14 * 0.25000000000000000 +
            u_2 * xi_17 * 0.16666666666666667 + xi_14 * -0.083333333333333333 +
            xi_17 * -0.083333333333333333;
        const double forceTerm_18 =
            omega_shear * u_0 * xi_14 * -0.083333333333333333 +
            omega_shear * u_0 * xi_17 * 0.12500000000000000 +
            omega_shear * u_1 * xi_8 * 0.041666666666666667 +
            omega_shear * u_2 * xi_14 * 0.12500000000000000 +
            omega_shear * u_2 * xi_17 * -0.083333333333333333 +
            rr_0 * xi_14 * -0.041666666666666667 +
            rr_0 * xi_17 * 0.041666666666666667 +
            u_0 * xi_14 * 0.16666666666666667 +
            u_0 * xi_17 * -0.25000000000000000 +
            u_1 * xi_8 * -0.083333333333333333 +
            u_2 * xi_14 * -0.25000000000000000 +
            u_2 * xi_17 * 0.16666666666666667 + xi_14 * 0.083333333333333333 +
            xi_17 * -0.083333333333333333;
        const double u0Mu1 = u_0 - u_1;
        const double u0Pu1 = u_0 + u_1;
        const double u1Pu2 = u_1 + u_2;
        const double u1Mu2 = u_1 - u_2;
        const double u0Mu2 = u_0 - u_2;
        const double u0Pu2 = u_0 + u_2;
        const double f_eq_common =
            -rho * (u_0 * u_0) - rho * (u_1 * u_1) - rho * (u_2 * u_2) + rho;
        _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_0 +
            omega_shear * (f_eq_common * 0.33333333333333333 - xi_13) + xi_13;
        _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_1 +
            omega_shear *
                (f_eq_common * 0.16666666666666667 +
                 rho *
                     ((u_1 * u_1) * 0.33333333333333333 - 0.11111111111111111) +
                 xi_5 * -0.50000000000000000 + xi_9 * -0.50000000000000000) +
            rr_0 * (rho * u_1 * 0.16666666666666667 +
                    xi_5 * -0.50000000000000000 + xi_9 * 0.50000000000000000) +
            xi_5 +
            ((-1 <= ctr_1 - grid_size)
                 ? (rho * v_s * (u_0 * 2.0 + v_s) * 0.16666666666666667)
                 : (0.0));
        _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_2 +
            omega_shear *
                (f_eq_common * 0.16666666666666667 +
                 rho *
                     ((u_1 * u_1) * 0.33333333333333333 - 0.11111111111111111) +
                 xi_5 * -0.50000000000000000 + xi_9 * -0.50000000000000000) +
            rr_0 * (rho * u_1 * -0.16666666666666667 +
                    xi_5 * 0.50000000000000000 + xi_9 * -0.50000000000000000) +
            xi_9 +
            ((ctr_1 <= 0)
                 ? (rho * v_s * (u_0 * -2.0 + v_s) * 0.16666666666666667)
                 : (0.0));
        _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_3 +
            omega_shear *
                (f_eq_common * 0.16666666666666667 +
                 rho *
                     ((u_0 * u_0) * 0.33333333333333333 - 0.11111111111111111) +
                 xi_12 * -0.50000000000000000 + xi_18 * -0.50000000000000000) +
            rr_0 *
                (rho * u_0 * -0.16666666666666667 +
                 xi_12 * -0.50000000000000000 + xi_18 * 0.50000000000000000) +
            xi_12;
        _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_4 +
            omega_shear *
                (f_eq_common * 0.16666666666666667 +
                 rho *
                     ((u_0 * u_0) * 0.33333333333333333 - 0.11111111111111111) +
                 xi_12 * -0.50000000000000000 + xi_18 * -0.50000000000000000) +
            rr_0 *
                (rho * u_0 * 0.16666666666666667 + xi_12 * 0.50000000000000000 +
                 xi_18 * -0.50000000000000000) +
            xi_18;
        _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_5 +
            omega_shear *
                (f_eq_common * 0.16666666666666667 +
                 rho *
                     ((u_2 * u_2) * 0.33333333333333333 - 0.11111111111111111) +
                 xi_20 * -0.50000000000000000 + xi_7 * -0.50000000000000000) +
            rr_0 * (rho * u_2 * 0.16666666666666667 +
                    xi_20 * -0.50000000000000000 + xi_7 * 0.50000000000000000) +
            xi_20;
        _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_6 +
            omega_shear *
                (f_eq_common * 0.16666666666666667 +
                 rho *
                     ((u_2 * u_2) * 0.33333333333333333 - 0.11111111111111111) +
                 xi_20 * -0.50000000000000000 + xi_7 * -0.50000000000000000) +
            rr_0 * (rho * u_2 * -0.16666666666666667 +
                    xi_20 * 0.50000000000000000 + xi_7 * -0.50000000000000000) +
            xi_7;
        _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_7 +
            omega_shear *
                (f_eq_common * 0.041666666666666667 +
                 rho * ((u0Mu1 * u0Mu1) * 0.12500000000000000 +
                        (u_2 * u_2) * 0.041666666666666667 -
                        0.013888888888888889) +
                 xi_11 * -0.50000000000000000 + xi_19 * -0.50000000000000000) +
            rr_0 *
                (rho * u0Mu1 * -0.083333333333333333 +
                 xi_11 * -0.50000000000000000 + xi_19 * 0.50000000000000000) +
            xi_11 +
            ((-1 <= ctr_1 - grid_size)
                 ? (rho * v_s * (u_0 * -2.0 + u_1 * 3.0 - v_s + 1.0) *
                    0.083333333333333333)
                 : (0.0));
        _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_8 +
            omega_shear *
                (f_eq_common * 0.041666666666666667 +
                 rho * ((u0Pu1 * u0Pu1) * 0.12500000000000000 +
                        (u_2 * u_2) * 0.041666666666666667 -
                        0.013888888888888889) +
                 xi_16 * -0.50000000000000000 + xi_23 * -0.50000000000000000) +
            rr_0 *
                (rho * u0Pu1 * 0.083333333333333333 +
                 xi_16 * -0.50000000000000000 + xi_23 * 0.50000000000000000) +
            xi_16 +
            ((-1 <= ctr_1 - grid_size)
                 ? (rho * v_s * (u_0 * 2.0 + u_1 * 3.0 + v_s + 1.0) *
                    -0.083333333333333333)
                 : (0.0));
        _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_9 +
            omega_shear *
                (f_eq_common * 0.041666666666666667 +
                 rho * ((u0Pu1 * u0Pu1) * 0.12500000000000000 +
                        (u_2 * u_2) * 0.041666666666666667 -
                        0.013888888888888889) +
                 xi_16 * -0.50000000000000000 + xi_23 * -0.50000000000000000) +
            rr_0 *
                (rho * u0Pu1 * -0.083333333333333333 +
                 xi_16 * 0.50000000000000000 + xi_23 * -0.50000000000000000) +
            xi_23 +
            ((ctr_1 <= 0) ? (rho * v_s * (u_0 * 2.0 + u_1 * 3.0 - v_s - 1.0) *
                             0.083333333333333333)
                          : (0.0));
        _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_10 +
            omega_shear *
                (f_eq_common * 0.041666666666666667 +
                 rho * ((u0Mu1 * u0Mu1) * 0.12500000000000000 +
                        (u_2 * u_2) * 0.041666666666666667 -
                        0.013888888888888889) +
                 xi_11 * -0.50000000000000000 + xi_19 * -0.50000000000000000) +
            rr_0 *
                (rho * u0Mu1 * 0.083333333333333333 +
                 xi_11 * 0.50000000000000000 + xi_19 * -0.50000000000000000) +
            xi_19 +
            ((ctr_1 <= 0) ? (rho * v_s * (u_0 * 2.0 + u_1 * -3.0 - v_s + 1.0) *
                             0.083333333333333333)
                          : (0.0));
        _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_11 +
            omega_shear *
                (f_eq_common * 0.041666666666666667 +
                 rho * ((u1Pu2 * u1Pu2) * 0.12500000000000000 +
                        (u_0 * u_0) * 0.041666666666666667 -
                        0.013888888888888889) +
                 xi_15 * -0.50000000000000000 + xi_21 * -0.50000000000000000) +
            rr_0 *
                (rho * u1Pu2 * 0.083333333333333333 +
                 xi_15 * -0.50000000000000000 + xi_21 * 0.50000000000000000) +
            xi_15;
        _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_12 +
            omega_shear *
                (f_eq_common * 0.041666666666666667 +
                 rho * ((u1Mu2 * u1Mu2) * 0.12500000000000000 +
                        (u_0 * u_0) * 0.041666666666666667 -
                        0.013888888888888889) +
                 xi_10 * -0.50000000000000000 + xi_24 * -0.50000000000000000) +
            rr_0 *
                (rho * u1Mu2 * -0.083333333333333333 +
                 xi_10 * 0.50000000000000000 + xi_24 * -0.50000000000000000) +
            xi_24;
        _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_13 +
            omega_shear *
                (f_eq_common * 0.041666666666666667 +
                 rho * ((u0Mu2 * u0Mu2) * 0.12500000000000000 +
                        (u_1 * u_1) * 0.041666666666666667 -
                        0.013888888888888889) +
                 xi_22 * -0.50000000000000000 + xi_3 * -0.50000000000000000) +
            rr_0 * (rho * u0Mu2 * -0.083333333333333333 +
                    xi_22 * 0.50000000000000000 + xi_3 * -0.50000000000000000) +
            xi_3;
        _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_14 +
            omega_shear *
                (f_eq_common * 0.041666666666666667 +
                 rho * ((u0Pu2 * u0Pu2) * 0.12500000000000000 +
                        (u_1 * u_1) * 0.041666666666666667 -
                        0.013888888888888889) +
                 xi_4 * -0.50000000000000000 + xi_6 * -0.50000000000000000) +
            rr_0 * (rho * u0Pu2 * 0.083333333333333333 +
                    xi_4 * 0.50000000000000000 + xi_6 * -0.50000000000000000) +
            xi_6;
        _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_15 +
            omega_shear *
                (f_eq_common * 0.041666666666666667 +
                 rho * ((u1Mu2 * u1Mu2) * 0.12500000000000000 +
                        (u_0 * u_0) * 0.041666666666666667 -
                        0.013888888888888889) +
                 xi_10 * -0.50000000000000000 + xi_24 * -0.50000000000000000) +
            rr_0 *
                (rho * u1Mu2 * 0.083333333333333333 +
                 xi_10 * -0.50000000000000000 + xi_24 * 0.50000000000000000) +
            xi_10;
        _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_16 +
            omega_shear *
                (f_eq_common * 0.041666666666666667 +
                 rho * ((u1Pu2 * u1Pu2) * 0.12500000000000000 +
                        (u_0 * u_0) * 0.041666666666666667 -
                        0.013888888888888889) +
                 xi_15 * -0.50000000000000000 + xi_21 * -0.50000000000000000) +
            rr_0 *
                (rho * u1Pu2 * -0.083333333333333333 +
                 xi_15 * 0.50000000000000000 + xi_21 * -0.50000000000000000) +
            xi_21;
        _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_17 +
            omega_shear *
                (f_eq_common * 0.041666666666666667 +
                 rho * ((u0Pu2 * u0Pu2) * 0.12500000000000000 +
                        (u_1 * u_1) * 0.041666666666666667 -
                        0.013888888888888889) +
                 xi_4 * -0.50000000000000000 + xi_6 * -0.50000000000000000) +
            rr_0 * (rho * u0Pu2 * -0.083333333333333333 +
                    xi_4 * -0.50000000000000000 + xi_6 * 0.50000000000000000) +
            xi_4;
        _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_18 +
            omega_shear *
                (f_eq_common * 0.041666666666666667 +
                 rho * ((u0Mu2 * u0Mu2) * 0.12500000000000000 +
                        (u_1 * u_1) * 0.041666666666666667 -
                        0.013888888888888889) +
                 xi_22 * -0.50000000000000000 + xi_3 * -0.50000000000000000) +
            rr_0 * (rho * u0Mu2 * 0.083333333333333333 +
                    xi_22 * -0.50000000000000000 + xi_3 * 0.50000000000000000) +
            xi_22;
      }
    }
  }
}
} // namespace internal_607d8a5c7ac58c25acf09ad94bb82cf4

void CollideSweepDoublePrecisionLeesEdwards::run(IBlock *block) {
  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);

  auto &omega_shear = this->omega_shear_;
  auto &grid_size = this->grid_size_;
  auto &v_s = this->v_s_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  double *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  double *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
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
  internal_607d8a5c7ac58c25acf09ad94bb82cf4::
      collidesweepdoubleprecisionleesedwards_collidesweepdoubleprecisionleesedwards(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3,
          _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3,
          grid_size, omega_shear, v_s);
}

void CollideSweepDoublePrecisionLeesEdwards::runOnCellInterval(
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

  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);

  auto &omega_shear = this->omega_shear_;
  auto &grid_size = this->grid_size_;
  auto &v_s = this->v_s_;
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
  internal_607d8a5c7ac58c25acf09ad94bb82cf4::
      collidesweepdoubleprecisionleesedwards_collidesweepdoubleprecisionleesedwards(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3,
          _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3,
          grid_size, omega_shear, v_s);
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