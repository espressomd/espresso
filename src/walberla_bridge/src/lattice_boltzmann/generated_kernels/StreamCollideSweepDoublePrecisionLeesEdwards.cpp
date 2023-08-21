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
//! \\file StreamCollideSweepDoublePrecisionLeesEdwards.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

// kernel generated with pystencils v1.2, lbmpy v1.3.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 065ce5f311850371a97ac4766f47dbb5ca8424ba

#include <cmath>

#include "StreamCollideSweepDoublePrecisionLeesEdwards.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#define FUNC_PREFIX

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
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

namespace internal_654880d20ef2a4dee16898b514dd6bf5 {
static FUNC_PREFIX void streamcollidesweepdoubleprecisionleesedwards_streamcollidesweepdoubleprecisionleesedwards(double *RESTRICT _data_dst, double *RESTRICT const _data_force, double *RESTRICT const _data_pdfs, int64_t const _size_dst_0, int64_t const _size_dst_1, int64_t const _size_dst_2, int64_t const _stride_dst_0, int64_t const _stride_dst_1, int64_t const _stride_dst_2, int64_t const _stride_dst_3, int64_t const _stride_force_0, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, double grid_size, double omega_shear, double v_s) {
  const double xi_0 = ((1.0) / (omega_shear * -0.25 + 2.0));
  const double rr_0 = xi_0 * (omega_shear * -2.0 + 4.0);
  for (int64_t ctr_2 = 1; ctr_2 < _size_dst_2 - 1; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_2m1_314 = _data_pdfs + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_318 = _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_311 = _data_pdfs + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_315 = _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_313 = _data_pdfs + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_312 = _data_pdfs + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_35 = _data_pdfs + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    double *RESTRICT _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_317 = _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_316 = _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_36 = _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    double *RESTRICT _data_force_20_31 = _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    double *RESTRICT _data_force_20_32 = _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    double *RESTRICT _data_dst_20_30 = _data_dst + _stride_dst_2 * ctr_2;
    double *RESTRICT _data_dst_20_31 = _data_dst + _stride_dst_2 * ctr_2 + _stride_dst_3;
    double *RESTRICT _data_dst_20_32 = _data_dst + _stride_dst_2 * ctr_2 + 2 * _stride_dst_3;
    double *RESTRICT _data_dst_20_33 = _data_dst + _stride_dst_2 * ctr_2 + 3 * _stride_dst_3;
    double *RESTRICT _data_dst_20_34 = _data_dst + _stride_dst_2 * ctr_2 + 4 * _stride_dst_3;
    double *RESTRICT _data_dst_20_35 = _data_dst + _stride_dst_2 * ctr_2 + 5 * _stride_dst_3;
    double *RESTRICT _data_dst_20_36 = _data_dst + _stride_dst_2 * ctr_2 + 6 * _stride_dst_3;
    double *RESTRICT _data_dst_20_37 = _data_dst + _stride_dst_2 * ctr_2 + 7 * _stride_dst_3;
    double *RESTRICT _data_dst_20_38 = _data_dst + _stride_dst_2 * ctr_2 + 8 * _stride_dst_3;
    double *RESTRICT _data_dst_20_39 = _data_dst + _stride_dst_2 * ctr_2 + 9 * _stride_dst_3;
    double *RESTRICT _data_dst_20_310 = _data_dst + _stride_dst_2 * ctr_2 + 10 * _stride_dst_3;
    double *RESTRICT _data_dst_20_311 = _data_dst + _stride_dst_2 * ctr_2 + 11 * _stride_dst_3;
    double *RESTRICT _data_dst_20_312 = _data_dst + _stride_dst_2 * ctr_2 + 12 * _stride_dst_3;
    double *RESTRICT _data_dst_20_313 = _data_dst + _stride_dst_2 * ctr_2 + 13 * _stride_dst_3;
    double *RESTRICT _data_dst_20_314 = _data_dst + _stride_dst_2 * ctr_2 + 14 * _stride_dst_3;
    double *RESTRICT _data_dst_20_315 = _data_dst + _stride_dst_2 * ctr_2 + 15 * _stride_dst_3;
    double *RESTRICT _data_dst_20_316 = _data_dst + _stride_dst_2 * ctr_2 + 16 * _stride_dst_3;
    double *RESTRICT _data_dst_20_317 = _data_dst + _stride_dst_2 * ctr_2 + 17 * _stride_dst_3;
    double *RESTRICT _data_dst_20_318 = _data_dst + _stride_dst_2 * ctr_2 + 18 * _stride_dst_3;
    for (int64_t ctr_1 = 1; ctr_1 < _size_dst_1 - 1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_2m1_314_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_314;
      double *RESTRICT _data_pdfs_20_310_11 = _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_310;
      double *RESTRICT _data_pdfs_20_38_1m1 = _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_38;
      double *RESTRICT _data_pdfs_21_318_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_21_318;
      double *RESTRICT _data_pdfs_20_34_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      double *RESTRICT _data_pdfs_2m1_311_1m1 = _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_311;
      double *RESTRICT _data_pdfs_20_37_1m1 = _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_37;
      double *RESTRICT _data_pdfs_20_31_1m1 = _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_31;
      double *RESTRICT _data_pdfs_21_315_1m1 = _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_21_315;
      double *RESTRICT _data_pdfs_2m1_313_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_313;
      double *RESTRICT _data_pdfs_2m1_312_11 = _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_312;
      double *RESTRICT _data_pdfs_2m1_35_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_35;
      double *RESTRICT _data_pdfs_20_30_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      double *RESTRICT _data_pdfs_20_33_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      double *RESTRICT _data_pdfs_20_39_11 = _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_39;
      double *RESTRICT _data_pdfs_20_32_11 = _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_32;
      double *RESTRICT _data_pdfs_21_317_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_21_317;
      double *RESTRICT _data_pdfs_21_316_11 = _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_21_316;
      double *RESTRICT _data_pdfs_21_36_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_21_36;
      double *RESTRICT _data_force_20_30_10 = _stride_force_1 * ctr_1 + _data_force_20_30;
      double *RESTRICT _data_force_20_31_10 = _stride_force_1 * ctr_1 + _data_force_20_31;
      double *RESTRICT _data_force_20_32_10 = _stride_force_1 * ctr_1 + _data_force_20_32;
      double *RESTRICT _data_dst_20_30_10 = _stride_dst_1 * ctr_1 + _data_dst_20_30;
      double *RESTRICT _data_dst_20_31_10 = _stride_dst_1 * ctr_1 + _data_dst_20_31;
      double *RESTRICT _data_dst_20_32_10 = _stride_dst_1 * ctr_1 + _data_dst_20_32;
      double *RESTRICT _data_dst_20_33_10 = _stride_dst_1 * ctr_1 + _data_dst_20_33;
      double *RESTRICT _data_dst_20_34_10 = _stride_dst_1 * ctr_1 + _data_dst_20_34;
      double *RESTRICT _data_dst_20_35_10 = _stride_dst_1 * ctr_1 + _data_dst_20_35;
      double *RESTRICT _data_dst_20_36_10 = _stride_dst_1 * ctr_1 + _data_dst_20_36;
      double *RESTRICT _data_dst_20_37_10 = _stride_dst_1 * ctr_1 + _data_dst_20_37;
      double *RESTRICT _data_dst_20_38_10 = _stride_dst_1 * ctr_1 + _data_dst_20_38;
      double *RESTRICT _data_dst_20_39_10 = _stride_dst_1 * ctr_1 + _data_dst_20_39;
      double *RESTRICT _data_dst_20_310_10 = _stride_dst_1 * ctr_1 + _data_dst_20_310;
      double *RESTRICT _data_dst_20_311_10 = _stride_dst_1 * ctr_1 + _data_dst_20_311;
      double *RESTRICT _data_dst_20_312_10 = _stride_dst_1 * ctr_1 + _data_dst_20_312;
      double *RESTRICT _data_dst_20_313_10 = _stride_dst_1 * ctr_1 + _data_dst_20_313;
      double *RESTRICT _data_dst_20_314_10 = _stride_dst_1 * ctr_1 + _data_dst_20_314;
      double *RESTRICT _data_dst_20_315_10 = _stride_dst_1 * ctr_1 + _data_dst_20_315;
      double *RESTRICT _data_dst_20_316_10 = _stride_dst_1 * ctr_1 + _data_dst_20_316;
      double *RESTRICT _data_dst_20_317_10 = _stride_dst_1 * ctr_1 + _data_dst_20_317;
      double *RESTRICT _data_dst_20_318_10 = _stride_dst_1 * ctr_1 + _data_dst_20_318;
      for (int64_t ctr_0 = 1; ctr_0 < _size_dst_0 - 1; ctr_0 += 1) {
        const double vel0Term = _data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] + _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] + _data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] + _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] + _data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const double vel1Term = _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0] + _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] + _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0] + _data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0];
        const double vel2Term = _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0] + _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] + _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0];
        const double rho = vel0Term + vel1Term + vel2Term + _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0] + _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0] + _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] + _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] + _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0] + _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] + _data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0];
        const double xi_1 = ((1.0) / (rho));
        const double u_0 = xi_1 * (vel0Term - 1.0 * _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] - 1.0 * _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] - 1.0 * _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] - 1.0 * _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] - 1.0 * _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0]) + xi_1 * 0.5 * _data_force_20_30_10[_stride_force_0 * ctr_0];
        const double u_1 = xi_1 * (vel1Term - 1.0 * _data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] - 1.0 * _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0] - 1.0 * _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] - 1.0 * _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0] - 1.0 * _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0] + _data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0]) + xi_1 * 0.5 * _data_force_20_31_10[_stride_force_0 * ctr_0];
        const double u_2 = xi_1 * (vel2Term - 1.0 * _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0] - 1.0 * _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0] - 1.0 * _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] - 1.0 * _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] - 1.0 * _data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0] + _data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0] + _data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0]) + xi_1 * 0.5 * _data_force_20_32_10[_stride_force_0 * ctr_0];
        const double forceTerm_0 = omega_shear * u_0 * 0.5 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * 0.5 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * 0.5 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_0 * -1.0 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_1 * -1.0 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_2 * -1.0 * _data_force_20_32_10[_stride_force_0 * ctr_0];
        const double forceTerm_1 = omega_shear * u_0 * 0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * -0.16666666666666666 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * 0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] + rr_0 * -0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_0 * -0.16666666666666666 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_1 * 0.33333333333333331 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_2 * -0.16666666666666666 * _data_force_20_32_10[_stride_force_0 * ctr_0] + 0.16666666666666666 * _data_force_20_31_10[_stride_force_0 * ctr_0];
        const double forceTerm_2 = omega_shear * u_0 * 0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * -0.16666666666666666 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * 0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] + rr_0 * 0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_0 * -0.16666666666666666 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_1 * 0.33333333333333331 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_2 * -0.16666666666666666 * _data_force_20_32_10[_stride_force_0 * ctr_0] - 0.16666666666666666 * _data_force_20_31_10[_stride_force_0 * ctr_0];
        const double forceTerm_3 = omega_shear * u_0 * -0.16666666666666666 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * 0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * 0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] + rr_0 * 0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_0 * 0.33333333333333331 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_1 * -0.16666666666666666 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_2 * -0.16666666666666666 * _data_force_20_32_10[_stride_force_0 * ctr_0] - 0.16666666666666666 * _data_force_20_30_10[_stride_force_0 * ctr_0];
        const double forceTerm_4 = omega_shear * u_0 * -0.16666666666666666 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * 0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * 0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] + rr_0 * -0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_0 * 0.33333333333333331 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_1 * -0.16666666666666666 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_2 * -0.16666666666666666 * _data_force_20_32_10[_stride_force_0 * ctr_0] + 0.16666666666666666 * _data_force_20_30_10[_stride_force_0 * ctr_0];
        const double forceTerm_5 = omega_shear * u_0 * 0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * 0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * -0.16666666666666666 * _data_force_20_32_10[_stride_force_0 * ctr_0] + rr_0 * -0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_0 * -0.16666666666666666 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_1 * -0.16666666666666666 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_2 * 0.33333333333333331 * _data_force_20_32_10[_stride_force_0 * ctr_0] + 0.16666666666666666 * _data_force_20_32_10[_stride_force_0 * ctr_0];
        const double forceTerm_6 = omega_shear * u_0 * 0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * 0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * -0.16666666666666666 * _data_force_20_32_10[_stride_force_0 * ctr_0] + rr_0 * 0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_0 * -0.16666666666666666 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_1 * -0.16666666666666666 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_2 * 0.33333333333333331 * _data_force_20_32_10[_stride_force_0 * ctr_0] - 0.16666666666666666 * _data_force_20_32_10[_stride_force_0 * ctr_0];
        const double forceTerm_7 = omega_shear * u_0 * -0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_0 * 0.125 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * -0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * 0.125 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * 0.041666666666666664 * _data_force_20_32_10[_stride_force_0 * ctr_0] + rr_0 * -0.041666666666666664 * _data_force_20_31_10[_stride_force_0 * ctr_0] + rr_0 * 0.041666666666666664 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_0 * -0.25 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_0 * 0.16666666666666666 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_1 * -0.25 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_1 * 0.16666666666666666 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_2 * -0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] - 0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + 0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0];
        const double forceTerm_8 = omega_shear * u_0 * -0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_0 * -0.125 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * -0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * -0.125 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * 0.041666666666666664 * _data_force_20_32_10[_stride_force_0 * ctr_0] + rr_0 * -0.041666666666666664 * _data_force_20_30_10[_stride_force_0 * ctr_0] + rr_0 * -0.041666666666666664 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_0 * 0.16666666666666666 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_0 * 0.25 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_1 * 0.16666666666666666 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_1 * 0.25 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_2 * -0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] + 0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + 0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0];
        const double forceTerm_9 = omega_shear * u_0 * -0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_0 * -0.125 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * -0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * -0.125 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * 0.041666666666666664 * _data_force_20_32_10[_stride_force_0 * ctr_0] + rr_0 * 0.041666666666666664 * _data_force_20_30_10[_stride_force_0 * ctr_0] + rr_0 * 0.041666666666666664 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_0 * 0.16666666666666666 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_0 * 0.25 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_1 * 0.16666666666666666 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_1 * 0.25 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_2 * -0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] - 0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] - 0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0];
        const double forceTerm_10 = omega_shear * u_0 * -0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_0 * 0.125 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * -0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * 0.125 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * 0.041666666666666664 * _data_force_20_32_10[_stride_force_0 * ctr_0] + rr_0 * -0.041666666666666664 * _data_force_20_30_10[_stride_force_0 * ctr_0] + rr_0 * 0.041666666666666664 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_0 * -0.25 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_0 * 0.16666666666666666 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_1 * -0.25 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_1 * 0.16666666666666666 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_2 * -0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] - 0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + 0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0];
        const double forceTerm_11 = omega_shear * u_0 * 0.041666666666666664 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * -0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * -0.125 * _data_force_20_32_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * -0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * -0.125 * _data_force_20_31_10[_stride_force_0 * ctr_0] + rr_0 * -0.041666666666666664 * _data_force_20_31_10[_stride_force_0 * ctr_0] + rr_0 * -0.041666666666666664 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_0 * -0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_1 * 0.16666666666666666 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_1 * 0.25 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_2 * 0.16666666666666666 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_2 * 0.25 * _data_force_20_31_10[_stride_force_0 * ctr_0] + 0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + 0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0];
        const double forceTerm_12 = omega_shear * u_0 * 0.041666666666666664 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * -0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * 0.125 * _data_force_20_32_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * -0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * 0.125 * _data_force_20_31_10[_stride_force_0 * ctr_0] + rr_0 * -0.041666666666666664 * _data_force_20_32_10[_stride_force_0 * ctr_0] + rr_0 * 0.041666666666666664 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_0 * -0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_1 * -0.25 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_1 * 0.16666666666666666 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_2 * -0.25 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_2 * 0.16666666666666666 * _data_force_20_32_10[_stride_force_0 * ctr_0] - 0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + 0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0];
        const double forceTerm_13 = omega_shear * u_0 * -0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_0 * 0.125 * _data_force_20_32_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * 0.041666666666666664 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * -0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * 0.125 * _data_force_20_30_10[_stride_force_0 * ctr_0] + rr_0 * -0.041666666666666664 * _data_force_20_32_10[_stride_force_0 * ctr_0] + rr_0 * 0.041666666666666664 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_0 * -0.25 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_0 * 0.16666666666666666 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_1 * -0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_2 * -0.25 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_2 * 0.16666666666666666 * _data_force_20_32_10[_stride_force_0 * ctr_0] - 0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + 0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0];
        const double forceTerm_14 = omega_shear * u_0 * -0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_0 * -0.125 * _data_force_20_32_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * 0.041666666666666664 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * -0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * -0.125 * _data_force_20_30_10[_stride_force_0 * ctr_0] + rr_0 * -0.041666666666666664 * _data_force_20_30_10[_stride_force_0 * ctr_0] + rr_0 * -0.041666666666666664 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_0 * 0.16666666666666666 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_0 * 0.25 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_1 * -0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_2 * 0.16666666666666666 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_2 * 0.25 * _data_force_20_30_10[_stride_force_0 * ctr_0] + 0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + 0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0];
        const double forceTerm_15 = omega_shear * u_0 * 0.041666666666666664 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * -0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * 0.125 * _data_force_20_32_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * -0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * 0.125 * _data_force_20_31_10[_stride_force_0 * ctr_0] + rr_0 * -0.041666666666666664 * _data_force_20_31_10[_stride_force_0 * ctr_0] + rr_0 * 0.041666666666666664 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_0 * -0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_1 * -0.25 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_1 * 0.16666666666666666 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_2 * -0.25 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_2 * 0.16666666666666666 * _data_force_20_32_10[_stride_force_0 * ctr_0] - 0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] + 0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0];
        const double forceTerm_16 = omega_shear * u_0 * 0.041666666666666664 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * -0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * -0.125 * _data_force_20_32_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * -0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * -0.125 * _data_force_20_31_10[_stride_force_0 * ctr_0] + rr_0 * 0.041666666666666664 * _data_force_20_31_10[_stride_force_0 * ctr_0] + rr_0 * 0.041666666666666664 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_0 * -0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_1 * 0.16666666666666666 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_1 * 0.25 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_2 * 0.16666666666666666 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_2 * 0.25 * _data_force_20_31_10[_stride_force_0 * ctr_0] - 0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] - 0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0];
        const double forceTerm_17 = omega_shear * u_0 * -0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_0 * -0.125 * _data_force_20_32_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * 0.041666666666666664 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * -0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * -0.125 * _data_force_20_30_10[_stride_force_0 * ctr_0] + rr_0 * 0.041666666666666664 * _data_force_20_30_10[_stride_force_0 * ctr_0] + rr_0 * 0.041666666666666664 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_0 * 0.16666666666666666 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_0 * 0.25 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_1 * -0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_2 * 0.16666666666666666 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_2 * 0.25 * _data_force_20_30_10[_stride_force_0 * ctr_0] - 0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] - 0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0];
        const double forceTerm_18 = omega_shear * u_0 * -0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0] + omega_shear * u_0 * 0.125 * _data_force_20_32_10[_stride_force_0 * ctr_0] + omega_shear * u_1 * 0.041666666666666664 * _data_force_20_31_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * -0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] + omega_shear * u_2 * 0.125 * _data_force_20_30_10[_stride_force_0 * ctr_0] + rr_0 * -0.041666666666666664 * _data_force_20_30_10[_stride_force_0 * ctr_0] + rr_0 * 0.041666666666666664 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_0 * -0.25 * _data_force_20_32_10[_stride_force_0 * ctr_0] + u_0 * 0.16666666666666666 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_1 * -0.083333333333333329 * _data_force_20_31_10[_stride_force_0 * ctr_0] + u_2 * -0.25 * _data_force_20_30_10[_stride_force_0 * ctr_0] + u_2 * 0.16666666666666666 * _data_force_20_32_10[_stride_force_0 * ctr_0] - 0.083333333333333329 * _data_force_20_32_10[_stride_force_0 * ctr_0] + 0.083333333333333329 * _data_force_20_30_10[_stride_force_0 * ctr_0];
        const double u0Mu1 = u_0 + u_1 * -1.0;
        const double u0Pu1 = u_0 + u_1;
        const double u1Pu2 = u_1 + u_2;
        const double u1Mu2 = u_1 + u_2 * -1.0;
        const double u0Mu2 = u_0 + u_2 * -1.0;
        const double u0Pu2 = u_0 + u_2;
        const double f_eq_common = rho * -1.0 * (u_0 * u_0) + rho * -1.0 * (u_1 * u_1) + rho * -1.0 * (u_2 * u_2) + rho;
        _data_dst_20_30_10[_stride_dst_0 * ctr_0] = forceTerm_0 + omega_shear * (f_eq_common * 0.33333333333333331 - 1.0 * _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0]) + _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        _data_dst_20_31_10[_stride_dst_0 * ctr_0] = forceTerm_1 + omega_shear * (f_eq_common * 0.16666666666666666 + rho * (-0.1111111111111111 + 0.33333333333333331 * (u_1 * u_1)) - 0.5 * _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0] - 0.5 * _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0]) + rr_0 * (rho * u_1 * 0.16666666666666666 - 0.5 * _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0] + 0.5 * _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0]) + ((-1.0 <= grid_size * -1.0 + ((double)(ctr_1))) ? (rho * v_s * (u_0 * 2.0 + v_s) * 0.16666666666666666) : (0.0)) + _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0];
        _data_dst_20_32_10[_stride_dst_0 * ctr_0] = forceTerm_2 + omega_shear * (f_eq_common * 0.16666666666666666 + rho * (-0.1111111111111111 + 0.33333333333333331 * (u_1 * u_1)) - 0.5 * _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0] - 0.5 * _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0]) + rr_0 * (rho * u_1 * -0.16666666666666666 - 0.5 * _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0] + 0.5 * _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0]) + ((0.0 >= ((double)(ctr_1))) ? (rho * v_s * (u_0 * -2.0 + v_s) * 0.16666666666666666) : (0.0)) + _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0];
        _data_dst_20_33_10[_stride_dst_0 * ctr_0] = forceTerm_3 + omega_shear * (f_eq_common * 0.16666666666666666 + rho * (-0.1111111111111111 + 0.33333333333333331 * (u_0 * u_0)) - 0.5 * _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] - 0.5 * _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0]) + rr_0 * (rho * u_0 * -0.16666666666666666 - 0.5 * _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] + 0.5 * _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0]) + _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_dst_20_34_10[_stride_dst_0 * ctr_0] = forceTerm_4 + omega_shear * (f_eq_common * 0.16666666666666666 + rho * (-0.1111111111111111 + 0.33333333333333331 * (u_0 * u_0)) - 0.5 * _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] - 0.5 * _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0]) + rr_0 * (rho * u_0 * 0.16666666666666666 - 0.5 * _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] + 0.5 * _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0]) + _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_dst_20_35_10[_stride_dst_0 * ctr_0] = forceTerm_5 + omega_shear * (f_eq_common * 0.16666666666666666 + rho * (-0.1111111111111111 + 0.33333333333333331 * (u_2 * u_2)) - 0.5 * _data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0] - 0.5 * _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0]) + rr_0 * (rho * u_2 * 0.16666666666666666 - 0.5 * _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0] + 0.5 * _data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0]) + _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0];
        _data_dst_20_36_10[_stride_dst_0 * ctr_0] = forceTerm_6 + omega_shear * (f_eq_common * 0.16666666666666666 + rho * (-0.1111111111111111 + 0.33333333333333331 * (u_2 * u_2)) - 0.5 * _data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0] - 0.5 * _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0]) + rr_0 * (rho * u_2 * -0.16666666666666666 - 0.5 * _data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0] + 0.5 * _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0]) + _data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0];
        _data_dst_20_37_10[_stride_dst_0 * ctr_0] = forceTerm_7 + omega_shear * (f_eq_common * 0.041666666666666664 + rho * (-0.013888888888888888 + 0.041666666666666664 * (u_2 * u_2) + 0.125 * (u0Mu1 * u0Mu1)) - 0.5 * _data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] - 0.5 * _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0]) + rr_0 * (rho * u0Mu1 * -0.083333333333333329 - 0.5 * _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] + 0.5 * _data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0]) + ((-1.0 <= grid_size * -1.0 + ((double)(ctr_1))) ? (rho * v_s * (u_0 * -2.0 + u_1 * 3.0 + v_s * -1.0 + 1.0) * 0.083333333333333329) : (0.0)) + _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_dst_20_38_10[_stride_dst_0 * ctr_0] = forceTerm_8 + omega_shear * (f_eq_common * 0.041666666666666664 + rho * (-0.013888888888888888 + 0.041666666666666664 * (u_2 * u_2) + 0.125 * (u0Pu1 * u0Pu1)) - 0.5 * _data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] - 0.5 * _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0]) + rr_0 * (rho * u0Pu1 * 0.083333333333333329 - 0.5 * _data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] + 0.5 * _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0]) + ((-1.0 <= grid_size * -1.0 + ((double)(ctr_1))) ? (rho * v_s * (u_0 * 2.0 + u_1 * 3.0 + v_s + 1.0) * -0.083333333333333329) : (0.0)) + _data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_dst_20_39_10[_stride_dst_0 * ctr_0] = forceTerm_9 + omega_shear * (f_eq_common * 0.041666666666666664 + rho * (-0.013888888888888888 + 0.041666666666666664 * (u_2 * u_2) + 0.125 * (u0Pu1 * u0Pu1)) - 0.5 * _data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] - 0.5 * _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0]) + rr_0 * (rho * u0Pu1 * -0.083333333333333329 - 0.5 * _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] + 0.5 * _data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0]) + ((0.0 >= ((double)(ctr_1))) ? (rho * v_s * (u_0 * 2.0 + u_1 * 3.0 + v_s * -1.0 - 1.0) * 0.083333333333333329) : (0.0)) + _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_dst_20_310_10[_stride_dst_0 * ctr_0] = forceTerm_10 + omega_shear * (f_eq_common * 0.041666666666666664 + rho * (-0.013888888888888888 + 0.041666666666666664 * (u_2 * u_2) + 0.125 * (u0Mu1 * u0Mu1)) - 0.5 * _data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] - 0.5 * _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0]) + rr_0 * (rho * u0Mu1 * 0.083333333333333329 - 0.5 * _data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] + 0.5 * _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0]) + ((0.0 >= ((double)(ctr_1))) ? (rho * v_s * (u_0 * 2.0 + u_1 * -3.0 + v_s * -1.0 + 1.0) * 0.083333333333333329) : (0.0)) + _data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_dst_20_311_10[_stride_dst_0 * ctr_0] = forceTerm_11 + omega_shear * (f_eq_common * 0.041666666666666664 + rho * (-0.013888888888888888 + 0.041666666666666664 * (u_0 * u_0) + 0.125 * (u1Pu2 * u1Pu2)) - 0.5 * _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0] - 0.5 * _data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0]) + rr_0 * (rho * u1Pu2 * 0.083333333333333329 - 0.5 * _data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0] + 0.5 * _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0]) + _data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0];
        _data_dst_20_312_10[_stride_dst_0 * ctr_0] = forceTerm_12 + omega_shear * (f_eq_common * 0.041666666666666664 + rho * (-0.013888888888888888 + 0.041666666666666664 * (u_0 * u_0) + 0.125 * (u1Mu2 * u1Mu2)) - 0.5 * _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0] - 0.5 * _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0]) + rr_0 * (rho * u1Mu2 * -0.083333333333333329 - 0.5 * _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0] + 0.5 * _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0]) + _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0];
        _data_dst_20_313_10[_stride_dst_0 * ctr_0] = forceTerm_13 + omega_shear * (f_eq_common * 0.041666666666666664 + rho * (-0.013888888888888888 + 0.041666666666666664 * (u_1 * u_1) + 0.125 * (u0Mu2 * u0Mu2)) - 0.5 * _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] - 0.5 * _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0]) + rr_0 * (rho * u0Mu2 * -0.083333333333333329 - 0.5 * _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] + 0.5 * _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0]) + _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_dst_20_314_10[_stride_dst_0 * ctr_0] = forceTerm_14 + omega_shear * (f_eq_common * 0.041666666666666664 + rho * (-0.013888888888888888 + 0.041666666666666664 * (u_1 * u_1) + 0.125 * (u0Pu2 * u0Pu2)) - 0.5 * _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] - 0.5 * _data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0]) + rr_0 * (rho * u0Pu2 * 0.083333333333333329 - 0.5 * _data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] + 0.5 * _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0]) + _data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_dst_20_315_10[_stride_dst_0 * ctr_0] = forceTerm_15 + omega_shear * (f_eq_common * 0.041666666666666664 + rho * (-0.013888888888888888 + 0.041666666666666664 * (u_0 * u_0) + 0.125 * (u1Mu2 * u1Mu2)) - 0.5 * _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0] - 0.5 * _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0]) + rr_0 * (rho * u1Mu2 * 0.083333333333333329 - 0.5 * _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0] + 0.5 * _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0]) + _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0];
        _data_dst_20_316_10[_stride_dst_0 * ctr_0] = forceTerm_16 + omega_shear * (f_eq_common * 0.041666666666666664 + rho * (-0.013888888888888888 + 0.041666666666666664 * (u_0 * u_0) + 0.125 * (u1Pu2 * u1Pu2)) - 0.5 * _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0] - 0.5 * _data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0]) + rr_0 * (rho * u1Pu2 * -0.083333333333333329 - 0.5 * _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0] + 0.5 * _data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0]) + _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0];
        _data_dst_20_317_10[_stride_dst_0 * ctr_0] = forceTerm_17 + omega_shear * (f_eq_common * 0.041666666666666664 + rho * (-0.013888888888888888 + 0.041666666666666664 * (u_1 * u_1) + 0.125 * (u0Pu2 * u0Pu2)) - 0.5 * _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] - 0.5 * _data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0]) + rr_0 * (rho * u0Pu2 * -0.083333333333333329 - 0.5 * _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] + 0.5 * _data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0]) + _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_dst_20_318_10[_stride_dst_0 * ctr_0] = forceTerm_18 + omega_shear * (f_eq_common * 0.041666666666666664 + rho * (-0.013888888888888888 + 0.041666666666666664 * (u_1 * u_1) + 0.125 * (u0Mu2 * u0Mu2)) - 0.5 * _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] - 0.5 * _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0]) + rr_0 * (rho * u0Mu2 * 0.083333333333333329 - 0.5 * _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] + 0.5 * _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0]) + _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
      }
    }
  }
}
} // namespace internal_654880d20ef2a4dee16898b514dd6bf5

void StreamCollideSweepDoublePrecisionLeesEdwards::run(IBlock *block) {
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);
  auto dst = block->getData<field::GhostLayerField<double, 19>>(dstID);

  auto &grid_size = this->grid_size_;
  auto &v_s = this->v_s_;
  auto &omega_shear = this->omega_shear_;
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(dst->nrOfGhostLayers()));
  double *RESTRICT _data_dst = dst->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(dst->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force->nrOfGhostLayers()));
  double *RESTRICT const _data_force = force->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs->nrOfGhostLayers()));
  double *RESTRICT const _data_pdfs = pdfs->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(dst->xSizeWithGhostLayer(), int64_t(cell_idx_c(dst->xSize()) + 2));
  const int64_t _size_dst_0 = int64_t(cell_idx_c(dst->xSize()) + 2);
  WALBERLA_ASSERT_EQUAL(dst->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(dst->ySizeWithGhostLayer(), int64_t(cell_idx_c(dst->ySize()) + 2));
  const int64_t _size_dst_1 = int64_t(cell_idx_c(dst->ySize()) + 2);
  WALBERLA_ASSERT_EQUAL(dst->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(dst->zSizeWithGhostLayer(), int64_t(cell_idx_c(dst->zSize()) + 2));
  const int64_t _size_dst_2 = int64_t(cell_idx_c(dst->zSize()) + 2);
  WALBERLA_ASSERT_EQUAL(dst->layout(), field::fzyx);
  const int64_t _stride_dst_0 = int64_t(dst->xStride());
  const int64_t _stride_dst_1 = int64_t(dst->yStride());
  const int64_t _stride_dst_2 = int64_t(dst->zStride());
  const int64_t _stride_dst_3 = int64_t(1 * int64_t(dst->fStride()));
  const int64_t _stride_force_0 = int64_t(force->xStride());
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  internal_654880d20ef2a4dee16898b514dd6bf5::streamcollidesweepdoubleprecisionleesedwards_streamcollidesweepdoubleprecisionleesedwards(_data_dst, _data_force, _data_pdfs, _size_dst_0, _size_dst_1, _size_dst_2, _stride_dst_0, _stride_dst_1, _stride_dst_2, _stride_dst_3, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, grid_size, omega_shear, v_s);
}

void StreamCollideSweepDoublePrecisionLeesEdwards::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block) {
  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);
  auto dst = block->getData<field::GhostLayerField<double, 19>>(dstID);

  auto &grid_size = this->grid_size_;
  auto &v_s = this->v_s_;
  auto &omega_shear = this->omega_shear_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(dst->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(dst->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(dst->nrOfGhostLayers()));
  double *RESTRICT _data_dst = dst->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(dst->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force->nrOfGhostLayers()));
  double *RESTRICT const _data_force = force->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs->nrOfGhostLayers()));
  double *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(dst->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 2));
  const int64_t _size_dst_0 = int64_t(cell_idx_c(ci.xSize()) + 2);
  WALBERLA_ASSERT_EQUAL(dst->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(dst->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 2));
  const int64_t _size_dst_1 = int64_t(cell_idx_c(ci.ySize()) + 2);
  WALBERLA_ASSERT_EQUAL(dst->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(dst->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 2));
  const int64_t _size_dst_2 = int64_t(cell_idx_c(ci.zSize()) + 2);
  WALBERLA_ASSERT_EQUAL(dst->layout(), field::fzyx);
  const int64_t _stride_dst_0 = int64_t(dst->xStride());
  const int64_t _stride_dst_1 = int64_t(dst->yStride());
  const int64_t _stride_dst_2 = int64_t(dst->zStride());
  const int64_t _stride_dst_3 = int64_t(1 * int64_t(dst->fStride()));
  const int64_t _stride_force_0 = int64_t(force->xStride());
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  internal_654880d20ef2a4dee16898b514dd6bf5::streamcollidesweepdoubleprecisionleesedwards_streamcollidesweepdoubleprecisionleesedwards(_data_dst, _data_force, _data_pdfs, _size_dst_0, _size_dst_1, _size_dst_2, _stride_dst_0, _stride_dst_1, _stride_dst_2, _stride_dst_3, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, grid_size, omega_shear, v_s);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif