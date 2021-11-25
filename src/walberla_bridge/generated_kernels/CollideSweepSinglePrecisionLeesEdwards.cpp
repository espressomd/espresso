// kernel generated with pystencils v0.4.3, lbmpy v0.4.3,
// lbmpy_walberla/pystencils_walberla from commit
// 88f85eb7a979f81d68e76009811aeed53ec3014e

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
//! \\file CollideSweepSinglePrecisionLeesEdwards.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "CollideSweepSinglePrecisionLeesEdwards.h"
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

namespace internal_ab1f3bc3368574afb482da84ccb58898 {
static FUNC_PREFIX void
collidesweepsingleprecisionleesedwards_collidesweepsingleprecisionleesedwards(
    float *RESTRICT const _data_force, float *RESTRICT _data_pdfs,
    float *RESTRICT const _data_velocity, int64_t const _size_force_0,
    int64_t const _size_force_1, int64_t const _size_force_2,
    int64_t const _stride_force_0, int64_t const _stride_force_1,
    int64_t const _stride_force_2, int64_t const _stride_force_3,
    int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1,
    int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3,
    int64_t const _stride_velocity_0, int64_t const _stride_velocity_1,
    int64_t const _stride_velocity_2, int64_t const _stride_velocity_3,
    float omega_bulk, bool points_down, bool points_up) {
  const float xi_1 = omega_bulk * -0.5f + 1.0f;
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    float *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    float *RESTRICT _data_velocity_20_30 =
        _data_velocity + _stride_velocity_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    float *RESTRICT _data_velocity_20_31 =
        _data_velocity + _stride_velocity_2 * ctr_2 + _stride_velocity_3;
    float *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    float *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    float *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    float *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    float *RESTRICT _data_velocity_20_32 =
        _data_velocity + _stride_velocity_2 * ctr_2 + 2 * _stride_velocity_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      float *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      float *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      float *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      float *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      float *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      float *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      float *RESTRICT _data_velocity_20_30_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_30;
      float *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      float *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      float *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      float *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      float *RESTRICT _data_velocity_20_31_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_31;
      float *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      float *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      float *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      float *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      float *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      float *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      float *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      float *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      float *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      float *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      float *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      float *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      float *RESTRICT _data_velocity_20_32_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_32;
      for (int64_t ctr_0 = 0; ctr_0 < _size_force_0; ctr_0 += 1) {
        const float xi_63 = _data_force_20_30_10[_stride_force_0 * ctr_0];
        const float xi_64 = _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0];
        const float xi_65 = _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0];
        const float xi_66 = _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0];
        const float xi_67 = _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0];
        const float xi_68 = _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0];
        const float xi_69 = _data_velocity_20_30_10[_stride_velocity_0 * ctr_0];
        const float xi_70 = _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0];
        const float xi_71 = _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0];
        const float xi_72 = _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const float xi_73 = _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0];
        const float xi_74 = _data_velocity_20_31_10[_stride_velocity_0 * ctr_0];
        const float xi_75 = _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0];
        const float xi_76 = _data_force_20_31_10[_stride_force_0 * ctr_0];
        const float xi_77 = _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0];
        const float xi_78 = _data_force_20_32_10[_stride_force_0 * ctr_0];
        const float xi_79 = _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0];
        const float xi_80 = _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0];
        const float xi_81 = _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0];
        const float xi_82 = _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0];
        const float xi_83 = _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0];
        const float xi_84 = _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0];
        const float xi_85 = _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0];
        const float xi_86 = _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0];
        const float xi_87 = _data_velocity_20_32_10[_stride_velocity_0 * ctr_0];
        const float xi_3 = xi_74 * xi_76;
        const float xi_4 = xi_78 * xi_87;
        const float xi_5 = xi_74 * 2.0f;
        const float xi_6 = xi_5 + 1.0f;
        const float xi_7 = xi_76 * 0.166666666666667f;
        const float xi_9 = xi_4 * -0.166666666666667f;
        const float xi_11 = xi_5 - 1.0f;
        const float xi_14 = xi_63 * 0.166666666666667f;
        const float xi_15 = xi_3 * -0.166666666666667f;
        const float xi_16 = xi_15 + xi_9;
        const float xi_18 = xi_87 * 2.0f;
        const float xi_19 = xi_18 + 1.0f;
        const float xi_20 = xi_78 * 0.166666666666667f;
        const float xi_22 = xi_18 - 1.0f;
        const float xi_23 = xi_4 * -0.0833333333333333f;
        const float xi_24 = xi_74 * 3.0f;
        const float xi_26 = xi_63 * 0.0833333333333333f;
        const float xi_29 = xi_76 * 0.0833333333333333f;
        const float xi_30 = -xi_24;
        const float xi_31 = -xi_5 + 1.0f;
        const float xi_33 = xi_87 * 3.0f;
        const float xi_34 = xi_78 * 0.0833333333333333f;
        const float xi_35 = xi_3 * -0.0833333333333333f;
        const float xi_36 = -xi_33;
        const float xi_37 = -xi_18 + 1.0f;
        const float xi_38 = -xi_87;
        const float xi_39 = (xi_74 * xi_74);
        const float xi_41 = (xi_87 * xi_87);
        const float xi_44 = xi_74 * 0.166666666666667f;
        const float xi_45 = xi_39 * 0.25f;
        const float xi_48 = xi_87 * 0.166666666666667f;
        const float xi_49 = xi_41 * 0.25f;
        const float rho = xi_64 + xi_65 + xi_66 + xi_67 + xi_68 + xi_70 +
                          xi_71 + xi_72 + xi_73 + xi_75 + xi_77 + xi_79 +
                          xi_80 + xi_81 + xi_82 + xi_83 + xi_84 + xi_85 + xi_86;
        const float xi_40 = rho * 1.5f;
        const float u_0 =
            xi_69 + ((points_down && ctr_1 <= 0)
                         ? (1.0f)
                         : ((points_up && ctr_1 >= 63) ? (-1.0f) : (0.0f))) *
                        0.050000000000000003f;
        const float xi_2 = u_0 * xi_63;
        const float xi_8 = xi_2 * -0.166666666666667f;
        const float xi_10 = xi_8 + xi_9;
        const float xi_12 = u_0 * 2.0f;
        const float xi_13 = xi_12 - 1.0f;
        const float xi_17 = xi_12 + 1.0f;
        const float xi_21 = xi_15 + xi_8;
        const float xi_25 = -xi_12 + 1.0f;
        const float xi_27 = u_0 * 3.0f;
        const float xi_28 = -xi_27;
        const float xi_32 = xi_2 * -0.0833333333333333f;
        const float xi_42 = (u_0 * u_0);
        const float xi_46 = u_0 * 0.166666666666667f;
        const float xi_47 = xi_42 * 0.25f;
        const float forceTerm_0 = xi_1 * (-xi_2 - xi_3 - xi_4);
        const float forceTerm_1 = xi_1 * (xi_10 + xi_6 * xi_7);
        const float forceTerm_2 = xi_1 * (xi_10 + xi_11 * xi_7);
        const float forceTerm_3 = xi_1 * (xi_13 * xi_14 + xi_16);
        const float forceTerm_4 = xi_1 * (xi_14 * xi_17 + xi_16);
        const float forceTerm_5 = xi_1 * (xi_19 * xi_20 + xi_21);
        const float forceTerm_6 = xi_1 * (xi_20 * xi_22 + xi_21);
        const float forceTerm_7 =
            xi_1 * (xi_23 - xi_26 * (xi_24 + xi_25) + xi_29 * (xi_28 + xi_6));
        const float forceTerm_8 =
            xi_1 * (xi_23 + xi_26 * (xi_17 + xi_24) + xi_29 * (xi_27 + xi_6));
        const float forceTerm_9 =
            xi_1 * (xi_23 + xi_26 * (xi_13 + xi_24) + xi_29 * (xi_11 + xi_27));
        const float forceTerm_10 =
            xi_1 * (xi_23 + xi_26 * (xi_17 + xi_30) - xi_29 * (xi_27 + xi_31));
        const float forceTerm_11 =
            xi_1 * (xi_29 * (xi_33 + xi_6) + xi_32 + xi_34 * (xi_19 + xi_24));
        const float forceTerm_12 =
            xi_1 * (-xi_29 * (xi_31 + xi_33) + xi_32 + xi_34 * (xi_19 + xi_30));
        const float forceTerm_13 =
            xi_1 * (-xi_26 * (xi_25 + xi_33) + xi_34 * (xi_19 + xi_28) + xi_35);
        const float forceTerm_14 =
            xi_1 * (xi_26 * (xi_17 + xi_33) + xi_34 * (xi_19 + xi_27) + xi_35);
        const float forceTerm_15 =
            xi_1 * (xi_29 * (xi_36 + xi_6) + xi_32 - xi_34 * (xi_24 + xi_37));
        const float forceTerm_16 =
            xi_1 * (xi_29 * (xi_11 + xi_33) + xi_32 + xi_34 * (xi_22 + xi_24));
        const float forceTerm_17 =
            xi_1 * (xi_26 * (xi_13 + xi_33) + xi_34 * (xi_22 + xi_27) + xi_35);
        const float forceTerm_18 =
            xi_1 * (xi_26 * (xi_17 + xi_36) - xi_34 * (xi_27 + xi_37) + xi_35);
        const float u0Mu1 = u_0 - xi_74;
        const float xi_51 = u0Mu1 * 0.0833333333333333f;
        const float xi_52 = (u0Mu1 * u0Mu1) * 0.125f;
        const float u0Pu1 = u_0 + xi_74;
        const float xi_53 = u0Pu1 * 0.0833333333333333f;
        const float xi_54 = (u0Pu1 * u0Pu1) * 0.125f;
        const float u1Pu2 = xi_74 + xi_87;
        const float xi_55 = u1Pu2 * 0.0833333333333333f;
        const float xi_56 = (u1Pu2 * u1Pu2) * 0.125f;
        const float u1Mu2 = xi_38 + xi_74;
        const float xi_57 = u1Mu2 * 0.0833333333333333f;
        const float xi_58 = (u1Mu2 * u1Mu2) * 0.125f;
        const float u0Mu2 = u_0 + xi_38;
        const float xi_59 = u0Mu2 * 0.0833333333333333f;
        const float xi_60 = (u0Mu2 * u0Mu2) * 0.125f;
        const float u0Pu2 = u_0 + xi_87;
        const float xi_61 = u0Pu2 * 0.0833333333333333f;
        const float xi_62 = (u0Pu2 * u0Pu2) * 0.125f;
        const float f_eq_common =
            rho - xi_39 * xi_40 - xi_40 * xi_41 - xi_40 * xi_42;
        const float xi_43 = f_eq_common * 0.0555555555555556f;
        const float xi_50 = f_eq_common * 0.0277777777777778f;
        _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_0 +
            omega_bulk * (f_eq_common * 0.333333333333333f - xi_72) + xi_72;
        _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_1 + omega_bulk * (rho * (xi_44 + xi_45) + xi_43 - xi_77) +
            xi_77;
        _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_2 +
            omega_bulk * (rho * (-xi_44 + xi_45) + xi_43 - xi_83) + xi_83;
        _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_3 +
            omega_bulk * (rho * (-xi_46 + xi_47) + xi_43 - xi_73) + xi_73;
        _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_4 + omega_bulk * (rho * (xi_46 + xi_47) + xi_43 - xi_80) +
            xi_80;
        _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_5 + omega_bulk * (rho * (xi_48 + xi_49) + xi_43 - xi_71) +
            xi_71;
        _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_6 +
            omega_bulk * (rho * (-xi_48 + xi_49) + xi_43 - xi_85) + xi_85;
        _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_7 +
            omega_bulk * (rho * (-xi_51 + xi_52) + xi_50 - xi_68) + xi_68;
        _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_8 + omega_bulk * (rho * (xi_53 + xi_54) + xi_50 - xi_79) +
            xi_79;
        _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_9 +
            omega_bulk * (rho * (-xi_53 + xi_54) + xi_50 - xi_82) + xi_82;
        _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_10 +
            omega_bulk * (rho * (xi_51 + xi_52) + xi_50 - xi_65) + xi_65;
        _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_11 +
            omega_bulk * (rho * (xi_55 + xi_56) + xi_50 - xi_66) + xi_66;
        _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_12 +
            omega_bulk * (rho * (-xi_57 + xi_58) + xi_50 - xi_67) + xi_67;
        _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_13 +
            omega_bulk * (rho * (-xi_59 + xi_60) + xi_50 - xi_70) + xi_70;
        _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_14 +
            omega_bulk * (rho * (xi_61 + xi_62) + xi_50 - xi_86) + xi_86;
        _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_15 +
            omega_bulk * (rho * (xi_57 + xi_58) + xi_50 - xi_64) + xi_64;
        _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_16 +
            omega_bulk * (rho * (-xi_55 + xi_56) + xi_50 - xi_81) + xi_81;
        _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_17 +
            omega_bulk * (rho * (-xi_61 + xi_62) + xi_50 - xi_84) + xi_84;
        _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_18 +
            omega_bulk * (rho * (xi_59 + xi_60) + xi_50 - xi_75) + xi_75;
      }
    }
  }
}
} // namespace internal_ab1f3bc3368574afb482da84ccb58898

void CollideSweepSinglePrecisionLeesEdwards::run(IBlock *block) {
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto velocity = block->getData<field::GhostLayerField<float, 3>>(velocityID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);

  auto &points_down = this->points_down_;
  auto &points_up = this->points_up_;
  auto &omega_bulk = this->omega_bulk_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(velocity->nrOfGhostLayers()));
  float *RESTRICT const _data_velocity = velocity->dataAt(0, 0, 0, 0);
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
  internal_ab1f3bc3368574afb482da84ccb58898::
      collidesweepsingleprecisionleesedwards_collidesweepsingleprecisionleesedwards(
          _data_force, _data_pdfs, _data_velocity, _size_force_0, _size_force_1,
          _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2,
          _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2,
          _stride_pdfs_3, _stride_velocity_0, _stride_velocity_1,
          _stride_velocity_2, _stride_velocity_3, omega_bulk, points_down,
          points_up);
}

void CollideSweepSinglePrecisionLeesEdwards::runOnCellInterval(
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

  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto velocity = block->getData<field::GhostLayerField<float, 3>>(velocityID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);

  auto &points_down = this->points_down_;
  auto &points_up = this->points_up_;
  auto &omega_bulk = this->omega_bulk_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force =
      force->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(velocity->nrOfGhostLayers()));
  float *RESTRICT const _data_velocity =
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
  internal_ab1f3bc3368574afb482da84ccb58898::
      collidesweepsingleprecisionleesedwards_collidesweepsingleprecisionleesedwards(
          _data_force, _data_pdfs, _data_velocity, _size_force_0, _size_force_1,
          _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2,
          _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2,
          _stride_pdfs_3, _stride_velocity_0, _stride_velocity_1,
          _stride_velocity_2, _stride_velocity_3, omega_bulk, points_down,
          points_up);
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