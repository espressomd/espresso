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
//! \\author Martin Bauer <martin.bauer@fau.de>
//======================================================================================================================

#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "lbm/field/PdfField.h"
#include "lbm/sweeps/Streaming.h"
#include "FluctuatingMRT_LatticeModel.h"

#ifdef _MSC_VER
#  pragma warning( disable : 4458 )
#endif

#define FUNC_PREFIX

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wfloat-equal"
#   pragma GCC diagnostic ignored "-Wshadow"
#   pragma GCC diagnostic ignored "-Wconversion"
#   pragma GCC diagnostic ignored "-Wunused-variable"
#   pragma GCC diagnostic ignored "-Wunused-parameter"
#endif


#include "philox_rand.h"



using namespace std;

namespace walberla {
namespace lbm {

namespace internal_kernel_streamCollide {
static FUNC_PREFIX void kernel_streamCollide(double * RESTRICT const _data_force, double * RESTRICT const _data_pdfs, double * RESTRICT _data_pdfs_tmp, int64_t const _size_force_0, int64_t const _size_force_1, int64_t const _size_force_2, int64_t const _stride_force_0, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_pdfs_tmp_0, int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2, int64_t const _stride_pdfs_tmp_3, uint32_t block_offset_0, uint32_t block_offset_1, uint32_t block_offset_2, double omega_bulk, double omega_even, double omega_odd, double omega_shear, uint32_t seed, double temperature, uint32_t time_step)
{
   const double xi_26 = omega_shear*-0.5 + 1.0;
   const double xi_64 = 2.4494897427831779;
   const double xi_89 = omega_odd*0.25;
   const double xi_105 = omega_odd*0.0833333333333333;
   const double xi_170 = omega_shear*0.25;
   const double xi_185 = omega_odd*0.0416666666666667;
   const double xi_187 = omega_odd*0.125;
   const int64_t rr_0 = 0.0;
   const double xi_94 = rr_0*0.166666666666667;
   const double xi_160 = rr_0*0.0833333333333333;
   for (int ctr_2 = 1; ctr_2 < _size_force_2 - 1; ctr_2 += 1)
   {
      double * RESTRICT _data_pdfs_2m1_314 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 14*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_318 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 18*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_311 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 11*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_315 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 15*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_312 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 12*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_35 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 5*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_316 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_317 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_21_36 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 6*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_2m1_313 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 13*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      double * RESTRICT _data_force_20_31 = _data_force + _stride_force_2*ctr_2 + _stride_force_3;
      double * RESTRICT _data_force_20_30 = _data_force + _stride_force_2*ctr_2;
      double * RESTRICT _data_force_20_32 = _data_force + _stride_force_2*ctr_2 + 2*_stride_force_3;
      double * RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2*ctr_2;
      double * RESTRICT _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_30 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2;
      double * RESTRICT _data_pdfs_tmp_20_31 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + _stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_32 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 2*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_33 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 3*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_34 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 4*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_35 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 5*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_36 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 6*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_37 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 7*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_38 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 8*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_39 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 9*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_310 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 10*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_311 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 11*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_312 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 12*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_313 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 13*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_314 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 14*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_315 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 15*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_316 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 16*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_317 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 17*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_tmp_20_318 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 18*_stride_pdfs_tmp_3;
      for (int ctr_1 = 1; ctr_1 < _size_force_1 - 1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_2m1_314_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_314;
         double * RESTRICT _data_pdfs_21_318_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_318;
         double * RESTRICT _data_pdfs_20_34_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_34;
         double * RESTRICT _data_pdfs_2m1_311_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_311;
         double * RESTRICT _data_pdfs_20_31_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_31;
         double * RESTRICT _data_pdfs_21_315_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_21_315;
         double * RESTRICT _data_pdfs_2m1_312_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_312;
         double * RESTRICT _data_pdfs_2m1_35_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_35;
         double * RESTRICT _data_pdfs_20_33_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_33;
         double * RESTRICT _data_pdfs_20_39_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_39;
         double * RESTRICT _data_pdfs_20_32_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_32;
         double * RESTRICT _data_pdfs_21_316_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_21_316;
         double * RESTRICT _data_pdfs_21_317_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_317;
         double * RESTRICT _data_pdfs_21_36_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_36;
         double * RESTRICT _data_pdfs_20_37_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_37;
         double * RESTRICT _data_pdfs_2m1_313_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_313;
         double * RESTRICT _data_pdfs_20_310_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_310;
         double * RESTRICT _data_force_20_31_10 = _stride_force_1*ctr_1 + _data_force_20_31;
         double * RESTRICT _data_force_20_30_10 = _stride_force_1*ctr_1 + _data_force_20_30;
         double * RESTRICT _data_force_20_32_10 = _stride_force_1*ctr_1 + _data_force_20_32;
         double * RESTRICT _data_pdfs_20_30_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_30;
         double * RESTRICT _data_pdfs_20_38_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_38;
         double * RESTRICT _data_pdfs_tmp_20_30_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_30;
         double * RESTRICT _data_pdfs_tmp_20_31_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_31;
         double * RESTRICT _data_pdfs_tmp_20_32_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_32;
         double * RESTRICT _data_pdfs_tmp_20_33_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_33;
         double * RESTRICT _data_pdfs_tmp_20_34_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_34;
         double * RESTRICT _data_pdfs_tmp_20_35_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_35;
         double * RESTRICT _data_pdfs_tmp_20_36_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_36;
         double * RESTRICT _data_pdfs_tmp_20_37_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_37;
         double * RESTRICT _data_pdfs_tmp_20_38_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_38;
         double * RESTRICT _data_pdfs_tmp_20_39_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_39;
         double * RESTRICT _data_pdfs_tmp_20_310_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_310;
         double * RESTRICT _data_pdfs_tmp_20_311_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_311;
         double * RESTRICT _data_pdfs_tmp_20_312_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_312;
         double * RESTRICT _data_pdfs_tmp_20_313_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_313;
         double * RESTRICT _data_pdfs_tmp_20_314_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_314;
         double * RESTRICT _data_pdfs_tmp_20_315_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_315;
         double * RESTRICT _data_pdfs_tmp_20_316_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_316;
         double * RESTRICT _data_pdfs_tmp_20_317_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_317;
         double * RESTRICT _data_pdfs_tmp_20_318_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_318;
         for (int ctr_0 = 1; ctr_0 < _size_force_0 - 1; ctr_0 += 1)
         {
            
            float Dummy_292;
            float Dummy_293;
            float Dummy_294;
            float Dummy_295;
            philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 3, seed, Dummy_292, Dummy_293, Dummy_294, Dummy_295);
            
            
            float Dummy_288;
            float Dummy_289;
            float Dummy_290;
            float Dummy_291;
            philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 2, seed, Dummy_288, Dummy_289, Dummy_290, Dummy_291);
            
            
            float Dummy_284;
            float Dummy_285;
            float Dummy_286;
            float Dummy_287;
            philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 1, seed, Dummy_284, Dummy_285, Dummy_286, Dummy_287);
            
            
            float Dummy_280;
            float Dummy_281;
            float Dummy_282;
            float Dummy_283;
            philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 0, seed, Dummy_280, Dummy_281, Dummy_282, Dummy_283);
            
            const double xi_0 = _data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double xi_1 = xi_0 + _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double xi_2 = _data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0] + _data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0] + _data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0];
            const double xi_3 = _data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0] + _data_pdfs_2m1_35_10[_stride_pdfs_0*ctr_0];
            const double xi_4 = _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double xi_5 = _data_pdfs_20_32_11[_stride_pdfs_0*ctr_0] + _data_pdfs_21_316_11[_stride_pdfs_0*ctr_0];
            const double xi_6 = _data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_21_36_10[_stride_pdfs_0*ctr_0];
            const double xi_9 = -_data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double xi_10 = xi_9 - _data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double xi_11 = -_data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double xi_12 = -_data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double xi_13 = -_data_pdfs_20_33_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double xi_14 = xi_11 + xi_12 + xi_13;
            const double xi_15 = -_data_pdfs_20_32_11[_stride_pdfs_0*ctr_0];
            const double xi_16 = -_data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double xi_17 = xi_15 + xi_16;
            const double xi_18 = -_data_pdfs_21_316_11[_stride_pdfs_0*ctr_0];
            const double xi_19 = -_data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0];
            const double xi_20 = xi_18 + xi_19;
            const double xi_21 = -_data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double xi_22 = xi_11 + xi_21;
            const double xi_23 = -_data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0];
            const double xi_24 = -_data_pdfs_21_36_10[_stride_pdfs_0*ctr_0];
            const double xi_25 = xi_18 + xi_23 + xi_24 + _data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0];
            const double xi_32 = 0.166666666666667*_data_force_20_31_10[_stride_force_0*ctr_0];
            const double xi_39 = 0.166666666666667*_data_force_20_30_10[_stride_force_0*ctr_0];
            const double xi_45 = 0.166666666666667*_data_force_20_32_10[_stride_force_0*ctr_0];
            const double xi_51 = 0.0833333333333333*_data_force_20_30_10[_stride_force_0*ctr_0];
            const double xi_54 = 0.0833333333333333*_data_force_20_31_10[_stride_force_0*ctr_0];
            const double xi_57 = 0.0833333333333333*_data_force_20_32_10[_stride_force_0*ctr_0];
            const double xi_67 = -_data_pdfs_20_30_10[_stride_pdfs_0*ctr_0];
            const double xi_68 = xi_67 + 3.0*_data_pdfs_21_36_10[_stride_pdfs_0*ctr_0] + 3.0*_data_pdfs_2m1_35_10[_stride_pdfs_0*ctr_0];
            const double xi_69 = omega_even*(xi_68 - 3.0*_data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0] - 3.0*_data_pdfs_21_316_11[_stride_pdfs_0*ctr_0] - 3.0*_data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0] - 3.0*_data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0] + 3.0*_data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0] + 3.0*_data_pdfs_20_32_11[_stride_pdfs_0*ctr_0]);
            const double xi_70 = 2.0*_data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0] + 2.0*_data_pdfs_21_316_11[_stride_pdfs_0*ctr_0] + 2.0*_data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0] + 2.0*_data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0];
            const double xi_71 = xi_70 + 5.0*_data_pdfs_20_33_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + 5.0*_data_pdfs_20_34_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double xi_72 = omega_even*(xi_68 + xi_71 - 2.0*_data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0] - 2.0*_data_pdfs_20_32_11[_stride_pdfs_0*ctr_0] - 5.0*_data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - 5.0*_data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] - 5.0*_data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - 5.0*_data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0]);
            const double xi_75 = -_data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0];
            const double xi_76 = xi_19 + xi_75;
            const double xi_77 = -_data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double xi_80 = -_data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double xi_81 = xi_12 + xi_16 + xi_22 + xi_80;
            const double xi_83 = 2.0*_data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double xi_84 = 2.0*_data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double xi_85 = 2.0*_data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + 2.0*_data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double xi_86 = omega_even*(xi_67 + xi_71 + xi_83 + xi_84 + xi_85 - 4.0*_data_pdfs_21_36_10[_stride_pdfs_0*ctr_0] - 4.0*_data_pdfs_2m1_35_10[_stride_pdfs_0*ctr_0] - 7.0*_data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] - 7.0*_data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - 7.0*_data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] - 7.0*_data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + 5.0*_data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0] + 5.0*_data_pdfs_20_32_11[_stride_pdfs_0*ctr_0]);
            const double xi_87 = xi_75 + _data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0];
            const double xi_88 = xi_15 + xi_23 + xi_87 + _data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0] + _data_pdfs_21_316_11[_stride_pdfs_0*ctr_0];
            const double xi_90 = xi_88*xi_89;
            const double xi_92 = xi_77 + _data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double xi_96 = Dummy_291 - 0.5;
            const double xi_101 = 2.0*_data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double xi_102 = 2.0*_data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double xi_103 = -2.0*_data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + 2.0*_data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double xi_104 = -xi_101 + xi_102 + xi_103 + xi_15 + xi_2 + xi_20;
            const double xi_106 = xi_104*xi_105;
            const double xi_107 = Dummy_286 - 0.5;
            const double xi_112 = Dummy_281 - 0.5;
            const double xi_116 = _data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double xi_130 = xi_80 + _data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double xi_131 = xi_13 + xi_130 + xi_21 + _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double xi_132 = xi_131*xi_89;
            const double xi_133 = Dummy_289 - 0.5;
            const double xi_135 = xi_1 + xi_101 - xi_102 + xi_103 + xi_14;
            const double xi_136 = xi_105*xi_135;
            const double xi_137 = Dummy_288 - 0.5;
            const double xi_142 = _data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0] + _data_pdfs_21_316_11[_stride_pdfs_0*ctr_0];
            const double xi_143 = xi_142 + xi_24 + xi_76 + _data_pdfs_2m1_35_10[_stride_pdfs_0*ctr_0];
            const double xi_144 = xi_143*xi_89;
            const double xi_147 = Dummy_290 - 0.5;
            const double xi_149 = xi_25 + xi_3 - xi_83 - xi_84 + xi_85;
            const double xi_150 = xi_105*xi_149;
            const double xi_151 = Dummy_287 - 0.5;
            const double xi_158 = xi_86*0.0138888888888889;
            const double xi_179 = xi_72*-0.00714285714285714;
            const double xi_181 = xi_69*0.025;
            const double xi_186 = xi_149*xi_185;
            const double xi_188 = xi_143*xi_187;
            const double xi_197 = xi_104*xi_185;
            const double xi_198 = xi_187*xi_88;
            const double xi_206 = xi_72*0.0178571428571429;
            const double xi_212 = xi_131*xi_187;
            const double xi_213 = xi_135*xi_185;
            const double vel0Term = xi_1 + _data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double vel1Term = xi_2 + _data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double vel2Term = xi_3 + _data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double rho = vel0Term + vel1Term + vel2Term + xi_4 + xi_5 + xi_6 + _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0];
            const double xi_7 = 1 / (rho);
            const double xi_8 = xi_7*0.5;
            const double xi_60 = rho*temperature;
            const double xi_61 = sqrt(xi_60*(-((-omega_even + 1.0)*(-omega_even + 1.0)) + 1.0));
            const double xi_62 = xi_61*(Dummy_292 - 0.5)*3.7416573867739413;
            const double xi_63 = xi_61*(Dummy_294 - 0.5)*5.4772255750516612;
            const double xi_65 = xi_64*sqrt(xi_60*(-((-omega_bulk + 1.0)*(-omega_bulk + 1.0)) + 1.0))*(Dummy_285 - 0.5);
            const double xi_66 = xi_61*(Dummy_293 - 0.5)*8.3666002653407556;
            const double xi_97 = sqrt(xi_60*(-((-omega_odd + 1.0)*(-omega_odd + 1.0)) + 1.0));
            const double xi_98 = xi_97*1.4142135623730951;
            const double xi_99 = xi_98*0.5;
            const double xi_100 = xi_96*xi_99;
            const double xi_108 = xi_64*xi_97;
            const double xi_109 = xi_108*0.166666666666667;
            const double xi_110 = xi_107*xi_109;
            const double xi_111 = -xi_106 - xi_110;
            const double xi_113 = sqrt(xi_60*(-((-omega_shear + 1.0)*(-omega_shear + 1.0)) + 1.0));
            const double xi_114 = xi_113*0.5;
            const double xi_115 = xi_112*xi_114;
            const double xi_120 = xi_62*-0.119047619047619 + xi_86*-0.0198412698412698;
            const double xi_122 = xi_113*(Dummy_280 - 0.5)*1.7320508075688772;
            const double xi_126 = xi_106 + xi_110;
            const double xi_134 = xi_133*xi_99;
            const double xi_138 = xi_109*xi_137;
            const double xi_139 = xi_136 + xi_138;
            const double xi_141 = -xi_136 - xi_138;
            const double xi_148 = xi_147*xi_99;
            const double xi_152 = xi_109*xi_151;
            const double xi_153 = -xi_150 - xi_152;
            const double xi_155 = xi_150 + xi_152;
            const double xi_156 = xi_112*xi_113*0.25;
            const double xi_159 = xi_62*0.0833333333333333;
            const double xi_169 = xi_114*(Dummy_282 - 0.5);
            const double xi_178 = xi_114*(Dummy_284 - 0.5);
            const double xi_182 = xi_66*-0.0142857142857143;
            const double xi_183 = xi_63*0.05;
            const double xi_189 = xi_108*0.0833333333333333;
            const double xi_190 = xi_151*xi_189;
            const double xi_191 = xi_98*0.25;
            const double xi_192 = xi_147*xi_191;
            const double xi_194 = xi_62*-0.0238095238095238 + xi_86*-0.00396825396825397;
            const double xi_199 = xi_107*xi_189;
            const double xi_200 = xi_191*xi_96;
            const double xi_204 = -xi_156;
            const double xi_207 = xi_66*0.0357142857142857;
            const double xi_209 = xi_114*(Dummy_283 - 0.5);
            const double xi_214 = xi_133*xi_191;
            const double xi_215 = xi_137*xi_189;
            const double u_0 = xi_7*(vel0Term + xi_10 + xi_14) + xi_8*_data_force_20_30_10[_stride_force_0*ctr_0];
            const double xi_27 = u_0*_data_force_20_30_10[_stride_force_0*ctr_0];
            const double xi_33 = xi_27*-0.166666666666667;
            const double xi_37 = u_0*2.0;
            const double xi_38 = xi_37 - 1.0;
            const double xi_42 = xi_37 + 1.0;
            const double xi_52 = u_0*3.0;
            const double xi_53 = -xi_52;
            const double xi_55 = xi_27*-0.0833333333333333;
            const double xi_73 = rho*(u_0*u_0);
            const double xi_127 = rho*u_0;
            const double xi_128 = -vel0Term + xi_116 + xi_127 + xi_4 + _data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            const double xi_129 = xi_128*xi_94;
            const double xi_165 = xi_128*xi_160;
            const double u_1 = xi_7*(vel1Term + xi_17 + xi_20 + xi_9 + _data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0]) + xi_8*_data_force_20_31_10[_stride_force_0*ctr_0];
            const double xi_28 = u_1*_data_force_20_31_10[_stride_force_0*ctr_0];
            const double xi_30 = u_1*2.0;
            const double xi_31 = xi_30 + 1.0;
            const double xi_36 = xi_30 - 1.0;
            const double xi_40 = xi_28*-0.166666666666667;
            const double xi_46 = xi_33 + xi_40;
            const double xi_49 = u_1*3.0;
            const double xi_50 = -xi_49;
            const double xi_59 = xi_28*-0.0833333333333333;
            const double xi_78 = rho*(u_1*u_1);
            const double xi_79 = xi_10 + xi_77 + xi_78;
            const double xi_91 = rho*u_1;
            const double xi_93 = -vel1Term + xi_5 + xi_91 + xi_92 + _data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] + _data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0];
            const double xi_95 = xi_93*xi_94;
            const double xi_161 = xi_160*xi_93;
            const double xi_171 = xi_170*(u_0*xi_91 + xi_9 + xi_92 + _data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0]);
            const double xi_172 = -xi_169 - xi_171;
            const double xi_173 = xi_169 + xi_171;
            const double u_2 = xi_7*(vel2Term + xi_22 + xi_25 + _data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0]) + xi_8*_data_force_20_32_10[_stride_force_0*ctr_0];
            const double xi_29 = u_2*_data_force_20_32_10[_stride_force_0*ctr_0];
            const double xi_34 = xi_29*-0.166666666666667;
            const double xi_35 = xi_33 + xi_34;
            const double xi_41 = xi_34 + xi_40;
            const double xi_43 = u_2*2.0;
            const double xi_44 = xi_43 + 1.0;
            const double xi_47 = xi_43 - 1.0;
            const double xi_48 = xi_29*-0.0833333333333333;
            const double xi_56 = u_2*3.0;
            const double xi_58 = -xi_56;
            const double xi_74 = rho*(u_2*u_2);
            const double xi_82 = omega_bulk*(xi_18 + xi_23 + xi_73 + xi_74 + xi_76 + xi_79 + xi_81 + _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0]);
            const double xi_117 = -xi_74 + _data_pdfs_21_36_10[_stride_pdfs_0*ctr_0] + _data_pdfs_2m1_35_10[_stride_pdfs_0*ctr_0];
            const double xi_118 = omega_shear*(xi_0 + xi_116 + xi_117 + xi_17 + xi_79 - _data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0]);
            const double xi_119 = xi_118*0.125;
            const double xi_121 = omega_shear*(xi_10 + xi_117 + xi_70 + xi_73*2.0 + xi_77 - xi_78 + xi_81 - 2.0*_data_pdfs_20_33_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0] - 2.0*_data_pdfs_20_34_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0] + _data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0] + _data_pdfs_20_32_11[_stride_pdfs_0*ctr_0]);
            const double xi_123 = xi_121*-0.0416666666666667 + xi_122*-0.166666666666667;
            const double xi_124 = xi_123 + xi_63*-0.1 + xi_69*-0.05;
            const double xi_125 = xi_115 + xi_119 + xi_120 + xi_124 + xi_66*0.0285714285714286 + xi_72*0.0142857142857143;
            const double xi_140 = xi_120 + xi_121*0.0833333333333333 + xi_122*0.333333333333333 + xi_66*-0.0714285714285714 + xi_72*-0.0357142857142857;
            const double xi_145 = rho*u_2 - vel2Term + xi_142 + xi_6 + xi_75 + xi_80 + _data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            const double xi_146 = xi_145*xi_94;
            const double xi_154 = -xi_115 - xi_119 + xi_124 + xi_62*0.0952380952380952 + xi_66*-0.0428571428571429 + xi_72*-0.0214285714285714 + xi_86*0.0158730158730159;
            const double xi_157 = xi_118*0.0625;
            const double xi_162 = xi_65*0.0833333333333333 + xi_82*0.0416666666666667;
            const double xi_163 = xi_161 + xi_162;
            const double xi_164 = xi_126 + xi_156 + xi_157 + xi_158 + xi_159 + xi_163;
            const double xi_166 = xi_121*0.0208333333333333 + xi_122*0.0833333333333333;
            const double xi_167 = -xi_165 + xi_166;
            const double xi_168 = xi_141 + xi_167;
            const double xi_174 = xi_165 + xi_166;
            const double xi_175 = xi_139 + xi_174;
            const double xi_176 = -xi_161 + xi_162;
            const double xi_177 = xi_111 + xi_156 + xi_157 + xi_158 + xi_159 + xi_176;
            const double xi_180 = xi_170*(u_2*xi_91 + xi_18 + xi_87 + _data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0]);
            const double xi_184 = xi_123 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183;
            const double xi_193 = xi_145*xi_160;
            const double xi_195 = xi_193 + xi_194;
            const double xi_196 = -xi_186 + xi_188 - xi_190 + xi_192 + xi_195;
            const double xi_201 = xi_163 - xi_197 + xi_198 - xi_199 + xi_200;
            const double xi_202 = xi_176 + xi_197 - xi_198 + xi_199 - xi_200;
            const double xi_203 = xi_123 - xi_178 + xi_179 - xi_180 + xi_181 + xi_182 + xi_183;
            const double xi_205 = -xi_157;
            const double xi_208 = xi_155 + xi_162 + xi_195 + xi_204 + xi_205 + xi_206 + xi_207;
            const double xi_210 = xi_170*(u_2*xi_127 + xi_11 + xi_130 + _data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0]);
            const double xi_211 = -xi_209 - xi_210;
            const double xi_216 = xi_167 - xi_212 + xi_213 - xi_214 + xi_215;
            const double xi_217 = xi_209 + xi_210;
            const double xi_218 = xi_174 + xi_212 - xi_213 + xi_214 - xi_215;
            const double xi_219 = -xi_193 + xi_194;
            const double xi_220 = xi_186 - xi_188 + xi_190 - xi_192 + xi_219;
            const double xi_221 = xi_153 + xi_162 + xi_204 + xi_205 + xi_206 + xi_207 + xi_219;
            const double forceTerm_0 = xi_26*(-xi_27 - xi_28 - xi_29);
            const double forceTerm_1 = xi_26*(xi_31*xi_32 + xi_35);
            const double forceTerm_2 = xi_26*(xi_32*xi_36 + xi_35);
            const double forceTerm_3 = xi_26*(xi_38*xi_39 + xi_41);
            const double forceTerm_4 = xi_26*(xi_39*xi_42 + xi_41);
            const double forceTerm_5 = xi_26*(xi_44*xi_45 + xi_46);
            const double forceTerm_6 = xi_26*(xi_45*xi_47 + xi_46);
            const double forceTerm_7 = xi_26*(xi_48 + xi_51*(xi_38 + xi_50) + xi_54*(xi_31 + xi_53));
            const double forceTerm_8 = xi_26*(xi_48 + xi_51*(xi_42 + xi_49) + xi_54*(xi_31 + xi_52));
            const double forceTerm_9 = xi_26*(xi_48 + xi_51*(xi_38 + xi_49) + xi_54*(xi_36 + xi_52));
            const double forceTerm_10 = xi_26*(xi_48 + xi_51*(xi_42 + xi_50) + xi_54*(xi_36 + xi_53));
            const double forceTerm_11 = xi_26*(xi_54*(xi_31 + xi_56) + xi_55 + xi_57*(xi_44 + xi_49));
            const double forceTerm_12 = xi_26*(xi_54*(xi_36 + xi_58) + xi_55 + xi_57*(xi_44 + xi_50));
            const double forceTerm_13 = xi_26*(xi_51*(xi_38 + xi_58) + xi_57*(xi_44 + xi_53) + xi_59);
            const double forceTerm_14 = xi_26*(xi_51*(xi_42 + xi_56) + xi_57*(xi_44 + xi_52) + xi_59);
            const double forceTerm_15 = xi_26*(xi_54*(xi_31 + xi_58) + xi_55 + xi_57*(xi_47 + xi_50));
            const double forceTerm_16 = xi_26*(xi_54*(xi_36 + xi_56) + xi_55 + xi_57*(xi_47 + xi_49));
            const double forceTerm_17 = xi_26*(xi_51*(xi_38 + xi_56) + xi_57*(xi_47 + xi_52) + xi_59);
            const double forceTerm_18 = xi_26*(xi_51*(xi_42 + xi_58) + xi_57*(xi_47 + xi_53) + xi_59);
            _data_pdfs_tmp_20_30_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_0 + xi_62*0.142857142857143 + xi_63*0.2 - xi_65 + xi_66*0.0857142857142857 + xi_69*0.1 + xi_72*0.0428571428571429 + xi_82*-0.5 + xi_86*0.0238095238095238 + _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_31_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_1 - xi_100 + xi_111 + xi_125 - xi_90 + xi_95 + _data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_32_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_2 + xi_100 + xi_125 + xi_126 + xi_90 - xi_95 + _data_pdfs_20_32_11[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_33_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_3 - xi_129 + xi_132 + xi_134 + xi_139 + xi_140 + _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_34_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_4 + xi_129 - xi_132 - xi_134 + xi_140 + xi_141 + _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_35_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_5 - xi_144 + xi_146 - xi_148 + xi_153 + xi_154 + _data_pdfs_2m1_35_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_36_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_6 + xi_144 - xi_146 + xi_148 + xi_154 + xi_155 + _data_pdfs_21_36_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_37_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_7 + xi_164 + xi_168 + xi_172 + _data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_38_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_8 + xi_164 + xi_173 + xi_175 + _data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_39_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_9 + xi_168 + xi_173 + xi_177 + _data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_310_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_10 + xi_172 + xi_175 + xi_177 + _data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_311_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_11 + xi_184 + xi_196 + xi_201 + _data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_312_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_12 + xi_196 + xi_202 + xi_203 + _data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_313_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_13 + xi_208 + xi_211 + xi_216 + _data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_314_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_14 + xi_208 + xi_217 + xi_218 + _data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_315_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_15 + xi_201 + xi_203 + xi_220 + _data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_316_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_16 + xi_184 + xi_202 + xi_220 + _data_pdfs_21_316_11[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_317_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_17 + xi_216 + xi_217 + xi_221 + _data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_318_10[_stride_pdfs_tmp_0*ctr_0] = forceTerm_18 + xi_211 + xi_218 + xi_221 + _data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
         }
      }
   }
}
}
namespace internal_kernel_collide {
static FUNC_PREFIX void kernel_collide(double * RESTRICT const _data_force, double * RESTRICT _data_pdfs, int64_t const _size_force_0, int64_t const _size_force_1, int64_t const _size_force_2, int64_t const _stride_force_0, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, uint32_t block_offset_0, uint32_t block_offset_1, uint32_t block_offset_2, double omega_bulk, double omega_even, double omega_odd, double omega_shear, uint32_t seed, double temperature, uint32_t time_step)
{
   const double xi_26 = omega_shear*-0.5 + 1.0;
   const double xi_64 = 2.4494897427831779;
   const double xi_89 = omega_odd*0.25;
   const double xi_105 = omega_odd*0.0833333333333333;
   const double xi_170 = omega_shear*0.25;
   const double xi_185 = omega_odd*0.0416666666666667;
   const double xi_187 = omega_odd*0.125;
   const int64_t rr_0 = 0.0;
   const double xi_94 = rr_0*0.166666666666667;
   const double xi_160 = rr_0*0.0833333333333333;
   for (int ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1)
   {
      double * RESTRICT _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_315 = _data_pdfs + _stride_pdfs_2*ctr_2 + 15*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_316 = _data_pdfs + _stride_pdfs_2*ctr_2 + 16*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2*ctr_2;
      double * RESTRICT _data_pdfs_20_312 = _data_pdfs + _stride_pdfs_2*ctr_2 + 12*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_35 = _data_pdfs + _stride_pdfs_2*ctr_2 + 5*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_311 = _data_pdfs + _stride_pdfs_2*ctr_2 + 11*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      double * RESTRICT _data_force_20_32 = _data_force + _stride_force_2*ctr_2 + 2*_stride_force_3;
      double * RESTRICT _data_pdfs_20_314 = _data_pdfs + _stride_pdfs_2*ctr_2 + 14*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_36 = _data_pdfs + _stride_pdfs_2*ctr_2 + 6*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_318 = _data_pdfs + _stride_pdfs_2*ctr_2 + 18*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_317 = _data_pdfs + _stride_pdfs_2*ctr_2 + 17*_stride_pdfs_3;
      double * RESTRICT _data_force_20_30 = _data_force + _stride_force_2*ctr_2;
      double * RESTRICT _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_20_313 = _data_pdfs + _stride_pdfs_2*ctr_2 + 13*_stride_pdfs_3;
      double * RESTRICT _data_force_20_31 = _data_force + _stride_force_2*ctr_2 + _stride_force_3;
      double * RESTRICT _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      for (int ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_20_34_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_34;
         double * RESTRICT _data_pdfs_20_37_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_37;
         double * RESTRICT _data_pdfs_20_32_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_32;
         double * RESTRICT _data_pdfs_20_315_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_315;
         double * RESTRICT _data_pdfs_20_316_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_316;
         double * RESTRICT _data_pdfs_20_30_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_30;
         double * RESTRICT _data_pdfs_20_312_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_312;
         double * RESTRICT _data_pdfs_20_35_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_35;
         double * RESTRICT _data_pdfs_20_310_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_310;
         double * RESTRICT _data_pdfs_20_33_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_33;
         double * RESTRICT _data_pdfs_20_311_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_311;
         double * RESTRICT _data_pdfs_20_38_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_38;
         double * RESTRICT _data_force_20_32_10 = _stride_force_1*ctr_1 + _data_force_20_32;
         double * RESTRICT _data_pdfs_20_314_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_314;
         double * RESTRICT _data_pdfs_20_36_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_36;
         double * RESTRICT _data_pdfs_20_318_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_318;
         double * RESTRICT _data_pdfs_20_317_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_317;
         double * RESTRICT _data_force_20_30_10 = _stride_force_1*ctr_1 + _data_force_20_30;
         double * RESTRICT _data_pdfs_20_39_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_39;
         double * RESTRICT _data_pdfs_20_313_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_313;
         double * RESTRICT _data_force_20_31_10 = _stride_force_1*ctr_1 + _data_force_20_31;
         double * RESTRICT _data_pdfs_20_31_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_31;
         for (int ctr_0 = 0; ctr_0 < _size_force_0; ctr_0 += 1)
         {
            const double xi_222 = _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0];
            const double xi_223 = _data_pdfs_20_37_10[_stride_pdfs_0*ctr_0];
            const double xi_224 = _data_pdfs_20_32_10[_stride_pdfs_0*ctr_0];
            const double xi_225 = _data_pdfs_20_315_10[_stride_pdfs_0*ctr_0];
            const double xi_226 = _data_pdfs_20_316_10[_stride_pdfs_0*ctr_0];
            const double xi_227 = _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0];
            const double xi_228 = _data_pdfs_20_312_10[_stride_pdfs_0*ctr_0];
            const double xi_229 = _data_pdfs_20_35_10[_stride_pdfs_0*ctr_0];
            const double xi_230 = _data_pdfs_20_310_10[_stride_pdfs_0*ctr_0];
            const double xi_231 = _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0];
            const double xi_232 = _data_pdfs_20_311_10[_stride_pdfs_0*ctr_0];
            const double xi_233 = _data_pdfs_20_38_10[_stride_pdfs_0*ctr_0];
            const double xi_234 = _data_force_20_32_10[_stride_force_0*ctr_0];
            const double xi_235 = _data_pdfs_20_314_10[_stride_pdfs_0*ctr_0];
            const double xi_236 = _data_pdfs_20_36_10[_stride_pdfs_0*ctr_0];
            const double xi_237 = _data_pdfs_20_318_10[_stride_pdfs_0*ctr_0];
            const double xi_238 = _data_pdfs_20_317_10[_stride_pdfs_0*ctr_0];
            const double xi_239 = _data_force_20_30_10[_stride_force_0*ctr_0];
            const double xi_240 = _data_pdfs_20_39_10[_stride_pdfs_0*ctr_0];
            const double xi_241 = _data_pdfs_20_313_10[_stride_pdfs_0*ctr_0];
            const double xi_242 = _data_force_20_31_10[_stride_force_0*ctr_0];
            const double xi_243 = _data_pdfs_20_31_10[_stride_pdfs_0*ctr_0];
            
            float Dummy_292;
            float Dummy_293;
            float Dummy_294;
            float Dummy_295;
            philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 3, seed, Dummy_292, Dummy_293, Dummy_294, Dummy_295);
            
            
            float Dummy_288;
            float Dummy_289;
            float Dummy_290;
            float Dummy_291;
            philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 2, seed, Dummy_288, Dummy_289, Dummy_290, Dummy_291);
            
            
            float Dummy_284;
            float Dummy_285;
            float Dummy_286;
            float Dummy_287;
            philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 1, seed, Dummy_284, Dummy_285, Dummy_286, Dummy_287);
            
            
            float Dummy_280;
            float Dummy_281;
            float Dummy_282;
            float Dummy_283;
            philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 0, seed, Dummy_280, Dummy_281, Dummy_282, Dummy_283);
            
            const double xi_0 = xi_235 + xi_237;
            const double xi_1 = xi_0 + xi_222;
            const double xi_2 = xi_225 + xi_232 + xi_243;
            const double xi_3 = xi_228 + xi_229;
            const double xi_4 = xi_231 + xi_240;
            const double xi_5 = xi_224 + xi_226;
            const double xi_6 = xi_236 + xi_238;
            const double xi_9 = -xi_240;
            const double xi_10 = -xi_223 + xi_9;
            const double xi_11 = -xi_238;
            const double xi_12 = -xi_241;
            const double xi_13 = -xi_231;
            const double xi_14 = xi_11 + xi_12 + xi_13;
            const double xi_15 = -xi_224;
            const double xi_16 = -xi_230;
            const double xi_17 = xi_15 + xi_16;
            const double xi_18 = -xi_226;
            const double xi_19 = -xi_228;
            const double xi_20 = xi_18 + xi_19;
            const double xi_21 = -xi_237;
            const double xi_22 = xi_11 + xi_21;
            const double xi_23 = -xi_225;
            const double xi_24 = -xi_236;
            const double xi_25 = xi_18 + xi_23 + xi_232 + xi_24;
            const double xi_32 = xi_242*0.166666666666667;
            const double xi_39 = xi_239*0.166666666666667;
            const double xi_45 = xi_234*0.166666666666667;
            const double xi_51 = xi_239*0.0833333333333333;
            const double xi_54 = xi_242*0.0833333333333333;
            const double xi_57 = xi_234*0.0833333333333333;
            const double xi_67 = -xi_227;
            const double xi_68 = xi_229*3.0 + xi_236*3.0 + xi_67;
            const double xi_69 = omega_even*(xi_224*3.0 + xi_225*-3.0 + xi_226*-3.0 + xi_228*-3.0 + xi_232*-3.0 + xi_243*3.0 + xi_68);
            const double xi_70 = xi_225*2.0 + xi_226*2.0 + xi_228*2.0 + xi_232*2.0;
            const double xi_71 = xi_222*5.0 + xi_231*5.0 + xi_70;
            const double xi_72 = omega_even*(xi_224*-2.0 + xi_235*-5.0 + xi_237*-5.0 + xi_238*-5.0 + xi_241*-5.0 + xi_243*-2.0 + xi_68 + xi_71);
            const double xi_75 = -xi_232;
            const double xi_76 = xi_19 + xi_75;
            const double xi_77 = -xi_233;
            const double xi_80 = -xi_235;
            const double xi_81 = xi_12 + xi_16 + xi_22 + xi_80;
            const double xi_83 = xi_241*2.0;
            const double xi_84 = xi_235*2.0;
            const double xi_85 = xi_237*2.0 + xi_238*2.0;
            const double xi_86 = omega_even*(xi_223*-7.0 + xi_224*5.0 + xi_229*-4.0 + xi_230*-7.0 + xi_233*-7.0 + xi_236*-4.0 + xi_240*-7.0 + xi_243*5.0 + xi_67 + xi_71 + xi_83 + xi_84 + xi_85);
            const double xi_87 = xi_228 + xi_75;
            const double xi_88 = xi_15 + xi_226 + xi_23 + xi_243 + xi_87;
            const double xi_90 = xi_88*xi_89;
            const double xi_92 = xi_230 + xi_77;
            const double xi_96 = Dummy_291 - 0.5;
            const double xi_101 = xi_223*2.0;
            const double xi_102 = xi_230*2.0;
            const double xi_103 = xi_233*-2.0 + xi_240*2.0;
            const double xi_104 = -xi_101 + xi_102 + xi_103 + xi_15 + xi_2 + xi_20;
            const double xi_106 = xi_104*xi_105;
            const double xi_107 = Dummy_286 - 0.5;
            const double xi_112 = Dummy_281 - 0.5;
            const double xi_116 = xi_238 + xi_241;
            const double xi_130 = xi_241 + xi_80;
            const double xi_131 = xi_13 + xi_130 + xi_21 + xi_222 + xi_238;
            const double xi_132 = xi_131*xi_89;
            const double xi_133 = Dummy_289 - 0.5;
            const double xi_135 = xi_1 + xi_101 - xi_102 + xi_103 + xi_14;
            const double xi_136 = xi_105*xi_135;
            const double xi_137 = Dummy_288 - 0.5;
            const double xi_142 = xi_225 + xi_226;
            const double xi_143 = xi_142 + xi_229 + xi_24 + xi_76;
            const double xi_144 = xi_143*xi_89;
            const double xi_147 = Dummy_290 - 0.5;
            const double xi_149 = xi_25 + xi_3 - xi_83 - xi_84 + xi_85;
            const double xi_150 = xi_105*xi_149;
            const double xi_151 = Dummy_287 - 0.5;
            const double xi_158 = xi_86*0.0138888888888889;
            const double xi_179 = xi_72*-0.00714285714285714;
            const double xi_181 = xi_69*0.025;
            const double xi_186 = xi_149*xi_185;
            const double xi_188 = xi_143*xi_187;
            const double xi_197 = xi_104*xi_185;
            const double xi_198 = xi_187*xi_88;
            const double xi_206 = xi_72*0.0178571428571429;
            const double xi_212 = xi_131*xi_187;
            const double xi_213 = xi_135*xi_185;
            const double vel0Term = xi_1 + xi_230 + xi_233;
            const double vel1Term = xi_2 + xi_223;
            const double vel2Term = xi_241 + xi_3;
            const double rho = vel0Term + vel1Term + vel2Term + xi_227 + xi_4 + xi_5 + xi_6;
            const double xi_7 = 1 / (rho);
            const double xi_8 = xi_7*0.5;
            const double xi_60 = rho*temperature;
            const double xi_61 = sqrt(xi_60*(-((-omega_even + 1.0)*(-omega_even + 1.0)) + 1.0));
            const double xi_62 = xi_61*(Dummy_292 - 0.5)*3.7416573867739413;
            const double xi_63 = xi_61*(Dummy_294 - 0.5)*5.4772255750516612;
            const double xi_65 = xi_64*sqrt(xi_60*(-((-omega_bulk + 1.0)*(-omega_bulk + 1.0)) + 1.0))*(Dummy_285 - 0.5);
            const double xi_66 = xi_61*(Dummy_293 - 0.5)*8.3666002653407556;
            const double xi_97 = sqrt(xi_60*(-((-omega_odd + 1.0)*(-omega_odd + 1.0)) + 1.0));
            const double xi_98 = xi_97*1.4142135623730951;
            const double xi_99 = xi_98*0.5;
            const double xi_100 = xi_96*xi_99;
            const double xi_108 = xi_64*xi_97;
            const double xi_109 = xi_108*0.166666666666667;
            const double xi_110 = xi_107*xi_109;
            const double xi_111 = -xi_106 - xi_110;
            const double xi_113 = sqrt(xi_60*(-((-omega_shear + 1.0)*(-omega_shear + 1.0)) + 1.0));
            const double xi_114 = xi_113*0.5;
            const double xi_115 = xi_112*xi_114;
            const double xi_120 = xi_62*-0.119047619047619 + xi_86*-0.0198412698412698;
            const double xi_122 = xi_113*(Dummy_280 - 0.5)*1.7320508075688772;
            const double xi_126 = xi_106 + xi_110;
            const double xi_134 = xi_133*xi_99;
            const double xi_138 = xi_109*xi_137;
            const double xi_139 = xi_136 + xi_138;
            const double xi_141 = -xi_136 - xi_138;
            const double xi_148 = xi_147*xi_99;
            const double xi_152 = xi_109*xi_151;
            const double xi_153 = -xi_150 - xi_152;
            const double xi_155 = xi_150 + xi_152;
            const double xi_156 = xi_112*xi_113*0.25;
            const double xi_159 = xi_62*0.0833333333333333;
            const double xi_169 = xi_114*(Dummy_282 - 0.5);
            const double xi_178 = xi_114*(Dummy_284 - 0.5);
            const double xi_182 = xi_66*-0.0142857142857143;
            const double xi_183 = xi_63*0.05;
            const double xi_189 = xi_108*0.0833333333333333;
            const double xi_190 = xi_151*xi_189;
            const double xi_191 = xi_98*0.25;
            const double xi_192 = xi_147*xi_191;
            const double xi_194 = xi_62*-0.0238095238095238 + xi_86*-0.00396825396825397;
            const double xi_199 = xi_107*xi_189;
            const double xi_200 = xi_191*xi_96;
            const double xi_204 = -xi_156;
            const double xi_207 = xi_66*0.0357142857142857;
            const double xi_209 = xi_114*(Dummy_283 - 0.5);
            const double xi_214 = xi_133*xi_191;
            const double xi_215 = xi_137*xi_189;
            const double u_0 = xi_239*xi_8 + xi_7*(vel0Term + xi_10 + xi_14);
            const double xi_27 = u_0*xi_239;
            const double xi_33 = xi_27*-0.166666666666667;
            const double xi_37 = u_0*2.0;
            const double xi_38 = xi_37 - 1.0;
            const double xi_42 = xi_37 + 1.0;
            const double xi_52 = u_0*3.0;
            const double xi_53 = -xi_52;
            const double xi_55 = xi_27*-0.0833333333333333;
            const double xi_73 = rho*(u_0*u_0);
            const double xi_127 = rho*u_0;
            const double xi_128 = -vel0Term + xi_116 + xi_127 + xi_223 + xi_4;
            const double xi_129 = xi_128*xi_94;
            const double xi_165 = xi_128*xi_160;
            const double u_1 = xi_242*xi_8 + xi_7*(vel1Term + xi_17 + xi_20 + xi_233 + xi_9);
            const double xi_28 = u_1*xi_242;
            const double xi_30 = u_1*2.0;
            const double xi_31 = xi_30 + 1.0;
            const double xi_36 = xi_30 - 1.0;
            const double xi_40 = xi_28*-0.166666666666667;
            const double xi_46 = xi_33 + xi_40;
            const double xi_49 = u_1*3.0;
            const double xi_50 = -xi_49;
            const double xi_59 = xi_28*-0.0833333333333333;
            const double xi_78 = rho*(u_1*u_1);
            const double xi_79 = xi_10 + xi_77 + xi_78;
            const double xi_91 = rho*u_1;
            const double xi_93 = -vel1Term + xi_228 + xi_240 + xi_5 + xi_91 + xi_92;
            const double xi_95 = xi_93*xi_94;
            const double xi_161 = xi_160*xi_93;
            const double xi_171 = xi_170*(u_0*xi_91 + xi_223 + xi_9 + xi_92);
            const double xi_172 = -xi_169 - xi_171;
            const double xi_173 = xi_169 + xi_171;
            const double u_2 = xi_234*xi_8 + xi_7*(vel2Term + xi_22 + xi_235 + xi_25);
            const double xi_29 = u_2*xi_234;
            const double xi_34 = xi_29*-0.166666666666667;
            const double xi_35 = xi_33 + xi_34;
            const double xi_41 = xi_34 + xi_40;
            const double xi_43 = u_2*2.0;
            const double xi_44 = xi_43 + 1.0;
            const double xi_47 = xi_43 - 1.0;
            const double xi_48 = xi_29*-0.0833333333333333;
            const double xi_56 = u_2*3.0;
            const double xi_58 = -xi_56;
            const double xi_74 = rho*(u_2*u_2);
            const double xi_82 = omega_bulk*(xi_18 + xi_227 + xi_23 + xi_73 + xi_74 + xi_76 + xi_79 + xi_81);
            const double xi_117 = xi_229 + xi_236 - xi_74;
            const double xi_118 = omega_shear*(xi_0 + xi_116 + xi_117 + xi_17 - xi_243 + xi_79);
            const double xi_119 = xi_118*0.125;
            const double xi_121 = omega_shear*(xi_10 + xi_117 + xi_222*-2.0 + xi_224 + xi_231*-2.0 + xi_243 + xi_70 + xi_73*2.0 + xi_77 - xi_78 + xi_81);
            const double xi_123 = xi_121*-0.0416666666666667 + xi_122*-0.166666666666667;
            const double xi_124 = xi_123 + xi_63*-0.1 + xi_69*-0.05;
            const double xi_125 = xi_115 + xi_119 + xi_120 + xi_124 + xi_66*0.0285714285714286 + xi_72*0.0142857142857143;
            const double xi_140 = xi_120 + xi_121*0.0833333333333333 + xi_122*0.333333333333333 + xi_66*-0.0714285714285714 + xi_72*-0.0357142857142857;
            const double xi_145 = rho*u_2 - vel2Term + xi_142 + xi_237 + xi_6 + xi_75 + xi_80;
            const double xi_146 = xi_145*xi_94;
            const double xi_154 = -xi_115 - xi_119 + xi_124 + xi_62*0.0952380952380952 + xi_66*-0.0428571428571429 + xi_72*-0.0214285714285714 + xi_86*0.0158730158730159;
            const double xi_157 = xi_118*0.0625;
            const double xi_162 = xi_65*0.0833333333333333 + xi_82*0.0416666666666667;
            const double xi_163 = xi_161 + xi_162;
            const double xi_164 = xi_126 + xi_156 + xi_157 + xi_158 + xi_159 + xi_163;
            const double xi_166 = xi_121*0.0208333333333333 + xi_122*0.0833333333333333;
            const double xi_167 = -xi_165 + xi_166;
            const double xi_168 = xi_141 + xi_167;
            const double xi_174 = xi_165 + xi_166;
            const double xi_175 = xi_139 + xi_174;
            const double xi_176 = -xi_161 + xi_162;
            const double xi_177 = xi_111 + xi_156 + xi_157 + xi_158 + xi_159 + xi_176;
            const double xi_180 = xi_170*(u_2*xi_91 + xi_18 + xi_225 + xi_87);
            const double xi_184 = xi_123 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183;
            const double xi_193 = xi_145*xi_160;
            const double xi_195 = xi_193 + xi_194;
            const double xi_196 = -xi_186 + xi_188 - xi_190 + xi_192 + xi_195;
            const double xi_201 = xi_163 - xi_197 + xi_198 - xi_199 + xi_200;
            const double xi_202 = xi_176 + xi_197 - xi_198 + xi_199 - xi_200;
            const double xi_203 = xi_123 - xi_178 + xi_179 - xi_180 + xi_181 + xi_182 + xi_183;
            const double xi_205 = -xi_157;
            const double xi_208 = xi_155 + xi_162 + xi_195 + xi_204 + xi_205 + xi_206 + xi_207;
            const double xi_210 = xi_170*(u_2*xi_127 + xi_11 + xi_130 + xi_237);
            const double xi_211 = -xi_209 - xi_210;
            const double xi_216 = xi_167 - xi_212 + xi_213 - xi_214 + xi_215;
            const double xi_217 = xi_209 + xi_210;
            const double xi_218 = xi_174 + xi_212 - xi_213 + xi_214 - xi_215;
            const double xi_219 = -xi_193 + xi_194;
            const double xi_220 = xi_186 - xi_188 + xi_190 - xi_192 + xi_219;
            const double xi_221 = xi_153 + xi_162 + xi_204 + xi_205 + xi_206 + xi_207 + xi_219;
            const double forceTerm_0 = xi_26*(-xi_27 - xi_28 - xi_29);
            const double forceTerm_1 = xi_26*(xi_31*xi_32 + xi_35);
            const double forceTerm_2 = xi_26*(xi_32*xi_36 + xi_35);
            const double forceTerm_3 = xi_26*(xi_38*xi_39 + xi_41);
            const double forceTerm_4 = xi_26*(xi_39*xi_42 + xi_41);
            const double forceTerm_5 = xi_26*(xi_44*xi_45 + xi_46);
            const double forceTerm_6 = xi_26*(xi_45*xi_47 + xi_46);
            const double forceTerm_7 = xi_26*(xi_48 + xi_51*(xi_38 + xi_50) + xi_54*(xi_31 + xi_53));
            const double forceTerm_8 = xi_26*(xi_48 + xi_51*(xi_42 + xi_49) + xi_54*(xi_31 + xi_52));
            const double forceTerm_9 = xi_26*(xi_48 + xi_51*(xi_38 + xi_49) + xi_54*(xi_36 + xi_52));
            const double forceTerm_10 = xi_26*(xi_48 + xi_51*(xi_42 + xi_50) + xi_54*(xi_36 + xi_53));
            const double forceTerm_11 = xi_26*(xi_54*(xi_31 + xi_56) + xi_55 + xi_57*(xi_44 + xi_49));
            const double forceTerm_12 = xi_26*(xi_54*(xi_36 + xi_58) + xi_55 + xi_57*(xi_44 + xi_50));
            const double forceTerm_13 = xi_26*(xi_51*(xi_38 + xi_58) + xi_57*(xi_44 + xi_53) + xi_59);
            const double forceTerm_14 = xi_26*(xi_51*(xi_42 + xi_56) + xi_57*(xi_44 + xi_52) + xi_59);
            const double forceTerm_15 = xi_26*(xi_54*(xi_31 + xi_58) + xi_55 + xi_57*(xi_47 + xi_50));
            const double forceTerm_16 = xi_26*(xi_54*(xi_36 + xi_56) + xi_55 + xi_57*(xi_47 + xi_49));
            const double forceTerm_17 = xi_26*(xi_51*(xi_38 + xi_56) + xi_57*(xi_47 + xi_52) + xi_59);
            const double forceTerm_18 = xi_26*(xi_51*(xi_42 + xi_58) + xi_57*(xi_47 + xi_53) + xi_59);
            _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0] = forceTerm_0 + xi_227 + xi_62*0.142857142857143 + xi_63*0.2 - xi_65 + xi_66*0.0857142857142857 + xi_69*0.1 + xi_72*0.0428571428571429 + xi_82*-0.5 + xi_86*0.0238095238095238;
            _data_pdfs_20_31_10[_stride_pdfs_0*ctr_0] = forceTerm_1 - xi_100 + xi_111 + xi_125 + xi_243 - xi_90 + xi_95;
            _data_pdfs_20_32_10[_stride_pdfs_0*ctr_0] = forceTerm_2 + xi_100 + xi_125 + xi_126 + xi_224 + xi_90 - xi_95;
            _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0] = forceTerm_3 - xi_129 + xi_132 + xi_134 + xi_139 + xi_140 + xi_231;
            _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0] = forceTerm_4 + xi_129 - xi_132 - xi_134 + xi_140 + xi_141 + xi_222;
            _data_pdfs_20_35_10[_stride_pdfs_0*ctr_0] = forceTerm_5 - xi_144 + xi_146 - xi_148 + xi_153 + xi_154 + xi_229;
            _data_pdfs_20_36_10[_stride_pdfs_0*ctr_0] = forceTerm_6 + xi_144 - xi_146 + xi_148 + xi_154 + xi_155 + xi_236;
            _data_pdfs_20_37_10[_stride_pdfs_0*ctr_0] = forceTerm_7 + xi_164 + xi_168 + xi_172 + xi_223;
            _data_pdfs_20_38_10[_stride_pdfs_0*ctr_0] = forceTerm_8 + xi_164 + xi_173 + xi_175 + xi_233;
            _data_pdfs_20_39_10[_stride_pdfs_0*ctr_0] = forceTerm_9 + xi_168 + xi_173 + xi_177 + xi_240;
            _data_pdfs_20_310_10[_stride_pdfs_0*ctr_0] = forceTerm_10 + xi_172 + xi_175 + xi_177 + xi_230;
            _data_pdfs_20_311_10[_stride_pdfs_0*ctr_0] = forceTerm_11 + xi_184 + xi_196 + xi_201 + xi_232;
            _data_pdfs_20_312_10[_stride_pdfs_0*ctr_0] = forceTerm_12 + xi_196 + xi_202 + xi_203 + xi_228;
            _data_pdfs_20_313_10[_stride_pdfs_0*ctr_0] = forceTerm_13 + xi_208 + xi_211 + xi_216 + xi_241;
            _data_pdfs_20_314_10[_stride_pdfs_0*ctr_0] = forceTerm_14 + xi_208 + xi_217 + xi_218 + xi_235;
            _data_pdfs_20_315_10[_stride_pdfs_0*ctr_0] = forceTerm_15 + xi_201 + xi_203 + xi_220 + xi_225;
            _data_pdfs_20_316_10[_stride_pdfs_0*ctr_0] = forceTerm_16 + xi_184 + xi_202 + xi_220 + xi_226;
            _data_pdfs_20_317_10[_stride_pdfs_0*ctr_0] = forceTerm_17 + xi_216 + xi_217 + xi_221 + xi_238;
            _data_pdfs_20_318_10[_stride_pdfs_0*ctr_0] = forceTerm_18 + xi_211 + xi_218 + xi_221 + xi_237;
         }
      }
   }
}
}
namespace internal_kernel_stream {
static FUNC_PREFIX void kernel_stream(double * RESTRICT const _data_pdfs, double * RESTRICT _data_pdfs_tmp, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_pdfs_tmp_0, int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2, int64_t const _stride_pdfs_tmp_3)
{
   for (int ctr_2 = 1; ctr_2 < _size_pdfs_2 - 1; ctr_2 += 1)
   {
      double * RESTRICT _data_pdfs_tmp_20_30 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2;
      double * RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2*ctr_2;
      double * RESTRICT _data_pdfs_tmp_20_31 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + _stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_32 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 2*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2*ctr_2 + 2*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_33 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 3*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2*ctr_2 + 3*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_34 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 4*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2*ctr_2 + 4*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_35 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 5*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_2m1_35 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 5*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_36 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 6*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_21_36 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 6*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_37 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 7*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2*ctr_2 + 7*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_38 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 8*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2*ctr_2 + 8*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_39 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 9*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2*ctr_2 + 9*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_310 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 10*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2*ctr_2 + 10*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_311 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 11*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_2m1_311 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 11*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_312 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 12*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_2m1_312 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 12*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_313 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 13*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_2m1_313 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 13*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_314 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 14*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_2m1_314 = _data_pdfs + _stride_pdfs_2*ctr_2 - _stride_pdfs_2 + 14*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_315 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 15*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_21_315 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 15*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_316 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 16*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_21_316 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 16*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_317 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 17*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_21_317 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 17*_stride_pdfs_3;
      double * RESTRICT _data_pdfs_tmp_20_318 = _data_pdfs_tmp + _stride_pdfs_tmp_2*ctr_2 + 18*_stride_pdfs_tmp_3;
      double * RESTRICT _data_pdfs_21_318 = _data_pdfs + _stride_pdfs_2*ctr_2 + _stride_pdfs_2 + 18*_stride_pdfs_3;
      for (int ctr_1 = 1; ctr_1 < _size_pdfs_1 - 1; ctr_1 += 1)
      {
         double * RESTRICT _data_pdfs_tmp_20_30_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_30;
         double * RESTRICT _data_pdfs_20_30_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_30;
         double * RESTRICT _data_pdfs_tmp_20_31_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_31;
         double * RESTRICT _data_pdfs_20_31_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_31;
         double * RESTRICT _data_pdfs_tmp_20_32_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_32;
         double * RESTRICT _data_pdfs_20_32_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_32;
         double * RESTRICT _data_pdfs_tmp_20_33_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_33;
         double * RESTRICT _data_pdfs_20_33_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_33;
         double * RESTRICT _data_pdfs_tmp_20_34_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_34;
         double * RESTRICT _data_pdfs_20_34_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_20_34;
         double * RESTRICT _data_pdfs_tmp_20_35_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_35;
         double * RESTRICT _data_pdfs_2m1_35_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_35;
         double * RESTRICT _data_pdfs_tmp_20_36_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_36;
         double * RESTRICT _data_pdfs_21_36_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_36;
         double * RESTRICT _data_pdfs_tmp_20_37_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_37;
         double * RESTRICT _data_pdfs_20_37_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_37;
         double * RESTRICT _data_pdfs_tmp_20_38_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_38;
         double * RESTRICT _data_pdfs_20_38_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_20_38;
         double * RESTRICT _data_pdfs_tmp_20_39_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_39;
         double * RESTRICT _data_pdfs_20_39_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_39;
         double * RESTRICT _data_pdfs_tmp_20_310_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_310;
         double * RESTRICT _data_pdfs_20_310_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_20_310;
         double * RESTRICT _data_pdfs_tmp_20_311_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_311;
         double * RESTRICT _data_pdfs_2m1_311_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_311;
         double * RESTRICT _data_pdfs_tmp_20_312_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_312;
         double * RESTRICT _data_pdfs_2m1_312_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_312;
         double * RESTRICT _data_pdfs_tmp_20_313_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_313;
         double * RESTRICT _data_pdfs_2m1_313_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_313;
         double * RESTRICT _data_pdfs_tmp_20_314_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_314;
         double * RESTRICT _data_pdfs_2m1_314_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_2m1_314;
         double * RESTRICT _data_pdfs_tmp_20_315_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_315;
         double * RESTRICT _data_pdfs_21_315_1m1 = _stride_pdfs_1*ctr_1 - _stride_pdfs_1 + _data_pdfs_21_315;
         double * RESTRICT _data_pdfs_tmp_20_316_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_316;
         double * RESTRICT _data_pdfs_21_316_11 = _stride_pdfs_1*ctr_1 + _stride_pdfs_1 + _data_pdfs_21_316;
         double * RESTRICT _data_pdfs_tmp_20_317_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_317;
         double * RESTRICT _data_pdfs_21_317_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_317;
         double * RESTRICT _data_pdfs_tmp_20_318_10 = _stride_pdfs_tmp_1*ctr_1 + _data_pdfs_tmp_20_318;
         double * RESTRICT _data_pdfs_21_318_10 = _stride_pdfs_1*ctr_1 + _data_pdfs_21_318;
         for (int ctr_0 = 1; ctr_0 < _size_pdfs_0 - 1; ctr_0 += 1)
         {
            _data_pdfs_tmp_20_30_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_30_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_31_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_31_1m1[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_32_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_32_11[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_33_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_33_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_34_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_34_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_35_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_2m1_35_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_36_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_21_36_10[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_37_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_37_1m1[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_38_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_38_1m1[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_39_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_39_11[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_310_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_20_310_11[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_311_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_2m1_311_1m1[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_312_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_2m1_312_11[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_313_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_2m1_313_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_314_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_2m1_314_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
            _data_pdfs_tmp_20_315_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_21_315_1m1[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_316_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_21_316_11[_stride_pdfs_0*ctr_0];
            _data_pdfs_tmp_20_317_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_21_317_10[_stride_pdfs_0*ctr_0 + _stride_pdfs_0];
            _data_pdfs_tmp_20_318_10[_stride_pdfs_tmp_0*ctr_0] = _data_pdfs_21_318_10[_stride_pdfs_0*ctr_0 - _stride_pdfs_0];
         }
      }
   }
}
}


const real_t FluctuatingMRT_LatticeModel::w[19] = { 0.333333333333333,0.0555555555555556,0.0555555555555556,0.0555555555555556,0.0555555555555556,0.0555555555555556,0.0555555555555556,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778,0.0277777777777778 };
const real_t FluctuatingMRT_LatticeModel::wInv[19] = { 3.00000000000000,18.0000000000000,18.0000000000000,18.0000000000000,18.0000000000000,18.0000000000000,18.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000,36.0000000000000 };

void FluctuatingMRT_LatticeModel::Sweep::streamCollide( IBlock * block, const uint_t numberOfGhostLayersToInclude )
{
    auto pdfs = block->getData< field::GhostLayerField<double, 19> >(pdfsID);
    field::GhostLayerField<double, 19> * pdfs_tmp;
    {
        // Getting temporary field pdfs_tmp
        auto it = cache_pdfs_.find( pdfs );
        if( it != cache_pdfs_.end() )
        {
            pdfs_tmp = *it;
        }
        else
        {
            pdfs_tmp = pdfs->cloneUninitialized();
            cache_pdfs_.insert(pdfs_tmp);
        }
    }


    auto & lm = dynamic_cast< lbm::PdfField<FluctuatingMRT_LatticeModel> * > (pdfs)->latticeModel();
    WALBERLA_ASSERT_EQUAL( *(lm.blockId_), block->getId() );

    auto & temperature = lm.temperature_;
    auto & block_offset_0 = lm.block_offset_0_;
    auto & force = lm.force_;
    auto & block_offset_2 = lm.block_offset_2_;
    auto & omega_odd = lm.omega_odd_;
    auto & block_offset_1 = lm.block_offset_1_;
    auto & seed = lm.seed_;
    auto & omega_shear = lm.omega_shear_;
    auto & time_step = lm.time_step_;
    auto & omega_even = lm.omega_even_;
    auto & omega_bulk = lm.omega_bulk_;
    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -int_c(force->nrOfGhostLayers()));
    double * RESTRICT const _data_force = force->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -int_c(pdfs->nrOfGhostLayers()));
    double * RESTRICT const _data_pdfs = pdfs->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -int_c(pdfs_tmp->nrOfGhostLayers()));
    double * RESTRICT _data_pdfs_tmp = pdfs_tmp->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(cell_idx_c(force->xSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_force_0 = int64_t(cell_idx_c(force->xSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(cell_idx_c(force->ySize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_force_1 = int64_t(cell_idx_c(force->ySize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(cell_idx_c(force->zSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_force_2 = int64_t(cell_idx_c(force->zSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    const int64_t _stride_force_0 = int64_t(force->xStride());
    const int64_t _stride_force_1 = int64_t(force->yStride());
    const int64_t _stride_force_2 = int64_t(force->zStride());
    const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    const int64_t _stride_pdfs_tmp_0 = int64_t(pdfs_tmp->xStride());
    const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
    const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
    const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
    internal_kernel_streamCollide::kernel_streamCollide(_data_force, _data_pdfs, _data_pdfs_tmp, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3, block_offset_0, block_offset_1, block_offset_2, omega_bulk, omega_even, omega_odd, omega_shear, seed, temperature, time_step);
    pdfs->swapDataPointers(pdfs_tmp);

}

void FluctuatingMRT_LatticeModel::Sweep::collide( IBlock * block, const uint_t numberOfGhostLayersToInclude )
{
   auto pdfs = block->getData< field::GhostLayerField<double, 19> >(pdfsID);


    auto & lm = dynamic_cast< lbm::PdfField<FluctuatingMRT_LatticeModel> * > (pdfs)->latticeModel();
    WALBERLA_ASSERT_EQUAL( *(lm.blockId_), block->getId() );

    auto & temperature = lm.temperature_;
    auto & block_offset_0 = lm.block_offset_0_;
    auto & force = lm.force_;
    auto & block_offset_2 = lm.block_offset_2_;
    auto & omega_odd = lm.omega_odd_;
    auto & block_offset_1 = lm.block_offset_1_;
    auto & seed = lm.seed_;
    auto & omega_shear = lm.omega_shear_;
    auto & time_step = lm.time_step_;
    auto & omega_even = lm.omega_even_;
    auto & omega_bulk = lm.omega_bulk_;
    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude), -int_c(force->nrOfGhostLayers()));
    double * RESTRICT const _data_force = force->dataAt(-cell_idx_c(numberOfGhostLayersToInclude), -cell_idx_c(numberOfGhostLayersToInclude), -cell_idx_c(numberOfGhostLayersToInclude), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude), -int_c(pdfs->nrOfGhostLayers()));
    double * RESTRICT _data_pdfs = pdfs->dataAt(-cell_idx_c(numberOfGhostLayersToInclude), -cell_idx_c(numberOfGhostLayersToInclude), -cell_idx_c(numberOfGhostLayersToInclude), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(cell_idx_c(force->xSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude)));
    const int64_t _size_force_0 = int64_t(cell_idx_c(force->xSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude));
    WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(cell_idx_c(force->ySize()) + 2*cell_idx_c(numberOfGhostLayersToInclude)));
    const int64_t _size_force_1 = int64_t(cell_idx_c(force->ySize()) + 2*cell_idx_c(numberOfGhostLayersToInclude));
    WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(cell_idx_c(force->zSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude)));
    const int64_t _size_force_2 = int64_t(cell_idx_c(force->zSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude));
    const int64_t _stride_force_0 = int64_t(force->xStride());
    const int64_t _stride_force_1 = int64_t(force->yStride());
    const int64_t _stride_force_2 = int64_t(force->zStride());
    const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_kernel_collide::kernel_collide(_data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, block_offset_0, block_offset_1, block_offset_2, omega_bulk, omega_even, omega_odd, omega_shear, seed, temperature, time_step);
}


void FluctuatingMRT_LatticeModel::Sweep::stream( IBlock * block, const uint_t numberOfGhostLayersToInclude )
{
    auto pdfs = block->getData< field::GhostLayerField<double, 19> >(pdfsID);
    field::GhostLayerField<double, 19> * pdfs_tmp;
    {
        // Getting temporary field pdfs_tmp
        auto it = cache_pdfs_.find( pdfs );
        if( it != cache_pdfs_.end() )
        {
            pdfs_tmp = *it;
        }
        else
        {
            pdfs_tmp = pdfs->cloneUninitialized();
            cache_pdfs_.insert(pdfs_tmp);
        }
    }


    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -int_c(pdfs->nrOfGhostLayers()));
    double * RESTRICT const _data_pdfs = pdfs->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -int_c(pdfs_tmp->nrOfGhostLayers()));
    double * RESTRICT _data_pdfs_tmp = pdfs_tmp->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(cell_idx_c(pdfs->xSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(pdfs->xSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(cell_idx_c(pdfs->ySize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(pdfs->ySize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(cell_idx_c(pdfs->zSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(pdfs->zSize()) + 2*cell_idx_c(numberOfGhostLayersToInclude) + 2);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    const int64_t _stride_pdfs_tmp_0 = int64_t(pdfs_tmp->xStride());
    const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
    const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
    const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
    internal_kernel_stream::kernel_stream(_data_pdfs, _data_pdfs_tmp, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3);
    
    pdfs->swapDataPointers(pdfs_tmp);

}


} // namespace lbm
} // namespace walberla




// Buffer Packing

namespace walberla {
namespace mpi {

mpi::SendBuffer & operator<< (mpi::SendBuffer & buf, const ::walberla::lbm::FluctuatingMRT_LatticeModel & lm)
{
    buf << lm.currentLevel;
    return buf;
}

mpi::RecvBuffer & operator>> (mpi::RecvBuffer & buf, ::walberla::lbm::FluctuatingMRT_LatticeModel & lm)
{
    buf >> lm.currentLevel;
    return buf;
}


} // namespace mpi
} // namespace walberla

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif