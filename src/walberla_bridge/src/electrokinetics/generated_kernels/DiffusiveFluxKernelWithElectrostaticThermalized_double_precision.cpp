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
//! \\file DiffusiveFluxKernelWithElectrostaticThermalized_double_precision.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.3, lbmpy v1.3.3, lbmpy_walberla/pystencils_walberla from waLBerla commit b0842e1a493ce19ef1bbb8d2cf382fc343970a7f

#include <cmath>

#include "DiffusiveFluxKernelWithElectrostaticThermalized_double_precision.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#include "philox_rand.h"

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

namespace internal_1ba80591fca01bd4852b82a8b7a2ca49 {
static FUNC_PREFIX void diffusivefluxkernelwithelectrostaticthermalized_double_precision_diffusivefluxkernelwithelectrostaticthermalized_double_precision(double D, double *RESTRICT const _data_j, double *RESTRICT const _data_phi, double *RESTRICT const _data_rho, int64_t const _size_j_0, int64_t const _size_j_1, int64_t const _size_j_2, int64_t const _stride_j_0, int64_t const _stride_j_1, int64_t const _stride_j_2, int64_t const _stride_j_3, int64_t const _stride_phi_0, int64_t const _stride_phi_1, int64_t const _stride_phi_2, int64_t const _stride_rho_0, int64_t const _stride_rho_1, int64_t const _stride_rho_2, uint32_t block_offset_0, uint32_t block_offset_1, uint32_t block_offset_2, double f_ext_0, double f_ext_1, double f_ext_2, uint32_t field_size_0, uint32_t field_size_1, uint32_t field_size_2, double kT, uint32_t seed, uint32_t time_step, double z) {
  {
    {
      {
        if (0 < _size_j_1 - 1 && 0 < _size_j_2 - 1) {

          double random_13_0;
          double random_13_1;
          philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

          double random_12_0;
          double random_12_1;
          philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

          double random_11_0;
          double random_11_1;
          philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

          double random_10_0;
          double random_10_1;
          philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

          double random_9_0;
          double random_9_1;
          philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

          double random_8_0;
          double random_8_1;
          philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

          double random_7_0;
          double random_7_1;
          philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

          _data_j[_stride_j_0 + 12 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0] + _data_rho[_stride_rho_1 + _stride_rho_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0] + _data_rho[_stride_rho_1 + _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0] + _data_rho[_stride_rho_1 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0] - _data_rho[_stride_rho_1 + _stride_rho_2]) * 4.0 - z * (_data_rho[_stride_rho_0] + _data_rho[_stride_rho_1 + _stride_rho_2]) * (-_data_phi[0] + _data_phi[_stride_phi_0 + _stride_phi_1 + _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_2] + _data_phi[_stride_phi_1]) + z * (_data_rho[_stride_rho_0] + _data_rho[_stride_rho_1 + _stride_rho_2]) * (_data_phi[0] - _data_phi[_stride_phi_0 + _stride_phi_1 + _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1] - _data_phi[_stride_phi_2]) + z * (_data_rho[_stride_rho_0] + _data_rho[_stride_rho_1 + _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1] + _data_phi[_stride_phi_0 + _stride_phi_2] - _data_phi[_stride_phi_1] - _data_phi[_stride_phi_2])) * 0.02351606505734748 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0] + 0.5 * _data_rho[_stride_rho_1 + _stride_rho_2]), 0.5) * (random_13_0 - 0.5) * 1.5025119784898082;
        }
        for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
          if (0 < _size_j_1 - 1 && 0 < _size_j_2 - 1) {

            double random_13_0;
            double random_13_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

            double random_12_0;
            double random_12_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

            double random_11_0;
            double random_11_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

            double random_10_0;
            double random_10_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

            double random_9_0;
            double random_9_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

            double random_8_0;
            double random_8_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

            double random_7_0;
            double random_7_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

            _data_j[_stride_j_0 * ctr_0 + 12 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0]) * 2.0 + kT * (-_data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0]) * (-_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0]) - z * (_data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0]) + z * (_data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_2])) * 0.02351606505734748 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0]), 0.5) * (random_13_0 - 0.5) * 1.5025119784898082;
          }
        }
        if (0 < _size_j_1 - 1 && 0 < _size_j_2 - 1) {

          double random_13_0;
          double random_13_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

          double random_12_0;
          double random_12_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

          double random_11_0;
          double random_11_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

          double random_10_0;
          double random_10_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

          double random_9_0;
          double random_9_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

          double random_8_0;
          double random_8_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

          double random_7_0;
          double random_7_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

          _data_j[_stride_j_0 * (_size_j_0 - 1) + 12 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1)]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1)]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1)]) * 2.0 + kT * (-_data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1)]) * 4.0 + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1)]) * (-_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0]) - z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1)]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 + _stride_phi_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0]) + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1)]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_2])) * 0.02351606505734748 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1)]), 0.5) * (random_13_0 - 0.5) * 1.5025119784898082;
        }
      }
      for (int64_t ctr_1 = 1; ctr_1 < _size_j_1 - 1; ctr_1 += 1) {
        {
          {
            if (ctr_1 > 0 && 0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + 6 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2]) * -2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] - _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2] + _data_phi[_stride_phi_1 * ctr_1]) + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2] - _data_phi[_stride_phi_1 * ctr_1])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + 0.5 * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2]), 0.5) * (random_10_0 - 0.5) * 1.6628028407278295;
            }
            if (ctr_1 > 0 && 0 < _size_j_2 - 1 && 1 < _size_j_0 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + 8 * _stride_j_3] = D * (f_ext_1 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1]) * -2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1]) * 2.0 + kT * (-_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1]) * (-_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1]) + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + 0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1]), 0.5) * (random_11_0 - 0.5) * 1.6628028407278295;
            }
            if (ctr_1 > 0 && 0 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + 10 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] - f_ext_0 * z * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] - f_ext_1 * z * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + f_ext_2 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + f_ext_2 * z * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1] * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1] * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] - z * _data_phi[_stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2] * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] - z * _data_phi[_stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2] * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + 0.5 * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2]), 0.5) * (random_12_0 - 0.5) * 1.5025119784898082;
            }
            if (0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + 12 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] - _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2] + _data_phi[_stride_phi_1 * ctr_1]) - z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] + _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_1] - _data_phi[_stride_phi_1 * ctr_1]) + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_1] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2])) * 0.02351606505734748 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + 0.5 * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]), 0.5) * (random_13_0 - 0.5) * 1.5025119784898082;
            }
          }
          for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
            if (ctr_1 > 0 && 0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + 6 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2]) * -2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1]) + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2]), 0.5) * (random_10_0 - 0.5) * 1.6628028407278295;
            }
            if (ctr_1 > 0 && 0 < _size_j_2 - 1 && ctr_0 < _size_j_0 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + 8 * _stride_j_3] = D * (f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1]) * -2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1]) * 2.0 + kT * (-_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1]) * (-_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1]) + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1]), 0.5) * (random_11_0 - 0.5) * 1.6628028407278295;
            }
            if (ctr_1 > 0 && 0 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + 10 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] - f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] - f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2]), 0.5) * (random_12_0 - 0.5) * 1.5025119784898082;
            }
            if (0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + 12 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_1] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1]) - z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1]) + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_1] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2])) * 0.02351606505734748 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]), 0.5) * (random_13_0 - 0.5) * 1.5025119784898082;
            }
          }
          {
            if (ctr_1 > 0 && 0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + 6 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2]) * -2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1]) + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2]), 0.5) * (random_10_0 - 0.5) * 1.6628028407278295;
            }
            if (ctr_1 > 0 && 0 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + 10 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] - f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] - f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2]), 0.5) * (random_12_0 - 0.5) * 1.5025119784898082;
            }
            if (0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + 12 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_1] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1]) - z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1]) + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_1] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2])) * 0.02351606505734748 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]), 0.5) * (random_13_0 - 0.5) * 1.5025119784898082;
            }
          }
        }
      }
      {
        {
          if (_size_j_1 - 1 > 0 && 0 < _size_j_2 - 1 && 1 < _size_j_0 - 1) {

            double random_13_0;
            double random_13_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

            double random_12_0;
            double random_12_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

            double random_11_0;
            double random_11_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

            double random_10_0;
            double random_10_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

            double random_9_0;
            double random_9_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

            double random_8_0;
            double random_8_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

            double random_7_0;
            double random_7_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

            _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + 8 * _stride_j_3] = D * (f_ext_1 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1)]) * -2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1)]) * 2.0 + kT * (-_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1)]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1)]) * (-_data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1] + _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1)]) + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1)]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1] + _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1)])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + 0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1)]), 0.5) * (random_11_0 - 0.5) * 1.6628028407278295;
          }
          if (_size_j_1 - 1 > 0 && 0 < _size_j_2 - 1) {

            double random_13_0;
            double random_13_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

            double random_12_0;
            double random_12_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

            double random_11_0;
            double random_11_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

            double random_10_0;
            double random_10_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

            double random_9_0;
            double random_9_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

            double random_8_0;
            double random_8_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

            double random_7_0;
            double random_7_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

            _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + 10 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1)] - f_ext_0 * z * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1)] - f_ext_1 * z * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + f_ext_2 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1)] + f_ext_2 * z * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1)] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1)] * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1)] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1)] * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] - z * _data_phi[_stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2] * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1)] - z * _data_phi[_stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2] * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1)] + 0.5 * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2]), 0.5) * (random_12_0 - 0.5) * 1.5025119784898082;
          }
        }
        for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
          if (_size_j_1 - 1 > 0 && 0 < _size_j_2 - 1 && ctr_0 < _size_j_0 - 1) {

            double random_13_0;
            double random_13_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

            double random_12_0;
            double random_12_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

            double random_11_0;
            double random_11_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

            double random_10_0;
            double random_10_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

            double random_9_0;
            double random_9_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

            double random_8_0;
            double random_8_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

            double random_7_0;
            double random_7_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

            _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + 8 * _stride_j_3] = D * (f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1)]) * -2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1)]) * 2.0 + kT * (-_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1)]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1)]) * (-_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1)]) + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1)]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1)])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1)]), 0.5) * (random_11_0 - 0.5) * 1.6628028407278295;
          }
          if (_size_j_1 - 1 > 0 && 0 < _size_j_2 - 1) {

            double random_13_0;
            double random_13_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

            double random_12_0;
            double random_12_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

            double random_11_0;
            double random_11_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

            double random_10_0;
            double random_10_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

            double random_9_0;
            double random_9_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

            double random_8_0;
            double random_8_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

            double random_7_0;
            double random_7_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

            _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + 10 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1)] - f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1)] - f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1)] + f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1)] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1)] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1)] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1)] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1)] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1)] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2]), 0.5) * (random_12_0 - 0.5) * 1.5025119784898082;
          }
        }
        if (_size_j_1 - 1 > 0 && 0 < _size_j_2 - 1) {

          double random_13_0;
          double random_13_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

          double random_12_0;
          double random_12_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

          double random_11_0;
          double random_11_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

          double random_10_0;
          double random_10_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

          double random_9_0;
          double random_9_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

          double random_8_0;
          double random_8_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

          double random_7_0;
          double random_7_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

          _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * (_size_j_1 - 1) + 10 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1)] - f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1)] - f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1)] + f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1)] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * (_size_j_1 - 1)] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1)] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * (_size_j_1 - 1)] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1)] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1)] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2]), 0.5) * (random_12_0 - 0.5) * 1.5025119784898082;
        }
      }
    }
    for (int64_t ctr_2 = 1; ctr_2 < _size_j_2 - 1; ctr_2 += 1) {
      {
        {
          {
            if (ctr_2 > 0 && 0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_2 * ctr_2 + 4 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2]) * (-_data_phi[_stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_9_0 - 0.5) * 1.6628028407278295;
            }
            if (ctr_2 > 0 && 0 < _size_j_1 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_2 * ctr_2 + 11 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + f_ext_1 * z * _data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + f_ext_1 * z * _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] - f_ext_2 * z * _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - z * _data_phi[_stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_12_1 - 0.5) * 1.5025119784898082;
            }
            if (0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_2 * ctr_2 + 12 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_2 * ctr_2]) - z * (_data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_2 * ctr_2 + _stride_phi_2])) * 0.02351606505734748 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_13_0 - 0.5) * 1.5025119784898082;
            }
          }
          for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
            if (ctr_2 > 0 && 0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_2 * ctr_2 + 4 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * (-_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_9_0 - 0.5) * 1.6628028407278295;
            }
            if (ctr_2 > 0 && 0 < _size_j_1 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_2 * ctr_2 + 11 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_12_1 - 0.5) * 1.5025119784898082;
            }
            if (0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_2 * ctr_2 + 12 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_2 * ctr_2]) - z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_2 * ctr_2 + _stride_phi_2])) * 0.02351606505734748 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_13_0 - 0.5) * 1.5025119784898082;
            }
          }
          {
            if (ctr_2 > 0 && 0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_2 * ctr_2 + 4 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * (-_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_9_0 - 0.5) * 1.6628028407278295;
            }
            if (ctr_2 > 0 && 0 < _size_j_1 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_2 * ctr_2 + 11 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_12_1 - 0.5) * 1.5025119784898082;
            }
            if (0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_2 * ctr_2 + 12 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_2 * ctr_2]) - z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_2 * ctr_2 + _stride_phi_2])) * 0.02351606505734748 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_13_0 - 0.5) * 1.5025119784898082;
            }
          }
        }
        for (int64_t ctr_1 = 1; ctr_1 < _size_j_1 - 1; ctr_1 += 1) {
          {
            {
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2] = D * (-f_ext_0 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]) + kT * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]) * 2.0 + z * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2])) * 0.081462038946841925 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_7_0 - 0.5) * 1.977416969040271;
              }
              if (ctr_1 > 0 && ctr_2 > 0 && 1 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3] = D * (-f_ext_1 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) + kT * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 2.0 + z * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2]) * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2])) * 0.081462038946841925 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_7_1 - 0.5) * 1.977416969040271;
              }
              if (ctr_1 > 0 && ctr_2 > 0 && 1 < _size_j_0 - 1 && ctr_1 < _size_j_1 - 1) {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3] = D * (f_ext_2 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]) + kT * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]) * 2.0 + z * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2])) * -0.081462038946841925 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + 0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_8_0 - 0.5) * 1.977416969040271;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + kT * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_8_1 - 0.5) * 1.6628028407278295;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * (-_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_9_0 - 0.5) * 1.6628028407278295;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 2.0 + kT * (-_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) - z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2])) * -0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_9_1 - 0.5) * 1.6628028407278295;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * -2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_10_0 - 0.5) * 1.6628028407278295;
              }
              if (ctr_1 > 0 && ctr_2 > 0 && 1 < _size_j_0 - 1) {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3] = D * (f_ext_1 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 2.0 + kT * (-_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2]) - z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2])) * -0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_10_1 - 0.5) * 1.6628028407278295;
              }
              if (ctr_1 > 0 && 1 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3] = D * (f_ext_1 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * -2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_11_0 - 0.5) * 1.6628028407278295;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_1 * z * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_2 * z * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - z * _data_phi[_stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_11_1 - 0.5) * 1.5025119784898082;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_1 * z * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + f_ext_2 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + f_ext_2 * z * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] - z * _data_phi[_stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_12_0 - 0.5) * 1.5025119784898082;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + f_ext_1 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + f_ext_1 * z * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_2 * z * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - z * _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_12_1 - 0.5) * 1.5025119784898082;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) - z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2])) * 0.02351606505734748 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_13_0 - 0.5) * 1.5025119784898082;
              }
            }
            for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2] = D * (-f_ext_0 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]) + kT * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]) * 2.0 + z * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2])) * 0.081462038946841925 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_7_0 - 0.5) * 1.977416969040271;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3] = D * (-f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) + kT * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 2.0 + z * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2]) * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2])) * 0.081462038946841925 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_7_1 - 0.5) * 1.977416969040271;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3] = D * (f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]) + kT * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]) * 2.0 + z * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2])) * -0.081462038946841925 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_8_0 - 0.5) * 1.977416969040271;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + kT * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_8_1 - 0.5) * 1.6628028407278295;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * (-_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_9_0 - 0.5) * 1.6628028407278295;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 2.0 + kT * (-_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) - z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2])) * -0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_9_1 - 0.5) * 1.6628028407278295;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * -2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_10_0 - 0.5) * 1.6628028407278295;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3] = D * (f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 2.0 + kT * (-_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2]) - z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2])) * -0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_10_1 - 0.5) * 1.6628028407278295;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3] = D * (f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * -2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_11_0 - 0.5) * 1.6628028407278295;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_11_1 - 0.5) * 1.5025119784898082;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_12_0 - 0.5) * 1.5025119784898082;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_12_1 - 0.5) * 1.5025119784898082;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) - z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2])) * 0.02351606505734748 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_13_0 - 0.5) * 1.5025119784898082;
              }
            }
            {
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2] = D * (-f_ext_0 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]) + kT * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]) * 2.0 + z * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2])) * 0.081462038946841925 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_7_0 - 0.5) * 1.977416969040271;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + kT * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_8_1 - 0.5) * 1.6628028407278295;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * (-_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_9_0 - 0.5) * 1.6628028407278295;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 2.0 + kT * (-_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) - z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2])) * -0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_9_1 - 0.5) * 1.6628028407278295;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * -2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_10_0 - 0.5) * 1.6628028407278295;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_11_1 - 0.5) * 1.5025119784898082;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_12_0 - 0.5) * 1.5025119784898082;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_12_1 - 0.5) * 1.5025119784898082;
              }
              {

                double random_13_0;
                double random_13_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

                double random_12_0;
                double random_12_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

                double random_11_0;
                double random_11_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

                double random_10_0;
                double random_10_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

                double random_9_0;
                double random_9_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

                double random_8_0;
                double random_8_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

                double random_7_0;
                double random_7_1;
                philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) - z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * ctr_2 + _stride_phi_2])) * 0.02351606505734748 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_13_0 - 0.5) * 1.5025119784898082;
              }
            }
          }
        }
        {
          {
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && 1 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + _stride_j_3] = D * (-f_ext_1 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) + kT * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 2.0 + z * (_data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2]) * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2])) * 0.081462038946841925 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_7_1 - 0.5) * 1.977416969040271;
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 3 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + kT * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_8_1 - 0.5) * 1.6628028407278295;
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && 1 < _size_j_0 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 7 * _stride_j_3] = D * (f_ext_1 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 2.0 + kT * (-_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2]) - z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2 - _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2])) * -0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_10_1 - 0.5) * 1.6628028407278295;
            }
            if (_size_j_1 - 1 > 0 && 1 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 8 * _stride_j_3] = D * (f_ext_1 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * -2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_11_0 - 0.5) * 1.6628028407278295;
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 9 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - f_ext_1 * z * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - f_ext_2 * z * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - z * _data_phi[_stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_11_1 - 0.5) * 1.5025119784898082;
            }
            if (_size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 10 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - f_ext_1 * z * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + f_ext_2 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + f_ext_2 * z * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] - z * _data_phi[_stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_12_0 - 0.5) * 1.5025119784898082;
            }
          }
          for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && ctr_0 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + _stride_j_3] = D * (-f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) + kT * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 2.0 + z * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2]) * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2])) * 0.081462038946841925 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_7_1 - 0.5) * 1.977416969040271;
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 3 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + kT * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_8_1 - 0.5) * 1.6628028407278295;
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && ctr_0 < _size_j_0 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 7 * _stride_j_3] = D * (f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 2.0 + kT * (-_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2]) - z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2 - _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2])) * -0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_10_1 - 0.5) * 1.6628028407278295;
            }
            if (_size_j_1 - 1 > 0 && ctr_0 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 8 * _stride_j_3] = D * (f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * -2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 2.0 + kT * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (-_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2 + _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_11_0 - 0.5) * 1.6628028407278295;
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 9 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_11_1 - 0.5) * 1.5025119784898082;
            }
            if (_size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 10 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_12_0 - 0.5) * 1.5025119784898082;
            }
          }
          {
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 3 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + f_ext_1 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * -2.0 + kT * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2]) + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2])) * 0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5) * (random_8_1 - 0.5) * 1.6628028407278295;
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 9 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 - _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5) * (random_11_1 - 0.5) * 1.5025119784898082;
            }
            if (_size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 10 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * ctr_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * ctr_2 + _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5) * (random_12_0 - 0.5) * 1.5025119784898082;
            }
          }
        }
      }
    }
    {
      {
        if (_size_j_2 - 1 > 0 && 0 < _size_j_1 - 1) {

          double random_13_0;
          double random_13_1;
          philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

          double random_12_0;
          double random_12_1;
          philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

          double random_11_0;
          double random_11_1;
          philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

          double random_10_0;
          double random_10_1;
          philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

          double random_9_0;
          double random_9_1;
          philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

          double random_8_0;
          double random_8_1;
          philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

          double random_7_0;
          double random_7_1;
          philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

          _data_j[_stride_j_0 + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_0 * z * _data_rho[_stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + f_ext_1 * z * _data_rho[_stride_rho_0 + _stride_rho_2 * (_size_j_2 - 1)] + f_ext_1 * z * _data_rho[_stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_2 * z * _data_rho[_stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - z * _data_phi[_stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 + _stride_rho_2 * (_size_j_2 - 1)] - z * _data_phi[_stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_12_1 - 0.5) * 1.5025119784898082;
        }
        for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
          if (_size_j_2 - 1 > 0 && 0 < _size_j_1 - 1) {

            double random_13_0;
            double random_13_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

            double random_12_0;
            double random_12_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

            double random_11_0;
            double random_11_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

            double random_10_0;
            double random_10_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

            double random_9_0;
            double random_9_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

            double random_8_0;
            double random_8_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

            double random_7_0;
            double random_7_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

            _data_j[_stride_j_0 * ctr_0 + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * (_size_j_2 - 1)] + f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * (_size_j_2 - 1)] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_12_1 - 0.5) * 1.5025119784898082;
          }
        }
        if (_size_j_2 - 1 > 0 && 0 < _size_j_1 - 1) {

          double random_13_0;
          double random_13_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

          double random_12_0;
          double random_12_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

          double random_11_0;
          double random_11_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

          double random_10_0;
          double random_10_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

          double random_9_0;
          double random_9_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

          double random_8_0;
          double random_8_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

          double random_7_0;
          double random_7_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

          _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_12_1 - 0.5) * 1.5025119784898082;
        }
      }
      for (int64_t ctr_1 = 1; ctr_1 < _size_j_1 - 1; ctr_1 += 1) {
        {
          {
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && 1 < _size_j_0 - 1 && ctr_1 < _size_j_1 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 2 * _stride_j_3] = D * (f_ext_2 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)]) + kT * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)]) * 2.0 + z * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)]) * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)])) * -0.081462038946841925 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + 0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)]), 0.5) * (random_8_0 - 0.5) * 1.977416969040271;
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 5 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 2.0 + kT * (-_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] + _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)]) - z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)])) * -0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_9_1 - 0.5) * 1.6628028407278295;
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && 1 < _size_j_0 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 7 * _stride_j_3] = D * (f_ext_1 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 2.0 + kT * (-_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1)]) - z * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1)])) * -0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_10_1 - 0.5) * 1.6628028407278295;
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_0 * z * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_1 * z * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_2 * z * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - z * _data_phi[_stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - z * _data_phi[_stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_11_1 - 0.5) * 1.5025119784898082;
            }
            if (_size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_0 * z * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + f_ext_1 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + f_ext_1 * z * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_2 * z * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - z * _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - z * _data_phi[_stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_12_1 - 0.5) * 1.5025119784898082;
            }
          }
          for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && ctr_0 < _size_j_0 - 1 && ctr_1 < _size_j_1 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 2 * _stride_j_3] = D * (f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)]) + kT * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)]) * 2.0 + z * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)]) * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)])) * -0.081462038946841925 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)]), 0.5) * (random_8_0 - 0.5) * 1.977416969040271;
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 5 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 2.0 + kT * (-_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] + _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)]) - z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)])) * -0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_9_1 - 0.5) * 1.6628028407278295;
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && ctr_0 < _size_j_0 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 7 * _stride_j_3] = D * (f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 2.0 + kT * (-_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1)]) - z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1)])) * -0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_10_1 - 0.5) * 1.6628028407278295;
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_11_1 - 0.5) * 1.5025119784898082;
            }
            if (_size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_12_1 - 0.5) * 1.5025119784898082;
            }
          }
          {
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 5 * _stride_j_3] = D * (f_ext_0 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 2.0 + kT * (-_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)]) - z * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * (_data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] + _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)])) * -0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_9_1 - 0.5) * 1.6628028407278295;
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_11_1 - 0.5) * 1.5025119784898082;
            }
            if (_size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {

              double random_13_0;
              double random_13_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

              double random_12_0;
              double random_12_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

              double random_11_0;
              double random_11_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

              double random_10_0;
              double random_10_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

              double random_9_0;
              double random_9_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

              double random_8_0;
              double random_8_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

              double random_7_0;
              double random_7_1;
              philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * ctr_1 + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * ctr_1 + _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_12_1 - 0.5) * 1.5025119784898082;
            }
          }
        }
      }
      {
        {
          if (_size_j_1 - 1 > 0 && _size_j_2 - 1 > 0 && 1 < _size_j_0 - 1) {

            double random_13_0;
            double random_13_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

            double random_12_0;
            double random_12_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

            double random_11_0;
            double random_11_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

            double random_10_0;
            double random_10_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

            double random_9_0;
            double random_9_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

            double random_8_0;
            double random_8_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

            double random_7_0;
            double random_7_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

            _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * (_size_j_2 - 1) + 7 * _stride_j_3] = D * (f_ext_1 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 2.0 + kT * (-_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * (_size_j_2 - 1)] + _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1)]) - z * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * (_data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] + _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * (_size_j_2 - 1)] - _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1)])) * -0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_10_1 - 0.5) * 1.6628028407278295;
          }
          if (_size_j_1 - 1 > 0 && _size_j_2 - 1 > 0) {

            double random_13_0;
            double random_13_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

            double random_12_0;
            double random_12_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

            double random_11_0;
            double random_11_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

            double random_10_0;
            double random_10_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

            double random_9_0;
            double random_9_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

            double random_8_0;
            double random_8_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

            double random_7_0;
            double random_7_1;
            philox_double2(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

            _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_0 * z * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_1 * z * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_2 * z * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - z * _data_phi[_stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - z * _data_phi[_stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_11_1 - 0.5) * 1.5025119784898082;
          }
        }
        for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
          if (_size_j_1 - 1 > 0 && _size_j_2 - 1 > 0 && ctr_0 < _size_j_0 - 1) {

            double random_13_0;
            double random_13_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

            double random_12_0;
            double random_12_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

            double random_11_0;
            double random_11_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

            double random_10_0;
            double random_10_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

            double random_9_0;
            double random_9_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

            double random_8_0;
            double random_8_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

            double random_7_0;
            double random_7_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

            _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * (_size_j_2 - 1) + 7 * _stride_j_3] = D * (f_ext_1 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 2.0 + f_ext_2 * z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 2.0 + kT * (-_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 4.0 + z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * (_size_j_2 - 1)] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1)]) - z * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * (_data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] + _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * (_size_j_2 - 1)] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] - _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1)])) * -0.028801180074297286 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_10_1 - 0.5) * 1.6628028407278295;
          }
          if (_size_j_1 - 1 > 0 && _size_j_2 - 1 > 0) {

            double random_13_0;
            double random_13_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

            double random_12_0;
            double random_12_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

            double random_11_0;
            double random_11_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

            double random_10_0;
            double random_10_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

            double random_9_0;
            double random_9_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

            double random_8_0;
            double random_8_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

            double random_7_0;
            double random_7_1;
            philox_double2(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

            _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_0 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_1 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_2 * z * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 * ctr_0 + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - z * _data_phi[_stride_phi_0 * ctr_0 - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_11_1 - 0.5) * 1.5025119784898082;
          }
        }
        if (_size_j_1 - 1 > 0 && _size_j_2 - 1 > 0) {

          double random_13_0;
          double random_13_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 6, seed, random_13_0, random_13_1);

          double random_12_0;
          double random_12_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 5, seed, random_12_0, random_12_1);

          double random_11_0;
          double random_11_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 4, seed, random_11_0, random_11_1);

          double random_10_0;
          double random_10_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_10_0, random_10_1);

          double random_9_0;
          double random_9_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_9_0, random_9_1);

          double random_8_0;
          double random_8_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_8_0, random_8_1);

          double random_7_0;
          double random_7_1;
          philox_double2(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_7_0, random_7_1);

          _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3] = D * (-f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_0 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_1 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - f_ext_2 * z * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * -2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + kT * 2.0 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) + _stride_phi_1 * (_size_j_1 - 1) + _stride_phi_2 * (_size_j_2 - 1)] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - z * _data_phi[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0 + _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2] * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.04703213011469496 * ((1.0) / (kT)) + pow(D * (0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + 0.5 * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5) * (random_11_1 - 0.5) * 1.5025119784898082;
        }
      }
    }
  }
}
} // namespace internal_1ba80591fca01bd4852b82a8b7a2ca49

void DiffusiveFluxKernelWithElectrostaticThermalized_double_precision::run(IBlock *block) {
  if (!this->configured_)
    WALBERLA_ABORT("This Sweep contains a configure function that needs to be called manually")

  auto j = block->getData<field::GhostLayerField<double, 13>>(jID);
  auto phi = block->getData<field::GhostLayerField<double, 1>>(phiID);
  auto rho = block->getData<field::GhostLayerField<double, 1>>(rhoID);

  auto &block_offset_1 = this->block_offset_1_;
  auto &f_ext_2 = this->f_ext_2_;
  auto &z = this->z_;
  auto &seed = this->seed_;
  auto &time_step = this->time_step_;
  auto &kT = this->kT_;
  auto &field_size_0 = this->field_size_0_;
  auto &f_ext_1 = this->f_ext_1_;
  auto &f_ext_0 = this->f_ext_0_;
  auto &D = this->D_;
  auto &field_size_1 = this->field_size_1_;
  auto &field_size_2 = this->field_size_2_;
  auto &block_offset_0 = this->block_offset_0_;
  auto &block_offset_2 = this->block_offset_2_;
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(j->nrOfGhostLayers()))
  double *RESTRICT const _data_j = j->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(phi->nrOfGhostLayers()))
  double *RESTRICT const _data_phi = phi->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho->nrOfGhostLayers()))
  double *RESTRICT const _data_rho = rho->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(), int64_t(int64_c(j->xSize()) + 2))
  const int64_t _size_j_0 = int64_t(int64_c(j->xSize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(), int64_t(int64_c(j->ySize()) + 2))
  const int64_t _size_j_1 = int64_t(int64_c(j->ySize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(), int64_t(int64_c(j->zSize()) + 2))
  const int64_t _size_j_2 = int64_t(int64_c(j->zSize()) + 2);
  const int64_t _stride_j_0 = int64_t(j->xStride());
  const int64_t _stride_j_1 = int64_t(j->yStride());
  const int64_t _stride_j_2 = int64_t(j->zStride());
  const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
  const int64_t _stride_phi_0 = int64_t(phi->xStride());
  const int64_t _stride_phi_1 = int64_t(phi->yStride());
  const int64_t _stride_phi_2 = int64_t(phi->zStride());
  const int64_t _stride_rho_0 = int64_t(rho->xStride());
  const int64_t _stride_rho_1 = int64_t(rho->yStride());
  const int64_t _stride_rho_2 = int64_t(rho->zStride());
  internal_1ba80591fca01bd4852b82a8b7a2ca49::diffusivefluxkernelwithelectrostaticthermalized_double_precision_diffusivefluxkernelwithelectrostaticthermalized_double_precision(D, _data_j, _data_phi, _data_rho, _size_j_0, _size_j_1, _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2, _stride_j_3, _stride_phi_0, _stride_phi_1, _stride_phi_2, _stride_rho_0, _stride_rho_1, _stride_rho_2, block_offset_0, block_offset_1, block_offset_2, f_ext_0, f_ext_1, f_ext_2, field_size_0, field_size_1, field_size_2, kT, seed, time_step, z);
}

void DiffusiveFluxKernelWithElectrostaticThermalized_double_precision::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block) {
  if (!this->configured_)
    WALBERLA_ABORT("This Sweep contains a configure function that needs to be called manually")

  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto j = block->getData<field::GhostLayerField<double, 13>>(jID);
  auto phi = block->getData<field::GhostLayerField<double, 1>>(phiID);
  auto rho = block->getData<field::GhostLayerField<double, 1>>(rhoID);

  auto &block_offset_1 = this->block_offset_1_;
  auto &f_ext_2 = this->f_ext_2_;
  auto &z = this->z_;
  auto &seed = this->seed_;
  auto &time_step = this->time_step_;
  auto &kT = this->kT_;
  auto &field_size_0 = this->field_size_0_;
  auto &f_ext_1 = this->f_ext_1_;
  auto &f_ext_0 = this->f_ext_0_;
  auto &D = this->D_;
  auto &field_size_1 = this->field_size_1_;
  auto &field_size_2 = this->field_size_2_;
  auto &block_offset_0 = this->block_offset_0_;
  auto &block_offset_2 = this->block_offset_2_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(j->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(j->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(j->nrOfGhostLayers()))
  double *RESTRICT const _data_j = j->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(phi->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(phi->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(phi->nrOfGhostLayers()))
  double *RESTRICT const _data_phi = phi->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho->nrOfGhostLayers()))
  double *RESTRICT const _data_rho = rho->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 2))
  const int64_t _size_j_0 = int64_t(int64_c(ci.xSize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 2))
  const int64_t _size_j_1 = int64_t(int64_c(ci.ySize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 2))
  const int64_t _size_j_2 = int64_t(int64_c(ci.zSize()) + 2);
  const int64_t _stride_j_0 = int64_t(j->xStride());
  const int64_t _stride_j_1 = int64_t(j->yStride());
  const int64_t _stride_j_2 = int64_t(j->zStride());
  const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
  const int64_t _stride_phi_0 = int64_t(phi->xStride());
  const int64_t _stride_phi_1 = int64_t(phi->yStride());
  const int64_t _stride_phi_2 = int64_t(phi->zStride());
  const int64_t _stride_rho_0 = int64_t(rho->xStride());
  const int64_t _stride_rho_1 = int64_t(rho->yStride());
  const int64_t _stride_rho_2 = int64_t(rho->zStride());
  internal_1ba80591fca01bd4852b82a8b7a2ca49::diffusivefluxkernelwithelectrostaticthermalized_double_precision_diffusivefluxkernelwithelectrostaticthermalized_double_precision(D, _data_j, _data_phi, _data_rho, _size_j_0, _size_j_1, _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2, _stride_j_3, _stride_phi_0, _stride_phi_1, _stride_phi_2, _stride_rho_0, _stride_rho_1, _stride_rho_2, block_offset_0, block_offset_1, block_offset_2, f_ext_0, f_ext_1, f_ext_2, field_size_0, field_size_1, field_size_2, kT, seed, time_step, z);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif
