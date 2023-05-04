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
//! \\file DiffusiveFluxKernelWithElectrostatic_single_precision.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

// kernel generated with pystencils v1.2, lbmpy v1.2, lbmpy_walberla/pystencils_walberla from waLBerla commit ref: a839fac6ef7d0c58e7710e4d50490e9dd7146b4a

#include <cmath>

#include "DiffusiveFluxKernelWithElectrostatic_single_precision.h"
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

namespace internal_823ab2463d465630661d5edc8f90930c {
static FUNC_PREFIX void diffusivefluxkernelwithelectrostatic_single_precision_diffusivefluxkernelwithelectrostatic_single_precision(float D, float *RESTRICT const _data_j, float *RESTRICT const _data_phi, float *RESTRICT const _data_rho, int64_t const _size_j_0, int64_t const _size_j_1, int64_t const _size_j_2, int64_t const _stride_j_0, int64_t const _stride_j_1, int64_t const _stride_j_2, int64_t const _stride_j_3, int64_t const _stride_phi_0, int64_t const _stride_phi_1, int64_t const _stride_phi_2, int64_t const _stride_rho_0, int64_t const _stride_rho_1, int64_t const _stride_rho_2, float f_ext_0, float f_ext_1, float f_ext_2, float kT, float z) {
  {
    {
      {
        if (0 < _size_j_1 - 1 && 0 < _size_j_2 - 1) {
          float *RESTRICT _data_j_20_312 = _data_j + 12 * _stride_j_3;
          float *RESTRICT _data_j_20_312_10 = _data_j_20_312;
          float *RESTRICT _data_phi_20 = _data_phi;
          float *RESTRICT _data_phi_20_10 = _data_phi_20;
          float *RESTRICT _data_rho_20 = _data_rho;
          float *RESTRICT _data_rho_20_10 = _data_rho_20;
          float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
          float *RESTRICT _data_rho_21_11 = _stride_rho_1 + _data_rho_21;
          float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
          float *RESTRICT _data_phi_21_11 = _stride_phi_1 + _data_phi_21;
          _data_j_20_312_10[_stride_j_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_0 * z * -1.0f * _data_rho_21_11[0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0] + f_ext_1 * z * _data_rho_21_11[0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0] + f_ext_2 * z * _data_rho_21_11[0] + kT * -2.0f * _data_rho_21_11[0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_21_11[0] * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_21_11[0] * _data_rho_21_11[0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_20_10[_stride_rho_0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_21_11[0]) * 0.04703213011469496f * ((1.0f) / (kT));
        }
        for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
          if (0 < _size_j_1 - 1 && 0 < _size_j_2 - 1) {
            float *RESTRICT _data_j_20_312 = _data_j + 12 * _stride_j_3;
            float *RESTRICT _data_j_20_312_10 = _data_j_20_312;
            float *RESTRICT _data_phi_20 = _data_phi;
            float *RESTRICT _data_phi_20_10 = _data_phi_20;
            float *RESTRICT _data_rho_20 = _data_rho;
            float *RESTRICT _data_rho_20_10 = _data_rho_20;
            float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
            float *RESTRICT _data_rho_21_11 = _stride_rho_1 + _data_rho_21;
            float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
            float *RESTRICT _data_phi_21_11 = _stride_phi_1 + _data_phi_21;
            _data_j_20_312_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_0 * z * -1.0f * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_1 * z * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_2 * z * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * -2.0f * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_21_11[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_21_11[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
          }
        }
        if (0 < _size_j_1 - 1 && 0 < _size_j_2 - 1) {
          float *RESTRICT _data_j_20_312 = _data_j + 12 * _stride_j_3;
          float *RESTRICT _data_j_20_312_10 = _data_j_20_312;
          float *RESTRICT _data_phi_20 = _data_phi;
          float *RESTRICT _data_phi_20_10 = _data_phi_20;
          float *RESTRICT _data_rho_20 = _data_rho;
          float *RESTRICT _data_rho_20_10 = _data_rho_20;
          float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
          float *RESTRICT _data_rho_21_11 = _stride_rho_1 + _data_rho_21;
          float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
          float *RESTRICT _data_phi_21_11 = _stride_phi_1 + _data_phi_21;
          _data_j_20_312_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_0 * z * -1.0f * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_1 * z * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_2 * z * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * -2.0f * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_21_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_21_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
        }
      }
      for (int64_t ctr_1 = 1; ctr_1 < _size_j_1 - 1; ctr_1 += 1) {
        {
          {
            if (ctr_1 > 0 && 0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {
              float *RESTRICT _data_j_20_36 = _data_j + 6 * _stride_j_3;
              float *RESTRICT _data_j_20_36_10 = _stride_j_1 * ctr_1 + _data_j_20_36;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
              float *RESTRICT _data_rho_21_10 = _stride_rho_1 * ctr_1 + _data_rho_21;
              float *RESTRICT _data_rho_20 = _data_rho;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
              float *RESTRICT _data_phi_21_10 = _stride_phi_1 * ctr_1 + _data_phi_21;
              float *RESTRICT _data_phi_20 = _data_phi;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              _data_j_20_36_10[_stride_j_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_10[0]) * -2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_10[0]) * 2.0f + kT * (-1.0f * _data_rho_21_10[0] + _data_rho_20_10[_stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_10[0]) * (-1.0f * _data_phi_20_10[0] - 1.0f * _data_phi_21_10[0] + _data_phi_20_10[_stride_phi_0] + _data_phi_21_10[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_10[0]) * (-1.0f * _data_phi_21_10[0] - 1.0f * _data_phi_21_10[_stride_phi_0] + _data_phi_20_10[0] + _data_phi_20_10[_stride_phi_0])) * 0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_1 > 0 && 0 < _size_j_2 - 1 && 1 < _size_j_0 - 1) {
              float *RESTRICT _data_j_20_38 = _data_j + 8 * _stride_j_3;
              float *RESTRICT _data_j_20_38_10 = _stride_j_1 * ctr_1 + _data_j_20_38;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
              float *RESTRICT _data_rho_21_1m1 = _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_21;
              float *RESTRICT _data_rho_20 = _data_rho;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_phi_20 = _data_phi;
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
              float *RESTRICT _data_phi_21_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_21;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_phi_21_10 = _stride_phi_1 * ctr_1 + _data_phi_21;
              _data_j_20_38_10[_stride_j_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_1m1[_stride_rho_0]) * -2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_1m1[_stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_21_1m1[_stride_rho_0] + _data_rho_20_10[_stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_1m1[_stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0] - 1.0f * _data_phi_21_1m1[_stride_phi_0] + _data_phi_20_10[_stride_phi_0] + _data_phi_21_10[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_1m1[_stride_rho_0]) * (-1.0f * _data_phi_21_10[_stride_phi_0] - 1.0f * _data_phi_21_1m1[_stride_phi_0] + _data_phi_20_10[_stride_phi_0] + _data_phi_20_1m1[_stride_phi_0])) * 0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_1 > 0 && 0 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_310 = _data_j + 10 * _stride_j_3;
              float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
              float *RESTRICT _data_phi_20 = _data_phi;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_rho_20 = _data_rho;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
              float *RESTRICT _data_rho_21_1m1 = _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_21;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
              float *RESTRICT _data_phi_21_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_21;
              _data_j_20_310_10[_stride_j_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_0 * z * -1.0f * _data_rho_21_1m1[0] + f_ext_1 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_1 * z * -1.0f * _data_rho_21_1m1[0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0] + f_ext_2 * z * _data_rho_21_1m1[0] + kT * -2.0f * _data_rho_21_1m1[0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_21_1m1[0] * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_21_1m1[0] * _data_rho_21_1m1[0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_20_10[_stride_rho_0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_21_1m1[0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
            if (0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {
              float *RESTRICT _data_j_20_312 = _data_j + 12 * _stride_j_3;
              float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
              float *RESTRICT _data_phi_20 = _data_phi;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_rho_20 = _data_rho;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
              float *RESTRICT _data_rho_21_11 = _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_21;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
              float *RESTRICT _data_phi_21_11 = _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_21;
              _data_j_20_312_10[_stride_j_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_0 * z * -1.0f * _data_rho_21_11[0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0] + f_ext_1 * z * _data_rho_21_11[0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0] + f_ext_2 * z * _data_rho_21_11[0] + kT * -2.0f * _data_rho_21_11[0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_21_11[0] * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_21_11[0] * _data_rho_21_11[0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_20_10[_stride_rho_0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_21_11[0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
          }
          for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
            if (ctr_1 > 0 && 0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {
              float *RESTRICT _data_j_20_36 = _data_j + 6 * _stride_j_3;
              float *RESTRICT _data_j_20_36_10 = _stride_j_1 * ctr_1 + _data_j_20_36;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
              float *RESTRICT _data_rho_21_10 = _stride_rho_1 * ctr_1 + _data_rho_21;
              float *RESTRICT _data_rho_20 = _data_rho;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
              float *RESTRICT _data_phi_21_10 = _stride_phi_1 * ctr_1 + _data_phi_21;
              float *RESTRICT _data_phi_20 = _data_phi;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              _data_j_20_36_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * -2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_21_10[_stride_rho_0 * ctr_0 - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * ctr_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_21_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_21_10[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_21_10[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_21_10[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * ctr_0])) * 0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_1 > 0 && 0 < _size_j_2 - 1 && ctr_0 < _size_j_0 - 1) {
              float *RESTRICT _data_j_20_38 = _data_j + 8 * _stride_j_3;
              float *RESTRICT _data_j_20_38_10 = _stride_j_1 * ctr_1 + _data_j_20_38;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
              float *RESTRICT _data_rho_21_1m1 = _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_21;
              float *RESTRICT _data_rho_20 = _data_rho;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_phi_20 = _data_phi;
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
              float *RESTRICT _data_phi_21_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_21;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_phi_21_10 = _stride_phi_1 * ctr_1 + _data_phi_21;
              _data_j_20_38_10[_stride_j_0 * ctr_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_1m1[_stride_rho_0 * ctr_0]) * -2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_1m1[_stride_rho_0 * ctr_0]) * 2.0f + kT * (-1.0f * _data_rho_21_1m1[_stride_rho_0 * ctr_0] + _data_rho_20_10[_stride_rho_0 * ctr_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_1m1[_stride_rho_0 * ctr_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] - 1.0f * _data_phi_21_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_21_10[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_1m1[_stride_rho_0 * ctr_0]) * (-1.0f * _data_phi_21_10[_stride_phi_0 * ctr_0] - 1.0f * _data_phi_21_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_20_1m1[_stride_phi_0 * ctr_0])) * 0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_1 > 0 && 0 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_310 = _data_j + 10 * _stride_j_3;
              float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
              float *RESTRICT _data_phi_20 = _data_phi;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_rho_20 = _data_rho;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
              float *RESTRICT _data_rho_21_1m1 = _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_21;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
              float *RESTRICT _data_phi_21_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_21;
              _data_j_20_310_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_0 * z * -1.0f * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_1 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_1 * z * -1.0f * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_2 * z * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * -2.0f * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_21_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_21_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
            if (0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {
              float *RESTRICT _data_j_20_312 = _data_j + 12 * _stride_j_3;
              float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
              float *RESTRICT _data_phi_20 = _data_phi;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_rho_20 = _data_rho;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
              float *RESTRICT _data_rho_21_11 = _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_21;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
              float *RESTRICT _data_phi_21_11 = _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_21;
              _data_j_20_312_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_0 * z * -1.0f * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_1 * z * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_2 * z * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * -2.0f * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_21_11[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_21_11[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
          }
          {
            if (ctr_1 > 0 && 0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {
              float *RESTRICT _data_j_20_36 = _data_j + 6 * _stride_j_3;
              float *RESTRICT _data_j_20_36_10 = _stride_j_1 * ctr_1 + _data_j_20_36;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
              float *RESTRICT _data_rho_21_10 = _stride_rho_1 * ctr_1 + _data_rho_21;
              float *RESTRICT _data_rho_20 = _data_rho;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
              float *RESTRICT _data_phi_21_10 = _stride_phi_1 * ctr_1 + _data_phi_21;
              float *RESTRICT _data_phi_20 = _data_phi;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              _data_j_20_36_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_21_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * -2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_21_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_21_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_21_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_21_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_21_10[_stride_phi_0 * (_size_j_0 - 1)]) + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_21_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_21_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_21_10[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)])) * 0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_1 > 0 && 0 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_310 = _data_j + 10 * _stride_j_3;
              float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
              float *RESTRICT _data_phi_20 = _data_phi;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_rho_20 = _data_rho;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
              float *RESTRICT _data_rho_21_1m1 = _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_21;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
              float *RESTRICT _data_phi_21_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_21;
              _data_j_20_310_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_0 * z * -1.0f * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_1 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_1 * z * -1.0f * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_2 * z * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * -2.0f * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_21_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_21_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
            if (0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {
              float *RESTRICT _data_j_20_312 = _data_j + 12 * _stride_j_3;
              float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
              float *RESTRICT _data_phi_20 = _data_phi;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_rho_20 = _data_rho;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
              float *RESTRICT _data_rho_21_11 = _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_21;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
              float *RESTRICT _data_phi_21_11 = _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_21;
              _data_j_20_312_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_0 * z * -1.0f * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_1 * z * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_2 * z * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * -2.0f * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_21_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_21_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
          }
        }
      }
      {
        {
          if (_size_j_1 - 1 > 0 && 0 < _size_j_2 - 1 && 1 < _size_j_0 - 1) {
            float *RESTRICT _data_j_20_38 = _data_j + 8 * _stride_j_3;
            float *RESTRICT _data_j_20_38_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_38;
            float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
            float *RESTRICT _data_rho_21_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_21;
            float *RESTRICT _data_rho_20 = _data_rho;
            float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
            float *RESTRICT _data_phi_20 = _data_phi;
            float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
            float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
            float *RESTRICT _data_phi_21_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_21;
            float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
            float *RESTRICT _data_phi_21_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_21;
            _data_j_20_38_10[_stride_j_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_1m1[_stride_rho_0]) * -2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_1m1[_stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_21_1m1[_stride_rho_0] + _data_rho_20_10[_stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_1m1[_stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0] - 1.0f * _data_phi_21_1m1[_stride_phi_0] + _data_phi_20_10[_stride_phi_0] + _data_phi_21_10[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_1m1[_stride_rho_0]) * (-1.0f * _data_phi_21_10[_stride_phi_0] - 1.0f * _data_phi_21_1m1[_stride_phi_0] + _data_phi_20_10[_stride_phi_0] + _data_phi_20_1m1[_stride_phi_0])) * 0.028801180074297286f * ((1.0f) / (kT));
          }
          if (_size_j_1 - 1 > 0 && 0 < _size_j_2 - 1) {
            float *RESTRICT _data_j_20_310 = _data_j + 10 * _stride_j_3;
            float *RESTRICT _data_j_20_310_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_310;
            float *RESTRICT _data_phi_20 = _data_phi;
            float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
            float *RESTRICT _data_rho_20 = _data_rho;
            float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
            float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
            float *RESTRICT _data_rho_21_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_21;
            float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
            float *RESTRICT _data_phi_21_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_21;
            _data_j_20_310_10[_stride_j_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_0 * z * -1.0f * _data_rho_21_1m1[0] + f_ext_1 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_1 * z * -1.0f * _data_rho_21_1m1[0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0] + f_ext_2 * z * _data_rho_21_1m1[0] + kT * -2.0f * _data_rho_21_1m1[0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_21_1m1[0] * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_21_1m1[0] * _data_rho_21_1m1[0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_20_10[_stride_rho_0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_21_1m1[0]) * 0.04703213011469496f * ((1.0f) / (kT));
          }
        }
        for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
          if (_size_j_1 - 1 > 0 && 0 < _size_j_2 - 1 && ctr_0 < _size_j_0 - 1) {
            float *RESTRICT _data_j_20_38 = _data_j + 8 * _stride_j_3;
            float *RESTRICT _data_j_20_38_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_38;
            float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
            float *RESTRICT _data_rho_21_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_21;
            float *RESTRICT _data_rho_20 = _data_rho;
            float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
            float *RESTRICT _data_phi_20 = _data_phi;
            float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
            float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
            float *RESTRICT _data_phi_21_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_21;
            float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
            float *RESTRICT _data_phi_21_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_21;
            _data_j_20_38_10[_stride_j_0 * ctr_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_1m1[_stride_rho_0 * ctr_0]) * -2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_1m1[_stride_rho_0 * ctr_0]) * 2.0f + kT * (-1.0f * _data_rho_21_1m1[_stride_rho_0 * ctr_0] + _data_rho_20_10[_stride_rho_0 * ctr_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_1m1[_stride_rho_0 * ctr_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] - 1.0f * _data_phi_21_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_21_10[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_1m1[_stride_rho_0 * ctr_0]) * (-1.0f * _data_phi_21_10[_stride_phi_0 * ctr_0] - 1.0f * _data_phi_21_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_20_1m1[_stride_phi_0 * ctr_0])) * 0.028801180074297286f * ((1.0f) / (kT));
          }
          if (_size_j_1 - 1 > 0 && 0 < _size_j_2 - 1) {
            float *RESTRICT _data_j_20_310 = _data_j + 10 * _stride_j_3;
            float *RESTRICT _data_j_20_310_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_310;
            float *RESTRICT _data_phi_20 = _data_phi;
            float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
            float *RESTRICT _data_rho_20 = _data_rho;
            float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
            float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
            float *RESTRICT _data_rho_21_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_21;
            float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
            float *RESTRICT _data_phi_21_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_21;
            _data_j_20_310_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_0 * z * -1.0f * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_1 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_1 * z * -1.0f * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_2 * z * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * -2.0f * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_21_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_21_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
          }
        }
        if (_size_j_1 - 1 > 0 && 0 < _size_j_2 - 1) {
          float *RESTRICT _data_j_20_310 = _data_j + 10 * _stride_j_3;
          float *RESTRICT _data_j_20_310_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_310;
          float *RESTRICT _data_phi_20 = _data_phi;
          float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
          float *RESTRICT _data_rho_20 = _data_rho;
          float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
          float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2;
          float *RESTRICT _data_rho_21_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_21;
          float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2;
          float *RESTRICT _data_phi_21_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_21;
          _data_j_20_310_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_0 * z * -1.0f * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_1 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_1 * z * -1.0f * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_2 * z * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * -2.0f * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_21_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_21_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
        }
      }
    }
    for (int64_t ctr_2 = 1; ctr_2 < _size_j_2 - 1; ctr_2 += 1) {
      float *RESTRICT _data_j_20_31 = _data_j + _stride_j_2 * ctr_2 + _stride_j_3;
      float *RESTRICT _data_j_20_32 = _data_j + _stride_j_2 * ctr_2 + 2 * _stride_j_3;
      float *RESTRICT _data_j_20_37 = _data_j + _stride_j_2 * ctr_2 + 7 * _stride_j_3;
      float *RESTRICT _data_j_20_38 = _data_j + _stride_j_2 * ctr_2 + 8 * _stride_j_3;
      float *RESTRICT _data_j_20_30 = _data_j + _stride_j_2 * ctr_2;
      float *RESTRICT _data_j_20_33 = _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
      float *RESTRICT _data_j_20_34 = _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
      float *RESTRICT _data_j_20_35 = _data_j + _stride_j_2 * ctr_2 + 5 * _stride_j_3;
      float *RESTRICT _data_j_20_36 = _data_j + _stride_j_2 * ctr_2 + 6 * _stride_j_3;
      float *RESTRICT _data_j_20_39 = _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
      float *RESTRICT _data_j_20_310 = _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
      float *RESTRICT _data_j_20_311 = _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
      float *RESTRICT _data_j_20_312 = _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
      {
        {
          {
            if (ctr_2 > 0 && 0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_34 = _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
              float *RESTRICT _data_j_20_34_10 = _data_j_20_34;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_11 = _stride_rho_1 + _data_rho_20;
              float *RESTRICT _data_rho_20_10 = _data_rho_20;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_11 = _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_20_10 = _data_phi_20;
              _data_j_20_34_10[_stride_j_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_11[0]) * -2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_11[0]) * 2.0f + kT * (-1.0f * _data_rho_20_11[0] + _data_rho_20_10[_stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_11[0]) * (-1.0f * _data_phi_20_10[0] - 1.0f * _data_phi_20_11[0] + _data_phi_20_10[_stride_phi_0] + _data_phi_20_11[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_11[0]) * (-1.0f * _data_phi_20_11[0] - 1.0f * _data_phi_20_11[_stride_phi_0] + _data_phi_20_10[0] + _data_phi_20_10[_stride_phi_0])) * 0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_2 > 0 && 0 < _size_j_1 - 1) {
              float *RESTRICT _data_j_20_311 = _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
              float *RESTRICT _data_j_20_311_10 = _data_j_20_311;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_10 = _data_phi_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_11 = _stride_rho_1 + _data_rho_2m1;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _data_rho_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_11 = _stride_phi_1 + _data_phi_2m1;
              _data_j_20_311_10[_stride_j_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_0 * z * -1.0f * _data_rho_2m1_11[0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0] + f_ext_1 * z * _data_rho_2m1_11[0] + f_ext_2 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_2 * z * -1.0f * _data_rho_2m1_11[0] + kT * -2.0f * _data_rho_2m1_11[0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_2m1_11[0] * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_2m1_11[0] * _data_rho_2m1_11[0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_20_10[_stride_rho_0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_2m1_11[0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
            if (0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_312 = _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
              float *RESTRICT _data_j_20_312_10 = _data_j_20_312;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_10 = _data_phi_20;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _data_rho_20;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
              float *RESTRICT _data_rho_21_11 = _stride_rho_1 + _data_rho_21;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
              float *RESTRICT _data_phi_21_11 = _stride_phi_1 + _data_phi_21;
              _data_j_20_312_10[_stride_j_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_0 * z * -1.0f * _data_rho_21_11[0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0] + f_ext_1 * z * _data_rho_21_11[0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0] + f_ext_2 * z * _data_rho_21_11[0] + kT * -2.0f * _data_rho_21_11[0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_21_11[0] * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_21_11[0] * _data_rho_21_11[0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_20_10[_stride_rho_0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_21_11[0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
          }
          for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
            if (ctr_2 > 0 && 0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_34 = _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
              float *RESTRICT _data_j_20_34_10 = _data_j_20_34;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_11 = _stride_rho_1 + _data_rho_20;
              float *RESTRICT _data_rho_20_10 = _data_rho_20;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_11 = _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_20_10 = _data_phi_20;
              _data_j_20_34_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_11[_stride_rho_0 * ctr_0 - _stride_rho_0]) * -2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_11[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * ctr_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_11[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_20_11[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_20_11[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_11[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_11[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_20_11[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * ctr_0])) * 0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_2 > 0 && 0 < _size_j_1 - 1) {
              float *RESTRICT _data_j_20_311 = _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
              float *RESTRICT _data_j_20_311_10 = _data_j_20_311;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_10 = _data_phi_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_11 = _stride_rho_1 + _data_rho_2m1;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _data_rho_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_11 = _stride_phi_1 + _data_phi_2m1;
              _data_j_20_311_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_0 * z * -1.0f * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_1 * z * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_2 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_2 * z * -1.0f * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * -2.0f * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_2m1_11[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_2m1_11[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
            if (0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_312 = _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
              float *RESTRICT _data_j_20_312_10 = _data_j_20_312;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_10 = _data_phi_20;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _data_rho_20;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
              float *RESTRICT _data_rho_21_11 = _stride_rho_1 + _data_rho_21;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
              float *RESTRICT _data_phi_21_11 = _stride_phi_1 + _data_phi_21;
              _data_j_20_312_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_0 * z * -1.0f * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_1 * z * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_2 * z * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * -2.0f * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_21_11[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_21_11[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
          }
          {
            if (ctr_2 > 0 && 0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_34 = _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
              float *RESTRICT _data_j_20_34_10 = _data_j_20_34;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_11 = _stride_rho_1 + _data_rho_20;
              float *RESTRICT _data_rho_20_10 = _data_rho_20;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_11 = _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_20_10 = _data_phi_20;
              _data_j_20_34_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_20_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * -2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_20_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_20_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_20_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_20_11[_stride_phi_0 * (_size_j_0 - 1)]) + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_20_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_20_11[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)])) * 0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_2 > 0 && 0 < _size_j_1 - 1) {
              float *RESTRICT _data_j_20_311 = _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
              float *RESTRICT _data_j_20_311_10 = _data_j_20_311;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_10 = _data_phi_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_11 = _stride_rho_1 + _data_rho_2m1;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _data_rho_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_11 = _stride_phi_1 + _data_phi_2m1;
              _data_j_20_311_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_0 * z * -1.0f * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_1 * z * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_2 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_2 * z * -1.0f * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * -2.0f * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_2m1_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_2m1_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
            if (0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_312 = _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
              float *RESTRICT _data_j_20_312_10 = _data_j_20_312;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_10 = _data_phi_20;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _data_rho_20;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
              float *RESTRICT _data_rho_21_11 = _stride_rho_1 + _data_rho_21;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
              float *RESTRICT _data_phi_21_11 = _stride_phi_1 + _data_phi_21;
              _data_j_20_312_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_0 * z * -1.0f * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_1 * z * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_2 * z * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * -2.0f * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_21_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_21_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
          }
        }
        for (int64_t ctr_1 = 1; ctr_1 < _size_j_1 - 1; ctr_1 += 1) {
          float *RESTRICT _data_j_20_31_10 = _stride_j_1 * ctr_1 + _data_j_20_31;
          float *RESTRICT _data_j_20_32_10 = _stride_j_1 * ctr_1 + _data_j_20_32;
          float *RESTRICT _data_j_20_37_10 = _stride_j_1 * ctr_1 + _data_j_20_37;
          float *RESTRICT _data_j_20_38_10 = _stride_j_1 * ctr_1 + _data_j_20_38;
          float *RESTRICT _data_j_20_30_10 = _stride_j_1 * ctr_1 + _data_j_20_30;
          float *RESTRICT _data_j_20_33_10 = _stride_j_1 * ctr_1 + _data_j_20_33;
          float *RESTRICT _data_j_20_34_10 = _stride_j_1 * ctr_1 + _data_j_20_34;
          float *RESTRICT _data_j_20_35_10 = _stride_j_1 * ctr_1 + _data_j_20_35;
          float *RESTRICT _data_j_20_36_10 = _stride_j_1 * ctr_1 + _data_j_20_36;
          float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
          float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
          float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
          float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
          {
            float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
            float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
            float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
            float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
            _data_j_20_30_10[_stride_j_0] = D * (f_ext_0 * z * (_data_rho_20_10[0] + _data_rho_20_10[_stride_rho_0]) * -1.0f + kT * (-1.0f * _data_rho_20_10[0] + _data_rho_20_10[_stride_rho_0]) * 2.0f + z * (-1.0f * _data_phi_20_10[0] + _data_phi_20_10[_stride_phi_0]) * (_data_rho_20_10[0] + _data_rho_20_10[_stride_rho_0])) * 0.081462038946841925f * ((1.0f) / (kT));
            float *RESTRICT _data_rho_20_1m1 = _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20;
            float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20;
            _data_j_20_33_10[_stride_j_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_1m1[0]) * -2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_1m1[0]) * -2.0f + kT * (-1.0f * _data_rho_20_1m1[0] + _data_rho_20_10[_stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_1m1[0]) * (-1.0f * _data_phi_20_10[0] - 1.0f * _data_phi_20_1m1[0] + _data_phi_20_10[_stride_phi_0] + _data_phi_20_1m1[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_1m1[0]) * (-1.0f * _data_phi_20_1m1[0] - 1.0f * _data_phi_20_1m1[_stride_phi_0] + _data_phi_20_10[0] + _data_phi_20_10[_stride_phi_0])) * 0.028801180074297286f * ((1.0f) / (kT));
            float *RESTRICT _data_rho_20_11 = _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_20;
            float *RESTRICT _data_phi_20_11 = _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_20;
            _data_j_20_34_10[_stride_j_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_11[0]) * -2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_11[0]) * 2.0f + kT * (-1.0f * _data_rho_20_11[0] + _data_rho_20_10[_stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_11[0]) * (-1.0f * _data_phi_20_10[0] - 1.0f * _data_phi_20_11[0] + _data_phi_20_10[_stride_phi_0] + _data_phi_20_11[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_11[0]) * (-1.0f * _data_phi_20_11[0] - 1.0f * _data_phi_20_11[_stride_phi_0] + _data_phi_20_10[0] + _data_phi_20_10[_stride_phi_0])) * 0.028801180074297286f * ((1.0f) / (kT));
            float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
            float *RESTRICT _data_rho_2m1_10 = _stride_rho_1 * ctr_1 + _data_rho_2m1;
            float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
            float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * ctr_1 + _data_phi_2m1;
            _data_j_20_35_10[_stride_j_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_10[0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_10[0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0] + _data_rho_2m1_10[0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_10[0]) * (-1.0f * _data_phi_20_10[0] - 1.0f * _data_phi_20_10[_stride_phi_0] + _data_phi_2m1_10[0] + _data_phi_2m1_10[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_10[0]) * (-1.0f * _data_phi_20_10[0] - 1.0f * _data_phi_2m1_10[0] + _data_phi_20_10[_stride_phi_0] + _data_phi_2m1_10[_stride_phi_0]) * -1.0f) * -0.028801180074297286f * ((1.0f) / (kT));
            float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
            float *RESTRICT _data_rho_21_10 = _stride_rho_1 * ctr_1 + _data_rho_21;
            float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
            float *RESTRICT _data_phi_21_10 = _stride_phi_1 * ctr_1 + _data_phi_21;
            _data_j_20_36_10[_stride_j_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_10[0]) * -2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_10[0]) * 2.0f + kT * (-1.0f * _data_rho_21_10[0] + _data_rho_20_10[_stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_10[0]) * (-1.0f * _data_phi_20_10[0] - 1.0f * _data_phi_21_10[0] + _data_phi_20_10[_stride_phi_0] + _data_phi_21_10[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_10[0]) * (-1.0f * _data_phi_21_10[0] - 1.0f * _data_phi_21_10[_stride_phi_0] + _data_phi_20_10[0] + _data_phi_20_10[_stride_phi_0])) * 0.028801180074297286f * ((1.0f) / (kT));
            float *RESTRICT _data_rho_2m1_1m1 = _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_2m1;
            float *RESTRICT _data_phi_2m1_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_2m1;
            _data_j_20_39_10[_stride_j_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * 2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * (-1.0f * _data_phi_20_10[0] - 1.0f * _data_phi_20_1m1[_stride_phi_0] + _data_phi_2m1_10[0] + _data_phi_2m1_1m1[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * (-1.0f * _data_phi_20_1m1[0] - 1.0f * _data_phi_2m1_10[0] + _data_phi_20_1m1[_stride_phi_0] + _data_phi_2m1_10[_stride_phi_0]) * -1.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * (-1.0f * _data_phi_20_1m1[0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0] + _data_phi_20_10[0] + _data_phi_2m1_10[_stride_phi_0]) * -1.0f) * -0.02351606505734748f * ((1.0f) / (kT));
            float *RESTRICT _data_rho_21_1m1 = _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_21;
            float *RESTRICT _data_phi_21_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_21;
            _data_j_20_310_10[_stride_j_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_0 * z * -1.0f * _data_rho_21_1m1[0] + f_ext_1 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_1 * z * -1.0f * _data_rho_21_1m1[0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0] + f_ext_2 * z * _data_rho_21_1m1[0] + kT * -2.0f * _data_rho_21_1m1[0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_21_1m1[0] * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_21_1m1[0] * _data_rho_21_1m1[0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_20_10[_stride_rho_0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_21_1m1[0]) * 0.04703213011469496f * ((1.0f) / (kT));
            float *RESTRICT _data_rho_2m1_11 = _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_2m1;
            float *RESTRICT _data_phi_2m1_11 = _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_2m1;
            _data_j_20_311_10[_stride_j_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_0 * z * -1.0f * _data_rho_2m1_11[0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0] + f_ext_1 * z * _data_rho_2m1_11[0] + f_ext_2 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_2 * z * -1.0f * _data_rho_2m1_11[0] + kT * -2.0f * _data_rho_2m1_11[0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_2m1_11[0] * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_2m1_11[0] * _data_rho_2m1_11[0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_20_10[_stride_rho_0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_2m1_11[0]) * 0.04703213011469496f * ((1.0f) / (kT));
            float *RESTRICT _data_rho_21_11 = _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_21;
            float *RESTRICT _data_phi_21_11 = _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_21;
            _data_j_20_312_10[_stride_j_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_0 * z * -1.0f * _data_rho_21_11[0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0] + f_ext_1 * z * _data_rho_21_11[0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0] + f_ext_2 * z * _data_rho_21_11[0] + kT * -2.0f * _data_rho_21_11[0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_21_11[0] * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_21_11[0] * _data_rho_21_11[0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_20_10[_stride_rho_0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_21_11[0]) * 0.04703213011469496f * ((1.0f) / (kT));
            {
              if (ctr_1 > 0 && ctr_2 > 0 && 1 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {
                float *RESTRICT _data_j_20_31 = _data_j + _stride_j_2 * ctr_2 + _stride_j_3;
                float *RESTRICT _data_j_20_31_10 = _stride_j_1 * ctr_1 + _data_j_20_31;
                float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
                float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20;
                float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
                float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
                float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
                float *RESTRICT _data_rho_20_1m1 = _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_20;
                _data_j_20_31_10[_stride_j_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_1m1[_stride_rho_0]) * -1.0f + kT * (-1.0f * _data_rho_20_1m1[_stride_rho_0] + _data_rho_20_10[_stride_rho_0]) * 2.0f + z * (-1.0f * _data_phi_20_1m1[_stride_phi_0] + _data_phi_20_10[_stride_phi_0]) * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_1m1[_stride_rho_0])) * 0.081462038946841925f * ((1.0f) / (kT));
              }
              if (ctr_1 > 0 && ctr_2 > 0 && 1 < _size_j_0 - 1 && ctr_1 < _size_j_1 - 1) {
                float *RESTRICT _data_j_20_32 = _data_j + _stride_j_2 * ctr_2 + 2 * _stride_j_3;
                float *RESTRICT _data_j_20_32_10 = _stride_j_1 * ctr_1 + _data_j_20_32;
                float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
                float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
                float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * ctr_1 + _data_phi_2m1;
                float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                float *RESTRICT _data_rho_2m1_10 = _stride_rho_1 * ctr_1 + _data_rho_2m1;
                float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
                float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
                _data_j_20_32_10[_stride_j_0] = D * (f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_10[_stride_rho_0]) + kT * (-1.0f * _data_rho_20_10[_stride_rho_0] + _data_rho_2m1_10[_stride_rho_0]) * 2.0f + z * (-1.0f * _data_phi_20_10[_stride_phi_0] + _data_phi_2m1_10[_stride_phi_0]) * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_10[_stride_rho_0])) * -0.081462038946841925f * ((1.0f) / (kT));
              }
              if (ctr_1 > 0 && ctr_2 > 0 && 1 < _size_j_0 - 1) {
                float *RESTRICT _data_j_20_37 = _data_j + _stride_j_2 * ctr_2 + 7 * _stride_j_3;
                float *RESTRICT _data_j_20_37_10 = _stride_j_1 * ctr_1 + _data_j_20_37;
                float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
                float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
                float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
                float *RESTRICT _data_rho_2m1_1m1 = _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_2m1;
                float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
                float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
                float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20;
                float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
                float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * ctr_1 + _data_phi_2m1;
                float *RESTRICT _data_phi_2m1_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_2m1;
                _data_j_20_37_10[_stride_j_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0] + _data_phi_2m1_10[_stride_phi_0] + _data_phi_2m1_1m1[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0] + _data_phi_20_10[_stride_phi_0] + _data_phi_2m1_10[_stride_phi_0]) * -1.0f) * -0.028801180074297286f * ((1.0f) / (kT));
              }
              if (ctr_1 > 0 && 1 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {
                float *RESTRICT _data_j_20_38 = _data_j + _stride_j_2 * ctr_2 + 8 * _stride_j_3;
                float *RESTRICT _data_j_20_38_10 = _stride_j_1 * ctr_1 + _data_j_20_38;
                float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
                float *RESTRICT _data_rho_21_1m1 = _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_21;
                float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
                float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
                float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
                float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20;
                float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
                float *RESTRICT _data_phi_21_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_21;
                float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
                float *RESTRICT _data_phi_21_10 = _stride_phi_1 * ctr_1 + _data_phi_21;
                _data_j_20_38_10[_stride_j_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_1m1[_stride_rho_0]) * -2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_1m1[_stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_21_1m1[_stride_rho_0] + _data_rho_20_10[_stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_1m1[_stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0] - 1.0f * _data_phi_21_1m1[_stride_phi_0] + _data_phi_20_10[_stride_phi_0] + _data_phi_21_10[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_1m1[_stride_rho_0]) * (-1.0f * _data_phi_21_10[_stride_phi_0] - 1.0f * _data_phi_21_1m1[_stride_phi_0] + _data_phi_20_10[_stride_phi_0] + _data_phi_20_1m1[_stride_phi_0])) * 0.028801180074297286f * ((1.0f) / (kT));
              }
            }
            for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
              _data_j_20_30_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0 - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * ctr_0]) * -1.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0 - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * ctr_0]) * 2.0f + z * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * ctr_0]) * (_data_rho_20_10[_stride_rho_0 * ctr_0 - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * ctr_0])) * 0.081462038946841925f * ((1.0f) / (kT));
              _data_j_20_31_10[_stride_j_0 * ctr_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_1m1[_stride_rho_0 * ctr_0]) * -1.0f + kT * (-1.0f * _data_rho_20_1m1[_stride_rho_0 * ctr_0] + _data_rho_20_10[_stride_rho_0 * ctr_0]) * 2.0f + z * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0]) * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_1m1[_stride_rho_0 * ctr_0])) * 0.081462038946841925f * ((1.0f) / (kT));
              _data_j_20_32_10[_stride_j_0 * ctr_0] = D * (f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_10[_stride_rho_0 * ctr_0]) + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_10[_stride_rho_0 * ctr_0]) * 2.0f + z * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_10[_stride_rho_0 * ctr_0])) * -0.081462038946841925f * ((1.0f) / (kT));
              _data_j_20_33_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * -2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * -2.0f + kT * (-1.0f * _data_rho_20_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * ctr_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_20_1m1[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * ctr_0])) * 0.028801180074297286f * ((1.0f) / (kT));
              _data_j_20_34_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_11[_stride_rho_0 * ctr_0 - _stride_rho_0]) * -2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_11[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * ctr_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_11[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_20_11[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_20_11[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_11[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_11[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_20_11[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * ctr_0])) * 0.028801180074297286f * ((1.0f) / (kT));
              _data_j_20_35_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_2m1_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) * -1.0f) * -0.028801180074297286f * ((1.0f) / (kT));
              _data_j_20_36_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * -2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_21_10[_stride_rho_0 * ctr_0 - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * ctr_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_21_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_21_10[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_21_10[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_21_10[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * ctr_0])) * 0.028801180074297286f * ((1.0f) / (kT));
              _data_j_20_37_10[_stride_j_0 * ctr_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0] + _data_phi_2m1_1m1[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) * -1.0f) * -0.028801180074297286f * ((1.0f) / (kT));
              _data_j_20_38_10[_stride_j_0 * ctr_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_1m1[_stride_rho_0 * ctr_0]) * -2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_1m1[_stride_rho_0 * ctr_0]) * 2.0f + kT * (-1.0f * _data_rho_21_1m1[_stride_rho_0 * ctr_0] + _data_rho_20_10[_stride_rho_0 * ctr_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_1m1[_stride_rho_0 * ctr_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] - 1.0f * _data_phi_21_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_21_10[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_1m1[_stride_rho_0 * ctr_0]) * (-1.0f * _data_phi_21_10[_stride_phi_0 * ctr_0] - 1.0f * _data_phi_21_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_20_1m1[_stride_phi_0 * ctr_0])) * 0.028801180074297286f * ((1.0f) / (kT));
              _data_j_20_39_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_2m1_1m1[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_2m1_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_1m1[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) * -1.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) * -1.0f) * -0.02351606505734748f * ((1.0f) / (kT));
              _data_j_20_310_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_0 * z * -1.0f * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_1 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_1 * z * -1.0f * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_2 * z * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * -2.0f * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_21_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_21_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
              _data_j_20_311_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_0 * z * -1.0f * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_1 * z * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_2 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_2 * z * -1.0f * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * -2.0f * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_2m1_11[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_2m1_11[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
              _data_j_20_312_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_0 * z * -1.0f * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_1 * z * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_2 * z * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * -2.0f * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_21_11[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_21_11[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_21_11[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
            _data_j_20_30_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)]) * -1.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)]) * 2.0f + z * (-1.0f * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)]) * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)])) * 0.081462038946841925f * ((1.0f) / (kT));
            _data_j_20_33_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_20_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * -2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_20_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * -2.0f + kT * (-1.0f * _data_rho_20_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_20_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1)]) + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_20_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)])) * 0.028801180074297286f * ((1.0f) / (kT));
            _data_j_20_34_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_20_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * -2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_20_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_20_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_20_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_20_11[_stride_phi_0 * (_size_j_0 - 1)]) + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_20_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_20_11[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)])) * 0.028801180074297286f * ((1.0f) / (kT));
            _data_j_20_35_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1)]) + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1)]) * -1.0f) * -0.028801180074297286f * ((1.0f) / (kT));
            _data_j_20_36_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_21_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * -2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_21_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_21_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_21_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_21_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_21_10[_stride_phi_0 * (_size_j_0 - 1)]) + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_21_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_21_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_21_10[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)])) * 0.028801180074297286f * ((1.0f) / (kT));
            _data_j_20_39_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_2m1_1m1[_stride_phi_0 * (_size_j_0 - 1)]) + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1)]) * -1.0f + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1)]) * -1.0f) * -0.02351606505734748f * ((1.0f) / (kT));
            _data_j_20_310_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_0 * z * -1.0f * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_1 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_1 * z * -1.0f * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_2 * z * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * -2.0f * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_21_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_21_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
            _data_j_20_311_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_0 * z * -1.0f * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_1 * z * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_2 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_2 * z * -1.0f * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * -2.0f * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_2m1_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_2m1_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
            _data_j_20_312_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_0 * z * -1.0f * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_1 * z * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_2 * z * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * -2.0f * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_21_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_21_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_21_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
            {
            }
          }
        }
        {
          {
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && 1 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_31 = _data_j + _stride_j_2 * ctr_2 + _stride_j_3;
              float *RESTRICT _data_j_20_31_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_31;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
              float *RESTRICT _data_rho_20_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_20;
              _data_j_20_31_10[_stride_j_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_1m1[_stride_rho_0]) * -1.0f + kT * (-1.0f * _data_rho_20_1m1[_stride_rho_0] + _data_rho_20_10[_stride_rho_0]) * 2.0f + z * (-1.0f * _data_phi_20_1m1[_stride_phi_0] + _data_phi_20_10[_stride_phi_0]) * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_1m1[_stride_rho_0])) * 0.081462038946841925f * ((1.0f) / (kT));
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_33 = _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
              float *RESTRICT _data_j_20_33_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_33;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_20;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
              _data_j_20_33_10[_stride_j_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_1m1[0]) * -2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_1m1[0]) * -2.0f + kT * (-1.0f * _data_rho_20_1m1[0] + _data_rho_20_10[_stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_1m1[0]) * (-1.0f * _data_phi_20_10[0] - 1.0f * _data_phi_20_1m1[0] + _data_phi_20_10[_stride_phi_0] + _data_phi_20_1m1[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_20_1m1[0]) * (-1.0f * _data_phi_20_1m1[0] - 1.0f * _data_phi_20_1m1[_stride_phi_0] + _data_phi_20_10[0] + _data_phi_20_10[_stride_phi_0])) * 0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && 1 < _size_j_0 - 1) {
              float *RESTRICT _data_j_20_37 = _data_j + _stride_j_2 * ctr_2 + 7 * _stride_j_3;
              float *RESTRICT _data_j_20_37_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_37;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_2m1;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1;
              float *RESTRICT _data_phi_2m1_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_2m1;
              _data_j_20_37_10[_stride_j_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0] + _data_phi_2m1_10[_stride_phi_0] + _data_phi_2m1_1m1[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0] + _data_phi_20_10[_stride_phi_0] + _data_phi_2m1_10[_stride_phi_0]) * -1.0f) * -0.028801180074297286f * ((1.0f) / (kT));
            }
            if (_size_j_1 - 1 > 0 && 1 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_38 = _data_j + _stride_j_2 * ctr_2 + 8 * _stride_j_3;
              float *RESTRICT _data_j_20_38_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_38;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
              float *RESTRICT _data_rho_21_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_21;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
              float *RESTRICT _data_phi_21_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_21;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
              float *RESTRICT _data_phi_21_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_21;
              _data_j_20_38_10[_stride_j_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_1m1[_stride_rho_0]) * -2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_1m1[_stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_21_1m1[_stride_rho_0] + _data_rho_20_10[_stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_1m1[_stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0] - 1.0f * _data_phi_21_1m1[_stride_phi_0] + _data_phi_20_10[_stride_phi_0] + _data_phi_21_10[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_21_1m1[_stride_rho_0]) * (-1.0f * _data_phi_21_10[_stride_phi_0] - 1.0f * _data_phi_21_1m1[_stride_phi_0] + _data_phi_20_10[_stride_phi_0] + _data_phi_20_1m1[_stride_phi_0])) * 0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0) {
              float *RESTRICT _data_j_20_39 = _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
              float *RESTRICT _data_j_20_39_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_39;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_2m1;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_2m1;
              float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1;
              _data_j_20_39_10[_stride_j_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * 2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * (-1.0f * _data_phi_20_10[0] - 1.0f * _data_phi_20_1m1[_stride_phi_0] + _data_phi_2m1_10[0] + _data_phi_2m1_1m1[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * (-1.0f * _data_phi_20_1m1[0] - 1.0f * _data_phi_2m1_10[0] + _data_phi_20_1m1[_stride_phi_0] + _data_phi_2m1_10[_stride_phi_0]) * -1.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * (-1.0f * _data_phi_20_1m1[0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0] + _data_phi_20_10[0] + _data_phi_2m1_10[_stride_phi_0]) * -1.0f) * -0.02351606505734748f * ((1.0f) / (kT));
            }
            if (_size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_310 = _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
              float *RESTRICT _data_j_20_310_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_310;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
              float *RESTRICT _data_rho_21_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_21;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
              float *RESTRICT _data_phi_21_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_21;
              _data_j_20_310_10[_stride_j_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_0 * z * -1.0f * _data_rho_21_1m1[0] + f_ext_1 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_1 * z * -1.0f * _data_rho_21_1m1[0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0] + f_ext_2 * z * _data_rho_21_1m1[0] + kT * -2.0f * _data_rho_21_1m1[0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_21_1m1[0] * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_21_1m1[0] * _data_rho_21_1m1[0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_20_10[_stride_rho_0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_21_1m1[0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
          }
          for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && ctr_0 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_31 = _data_j + _stride_j_2 * ctr_2 + _stride_j_3;
              float *RESTRICT _data_j_20_31_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_31;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
              float *RESTRICT _data_rho_20_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_20;
              _data_j_20_31_10[_stride_j_0 * ctr_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_1m1[_stride_rho_0 * ctr_0]) * -1.0f + kT * (-1.0f * _data_rho_20_1m1[_stride_rho_0 * ctr_0] + _data_rho_20_10[_stride_rho_0 * ctr_0]) * 2.0f + z * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0]) * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_1m1[_stride_rho_0 * ctr_0])) * 0.081462038946841925f * ((1.0f) / (kT));
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_33 = _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
              float *RESTRICT _data_j_20_33_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_33;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_20;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
              _data_j_20_33_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * -2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * -2.0f + kT * (-1.0f * _data_rho_20_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * ctr_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_20_1m1[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_20_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * ctr_0])) * 0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && ctr_0 < _size_j_0 - 1) {
              float *RESTRICT _data_j_20_37 = _data_j + _stride_j_2 * ctr_2 + 7 * _stride_j_3;
              float *RESTRICT _data_j_20_37_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_37;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_2m1;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1;
              float *RESTRICT _data_phi_2m1_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_2m1;
              _data_j_20_37_10[_stride_j_0 * ctr_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0] + _data_phi_2m1_1m1[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) * -1.0f) * -0.028801180074297286f * ((1.0f) / (kT));
            }
            if (_size_j_1 - 1 > 0 && ctr_0 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_38 = _data_j + _stride_j_2 * ctr_2 + 8 * _stride_j_3;
              float *RESTRICT _data_j_20_38_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_38;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
              float *RESTRICT _data_rho_21_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_21;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
              float *RESTRICT _data_phi_21_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_21;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
              float *RESTRICT _data_phi_21_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_21;
              _data_j_20_38_10[_stride_j_0 * ctr_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_1m1[_stride_rho_0 * ctr_0]) * -2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_1m1[_stride_rho_0 * ctr_0]) * 2.0f + kT * (-1.0f * _data_rho_21_1m1[_stride_rho_0 * ctr_0] + _data_rho_20_10[_stride_rho_0 * ctr_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_1m1[_stride_rho_0 * ctr_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] - 1.0f * _data_phi_21_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_21_10[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_21_1m1[_stride_rho_0 * ctr_0]) * (-1.0f * _data_phi_21_10[_stride_phi_0 * ctr_0] - 1.0f * _data_phi_21_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_20_1m1[_stride_phi_0 * ctr_0])) * 0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0) {
              float *RESTRICT _data_j_20_39 = _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
              float *RESTRICT _data_j_20_39_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_39;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_2m1;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_2m1;
              float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1;
              _data_j_20_39_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_2m1_1m1[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_2m1_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_1m1[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) * -1.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) * -1.0f) * -0.02351606505734748f * ((1.0f) / (kT));
            }
            if (_size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_310 = _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
              float *RESTRICT _data_j_20_310_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_310;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
              float *RESTRICT _data_rho_21_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_21;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
              float *RESTRICT _data_phi_21_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_21;
              _data_j_20_310_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_0 * z * -1.0f * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_1 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_1 * z * -1.0f * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_2 * z * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * -2.0f * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_21_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_21_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_21_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
          }
          {
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_33 = _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
              float *RESTRICT _data_j_20_33_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_33;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_20;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
              _data_j_20_33_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_20_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * -2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_20_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * -2.0f + kT * (-1.0f * _data_rho_20_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_20_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1)]) + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_20_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)])) * 0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0) {
              float *RESTRICT _data_j_20_39 = _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
              float *RESTRICT _data_j_20_39_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_39;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * ctr_2 - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_2m1;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * ctr_2 - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_2m1;
              float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1;
              _data_j_20_39_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_2m1_1m1[_stride_phi_0 * (_size_j_0 - 1)]) + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1)]) * -1.0f + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1)]) * -1.0f) * -0.02351606505734748f * ((1.0f) / (kT));
            }
            if (_size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {
              float *RESTRICT _data_j_20_310 = _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
              float *RESTRICT _data_j_20_310_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_310;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * ctr_2;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * ctr_2;
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
              float *RESTRICT _data_rho_21 = _data_rho + _stride_rho_2 * ctr_2 + _stride_rho_2;
              float *RESTRICT _data_rho_21_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_21;
              float *RESTRICT _data_phi_21 = _data_phi + _stride_phi_2 * ctr_2 + _stride_phi_2;
              float *RESTRICT _data_phi_21_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_21;
              _data_j_20_310_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_0 * z * -1.0f * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_1 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_1 * z * -1.0f * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_2 * z * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_2 * z * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * -2.0f * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_21_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_21_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_21_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
          }
        }
      }
    }
    {
      {
        if (_size_j_2 - 1 > 0 && 0 < _size_j_1 - 1) {
          float *RESTRICT _data_j_20_311 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3;
          float *RESTRICT _data_j_20_311_10 = _data_j_20_311;
          float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
          float *RESTRICT _data_phi_20_10 = _data_phi_20;
          float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
          float *RESTRICT _data_rho_2m1_11 = _stride_rho_1 + _data_rho_2m1;
          float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
          float *RESTRICT _data_rho_20_10 = _data_rho_20;
          float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
          float *RESTRICT _data_phi_2m1_11 = _stride_phi_1 + _data_phi_2m1;
          _data_j_20_311_10[_stride_j_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_0 * z * -1.0f * _data_rho_2m1_11[0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0] + f_ext_1 * z * _data_rho_2m1_11[0] + f_ext_2 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_2 * z * -1.0f * _data_rho_2m1_11[0] + kT * -2.0f * _data_rho_2m1_11[0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_2m1_11[0] * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_2m1_11[0] * _data_rho_2m1_11[0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_20_10[_stride_rho_0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_2m1_11[0]) * 0.04703213011469496f * ((1.0f) / (kT));
        }
        for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
          if (_size_j_2 - 1 > 0 && 0 < _size_j_1 - 1) {
            float *RESTRICT _data_j_20_311 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3;
            float *RESTRICT _data_j_20_311_10 = _data_j_20_311;
            float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
            float *RESTRICT _data_phi_20_10 = _data_phi_20;
            float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
            float *RESTRICT _data_rho_2m1_11 = _stride_rho_1 + _data_rho_2m1;
            float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
            float *RESTRICT _data_rho_20_10 = _data_rho_20;
            float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
            float *RESTRICT _data_phi_2m1_11 = _stride_phi_1 + _data_phi_2m1;
            _data_j_20_311_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_0 * z * -1.0f * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_1 * z * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_2 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_2 * z * -1.0f * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * -2.0f * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_2m1_11[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_2m1_11[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
          }
        }
        if (_size_j_2 - 1 > 0 && 0 < _size_j_1 - 1) {
          float *RESTRICT _data_j_20_311 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3;
          float *RESTRICT _data_j_20_311_10 = _data_j_20_311;
          float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
          float *RESTRICT _data_phi_20_10 = _data_phi_20;
          float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
          float *RESTRICT _data_rho_2m1_11 = _stride_rho_1 + _data_rho_2m1;
          float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
          float *RESTRICT _data_rho_20_10 = _data_rho_20;
          float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
          float *RESTRICT _data_phi_2m1_11 = _stride_phi_1 + _data_phi_2m1;
          _data_j_20_311_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_0 * z * -1.0f * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_1 * z * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_2 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_2 * z * -1.0f * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * -2.0f * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_2m1_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_2m1_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
        }
      }
      for (int64_t ctr_1 = 1; ctr_1 < _size_j_1 - 1; ctr_1 += 1) {
        {
          {
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && 1 < _size_j_0 - 1 && ctr_1 < _size_j_1 - 1) {
              float *RESTRICT _data_j_20_32 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 2 * _stride_j_3;
              float *RESTRICT _data_j_20_32_10 = _stride_j_1 * ctr_1 + _data_j_20_32;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * ctr_1 + _data_phi_2m1;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_10 = _stride_rho_1 * ctr_1 + _data_rho_2m1;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              _data_j_20_32_10[_stride_j_0] = D * (f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_10[_stride_rho_0]) + kT * (-1.0f * _data_rho_20_10[_stride_rho_0] + _data_rho_2m1_10[_stride_rho_0]) * 2.0f + z * (-1.0f * _data_phi_20_10[_stride_phi_0] + _data_phi_2m1_10[_stride_phi_0]) * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_10[_stride_rho_0])) * -0.081462038946841925f * ((1.0f) / (kT));
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {
              float *RESTRICT _data_j_20_35 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 5 * _stride_j_3;
              float *RESTRICT _data_j_20_35_10 = _stride_j_1 * ctr_1 + _data_j_20_35;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_10 = _stride_rho_1 * ctr_1 + _data_rho_2m1;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * ctr_1 + _data_phi_2m1;
              _data_j_20_35_10[_stride_j_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_10[0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_10[0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0] + _data_rho_2m1_10[0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_10[0]) * (-1.0f * _data_phi_20_10[0] - 1.0f * _data_phi_20_10[_stride_phi_0] + _data_phi_2m1_10[0] + _data_phi_2m1_10[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_10[0]) * (-1.0f * _data_phi_20_10[0] - 1.0f * _data_phi_2m1_10[0] + _data_phi_20_10[_stride_phi_0] + _data_phi_2m1_10[_stride_phi_0]) * -1.0f) * -0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && 1 < _size_j_0 - 1) {
              float *RESTRICT _data_j_20_37 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 7 * _stride_j_3;
              float *RESTRICT _data_j_20_37_10 = _stride_j_1 * ctr_1 + _data_j_20_37;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_1m1 = _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_2m1;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * ctr_1 + _data_phi_2m1;
              float *RESTRICT _data_phi_2m1_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_2m1;
              _data_j_20_37_10[_stride_j_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0] + _data_phi_2m1_10[_stride_phi_0] + _data_phi_2m1_1m1[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0] + _data_phi_20_10[_stride_phi_0] + _data_phi_2m1_10[_stride_phi_0]) * -1.0f) * -0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0) {
              float *RESTRICT _data_j_20_39 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3;
              float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_1m1 = _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_2m1;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_2m1;
              float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * ctr_1 + _data_phi_2m1;
              _data_j_20_39_10[_stride_j_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * 2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * (-1.0f * _data_phi_20_10[0] - 1.0f * _data_phi_20_1m1[_stride_phi_0] + _data_phi_2m1_10[0] + _data_phi_2m1_1m1[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * (-1.0f * _data_phi_20_1m1[0] - 1.0f * _data_phi_2m1_10[0] + _data_phi_20_1m1[_stride_phi_0] + _data_phi_2m1_10[_stride_phi_0]) * -1.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * (-1.0f * _data_phi_20_1m1[0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0] + _data_phi_20_10[0] + _data_phi_2m1_10[_stride_phi_0]) * -1.0f) * -0.02351606505734748f * ((1.0f) / (kT));
            }
            if (_size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {
              float *RESTRICT _data_j_20_311 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3;
              float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_11 = _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_2m1;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_11 = _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_2m1;
              _data_j_20_311_10[_stride_j_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_0 * z * -1.0f * _data_rho_2m1_11[0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0] + f_ext_1 * z * _data_rho_2m1_11[0] + f_ext_2 * z * -1.0f * _data_rho_20_10[_stride_rho_0] + f_ext_2 * z * -1.0f * _data_rho_2m1_11[0] + kT * -2.0f * _data_rho_2m1_11[0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_2m1_11[0] * _data_rho_20_10[_stride_rho_0] + z * -1.0f * _data_phi_2m1_11[0] * _data_rho_2m1_11[0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_20_10[_stride_rho_0] + z * _data_phi_20_10[_stride_phi_0] * _data_rho_2m1_11[0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
          }
          for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && ctr_0 < _size_j_0 - 1 && ctr_1 < _size_j_1 - 1) {
              float *RESTRICT _data_j_20_32 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 2 * _stride_j_3;
              float *RESTRICT _data_j_20_32_10 = _stride_j_1 * ctr_1 + _data_j_20_32;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * ctr_1 + _data_phi_2m1;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_10 = _stride_rho_1 * ctr_1 + _data_rho_2m1;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              _data_j_20_32_10[_stride_j_0 * ctr_0] = D * (f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_10[_stride_rho_0 * ctr_0]) + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_10[_stride_rho_0 * ctr_0]) * 2.0f + z * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_10[_stride_rho_0 * ctr_0])) * -0.081462038946841925f * ((1.0f) / (kT));
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {
              float *RESTRICT _data_j_20_35 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 5 * _stride_j_3;
              float *RESTRICT _data_j_20_35_10 = _stride_j_1 * ctr_1 + _data_j_20_35;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_10 = _stride_rho_1 * ctr_1 + _data_rho_2m1;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * ctr_1 + _data_phi_2m1;
              _data_j_20_35_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_10[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_2m1_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) * -1.0f) * -0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && ctr_0 < _size_j_0 - 1) {
              float *RESTRICT _data_j_20_37 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 7 * _stride_j_3;
              float *RESTRICT _data_j_20_37_10 = _stride_j_1 * ctr_1 + _data_j_20_37;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_1m1 = _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_2m1;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * ctr_1 + _data_phi_2m1;
              float *RESTRICT _data_phi_2m1_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_2m1;
              _data_j_20_37_10[_stride_j_0 * ctr_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0] + _data_phi_2m1_1m1[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) * -1.0f) * -0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0) {
              float *RESTRICT _data_j_20_39 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3;
              float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_1m1 = _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_2m1;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_2m1;
              float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * ctr_1 + _data_phi_2m1;
              _data_j_20_39_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_2m1_1m1[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_2m1_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_1m1[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) * -1.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) * -1.0f) * -0.02351606505734748f * ((1.0f) / (kT));
            }
            if (_size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {
              float *RESTRICT _data_j_20_311 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3;
              float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_11 = _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_2m1;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_11 = _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_2m1;
              _data_j_20_311_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_0 * z * -1.0f * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_1 * z * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + f_ext_2 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + f_ext_2 * z * -1.0f * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * -2.0f * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_2m1_11[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * -1.0f * _data_phi_2m1_11[_stride_phi_0 * ctr_0 - _stride_phi_0] * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_20_10[_stride_rho_0 * ctr_0] + z * _data_phi_20_10[_stride_phi_0 * ctr_0] * _data_rho_2m1_11[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
          }
          {
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {
              float *RESTRICT _data_j_20_35 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 5 * _stride_j_3;
              float *RESTRICT _data_j_20_35_10 = _stride_j_1 * ctr_1 + _data_j_20_35;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_10 = _stride_rho_1 * ctr_1 + _data_rho_2m1;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * ctr_1 + _data_phi_2m1;
              _data_j_20_35_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1)]) + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_10[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1)]) * -1.0f) * -0.028801180074297286f * ((1.0f) / (kT));
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0) {
              float *RESTRICT _data_j_20_39 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3;
              float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_1m1 = _stride_rho_1 * ctr_1 - _stride_rho_1 + _data_rho_2m1;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_20;
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_1m1 = _stride_phi_1 * ctr_1 - _stride_phi_1 + _data_phi_2m1;
              float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * ctr_1 + _data_phi_2m1;
              _data_j_20_39_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_2m1_1m1[_stride_phi_0 * (_size_j_0 - 1)]) + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1)]) * -1.0f + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1)]) * -1.0f) * -0.02351606505734748f * ((1.0f) / (kT));
            }
            if (_size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {
              float *RESTRICT _data_j_20_311 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3;
              float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
              float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
              float *RESTRICT _data_phi_20_10 = _stride_phi_1 * ctr_1 + _data_phi_20;
              float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
              float *RESTRICT _data_rho_2m1_11 = _stride_rho_1 * ctr_1 + _stride_rho_1 + _data_rho_2m1;
              float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
              float *RESTRICT _data_rho_20_10 = _stride_rho_1 * ctr_1 + _data_rho_20;
              float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
              float *RESTRICT _data_phi_2m1_11 = _stride_phi_1 * ctr_1 + _stride_phi_1 + _data_phi_2m1;
              _data_j_20_311_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_0 * z * -1.0f * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_1 * z * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_1 * z * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + f_ext_2 * z * -1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + f_ext_2 * z * -1.0f * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * -2.0f * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + kT * 2.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_2m1_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * -1.0f * _data_phi_2m1_11[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + z * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1)] * _data_rho_2m1_11[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 0.04703213011469496f * ((1.0f) / (kT));
            }
          }
        }
      }
      {
        {
          if (_size_j_1 - 1 > 0 && _size_j_2 - 1 > 0 && 1 < _size_j_0 - 1) {
            float *RESTRICT _data_j_20_37 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 7 * _stride_j_3;
            float *RESTRICT _data_j_20_37_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_37;
            float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
            float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
            float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
            float *RESTRICT _data_rho_2m1_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_2m1;
            float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
            float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
            float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
            float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
            float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1;
            float *RESTRICT _data_phi_2m1_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_2m1;
            _data_j_20_37_10[_stride_j_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0] + _data_phi_2m1_10[_stride_phi_0] + _data_phi_2m1_1m1[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[_stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0] + _data_phi_20_10[_stride_phi_0] + _data_phi_2m1_10[_stride_phi_0]) * -1.0f) * -0.028801180074297286f * ((1.0f) / (kT));
          }
          if (_size_j_1 - 1 > 0 && _size_j_2 - 1 > 0) {
            float *RESTRICT _data_j_20_39 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3;
            float *RESTRICT _data_j_20_39_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_39;
            float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
            float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
            float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
            float *RESTRICT _data_rho_2m1_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_2m1;
            float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
            float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
            float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
            float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
            float *RESTRICT _data_phi_2m1_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_2m1;
            float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1;
            _data_j_20_39_10[_stride_j_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * 2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * (-1.0f * _data_phi_20_10[0] - 1.0f * _data_phi_20_1m1[_stride_phi_0] + _data_phi_2m1_10[0] + _data_phi_2m1_1m1[_stride_phi_0]) + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * (-1.0f * _data_phi_20_1m1[0] - 1.0f * _data_phi_2m1_10[0] + _data_phi_20_1m1[_stride_phi_0] + _data_phi_2m1_10[_stride_phi_0]) * -1.0f + z * (_data_rho_20_10[_stride_rho_0] + _data_rho_2m1_1m1[0]) * (-1.0f * _data_phi_20_1m1[0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0] + _data_phi_20_10[0] + _data_phi_2m1_10[_stride_phi_0]) * -1.0f) * -0.02351606505734748f * ((1.0f) / (kT));
          }
        }
        for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
          if (_size_j_1 - 1 > 0 && _size_j_2 - 1 > 0 && ctr_0 < _size_j_0 - 1) {
            float *RESTRICT _data_j_20_37 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 7 * _stride_j_3;
            float *RESTRICT _data_j_20_37_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_37;
            float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
            float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
            float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
            float *RESTRICT _data_rho_2m1_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_2m1;
            float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
            float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
            float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
            float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
            float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1;
            float *RESTRICT _data_phi_2m1_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_2m1;
            _data_j_20_37_10[_stride_j_0 * ctr_0] = D * (f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0] + _data_phi_2m1_1m1[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) * -1.0f) * -0.028801180074297286f * ((1.0f) / (kT));
          }
          if (_size_j_1 - 1 > 0 && _size_j_2 - 1 > 0) {
            float *RESTRICT _data_j_20_39 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3;
            float *RESTRICT _data_j_20_39_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_39;
            float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
            float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
            float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
            float *RESTRICT _data_rho_2m1_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_2m1;
            float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
            float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
            float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
            float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
            float *RESTRICT _data_phi_2m1_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_2m1;
            float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1;
            _data_j_20_39_10[_stride_j_0 * ctr_0] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_2m1_1m1[_stride_phi_0 * ctr_0]) + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_2m1_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_20_1m1[_stride_phi_0 * ctr_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) * -1.0f + z * (_data_rho_20_10[_stride_rho_0 * ctr_0] + _data_rho_2m1_1m1[_stride_rho_0 * ctr_0 - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * ctr_0 - _stride_phi_0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0 * ctr_0] + _data_phi_20_10[_stride_phi_0 * ctr_0 - _stride_phi_0] + _data_phi_2m1_10[_stride_phi_0 * ctr_0]) * -1.0f) * -0.02351606505734748f * ((1.0f) / (kT));
          }
        }
        if (_size_j_1 - 1 > 0 && _size_j_2 - 1 > 0) {
          float *RESTRICT _data_j_20_39 = _data_j + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3;
          float *RESTRICT _data_j_20_39_10 = _stride_j_1 * (_size_j_1 - 1) + _data_j_20_39;
          float *RESTRICT _data_rho_20 = _data_rho + _stride_rho_2 * (_size_j_2 - 1);
          float *RESTRICT _data_rho_20_10 = _stride_rho_1 * (_size_j_1 - 1) + _data_rho_20;
          float *RESTRICT _data_rho_2m1 = _data_rho + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2;
          float *RESTRICT _data_rho_2m1_1m1 = _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _data_rho_2m1;
          float *RESTRICT _data_phi_20 = _data_phi + _stride_phi_2 * (_size_j_2 - 1);
          float *RESTRICT _data_phi_20_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_20;
          float *RESTRICT _data_phi_20_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_20;
          float *RESTRICT _data_phi_2m1 = _data_phi + _stride_phi_2 * (_size_j_2 - 1) - _stride_phi_2;
          float *RESTRICT _data_phi_2m1_1m1 = _stride_phi_1 * (_size_j_1 - 1) - _stride_phi_1 + _data_phi_2m1;
          float *RESTRICT _data_phi_2m1_10 = _stride_phi_1 * (_size_j_1 - 1) + _data_phi_2m1;
          _data_j_20_39_10[_stride_j_0 * (_size_j_0 - 1)] = D * (f_ext_0 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + f_ext_1 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + f_ext_2 * z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 2.0f + kT * (-1.0f * _data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * 4.0f + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_2m1_1m1[_stride_phi_0 * (_size_j_0 - 1)]) + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1)]) * -1.0f + z * (_data_rho_20_10[_stride_rho_0 * (_size_j_0 - 1)] + _data_rho_2m1_1m1[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0]) * (-1.0f * _data_phi_20_1m1[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] - 1.0f * _data_phi_2m1_1m1[_stride_phi_0 * (_size_j_0 - 1)] + _data_phi_20_10[_stride_phi_0 * (_size_j_0 - 1) - _stride_phi_0] + _data_phi_2m1_10[_stride_phi_0 * (_size_j_0 - 1)]) * -1.0f) * -0.02351606505734748f * ((1.0f) / (kT));
        }
      }
    }
  }
}
} // namespace internal_823ab2463d465630661d5edc8f90930c

void DiffusiveFluxKernelWithElectrostatic_single_precision::run(IBlock *block) {
  auto phi = block->getData<field::GhostLayerField<float, 1>>(phiID);
  auto j = block->getData<field::GhostLayerField<float, 13>>(jID);
  auto rho = block->getData<field::GhostLayerField<float, 1>>(rhoID);

  auto &kT = this->kT_;
  auto &f_ext_2 = this->f_ext_2_;
  auto &f_ext_1 = this->f_ext_1_;
  auto &D = this->D_;
  auto &z = this->z_;
  auto &f_ext_0 = this->f_ext_0_;
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(j->nrOfGhostLayers()));
  float *RESTRICT const _data_j = j->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(phi->nrOfGhostLayers()));
  float *RESTRICT const _data_phi = phi->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho->nrOfGhostLayers()));
  float *RESTRICT const _data_rho = rho->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(), int64_t(cell_idx_c(j->xSize()) + 2));
  const int64_t _size_j_0 = int64_t(cell_idx_c(j->xSize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(), int64_t(cell_idx_c(j->ySize()) + 2));
  const int64_t _size_j_1 = int64_t(cell_idx_c(j->ySize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(), int64_t(cell_idx_c(j->zSize()) + 2));
  const int64_t _size_j_2 = int64_t(cell_idx_c(j->zSize()) + 2);
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
  internal_823ab2463d465630661d5edc8f90930c::diffusivefluxkernelwithelectrostatic_single_precision_diffusivefluxkernelwithelectrostatic_single_precision(D, _data_j, _data_phi, _data_rho, _size_j_0, _size_j_1, _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2, _stride_j_3, _stride_phi_0, _stride_phi_1, _stride_phi_2, _stride_rho_0, _stride_rho_1, _stride_rho_2, f_ext_0, f_ext_1, f_ext_2, kT, z);
}

void DiffusiveFluxKernelWithElectrostatic_single_precision::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block) {
  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto phi = block->getData<field::GhostLayerField<float, 1>>(phiID);
  auto j = block->getData<field::GhostLayerField<float, 13>>(jID);
  auto rho = block->getData<field::GhostLayerField<float, 1>>(rhoID);

  auto &kT = this->kT_;
  auto &f_ext_2 = this->f_ext_2_;
  auto &f_ext_1 = this->f_ext_1_;
  auto &D = this->D_;
  auto &z = this->z_;
  auto &f_ext_0 = this->f_ext_0_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(j->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(j->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(j->nrOfGhostLayers()));
  float *RESTRICT const _data_j = j->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(phi->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(phi->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(phi->nrOfGhostLayers()));
  float *RESTRICT const _data_phi = phi->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho->nrOfGhostLayers()));
  float *RESTRICT const _data_rho = rho->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 2));
  const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 2));
  const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 2));
  const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 2);
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
  internal_823ab2463d465630661d5edc8f90930c::diffusivefluxkernelwithelectrostatic_single_precision_diffusivefluxkernelwithelectrostatic_single_precision(D, _data_j, _data_phi, _data_rho, _size_j_0, _size_j_1, _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2, _stride_j_3, _stride_phi_0, _stride_phi_1, _stride_phi_2, _stride_rho_0, _stride_rho_1, _stride_rho_2, f_ext_0, f_ext_1, f_ext_2, kT, z);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif