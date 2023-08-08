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
//! \\file DiffusiveFluxKernelWithElectrostaticThermalized_double_precision.h
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.2, lbmpy v1.2,
// lbmpy_walberla/pystencils_walberla from waLBerla commit
// 065ce5f311850371a97ac4766f47dbb5ca8424ba

#pragma once
#include "core/DataTypes.h"

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/GhostLayerField.h"
#include "field/SwapableCompare.h"
#include <set>

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wreorder"
#endif

namespace walberla {
namespace pystencils {

class DiffusiveFluxKernelWithElectrostaticThermalized_double_precision {
public:
  DiffusiveFluxKernelWithElectrostaticThermalized_double_precision(
      BlockDataID jID_, BlockDataID phiID_, BlockDataID rhoID_, double D,
      uint32_t block_offset_0, uint32_t block_offset_1, uint32_t block_offset_2,
      double f_ext_0, double f_ext_1, double f_ext_2, uint32_t field_size_0,
      uint32_t field_size_1, uint32_t field_size_2, double kT, uint32_t seed,
      uint32_t time_step, double z)
      : jID(jID_), phiID(phiID_), rhoID(rhoID_), D_(D),
        block_offset_0_(block_offset_0), block_offset_1_(block_offset_1),
        block_offset_2_(block_offset_2), f_ext_0_(f_ext_0), f_ext_1_(f_ext_1),
        f_ext_2_(f_ext_2), field_size_0_(field_size_0),
        field_size_1_(field_size_1), field_size_2_(field_size_2), kT_(kT),
        seed_(seed), time_step_(time_step), z_(z){};

  void run(IBlock *block);

  void runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks,
                         const CellInterval &globalCellInterval,
                         cell_idx_t ghostLayers, IBlock *block);

  void operator()(IBlock *block) { run(block); }

  static std::function<void(IBlock *)>
  getSweep(const shared_ptr<
           DiffusiveFluxKernelWithElectrostaticThermalized_double_precision>
               &kernel) {
    return [kernel](IBlock *b) { kernel->run(b); };
  }

  static std::function<void(IBlock *)> getSweepOnCellInterval(
      const shared_ptr<
          DiffusiveFluxKernelWithElectrostaticThermalized_double_precision>
          &kernel,
      const shared_ptr<StructuredBlockStorage> &blocks,
      const CellInterval &globalCellInterval, cell_idx_t ghostLayers = 1) {
    return [kernel, blocks, globalCellInterval, ghostLayers](IBlock *b) {
      kernel->runOnCellInterval(blocks, globalCellInterval, ghostLayers, b);
    };
  }

  std::function<void(IBlock *)> getSweep() {
    return [this](IBlock *b) { this->run(b); };
  }

  std::function<void(IBlock *)>
  getSweepOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks,
                         const CellInterval &globalCellInterval,
                         cell_idx_t ghostLayers = 1) {
    return [this, blocks, globalCellInterval, ghostLayers](IBlock *b) {
      this->runOnCellInterval(blocks, globalCellInterval, ghostLayers, b);
    };
  }

  BlockDataID jID;
  BlockDataID phiID;
  BlockDataID rhoID;
  double D_;
  uint32_t block_offset_0_;
  uint32_t block_offset_1_;
  uint32_t block_offset_2_;
  double f_ext_0_;
  double f_ext_1_;
  double f_ext_2_;
  uint32_t field_size_0_;
  uint32_t field_size_1_;
  uint32_t field_size_2_;
  double kT_;
  uint32_t seed_;
  uint32_t time_step_;
  double z_;
  std::function<void(IBlock *, uint32_t &, uint32_t &, uint32_t &)>
      block_offset_generator =
          [](IBlock *const, uint32_t &, uint32_t &, uint32_t &) {};
};

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif