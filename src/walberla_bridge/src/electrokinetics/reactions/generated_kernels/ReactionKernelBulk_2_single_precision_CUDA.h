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
//! \\file ReactionKernelBulk_2_single_precision_CUDA.h
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.7+13.gdfd203a, lbmpy
// v1.3.7+10.gd3f6236, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from
// waLBerla commit c69cb11d6a95d32b2280544d3d9abde1fe5fdbb5

#pragma once
#include "core/DataTypes.h"
#include "core/logging/Logging.h"

#include "gpu/GPUField.h"
#include "gpu/GPUWrapper.h"

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/SwapableCompare.h"

#include <functional>
#include <unordered_map>

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

class ReactionKernelBulk_2_single_precision_CUDA {
public:
  ReactionKernelBulk_2_single_precision_CUDA(BlockDataID rho_0ID_,
                                             BlockDataID rho_1ID_,
                                             float order_0, float order_1,
                                             float rate_coefficient,
                                             float stoech_0, float stoech_1)
      : rho_0ID(rho_0ID_), rho_1ID(rho_1ID_), order_0_(order_0),
        order_1_(order_1), rate_coefficient_(rate_coefficient),
        stoech_0_(stoech_0), stoech_1_(stoech_1) {}

  void run(IBlock *block, gpuStream_t stream = nullptr);

  void runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks,
                         const CellInterval &globalCellInterval,
                         cell_idx_t ghostLayers, IBlock *block,
                         gpuStream_t stream = nullptr);

  void operator()(IBlock *block, gpuStream_t stream = nullptr) {
    run(block, stream);
  }

  static std::function<void(IBlock *)> getSweep(
      const shared_ptr<ReactionKernelBulk_2_single_precision_CUDA> &kernel) {
    return [kernel](IBlock *b) { kernel->run(b); };
  }

  static std::function<void(IBlock *, gpuStream_t)> getSweepOnCellInterval(
      const shared_ptr<ReactionKernelBulk_2_single_precision_CUDA> &kernel,
      const shared_ptr<StructuredBlockStorage> &blocks,
      const CellInterval &globalCellInterval, cell_idx_t ghostLayers = 1) {
    return [kernel, blocks, globalCellInterval,
            ghostLayers](IBlock *b, gpuStream_t stream = nullptr) {
      kernel->runOnCellInterval(blocks, globalCellInterval, ghostLayers, b,
                                stream);
    };
  }

  std::function<void(IBlock *)> getSweep(gpuStream_t stream = nullptr) {
    return [this, stream](IBlock *b) { this->run(b, stream); };
  }

  std::function<void(IBlock *)>
  getSweepOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks,
                         const CellInterval &globalCellInterval,
                         cell_idx_t ghostLayers = 1,
                         gpuStream_t stream = nullptr) {
    return [this, blocks, globalCellInterval, ghostLayers, stream](IBlock *b) {
      this->runOnCellInterval(blocks, globalCellInterval, ghostLayers, b,
                              stream);
    };
  }

  void configure(const shared_ptr<StructuredBlockStorage> & /*blocks*/,
                 IBlock * /*block*/) {}

  inline float getOrder_0() const { return order_0_; }
  inline float getOrder_1() const { return order_1_; }
  inline float getRate_coefficient() const { return rate_coefficient_; }
  inline float getStoech_0() const { return stoech_0_; }
  inline float getStoech_1() const { return stoech_1_; }
  inline void setOrder_0(const float value) { order_0_ = value; }
  inline void setOrder_1(const float value) { order_1_ = value; }
  inline void setRate_coefficient(const float value) {
    rate_coefficient_ = value;
  }
  inline void setStoech_0(const float value) { stoech_0_ = value; }
  inline void setStoech_1(const float value) { stoech_1_ = value; }

private:
  BlockDataID rho_0ID;
  BlockDataID rho_1ID;
  float order_0_;
  float order_1_;
  float rate_coefficient_;
  float stoech_0_;
  float stoech_1_;
};

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif
