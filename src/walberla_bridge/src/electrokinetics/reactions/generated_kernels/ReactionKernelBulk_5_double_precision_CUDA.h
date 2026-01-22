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
//! \\file ReactionKernelBulk_5_double_precision_CUDA.h
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34,
// sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit
// 007e77e077ad9d22b5eed6f3d3118240993e553c

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

class ReactionKernelBulk_5_double_precision_CUDA {
public:
  ReactionKernelBulk_5_double_precision_CUDA(
      BlockDataID rho_0ID_, BlockDataID rho_1ID_, BlockDataID rho_2ID_,
      BlockDataID rho_3ID_, BlockDataID rho_4ID_, double order_0,
      double order_1, double order_2, double order_3, double order_4,
      double rate_coefficient, double stoech_0, double stoech_1,
      double stoech_2, double stoech_3, double stoech_4)
      : rho_0ID(rho_0ID_), rho_1ID(rho_1ID_), rho_2ID(rho_2ID_),
        rho_3ID(rho_3ID_), rho_4ID(rho_4ID_), order_0_(order_0),
        order_1_(order_1), order_2_(order_2), order_3_(order_3),
        order_4_(order_4), rate_coefficient_(rate_coefficient),
        stoech_0_(stoech_0), stoech_1_(stoech_1), stoech_2_(stoech_2),
        stoech_3_(stoech_3), stoech_4_(stoech_4) {}

  void run(IBlock *block, gpuStream_t stream = nullptr);

  void runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks,
                         const CellInterval &globalCellInterval,
                         cell_idx_t ghostLayers, IBlock *block,
                         gpuStream_t stream = nullptr);

  void operator()(IBlock *block, gpuStream_t stream = nullptr) {
    run(block, stream);
  }

  static std::function<void(IBlock *)> getSweep(
      const shared_ptr<ReactionKernelBulk_5_double_precision_CUDA> &kernel) {
    return [kernel](IBlock *b) { kernel->run(b); };
  }

  static std::function<void(IBlock *, gpuStream_t)> getSweepOnCellInterval(
      const shared_ptr<ReactionKernelBulk_5_double_precision_CUDA> &kernel,
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

  inline double getOrder_0() const { return order_0_; }
  inline double getOrder_1() const { return order_1_; }
  inline double getOrder_2() const { return order_2_; }
  inline double getOrder_3() const { return order_3_; }
  inline double getOrder_4() const { return order_4_; }
  inline double getRate_coefficient() const { return rate_coefficient_; }
  inline double getStoech_0() const { return stoech_0_; }
  inline double getStoech_1() const { return stoech_1_; }
  inline double getStoech_2() const { return stoech_2_; }
  inline double getStoech_3() const { return stoech_3_; }
  inline double getStoech_4() const { return stoech_4_; }
  inline void setOrder_0(const double value) { order_0_ = value; }
  inline void setOrder_1(const double value) { order_1_ = value; }
  inline void setOrder_2(const double value) { order_2_ = value; }
  inline void setOrder_3(const double value) { order_3_ = value; }
  inline void setOrder_4(const double value) { order_4_ = value; }
  inline void setRate_coefficient(const double value) {
    rate_coefficient_ = value;
  }
  inline void setStoech_0(const double value) { stoech_0_ = value; }
  inline void setStoech_1(const double value) { stoech_1_ = value; }
  inline void setStoech_2(const double value) { stoech_2_ = value; }
  inline void setStoech_3(const double value) { stoech_3_ = value; }
  inline void setStoech_4(const double value) { stoech_4_ = value; }

private:
  BlockDataID rho_0ID;
  BlockDataID rho_1ID;
  BlockDataID rho_2ID;
  BlockDataID rho_3ID;
  BlockDataID rho_4ID;
  double order_0_;
  double order_1_;
  double order_2_;
  double order_3_;
  double order_4_;
  double rate_coefficient_;
  double stoech_0_;
  double stoech_1_;
  double stoech_2_;
  double stoech_3_;
  double stoech_4_;
};

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif
