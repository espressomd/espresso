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
//! \\file FrictionCouplingKernel_double_precision.h
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.7+13.gdfd203a, lbmpy
// v1.3.7+10.gd3f6236, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from
// waLBerla commit c69cb11d6a95d32b2280544d3d9abde1fe5fdbb5

#pragma once
#include "core/DataTypes.h"
#include "core/logging/Logging.h"

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/GhostLayerField.h"
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

class FrictionCouplingKernel_double_precision {
public:
  FrictionCouplingKernel_double_precision(BlockDataID fID_, BlockDataID jID_,
                                          double D, double kT, double rho_lb)
      : fID(fID_), jID(jID_), D_(D), kT_(kT), rho_lb_(rho_lb) {}

  void run(IBlock *block);

  void runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks,
                         const CellInterval &globalCellInterval,
                         cell_idx_t ghostLayers, IBlock *block);

  void operator()(IBlock *block) { run(block); }

  static std::function<void(IBlock *)>
  getSweep(const shared_ptr<FrictionCouplingKernel_double_precision> &kernel) {
    return [kernel](IBlock *b) { kernel->run(b); };
  }

  static std::function<void(IBlock *)> getSweepOnCellInterval(
      const shared_ptr<FrictionCouplingKernel_double_precision> &kernel,
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

  void configure(const shared_ptr<StructuredBlockStorage> & /*blocks*/,
                 IBlock * /*block*/) {}

  inline double getD() const { return D_; }
  inline double getKt() const { return kT_; }
  inline double getRho_lb() const { return rho_lb_; }
  inline void setD(const double value) { D_ = value; }
  inline void setKt(const double value) { kT_ = value; }
  inline void setRho_lb(const double value) { rho_lb_ = value; }

private:
  BlockDataID fID;
  BlockDataID jID;
  double D_;
  double kT_;
  double rho_lb_;
};

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif
