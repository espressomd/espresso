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
//! \\file StreamCollideSweepLeesEdwardsDoublePrecisionCUDA.h
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34,
// sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit
// 272d4a09ec35da50685afc9586645e1b9984b423

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

class StreamCollideSweepLeesEdwardsDoublePrecisionCUDA {
public:
  StreamCollideSweepLeesEdwardsDoublePrecisionCUDA(
      BlockDataID forceID_, BlockDataID pdfsID_, int64_t lebc_bot_index,
      int64_t lebc_top_index, double omega_bulk, double omega_even,
      double omega_odd, double omega_shear, double v_s)
      : forceID(forceID_), pdfsID(pdfsID_), lebc_bot_index_(lebc_bot_index),
        lebc_top_index_(lebc_top_index), omega_bulk_(omega_bulk),
        omega_even_(omega_even), omega_odd_(omega_odd),
        omega_shear_(omega_shear), v_s_(v_s) {}

  ~StreamCollideSweepLeesEdwardsDoublePrecisionCUDA() {
    for (auto p : cache_pdfs_) {
      delete p.second;
    }
  }

  void run(IBlock *block, gpuStream_t stream = nullptr);

  void runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks,
                         const CellInterval &globalCellInterval,
                         cell_idx_t ghostLayers, IBlock *block,
                         gpuStream_t stream = nullptr);

  void operator()(IBlock *block, gpuStream_t stream = nullptr) {
    run(block, stream);
  }

  static std::function<void(IBlock *)>
  getSweep(const shared_ptr<StreamCollideSweepLeesEdwardsDoublePrecisionCUDA>
               &kernel) {
    return [kernel](IBlock *b) { kernel->run(b); };
  }

  static std::function<void(IBlock *, gpuStream_t)> getSweepOnCellInterval(
      const shared_ptr<StreamCollideSweepLeesEdwardsDoublePrecisionCUDA>
          &kernel,
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

  inline int64_t getLebc_bot_index() const { return lebc_bot_index_; }
  inline int64_t getLebc_top_index() const { return lebc_top_index_; }
  inline double getOmega_bulk() const { return omega_bulk_; }
  inline double getOmega_even() const { return omega_even_; }
  inline double getOmega_odd() const { return omega_odd_; }
  inline double getOmega_shear() const { return omega_shear_; }
  inline double getV_s() const { return v_s_; }
  inline void setLebc_bot_index(const int64_t value) {
    lebc_bot_index_ = value;
  }
  inline void setLebc_top_index(const int64_t value) {
    lebc_top_index_ = value;
  }
  inline void setOmega_bulk(const double value) { omega_bulk_ = value; }
  inline void setOmega_even(const double value) { omega_even_ = value; }
  inline void setOmega_odd(const double value) { omega_odd_ = value; }
  inline void setOmega_shear(const double value) { omega_shear_ = value; }
  inline void setV_s(const double value) { v_s_ = value; }

private:
  BlockDataID forceID;
  BlockDataID pdfsID;
  int64_t lebc_bot_index_;
  int64_t lebc_top_index_;
  double omega_bulk_;
  double omega_even_;
  double omega_odd_;
  double omega_shear_;
  double v_s_;
  std::unordered_map<IBlock *, gpu::GPUField<double> *> cache_pdfs_;
};

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif
