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
//! \\file StreamCollideSweepThermalizedDoublePrecision.h
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

class StreamCollideSweepThermalizedDoublePrecision {
public:
  StreamCollideSweepThermalizedDoublePrecision(
      BlockDataID forceID_, BlockDataID pdfsID_, double kT, double omega_bulk,
      double omega_even, double omega_odd, double omega_shear, uint32_t seed,
      uint32_t time_step)
      : forceID(forceID_), pdfsID(pdfsID_), kT_(kT), omega_bulk_(omega_bulk),
        omega_even_(omega_even), omega_odd_(omega_odd),
        omega_shear_(omega_shear), seed_(seed), time_step_(time_step),
        block_offset_0_(uint32_t(0)), block_offset_1_(uint32_t(0)),
        block_offset_2_(uint32_t(0)), configured_(false) {}

  ~StreamCollideSweepThermalizedDoublePrecision() {
    for (auto p : cache_pdfs_) {
      delete p.second;
    }
  }

  void run(IBlock *block);

  void runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks,
                         const CellInterval &globalCellInterval,
                         cell_idx_t ghostLayers, IBlock *block);

  void operator()(IBlock *block) { run(block); }

  static std::function<void(IBlock *)> getSweep(
      const shared_ptr<StreamCollideSweepThermalizedDoublePrecision> &kernel) {
    return [kernel](IBlock *b) { kernel->run(b); };
  }

  static std::function<void(IBlock *)> getSweepOnCellInterval(
      const shared_ptr<StreamCollideSweepThermalizedDoublePrecision> &kernel,
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

  void configure(const shared_ptr<StructuredBlockStorage> &blocks,
                 IBlock *block) {
    Cell BlockCellBB = blocks->getBlockCellBB(*block).min();
    block_offset_0_ = uint32_t(BlockCellBB[0]);
    block_offset_1_ = uint32_t(BlockCellBB[1]);
    block_offset_2_ = uint32_t(BlockCellBB[2]);
    configured_ = true;
  }

  inline uint32_t getBlock_offset_0() const { return block_offset_0_; }
  inline uint32_t getBlock_offset_1() const { return block_offset_1_; }
  inline uint32_t getBlock_offset_2() const { return block_offset_2_; }
  inline double getKt() const { return kT_; }
  inline double getOmega_bulk() const { return omega_bulk_; }
  inline double getOmega_even() const { return omega_even_; }
  inline double getOmega_odd() const { return omega_odd_; }
  inline double getOmega_shear() const { return omega_shear_; }
  inline uint32_t getSeed() const { return seed_; }
  inline uint32_t getTime_step() const { return time_step_; }
  inline void setBlock_offset_0(const uint32_t value) {
    block_offset_0_ = value;
  }
  inline void setBlock_offset_1(const uint32_t value) {
    block_offset_1_ = value;
  }
  inline void setBlock_offset_2(const uint32_t value) {
    block_offset_2_ = value;
  }
  inline void setKt(const double value) { kT_ = value; }
  inline void setOmega_bulk(const double value) { omega_bulk_ = value; }
  inline void setOmega_even(const double value) { omega_even_ = value; }
  inline void setOmega_odd(const double value) { omega_odd_ = value; }
  inline void setOmega_shear(const double value) { omega_shear_ = value; }
  inline void setSeed(const uint32_t value) { seed_ = value; }
  inline void setTime_step(const uint32_t value) { time_step_ = value; }

private:
  BlockDataID forceID;
  BlockDataID pdfsID;
  uint32_t block_offset_0_;
  uint32_t block_offset_1_;
  uint32_t block_offset_2_;
  double kT_;
  double omega_bulk_;
  double omega_even_;
  double omega_odd_;
  double omega_shear_;
  uint32_t seed_;
  uint32_t time_step_;
  std::unordered_map<IBlock *, field::GhostLayerField<double, 19> *>
      cache_pdfs_;

  bool configured_;
};

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif
