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
//! \\file
//! DiffusiveFluxKernelWithElectrostaticThermalized_single_precision_CUDA.h
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

class DiffusiveFluxKernelWithElectrostaticThermalized_single_precision_CUDA {
public:
  DiffusiveFluxKernelWithElectrostaticThermalized_single_precision_CUDA(
      BlockDataID jID_, BlockDataID phiID_, BlockDataID rhoID_, float D,
      float f_ext_0, float f_ext_1, float f_ext_2, uint32_t field_size_0,
      uint32_t field_size_1, uint32_t field_size_2, float kT, uint32_t seed,
      uint32_t time_step, float z)
      : jID(jID_), phiID(phiID_), rhoID(rhoID_), D_(D), f_ext_0_(f_ext_0),
        f_ext_1_(f_ext_1), f_ext_2_(f_ext_2), field_size_0_(field_size_0),
        field_size_1_(field_size_1), field_size_2_(field_size_2), kT_(kT),
        seed_(seed), time_step_(time_step), z_(z), block_offset_0_(uint32_t(0)),
        block_offset_1_(uint32_t(0)), block_offset_2_(uint32_t(0)),
        configured_(false) {}

  void run(IBlock *block, gpuStream_t stream = nullptr);

  void runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks,
                         const CellInterval &globalCellInterval,
                         cell_idx_t ghostLayers, IBlock *block,
                         gpuStream_t stream = nullptr);

  void operator()(IBlock *block, gpuStream_t stream = nullptr) {
    run(block, stream);
  }

  static std::function<void(IBlock *)> getSweep(
      const shared_ptr<
          DiffusiveFluxKernelWithElectrostaticThermalized_single_precision_CUDA>
          &kernel) {
    return [kernel](IBlock *b) { kernel->run(b); };
  }

  static std::function<void(IBlock *, gpuStream_t)> getSweepOnCellInterval(
      const shared_ptr<
          DiffusiveFluxKernelWithElectrostaticThermalized_single_precision_CUDA>
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

  void configure(const shared_ptr<StructuredBlockStorage> &blocks,
                 IBlock *block) {
    Cell BlockCellBB = blocks->getBlockCellBB(*block).min();
    block_offset_0_ = uint32_t(BlockCellBB[0]);
    block_offset_1_ = uint32_t(BlockCellBB[1]);
    block_offset_2_ = uint32_t(BlockCellBB[2]);
    configured_ = true;
  }

  inline float getD() const { return D_; }
  inline uint32_t getBlock_offset_0() const { return block_offset_0_; }
  inline uint32_t getBlock_offset_1() const { return block_offset_1_; }
  inline uint32_t getBlock_offset_2() const { return block_offset_2_; }
  inline float getF_ext_0() const { return f_ext_0_; }
  inline float getF_ext_1() const { return f_ext_1_; }
  inline float getF_ext_2() const { return f_ext_2_; }
  inline uint32_t getField_size_0() const { return field_size_0_; }
  inline uint32_t getField_size_1() const { return field_size_1_; }
  inline uint32_t getField_size_2() const { return field_size_2_; }
  inline float getKt() const { return kT_; }
  inline uint32_t getSeed() const { return seed_; }
  inline uint32_t getTime_step() const { return time_step_; }
  inline float getZ() const { return z_; }
  inline void setD(const float value) { D_ = value; }
  inline void setBlock_offset_0(const uint32_t value) {
    block_offset_0_ = value;
  }
  inline void setBlock_offset_1(const uint32_t value) {
    block_offset_1_ = value;
  }
  inline void setBlock_offset_2(const uint32_t value) {
    block_offset_2_ = value;
  }
  inline void setF_ext_0(const float value) { f_ext_0_ = value; }
  inline void setF_ext_1(const float value) { f_ext_1_ = value; }
  inline void setF_ext_2(const float value) { f_ext_2_ = value; }
  inline void setField_size_0(const uint32_t value) { field_size_0_ = value; }
  inline void setField_size_1(const uint32_t value) { field_size_1_ = value; }
  inline void setField_size_2(const uint32_t value) { field_size_2_ = value; }
  inline void setKt(const float value) { kT_ = value; }
  inline void setSeed(const uint32_t value) { seed_ = value; }
  inline void setTime_step(const uint32_t value) { time_step_ = value; }
  inline void setZ(const float value) { z_ = value; }

private:
  BlockDataID jID;
  BlockDataID phiID;

public:
  inline void setPhiID(BlockDataID phiID_) { phiID = phiID_; }

private:
  BlockDataID rhoID;
  float D_;
  uint32_t block_offset_0_;
  uint32_t block_offset_1_;
  uint32_t block_offset_2_;
  float f_ext_0_;
  float f_ext_1_;
  float f_ext_2_;
  uint32_t field_size_0_;
  uint32_t field_size_1_;
  uint32_t field_size_2_;
  float kT_;
  uint32_t seed_;
  uint32_t time_step_;
  float z_;

  bool configured_;
};

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif
