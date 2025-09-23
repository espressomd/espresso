/*
 * Copyright (C) 2022-2023 The ESPResSo project
 * Copyright (C) 2020-2023 The waLBerla project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// kernel generated with pystencils v1.3.7+13.gdfd203a, lbmpy
// v1.3.7+10.gd3f6236, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from
// waLBerla commit e12db9965373887d86aab4aaaf4dd7b38fa588e8

/*
 * Boundary class.
 * Adapted from the waLBerla source file
 * https://i10git.cs.fau.de/walberla/walberla/-/blob/e12db9965373887d86aab4aaaf4dd7b38fa588e8/python/pystencils_walberla/templates/Boundary.tmpl.h
 */

#pragma once

#include <core/DataTypes.h>

#include <blockforest/StructuredBlockForest.h>
#include <core/debug/Debug.h>
#include <domain_decomposition/BlockDataID.h>
#include <domain_decomposition/IBlock.h>
#include <field/FlagField.h>
#include <gpu/FieldCopy.h>
#include <gpu/GPUField.h>
#include <gpu/GPUWrapper.h>

#include <cassert>
#include <functional>
#include <memory>
#include <vector>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#elif defined(__GNUC__) or defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
using walberla::half;
#endif

namespace walberla {
namespace pystencils {

class ReactionKernelIndexed_1_double_precision_CUDA {
public:
  struct IndexInfo {
    int32_t x;
    int32_t y;
    int32_t z;
    IndexInfo(int32_t x_, int32_t y_, int32_t z_) : x(x_), y(y_), z(z_) {}
    bool operator==(const IndexInfo &o) const {
      return x == o.x && y == o.y && z == o.z;
    }
  };

  class IndexVectors {
  public:
    using CpuIndexVector = std::vector<IndexInfo>;

    enum Type { ALL = 0, INNER = 1, OUTER = 2, NUM_TYPES = 3 };

    IndexVectors() = default;
    bool operator==(IndexVectors const &other) const {
      return other.cpuVectors_ == cpuVectors_;
    }

    ~IndexVectors() {
      for (auto &gpuVec : gpuVectors_) {
        if (gpuVec) {
          WALBERLA_GPU_CHECK(gpuFree(gpuVec));
        }
      }
    }
    CpuIndexVector &indexVector(Type t) { return cpuVectors_[t]; }
    IndexInfo *pointerCpu(Type t) {
      return cpuVectors_[t].empty() ? nullptr : cpuVectors_[t].data();
    }

    IndexInfo *pointerGpu(Type t) { return gpuVectors_[t]; }
    void syncGPU() {
      for (auto &gpuVec : gpuVectors_)
        WALBERLA_GPU_CHECK(gpuFree(gpuVec));
      gpuVectors_.resize(cpuVectors_.size());

      WALBERLA_ASSERT_EQUAL(cpuVectors_.size(), NUM_TYPES);
      for (size_t i = 0; i < cpuVectors_.size(); ++i) {
        auto &gpuVec = gpuVectors_[i];
        auto &cpuVec = cpuVectors_[i];
        if (cpuVec.empty()) {
          continue;
        }
        WALBERLA_GPU_CHECK(
            gpuMalloc(&gpuVec, sizeof(IndexInfo) * cpuVec.size()));
        WALBERLA_GPU_CHECK(gpuMemcpy(gpuVec, cpuVec.data(),
                                     sizeof(IndexInfo) * cpuVec.size(),
                                     gpuMemcpyHostToDevice));
      }
    }

  private:
    std::vector<CpuIndexVector> cpuVectors_{NUM_TYPES};

    using GpuIndexVector = IndexInfo *;
    std::vector<GpuIndexVector> gpuVectors_;
  };

  ReactionKernelIndexed_1_double_precision_CUDA(
      const std::shared_ptr<StructuredBlockForest> &blocks,
      BlockDataID rho_0ID_, double order_0, double rate_coefficient,
      double stoech_0)
      : rho_0ID(rho_0ID_), order_0_(order_0),
        rate_coefficient_(rate_coefficient), stoech_0_(stoech_0) {
    auto createIdxVector = [](IBlock *const, StructuredBlockStorage *const) {
      return new IndexVectors();
    };
    indexVectorID = blocks->addStructuredBlockData<IndexVectors>(
        createIdxVector,
        "IndexField_ReactionKernelIndexed_1_double_precision_CUDA");
  }

  ReactionKernelIndexed_1_double_precision_CUDA(BlockDataID indexVectorID_,
                                                BlockDataID rho_0ID_,
                                                double order_0,
                                                double rate_coefficient,
                                                double stoech_0)
      : indexVectorID(indexVectorID_), rho_0ID(rho_0ID_), order_0_(order_0),
        rate_coefficient_(rate_coefficient), stoech_0_(stoech_0) {}

  void run(IBlock *block, gpuStream_t stream = nullptr);

  void operator()(IBlock *block, gpuStream_t stream = nullptr) {
    run(block, stream);
  }

  void inner(IBlock *block, gpuStream_t stream = nullptr);

  void outer(IBlock *block, gpuStream_t stream = nullptr);

  Vector3<real_t> getForce(IBlock * /*block*/) {

    WALBERLA_ABORT(
        "Boundary condition was not generated including force calculation.")
    return Vector3<real_t>(real_c(0.0));
  }

  std::function<void(IBlock *)> getSweep(gpuStream_t stream = nullptr) {
    return [this, stream](IBlock *b) { this->run(b, stream); };
  }

  std::function<void(IBlock *)> getInnerSweep(gpuStream_t stream = nullptr) {
    return [this, stream](IBlock *b) { this->inner(b, stream); };
  }

  std::function<void(IBlock *)> getOuterSweep(gpuStream_t stream = nullptr) {
    return [this, stream](IBlock *b) { this->outer(b, stream); };
  }

  template <typename FlagField_T>
  void fillFromFlagField(const std::shared_ptr<StructuredBlockForest> &blocks,
                         ConstBlockDataID flagFieldID, FlagUID boundaryFlagUID,
                         FlagUID domainFlagUID) {
    for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
      fillFromFlagField<FlagField_T>(&*blockIt, flagFieldID, boundaryFlagUID,
                                     domainFlagUID);
  }

  template <typename FlagField_T>
  void fillFromFlagField(IBlock *block, ConstBlockDataID flagFieldID,
                         FlagUID boundaryFlagUID, FlagUID domainFlagUID) {
    auto *indexVectors = block->getData<IndexVectors>(indexVectorID);
    auto &indexVectorAll = indexVectors->indexVector(IndexVectors::ALL);
    auto &indexVectorInner = indexVectors->indexVector(IndexVectors::INNER);
    auto &indexVectorOuter = indexVectors->indexVector(IndexVectors::OUTER);

    auto *flagField = block->getData<FlagField_T>(flagFieldID);

    if (!(flagField->flagExists(boundaryFlagUID) and
          flagField->flagExists(domainFlagUID)))
      return;

    auto boundaryFlag = flagField->getFlag(boundaryFlagUID);
    auto domainFlag = flagField->getFlag(domainFlagUID);

    auto inner = flagField->xyzSize();
    inner.expand(cell_idx_t(-1));

    indexVectorAll.clear();
    indexVectorInner.clear();
    indexVectorOuter.clear();

    auto flagWithGLayers = flagField->xyzSizeWithGhostLayer();
    for (auto it = flagField->beginWithGhostLayerXYZ(); it != flagField->end();
         ++it) {

      if (!isFlagSet(it, boundaryFlag))
        continue;
      if (flagWithGLayers.contains(it.x() + cell_idx_c(0),
                                   it.y() + cell_idx_c(0),
                                   it.z() + cell_idx_c(0)) &&
          isFlagSet(it.neighbor(0, 0, 0, 0), domainFlag)) {

        auto element = IndexInfo(it.x(), it.y(), it.z(), 0);

        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    indexVectors->syncGPU();
  }

private:
  void run_impl(IBlock *block, IndexVectors::Type type,
                gpuStream_t stream = nullptr);

  BlockDataID indexVectorID;

public:
  BlockDataID rho_0ID;
  double order_0_;
  double rate_coefficient_;
  double stoech_0_;
};

} // namespace pystencils
} // namespace walberla
