/*
 * Copyright (C) 2022-2026 The ESPResSo project
 * Copyright (C) 2020-2025 The waLBerla project
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

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34,
// sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit
// 007e77e077ad9d22b5eed6f3d3118240993e553c

/*
 * Boundary class.
 * Adapted from the waLBerla source file
 * https://i10git.cs.fau.de/walberla/walberla/-/blob/3e54d4f2336e47168ad87e3caaf7b3b082d86ca7/python/pystencils_walberla/templates/Boundary.tmpl.h
 */

#pragma once

#include <core/DataTypes.h>

#include <blockforest/StructuredBlockForest.h>
#include <core/debug/Debug.h>
#include <domain_decomposition/BlockDataID.h>
#include <domain_decomposition/IBlock.h>
#include <field/FlagField.h>
#include <field/GhostLayerField.h>

#include <array>
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
namespace lbm {

class DynamicUBBDoublePrecision {
public:
  struct IndexInfo {
    int32_t x;
    int32_t y;
    int32_t z;
    int32_t dir;
    double vel_0;
    double vel_1;
    double vel_2;
    IndexInfo(int32_t x_, int32_t y_, int32_t z_, int32_t dir_)
        : x(x_), y(y_), z(z_), dir(dir_), vel_0(), vel_1(), vel_2() {}
    bool operator==(const IndexInfo &o) const {
      return x == o.x && y == o.y && z == o.z && dir == o.dir &&
             floatIsEqual(vel_0, o.vel_0) && floatIsEqual(vel_1, o.vel_1) &&
             floatIsEqual(vel_2, o.vel_2);
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

    auto &indexVector(Type t) { return cpuVectors_[t]; }
    auto const &indexVector(Type t) const { return cpuVectors_[t]; }
    IndexInfo *pointerCpu(Type t) {
      return cpuVectors_[t].empty() ? nullptr : cpuVectors_[t].data();
    }

    void syncGPU() {}

  private:
    std::vector<CpuIndexVector> cpuVectors_{NUM_TYPES};
  };

  struct ForceStruct {
    double F_0;
    double F_1;
    double F_2;
    ForceStruct()
        : F_0(double_c(0.0)), F_1(double_c(0.0)), F_2(double_c(0.0)) {}
    bool operator==(const ForceStruct &o) const {
      return floatIsEqual(F_0, o.F_0) && floatIsEqual(F_1, o.F_1) &&
             floatIsEqual(F_2, o.F_2);
    }
  };

  class ForceVector {
  public:
    ForceVector() = default;
    bool operator==(ForceVector const &other) const {
      return other.cpuVector_ == cpuVector_;
    }

    auto &forceVector() { return cpuVector_; }
    auto const &forceVector() const { return cpuVector_; }
    ForceStruct *pointerCpu() {
      return cpuVector_.empty() ? nullptr : cpuVector_.data();
    }
    bool empty() const { return cpuVector_.empty(); }

    Vector3<double> getForce() {
      syncCPU();
      Vector3<double> result(double_c(0.0));
      for (auto const &force : cpuVector_) {
        result[0] += force.F_0;
        result[1] += force.F_1;
        result[2] += force.F_2;
      }
      return result;
    }

    void syncGPU() {}

    void syncCPU() {}

  private:
    std::vector<ForceStruct> cpuVector_;
  };

  DynamicUBBDoublePrecision(
      const std::shared_ptr<StructuredBlockForest> &blocks, BlockDataID pdfsID_,
      std::function<Vector3<double>(
          const Cell &, const shared_ptr<StructuredBlockForest> &, IBlock &)>
          &velocityCallbackDynamicUBBDoublePrecision)
      : elementInitialiser(velocityCallbackDynamicUBBDoublePrecision),
        pdfsID(pdfsID_) {
    auto createIdxVector = [](IBlock *const, StructuredBlockStorage *const) {
      return new IndexVectors();
    };
    indexVectorID = blocks->addStructuredBlockData<IndexVectors>(
        createIdxVector, "IndexField_DynamicUBBDoublePrecision");
    auto createForceVector = [](IBlock *const, StructuredBlockStorage *const) {
      return new ForceVector();
    };
    forceVectorID = blocks->addStructuredBlockData<ForceVector>(
        createForceVector, "forceVector_DynamicUBBDoublePrecision");
  }

  void run(IBlock *block);

  void operator()(IBlock *block) { run(block); }

  void inner(IBlock *block);

  void outer(IBlock *block);

  Vector3<double> getForce(IBlock *block) {
    auto *forceVector = block->getData<ForceVector>(forceVectorID);
    if (forceVector->empty())
      return Vector3<double>(double_c(0.0));
    return forceVector->getForce();
  }

  std::function<void(IBlock *)> getSweep() {
    return [this](IBlock *b) { this->run(b); };
  }

  std::function<void(IBlock *)> getInnerSweep() {
    return [this](IBlock *b) { this->inner(b); };
  }

  std::function<void(IBlock *)> getOuterSweep() {
    return [this](IBlock *b) { this->outer(b); };
  }

  template <typename FlagField_T>
  void fillFromFlagField(const std::shared_ptr<StructuredBlockForest> &blocks,
                         ConstBlockDataID flagFieldID, FlagUID boundaryFlagUID,
                         FlagUID domainFlagUID) {
    for (auto &block : *blocks)
      fillFromFlagField<FlagField_T>(blocks, &block, flagFieldID,
                                     boundaryFlagUID, domainFlagUID);
  }

  template <typename FlagField_T>
  void fillFromFlagField(const shared_ptr<StructuredBlockForest> &blocks,
                         IBlock *block, ConstBlockDataID flagFieldID,
                         FlagUID boundaryFlagUID, FlagUID domainFlagUID) {
    auto *indexVectors = block->getData<IndexVectors>(indexVectorID);
    auto &indexVectorAll = indexVectors->indexVector(IndexVectors::ALL);
    auto &indexVectorInner = indexVectors->indexVector(IndexVectors::INNER);
    auto &indexVectorOuter = indexVectors->indexVector(IndexVectors::OUTER);
    auto *forceVector = block->getData<ForceVector>(forceVectorID);

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

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(0, 0, 0, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 0);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() + 0, it.y() + 0, it.z() + 0), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(0, 1, 0, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 1);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() + 0, it.y() + 1, it.z() + 0), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(0, -1, 0, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 2);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() + 0, it.y() - 1, it.z() + 0), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(-1, 0, 0, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 3);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() - 1, it.y() + 0, it.z() + 0), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(1, 0, 0, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 4);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() + 1, it.y() + 0, it.z() + 0), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(0, 0, 1, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 5);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() + 0, it.y() + 0, it.z() + 1), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(0, 0, -1, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 6);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() + 0, it.y() + 0, it.z() - 1), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(-1, 1, 0, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 7);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() - 1, it.y() + 1, it.z() + 0), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(1, 1, 0, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 8);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() + 1, it.y() + 1, it.z() + 0), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(-1, -1, 0, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 9);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() - 1, it.y() - 1, it.z() + 0), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(1, -1, 0, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 10);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() + 1, it.y() - 1, it.z() + 0), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(0, 1, 1, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 11);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() + 0, it.y() + 1, it.z() + 1), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(0, -1, 1, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 12);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() + 0, it.y() - 1, it.z() + 1), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(-1, 0, 1, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 13);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() - 1, it.y() + 0, it.z() + 1), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(1, 0, 1, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 14);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() + 1, it.y() + 0, it.z() + 1), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(0, 1, -1, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 15);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() + 0, it.y() + 1, it.z() - 1), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(0, -1, -1, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 16);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() + 0, it.y() - 1, it.z() - 1), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(-1, 0, -1, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 17);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() - 1, it.y() + 0, it.z() - 1), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    for (auto it = flagField->beginWithGhostLayerXYZ(
             cell_idx_c(flagField->nrOfGhostLayers() - 1));
         it != flagField->end(); ++it) {
      if (!isFlagSet(it, domainFlag) || isFlagSet(it, boundaryFlag))
        continue;

      if (isFlagSet(it.neighbor(1, 0, -1, 0), boundaryFlag)) {
        auto element = IndexInfo(it.x(), it.y(), it.z(), 18);
        auto const InitialisationAdditionalData = elementInitialiser(
            Cell(it.x() + 1, it.y() + 0, it.z() - 1), blocks, *block);
        element.vel_0 = InitialisationAdditionalData[0];
        element.vel_1 = InitialisationAdditionalData[1];
        element.vel_2 = InitialisationAdditionalData[2];
        indexVectorAll.emplace_back(element);
        if (inner.contains(it.x(), it.y(), it.z()))
          indexVectorInner.emplace_back(element);
        else
          indexVectorOuter.emplace_back(element);
      }
    }

    indexVectors->syncGPU();
    forceVector->forceVector().resize(indexVectorAll.size());
    forceVector->syncGPU();
  }

private:
  void run_impl(IBlock *block, IndexVectors::Type type);

  BlockDataID indexVectorID;
  BlockDataID forceVectorID;
  std::function<Vector3<double>(
      const Cell &, const shared_ptr<StructuredBlockForest> &, IBlock &)>
      elementInitialiser;

public:
  static constexpr std::array<std::array<int, 19u>, 3u> neighborOffset = {{
      {0, 0, 0, -1, 1, 0, 0, -1, 1, -1, 1, 0, 0, -1, 1, 0, 0, -1, 1},
      {0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, -1, 0, 0, 1, -1, 0, 0},
      {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1},
  }};

  auto const &getForceVector(IBlock const *block) const {
    auto const *forceVector = block->getData<ForceVector>(forceVectorID);
    return forceVector->forceVector();
  }

  auto const &getIndexVector(IBlock const *block) const {
    auto const *indexVectors = block->getData<IndexVectors>(indexVectorID);
    return indexVectors->indexVector(IndexVectors::ALL);
  }

  BlockDataID getIndexVectorID() const { return indexVectorID; }
  BlockDataID getForceVectorID() const { return forceVectorID; }

public:
  BlockDataID pdfsID;
};

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__) or defined(__GNUG__)
#pragma GCC diagnostic pop
#endif

} // namespace lbm
} // namespace walberla
