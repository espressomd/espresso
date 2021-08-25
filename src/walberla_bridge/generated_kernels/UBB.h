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
//! \file UBB.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/field/DensityAndVelocity.h"
#include "lbm/field/PdfField.h"

#include "boundary/Boundary.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/config/Config.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/math/Vector3.h"

#include "field/FlagUID.h"

#include "stencil/Directions.h"

#include <vector>



namespace walberla {
namespace lbm {
namespace espresso {

template <typename LatticeModel_T, typename flag_t,
          bool AdaptVelocityToExternalForce = false, bool StoreForce = false>
class UBB : public Boundary<flag_t> {
  using PDFField = GhostLayerField<real_t, LatticeModel_T::Stencil::Size>;
  using Stencil = typename LatticeModel_T::Stencil;

  using VelField = GhostLayerField<Vector3<real_t>, 1>;

  using ForceField = GhostLayerField<Vector3<real_t>, 1>;

public:
  static const bool threadsafe = true;

  class Velocity : public BoundaryConfiguration {
  public:
    Velocity(const Vector3<real_t> &_velocity) : velocity_(_velocity) {}
    Velocity(const real_t _x, const real_t _y, const real_t _z)
        : velocity_(_x, _y, _z) {}
    inline Velocity(const Config::BlockHandle &config);

    const Vector3<real_t> &velocity() const { return velocity_; }

    const real_t &x() const { return velocity_[0]; }
    const real_t &y() const { return velocity_[1]; }
    const real_t &z() const { return velocity_[2]; }

    Vector3<real_t> &velocity() { return velocity_; }

    real_t &x() { return velocity_[0]; }
    real_t &y() { return velocity_[1]; }
    real_t &z() { return velocity_[2]; }

  private:
    Vector3<real_t> velocity_;
  };

  static shared_ptr<Velocity>
  createConfiguration(const Config::BlockHandle &config) {
    return make_shared<Velocity>(config);
  }

  inline UBB(const BoundaryUID &boundaryUID, const FlagUID &uid,
             const LatticeModel_T * lm, PDFField *const pdfField, GhostLayerField<real_t, 3> *force_field,
             FlagField<flag_t> *const flagField = nullptr);

  void pushFlags(std::vector<FlagUID> &uids) const { uids.push_back(uid_); }

  inline void beforeBoundaryTreatment() const;
  void afterBoundaryTreatment() const {}

  template <typename Buffer_T>
  inline void packCell(Buffer_T &buffer, const cell_idx_t x, const cell_idx_t y,
                       const cell_idx_t z) const;

  template <typename Buffer_T>
  inline void registerCell(Buffer_T &buffer, const flag_t, const cell_idx_t x,
                           const cell_idx_t y, const cell_idx_t z);

  inline void registerCell(const flag_t, const cell_idx_t x, const cell_idx_t y,
                           const cell_idx_t z,
                           const BoundaryConfiguration &velocity);
  inline void registerCells(const flag_t, const CellInterval &cells,
                            const BoundaryConfiguration &velocity);
  template <typename CellIterator>
  inline void registerCells(const flag_t, const CellIterator &begin,
                            const CellIterator &end,
                            const BoundaryConfiguration &velocity);

  void unregisterCell(const flag_t, const cell_idx_t, const cell_idx_t,
                      const cell_idx_t) const {}

  inline void treatDirection(const cell_idx_t x, const cell_idx_t y,
                             const cell_idx_t z, const stencil::Direction dir,
                             const cell_idx_t nx, const cell_idx_t ny,
                             const cell_idx_t nz, const flag_t mask);

  inline const Vector3<real_t> &getValue(const cell_idx_t x, cell_idx_t y,
                                         cell_idx_t z) const {
    return vel_->get(x, y, z);
  }

  const typename ForceField::value_type &
  getForce(const cell_idx_t x, cell_idx_t y, cell_idx_t z) const {
    static_assert(StoreForce, "this member function is only available if the "
                              "fourth template argument on the class is true");
    return force_->get(x, y, z);
  }

private:
  const FlagUID uid_;

  PDFField *const pdfField_;
  GhostLayerField<real_t, 3> *force_field_;
  shared_ptr<VelField> vel_;
  shared_ptr<ForceField> force_;
  const LatticeModel_T *const lm_;

}; // class UBB

template <typename LatticeModel_T, typename flag_t,
          bool AdaptVelocityToExternalForce, bool StoreForce>
inline UBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce,
           StoreForce>::Velocity::Velocity(const Config::BlockHandle &config) {
  velocity_[0] = (config && config.isDefined("x"))
                     ? config.getParameter<real_t>("x")
                     : real_c(0.0);
  velocity_[1] = (config && config.isDefined("y"))
                     ? config.getParameter<real_t>("y")
                     : real_c(0.0);
  velocity_[2] = (config && config.isDefined("z"))
                     ? config.getParameter<real_t>("z")
                     : real_c(0.0);
}

template <typename LatticeModel_T, typename flag_t,
          bool AdaptVelocityToExternalForce, bool StoreForce>
inline UBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce,
           StoreForce>::UBB(const BoundaryUID &boundaryUID, const FlagUID &uid,
                            const LatticeModel_T * lm, PDFField *const pdfField,
                            GhostLayerField<real_t, 3> *force_field, FlagField<flag_t> *const flagField)
    :

    Boundary<flag_t>(boundaryUID), uid_(uid), pdfField_(pdfField), force_field_(force_field), lm_(lm) {
  WALBERLA_ASSERT_NOT_NULLPTR(pdfField_);
  if (flagField != nullptr)
    vel_ = make_shared<VelField>(pdfField_->xSize(), pdfField_->ySize(),
                                 pdfField_->zSize(),
                                 flagField->nrOfGhostLayers(), field::zyxf);
  else
    vel_ = make_shared<VelField>(pdfField_->xSize(), pdfField_->ySize(),
                                 pdfField_->zSize(),
                                 pdfField_->nrOfGhostLayers(), field::zyxf);

  if (StoreForce)
    force_ = make_shared<ForceField>(pdfField_->xSize(), pdfField_->ySize(),
                                     pdfField_->zSize(),
                                     pdfField_->nrOfGhostLayers(), field::zyxf);
}

template <typename LatticeModel_T, typename flag_t,
          bool AdaptVelocityToExternalForce, bool StoreForce>
inline void UBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce,
                StoreForce>::beforeBoundaryTreatment() const {
  if (StoreForce)
    force_->setWithGhostLayer(Vector3<real_t>());
}

template <typename LatticeModel_T, typename flag_t,
          bool AdaptVelocityToExternalForce, bool StoreForce>
template <typename Buffer_T>
inline void
UBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce>::packCell(
    Buffer_T &buffer, const cell_idx_t x, const cell_idx_t y,
    const cell_idx_t z) const {
  buffer << vel_->get(x, y, z);
}

template <typename LatticeModel_T, typename flag_t,
          bool AdaptVelocityToExternalForce, bool StoreForce>
template <typename Buffer_T>
inline void UBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce,
                StoreForce>::registerCell(Buffer_T &buffer, const flag_t,
                                          const cell_idx_t x,
                                          const cell_idx_t y,
                                          const cell_idx_t z) {
  buffer >> vel_->get(x, y, z);
}

template <typename LatticeModel_T, typename flag_t,
          bool AdaptVelocityToExternalForce, bool StoreForce>
inline void
UBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce,
    StoreForce>::registerCell(const flag_t, const cell_idx_t x,
                              const cell_idx_t y, const cell_idx_t z,
                              const BoundaryConfiguration &velocity) {
  WALBERLA_ASSERT_EQUAL(dynamic_cast<const Velocity *>(&velocity), &velocity);
  WALBERLA_ASSERT_NOT_NULLPTR(vel_);

  const Velocity &vel = dynamic_cast<const Velocity &>(velocity);

  vel_->get(x, y, z) = vel.velocity();
}

template <typename LatticeModel_T, typename flag_t,
          bool AdaptVelocityToExternalForce, bool StoreForce>
inline void
UBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce,
    StoreForce>::registerCells(const flag_t, const CellInterval &cells,
                               const BoundaryConfiguration &velocity) {
  WALBERLA_ASSERT_EQUAL(dynamic_cast<const Velocity *>(&velocity), &velocity);
  WALBERLA_ASSERT_NOT_NULLPTR(vel_);

  const Velocity &vel = dynamic_cast<const Velocity &>(velocity);

  for (auto cell = vel_->beginSliceXYZ(cells); cell != vel_->end(); ++cell)
    *cell = vel.velocity();
}

template <typename LatticeModel_T, typename flag_t,
          bool AdaptVelocityToExternalForce, bool StoreForce>
template <typename CellIterator>
inline void
UBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce,
    StoreForce>::registerCells(const flag_t, const CellIterator &begin,
                               const CellIterator &end,
                               const BoundaryConfiguration &velocity) {
  WALBERLA_ASSERT_EQUAL(dynamic_cast<const Velocity *>(&velocity), &velocity);
  WALBERLA_ASSERT_NOT_NULLPTR(vel_);

  const Velocity &vel = dynamic_cast<const Velocity &>(velocity);

  for (auto cell = begin; cell != end; ++cell)
    vel_->get(cell->x(), cell->y(), cell->z()) = vel.velocity();
}

template <typename LatticeModel_T, typename flag_t,
          bool AdaptVelocityToExternalForce, bool StoreForce>
#ifndef NDEBUG
inline void
UBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce,
    StoreForce>::treatDirection(const cell_idx_t x, const cell_idx_t y,
                                const cell_idx_t z,
                                const stencil::Direction dir,
                                const cell_idx_t nx, const cell_idx_t ny,
                                const cell_idx_t nz, const flag_t mask)
#else
inline void
UBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce,
    StoreForce>::treatDirection(const cell_idx_t x, const cell_idx_t y,
                                const cell_idx_t z,
                                const stencil::Direction dir,
                                const cell_idx_t nx, const cell_idx_t ny,
                                const cell_idx_t nz, const flag_t /*mask*/)
#endif
{
  WALBERLA_ASSERT_EQUAL(nx, x + cell_idx_c(stencil::cx[dir]));
  WALBERLA_ASSERT_EQUAL(ny, y + cell_idx_c(stencil::cy[dir]));
  WALBERLA_ASSERT_EQUAL(nz, z + cell_idx_c(stencil::cz[dir]));
  WALBERLA_ASSERT_UNEQUAL(mask & this->mask_, numeric_cast<flag_t>(0));
  WALBERLA_ASSERT_EQUAL(
      mask & this->mask_,
      this->mask_); // only true if "this->mask_" only contains one single flag,
                    // which is the case for the current implementation of this
                    // boundary condition (UBB)

  const real_t pdf_old = pdfField_->get(x, y, z, Stencil::idx[dir]);

  if constexpr (LatticeModel_T::compressible) {
    const auto density = Density<LatticeModel_T>::get(*lm_, *pdfField_, x, y, z);
    const auto velocity =
        AdaptVelocityToExternalForce
            ? internal::AdaptVelocityToForce<LatticeModel_T>::get(
                  x, y, z, force_field_, vel_->get(nx, ny, nz), density)
            : vel_->get(nx, ny, nz);

    pdfField_->get(nx, ny, nz, Stencil::invDirIdx(dir)) =
        pdfField_->get(x, y, z, Stencil::idx[dir]) -
        (real_c(6) * density * real_c(LatticeModel_T::w[Stencil::idx[dir]]) *
         (real_c(stencil::cx[dir]) * velocity[0] +
          real_c(stencil::cy[dir]) * velocity[1] +
          real_c(stencil::cz[dir]) * velocity[2]));
  } else {
    const auto velocity =
        AdaptVelocityToExternalForce
            ? internal::AdaptVelocityToForce<LatticeModel_T>::get(
                x, y, z, *lm_, vel_->get(nx, ny, nz), real_t(1))
            : vel_->get(nx, ny, nz);

    pdfField_->get(nx, ny, nz, Stencil::invDirIdx(dir)) =
        pdfField_->get(x, y, z, Stencil::idx[dir]) -
        (real_c(6) * real_c(LatticeModel_T::w[Stencil::idx[dir]]) *
         (real_c(stencil::cx[dir]) * velocity[0] +
          real_c(stencil::cy[dir]) * velocity[1] +
          real_c(stencil::cz[dir]) * velocity[2]));
  }

  if (StoreForce && pdfField_->isInInnerPart(Cell(x, y, z))) {
    const real_t forceMEM =
        pdf_old + pdfField_->get(nx, ny, nz, Stencil::invDirIdx(dir));
    Vector3<real_t> force(real_c(stencil::cx[dir]) * forceMEM,
                          real_c(stencil::cy[dir]) * forceMEM,
                          real_c(stencil::cz[dir]) * forceMEM);
    force_->get(nx, ny, nz) += force;
  }
}

} // namespace espresso
} // namespace lbm
} // namespace walberla
