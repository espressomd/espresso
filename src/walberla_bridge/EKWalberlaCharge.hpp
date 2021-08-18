#ifndef ESPRESSO_EKWALBERLACHARGE_HPP
#define ESPRESSO_EKWALBERLACHARGE_HPP

#include "PoissonSolver/PoissonSolver.hpp"
#include "WalberlaBlockForest.hpp"

#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"

namespace walberla {
template <typename FloatType = double> class EKWalberlaCharge {
private:
  const WalberlaBlockForest *m_blockforest;
  const PoissonSolver<FloatType> *m_poissonsolver;

  // Type definitions
  // in principle the charge field does not need ghost-layers. This causes some
  // trouble when traversing field with and without ghost layers.
  using ChargeField = GhostLayerField<FloatType, 1>;

  // Block data access handles
  BlockDataID m_charge_field_id;
  BlockDataID m_charge_field_flattened_id;

public:
  EKWalberlaCharge(const WalberlaBlockForest *blockforest,
                   const PoissonSolver<FloatType> *poissonsolver)
      : m_blockforest{blockforest}, m_poissonsolver{poissonsolver} {
    m_charge_field_id = field::addToStorage<ChargeField>(
        get_blockforest()->get_blocks(), "charge field", 0.0, field::fzyx,
        get_blockforest()->get_ghost_layers());
    m_charge_field_flattened_id =
        field::addFlattenedShallowCopyToStorage<ChargeField>(
            get_blockforest()->get_blocks(), m_charge_field_id,
            "flattened charge field");
  }

  [[nodiscard]] const walberla::WalberlaBlockForest *get_blockforest() const {
    return m_blockforest;
  };

  //  [[nodiscard]] BlockDataID get_potential_id() const {
  //    return m_potential_field_id;
  //  }

  void reset_charge_field() const {
    for (auto &block : *get_blockforest()->get_blocks()) {
      auto field = block.template getData<ChargeField>(m_charge_field_id);
      WALBERLA_FOR_ALL_CELLS_XYZ(field, field->get(x, y, z) = 0.;)
    }
  }

  void add_charge_to_field(const BlockDataID &id, FloatType valency) {
    for (auto &block : *get_blockforest()->get_blocks()) {
      auto charge_field =
          block.template getData<ChargeField>(m_charge_field_id);
      // TODO: how do i know which field type is used?
      auto density_field = block.template getData<ChargeField>(id);
      WALBERLA_FOR_ALL_CELLS_XYZ(charge_field,
                                 charge_field->get(x, y, z) +=
                                 valency * density_field->get(x, y, z);)
    }
  }

  void solve_potential() const { m_poissonsolver->solve(); }

  // - methods for solving potential using charge
  //   -> maybe a own sub-class for that?
  // - boundary handling?
  // - needs to be able to calculate the charge using a kernel?
  //   -> what about multiple species?
};

} // namespace walberla
#endif // ESPRESSO_EKWALBERLACHARGE_HPP
