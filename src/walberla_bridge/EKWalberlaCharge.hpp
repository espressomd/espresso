#ifndef ESPRESSO_EKWALBERLACHARGE_HPP
#define ESPRESSO_EKWALBERLACHARGE_HPP

#include "WalberlaBlockForest.hpp"

#include "blockforest/communication/UniformBufferedScheme.h"
#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"
#include "field/communication/PackInfo.h"
#include "stencil/D3Q27.h"
#include "timeloop/SweepTimeloop.h"

namespace walberla {
template <typename FloatType = double> class EKWalberlaCharge {
private:
  const WalberlaBlockForest *m_blockforest;

  // Type definitions
  // in principle the charge field does not need ghost-layers. This causes some
  // trouble when traversing field with and without ghost layers.
  using ChargeField = GhostLayerField<FloatType, 1>;
  using PotentialField = GhostLayerField<FloatType, 1>;

  // Block data access handles
  BlockDataID m_charge_field_id;
  BlockDataID m_charge_field_flattened_id;
  BlockDataID m_potential_field_id;
  // TODO: check that this is necessary
  BlockDataID m_potential_field_flattened_id;

  // TODO: check that this is the right approach
  std::shared_ptr<timeloop::SweepTimeloop> m_time_loop;

  using FullCommunicator = blockforest::communication::UniformBufferedScheme<
      typename stencil::D3Q27>;
  std::shared_ptr<FullCommunicator> m_full_communication;

public:
  EKWalberlaCharge(const WalberlaBlockForest *blockforest)
      : m_blockforest{blockforest} {
    m_charge_field_id = field::addToStorage<ChargeField>(
        get_blockforest()->get_blocks(), "charge field", 0.0, field::fzyx,
        get_blockforest()->get_ghost_layers());
    m_charge_field_flattened_id =
        field::addFlattenedShallowCopyToStorage<ChargeField>(
            get_blockforest()->get_blocks(), m_charge_field_id,
            "flattened charge field");
    m_potential_field_id = field::addToStorage<PotentialField>(
        get_blockforest()->get_blocks(), "potential field", 0.0, field::fzyx,
        get_blockforest()->get_ghost_layers());
    m_potential_field_flattened_id =
        field::addFlattenedShallowCopyToStorage<PotentialField>(
            get_blockforest()->get_blocks(), m_potential_field_id,
            "flattened potential field");

    m_full_communication =
        std::make_shared<FullCommunicator>(get_blockforest()->get_blocks());
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<PotentialField>>(
            m_potential_field_id));
  }

  [[nodiscard]] const walberla::WalberlaBlockForest *get_blockforest() const {
    return m_blockforest;
  };

  [[nodiscard]] BlockDataID get_potential_id() const {
    return m_potential_field_id;
  }

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

  void solve_potential() const {
    // TODO: call poisson-solver to solve the potential
    // TODO: make sure to sync the field!
  }

  // - methods for solving potential using charge
  //   -> maybe a own sub-class for that?
  // - boundary handling?
  // - needs to be able to calculate the charge using a kernel?
  //   -> what about multiple species?
};

} // namespace walberla
#endif // ESPRESSO_EKWALBERLACHARGE_HPP
