#ifndef ESPRESSO_POISSONSOLVER_HPP
#define ESPRESSO_POISSONSOLVER_HPP

#include <utility>

#include "WalberlaBlockForest.hpp"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"
#include "field/communication/PackInfo.h"
#include "stencil/D3Q27.h"

namespace walberla {
template <typename FloatType = double> class PoissonSolver {
private:
  const std::shared_ptr<WalberlaBlockForest> m_blockforest;
  FloatType m_permittivity;

protected:
  BlockDataID m_potential_field_id;
  // TODO: check that this is necessary
  BlockDataID m_potential_field_flattened_id;

  using PotentialField = GhostLayerField<FloatType, 1>;
  using ChargeField = GhostLayerField<FloatType, 1>;

  using FullCommunicator = blockforest::communication::UniformBufferedScheme<
      typename stencil::D3Q27>;
  std::shared_ptr<FullCommunicator> m_full_communication;

public:
  PoissonSolver(std::shared_ptr<WalberlaBlockForest> blockforest,
                FloatType permittivity)
      : m_blockforest{std::move(blockforest)}, m_permittivity{permittivity} {
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

  virtual void reset_charge_field() = 0;
  virtual void add_charge_to_field(const BlockDataID &id,
                                   FloatType valency) = 0;
  virtual BlockDataID get_potential_field_id() = 0;

  [[nodiscard]] const std::shared_ptr<WalberlaBlockForest> &
  get_blockforest() const {
    return m_blockforest;
  };
  void set_permittivity(FloatType permittivity) {
    m_permittivity = permittivity;
  }
  [[nodiscard]] FloatType get_permittivity() const { return m_permittivity; }

  virtual void solve() = 0;
  void ghost_communication() { (*m_full_communication)(); };
};
} // namespace walberla

#endif // ESPRESSO_POISSONSOLVER_HPP
