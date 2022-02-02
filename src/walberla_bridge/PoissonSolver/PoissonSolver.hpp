#ifndef ESPRESSO_POISSONSOLVER_HPP
#define ESPRESSO_POISSONSOLVER_HPP

#include <utility>

#include "LatticeWalberla.hpp"

#include "blockforest/communication/UniformBufferedScheme.h"
#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"
#include "field/communication/PackInfo.h"
#include "stencil/D3Q27.h"

namespace walberla {
template <typename FloatType = double> class PoissonSolver {
protected:
  std::shared_ptr<LatticeWalberla> m_lattice;

private:
  FloatType m_permittivity;
  BlockDataID m_potential_field_id;

protected:
  using PotentialField = GhostLayerField<FloatType, 1>;
  using ChargeField = GhostLayerField<FloatType, 1>;

  using FullCommunicator = blockforest::communication::UniformBufferedScheme<
      typename stencil::D3Q27>;
  std::shared_ptr<FullCommunicator> m_full_communication;

public:
  PoissonSolver(std::shared_ptr<LatticeWalberla> lattice,
                FloatType permittivity)
      : m_lattice{std::move(lattice)}, m_permittivity{permittivity} {
    m_potential_field_id = field::addToStorage<PotentialField>(
        m_lattice->get_blocks(), "potential field", 0.0, field::fzyx,
        m_lattice->get_ghost_layers());

    m_full_communication =
        std::make_shared<FullCommunicator>(m_lattice->get_blocks());
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<PotentialField>>(
            m_potential_field_id));
  }

  virtual void reset_charge_field() = 0;
  virtual void add_charge_to_field(const std::size_t &id,
                                   FloatType valency) = 0;
  [[nodiscard]] std::size_t get_potential_field_id() {
    return m_potential_field_id;
  }

  void set_permittivity(FloatType permittivity) {
    m_permittivity = permittivity;
  }
  [[nodiscard]] FloatType get_permittivity() const { return m_permittivity; }

  virtual void solve() = 0;
  void ghost_communication() { (*m_full_communication)(); };
};
} // namespace walberla

#endif // ESPRESSO_POISSONSOLVER_HPP
