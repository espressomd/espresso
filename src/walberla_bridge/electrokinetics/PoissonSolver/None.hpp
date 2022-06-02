#ifndef ESPRESSO_NONE_HPP
#define ESPRESSO_NONE_HPP

#include "LatticeWalberla.hpp"
#include "PoissonSolver.hpp"

namespace walberla {
template <typename FloatType> class None : public PoissonSolver {
private:
  BlockDataID m_potential_field_id;

  using PotentialField = GhostLayerField<FloatType, 1>;

public:
  explicit None(std::shared_ptr<LatticeWalberla> lattice)
      : PoissonSolver(std::move(lattice), 0.0) {
    m_potential_field_id = field::addToStorage<PotentialField>(
        get_lattice()->get_blocks(), "potential field", 0.0, field::fzyx,
        get_lattice()->get_ghost_layers());
  }

  void reset_charge_field() override {}
  void add_charge_to_field(const BlockDataID &, double, bool) override {}

  [[nodiscard]] domain_decomposition::BlockDataID get_potential_field_id() const
      noexcept override {
    return m_potential_field_id;
  }

  void solve() override {}
  void ghost_communication() override {}
};
} // namespace walberla

#endif // ESPRESSO_NONE_HPP
