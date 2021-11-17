#ifndef ESPRESSO_NONE_HPP
#define ESPRESSO_NONE_HPP

#include "LatticeWalberla.hpp.hpp"
#include "PoissonSolver.hpp"

#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"

namespace walberla {
template <typename FloatType = double>
class None : public PoissonSolver<FloatType> {
private:
  using PS = PoissonSolver<FloatType>;
  using typename PS::PotentialField;
  using PoissonSolver<FloatType>::m_potential_field_id;

public:
  explicit None(std::shared_ptr<LatticeWalberla> lattice)
      : PS(std::move(lattice), 0.0) {}

  void reset_charge_field() override {}

  void add_charge_to_field(const BlockDataID &id, FloatType valency) override {}

  [[nodiscard]] BlockDataID get_potential_field_id() override {
    return m_potential_field_id;
  }

  void solve() override {}
};
} // namespace walberla

#endif // ESPRESSO_NONE_HPP
