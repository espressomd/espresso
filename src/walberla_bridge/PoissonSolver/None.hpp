#ifndef ESPRESSO_NONE_HPP
#define ESPRESSO_NONE_HPP

#include "LatticeWalberla.hpp"
#include "PoissonSolver.hpp"

namespace walberla {
template <typename FloatType = double>
class None : public PoissonSolver<FloatType> {
private:
  using PS = PoissonSolver<FloatType>;

public:
  explicit None(std::shared_ptr<LatticeWalberla> lattice)
      : PS(std::move(lattice), 0.0) {}

  void reset_charge_field() override {}

  void add_charge_to_field(const std::size_t &id, FloatType valency) override {}

  void solve() override {}
};
} // namespace walberla

#endif // ESPRESSO_NONE_HPP
