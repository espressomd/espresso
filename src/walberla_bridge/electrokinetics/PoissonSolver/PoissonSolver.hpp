#ifndef ESPRESSO_POISSONSOLVER_HPP
#define ESPRESSO_POISSONSOLVER_HPP

#include <memory>

#include "LatticeWalberla.hpp"

// TODO: what about this here?
#include <domain_decomposition/BlockDataID.h>

namespace walberla {
class PoissonSolver {
private:
  std::shared_ptr<LatticeWalberla> m_lattice;
  double m_permittivity;

public:
  PoissonSolver(std::shared_ptr<LatticeWalberla> lattice, double permittivity)
      : m_lattice(std::move(lattice)), m_permittivity(permittivity){};

  virtual void reset_charge_field() = 0;

  virtual void add_charge_to_field(const domain_decomposition::BlockDataID &id,
                                   double valency,
                                   bool is_double_precision) = 0;
  [[nodiscard]] virtual domain_decomposition::BlockDataID
  get_potential_field_id() const noexcept = 0;

  void set_permittivity(double permittivity) noexcept {
    m_permittivity = permittivity;
  }
  [[nodiscard]] double get_permittivity() const noexcept {
    return m_permittivity;
  }

  [[nodiscard]] auto get_lattice() const noexcept { return m_lattice; }

  virtual void solve() = 0;
  virtual void ghost_communication() = 0;
};
} // namespace walberla

#endif // ESPRESSO_POISSONSOLVER_HPP
