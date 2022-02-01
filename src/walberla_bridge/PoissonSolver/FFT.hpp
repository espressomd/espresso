#ifndef ESPRESSO_FFT_HPP
#define ESPRESSO_FFT_HPP

#include "PoissonSolver.hpp"

#include <memory>

namespace walberla::blockforest {
// forward declare
class StructuredBlockForest;
} // namespace walberla::blockforest

namespace walberla::fft {
// forward declare
template <typename Field_T> class FourierTransform;
} // namespace walberla::fft

namespace walberla {
template <typename FloatType = double>
class FFT : public PoissonSolver<FloatType> {
private:
  using PS = PoissonSolver<FloatType>;
  using PS::ghost_communication;
  using PS::m_lattice;
  using PS::m_potential_field_id;
  using typename PS::ChargeField;
  using typename PS::PotentialField;

  std::unique_ptr<fft::FourierTransform<PotentialField>> m_ft;

  // TODO: figure out how to pass this to the fourier-transformation without
  //       taking partial-ownership of the BlockForest...
  std::shared_ptr<blockforest::StructuredBlockForest> m_blocks;

public:
  FFT(std::shared_ptr<LatticeWalberla> lattice, FloatType permittivity);

  using PS::get_permittivity;
  using PS::set_permittivity;

  void reset_charge_field() override;
  void add_charge_to_field(const BlockDataID &id, FloatType valency) override;

  [[nodiscard]] BlockDataID get_potential_field_id();

  void solve();
};
} // namespace walberla

#endif // ESPRESSO_FFT_HPP
