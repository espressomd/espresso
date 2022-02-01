#include "FFT.hpp"

#include "LatticeWalberla.hpp"
#include "PoissonSolver.hpp"

#include "utils/constants.hpp"

#include "field/GhostLayerField.h"

#include "blockforest/StructuredBlockForest.h"
#include "fft/Fft.h"

namespace walberla {

template <typename FloatType>
FFT<FloatType>::FFT(std::shared_ptr<LatticeWalberla> lattice,
                    FloatType permittivity)
    : PS(std::move(lattice), permittivity) {
  m_blocks = m_lattice->get_blocks();

  Vector3<uint_t> dim(m_blocks->getNumberOfXCells(),
                      m_blocks->getNumberOfYCells(),
                      m_blocks->getNumberOfZCells());
  const auto greens = [dim](uint_t x, uint_t y, uint_t z) -> real_t {
    if (x == 0 && y == 0 && z == 0)
      return 0;
    return -0.5 /
           (std::cos(2 * Utils::pi() * real_c(x) / real_c(dim[0])) +
            std::cos(2 * Utils::pi() * real_c(y) / real_c(dim[1])) +
            std::cos(2 * Utils::pi() * real_c(z) / real_c(dim[2])) - 3) /
           real_c(dim[0] * dim[1] * dim[2]);
  };

  m_ft = std::make_unique<fft::FourierTransform<PotentialField>>(
      m_blocks, m_potential_field_id, greens);
}

template <typename FloatType> void FFT<FloatType>::reset_charge_field() {
  // the FFT-solver re-uses the potential field for the charge
  for (auto &block : *m_lattice->get_blocks()) {
    auto field = block.template getData<PotentialField>(m_potential_field_id);
    WALBERLA_FOR_ALL_CELLS_XYZ(field, field->get(x, y, z) = 0.;)
  }
}

template <typename FloatType>
void FFT<FloatType>::add_charge_to_field(const BlockDataID &id,
                                         FloatType valency) {
  auto const factor = valency / get_permittivity();
  // the FFT-solver re-uses the potential field for the charge
  for (auto &block : *m_lattice->get_blocks()) {
    auto charge_field =
        block.template getData<PotentialField>(m_potential_field_id);
    auto density_field = block.template getData<ChargeField>(id);
    WALBERLA_FOR_ALL_CELLS_XYZ(charge_field,
                               charge_field->get(x, y, z) +=
                               factor * density_field->get(x, y, z);)
  }
}

template <typename FloatType>
BlockDataID FFT<FloatType>::get_potential_field_id() {
  return PS::m_potential_field_id;
}

template <typename FloatType> void FFT<FloatType>::solve() {
  (*m_ft)();
  ghost_communication();
}

// explicit template instatiation
template class FFT<float>;
template class FFT<double>;
} // namespace walberla
