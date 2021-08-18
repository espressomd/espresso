#ifndef ESPRESSO_FFT_HPP
#define ESPRESSO_FFT_HPP

#include "PoissonSolver.hpp"
#include "WalberlaBlockForest.hpp"

#include "utils/constants.hpp"

#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"

#include "fft/Fft.h"

namespace walberla {
template <typename FloatType = double>
class FFT : public PoissonSolver<FloatType> {
private:
  BlockDataID m_greens_field;
  BlockDataID m_fourier_field;

  using PS = PoissonSolver<FloatType>;
  using PotentialField = typename PS::PotentialField;

  std::unique_ptr<fft::FourierTransform<PotentialField>> m_ft;

public:
  explicit FFT(const WalberlaBlockForest *blockforest)
      : PoissonSolver<FloatType>(blockforest) {
    auto &blocks = PS::get_blockforest()->get_blocks();

    Vector3<uint_t> dim(blocks->getNumberOfXCells(),
                        blocks->getNumberOfYCells(),
                        blocks->getNumberOfZCells());
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
        blocks, PS::m_potential_field_id, greens);
  }

  void solve() override {
    (*m_ft)();
    PS::ghost_communication();
  }
};
} // namespace walberla

#endif // ESPRESSO_FFT_HPP
