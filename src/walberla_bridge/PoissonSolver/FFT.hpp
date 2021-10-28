#ifndef ESPRESSO_FFT_HPP
#define ESPRESSO_FFT_HPP

#include "PoissonSolver.hpp"
#include "WalberlaBlockForest.hpp"

#include "utils/constants.hpp"

#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"

#include "blockforest/StructuredBlockForest.h"
#include "fft/Fft.h"

namespace walberla {
template <typename FloatType = double>
class FFT : public PoissonSolver<FloatType> {
private:
  using PS = PoissonSolver<FloatType>;
  using PS::get_blockforest;
  using PS::ghost_communication;
  using PS::m_potential_field_id;
  using typename PS::ChargeField;
  using typename PS::PotentialField;

  std::unique_ptr<fft::FourierTransform<PotentialField>> m_ft;

  // TODO: figure out how to pass this to the fourier-transformation without
  //       taking partial-ownership of the BlockForest...
  std::shared_ptr<blockforest::StructuredBlockForest> m_blocks;

public:
  FFT(std::shared_ptr<WalberlaBlockForest> blockforest, FloatType permittivity)
      : PS(std::move(blockforest), permittivity) {
    m_blocks = get_blockforest()->get_blocks();

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

  using PS::get_permittivity;
  using PS::set_permittivity;

  void reset_charge_field() override {
    // the FFT-solver re-uses the potential field for the charge
    for (auto &block : *get_blockforest()->get_blocks()) {
      auto field = block.template getData<PotentialField>(m_potential_field_id);
      WALBERLA_FOR_ALL_CELLS_XYZ(field, field->get(x, y, z) = 0.;)
    }
  }

  void add_charge_to_field(const BlockDataID &id, FloatType valency) override {
    auto const factor = valency / get_permittivity();
    // the FFT-solver re-uses the potential field for the charge
    for (auto &block : *get_blockforest()->get_blocks()) {
      auto charge_field =
          block.template getData<PotentialField>(m_potential_field_id);
      auto density_field = block.template getData<ChargeField>(id);
      WALBERLA_FOR_ALL_CELLS_XYZ(charge_field,
                                 charge_field->get(x, y, z) +=
                                 factor * density_field->get(x, y, z);)
    }
  }

  [[nodiscard]] BlockDataID get_potential_field_id() override {
    return PS::m_potential_field_id;
  }

  void solve() override {
    (*m_ft)();
    ghost_communication();
  }
};
} // namespace walberla

#endif // ESPRESSO_FFT_HPP
