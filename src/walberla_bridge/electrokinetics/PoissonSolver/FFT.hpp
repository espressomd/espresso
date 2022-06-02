#ifndef ESPRESSO_FFT_HPP
#define ESPRESSO_FFT_HPP

#include "PoissonSolver.hpp"

#include <blockforest/communication/UniformBufferedScheme.h>
#include <fft/Fft.h>
#include <field/AddToStorage.h>
#include <field/GhostLayerField.h>
#include <field/communication/PackInfo.h>
#include <stencil/D3Q27.h>

#include <utils/constants.hpp>

#include <memory>

namespace walberla {
template <typename FloatType> class FFT : public PoissonSolver {
private:
  template <typename T> inline FloatType FloatType_c(T t) {
    return numeric_cast<FloatType>(t);
  }

  BlockDataID m_potential_field_id;

  using PoissonSolver::get_lattice;
  using PoissonSolver::ghost_communication;
  using PotentialField = GhostLayerField<FloatType, 1>;

  std::shared_ptr<fft::FourierTransform<PotentialField>> m_ft;
  std::shared_ptr<blockforest::StructuredBlockForest> m_blocks;

  using FullCommunicator = blockforest::communication::UniformBufferedScheme<
      typename stencil::D3Q27>;
  std::shared_ptr<FullCommunicator> m_full_communication;

public:
  FFT(std::shared_ptr<LatticeWalberla> lattice, double permittivity)
      : PoissonSolver(std::move(lattice), permittivity) {
    m_blocks = get_lattice()->get_blocks();

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

    m_potential_field_id = field::addToStorage<PotentialField>(
        get_lattice()->get_blocks(), "potential field", 0.0, field::fzyx,
        get_lattice()->get_ghost_layers());

    m_ft = std::make_shared<fft::FourierTransform<PotentialField>>(
        m_blocks, m_potential_field_id, greens);

    m_full_communication =
        std::make_shared<FullCommunicator>(get_lattice()->get_blocks());
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<PotentialField>>(
            m_potential_field_id));
  }

  using PoissonSolver::get_permittivity;
  using PoissonSolver::set_permittivity;

  void reset_charge_field() override {
    // the FFT-solver re-uses the potential field for the charge
    const auto potential_id = get_potential_field_id();

    for (auto &block : *get_lattice()->get_blocks()) {
      auto field = block.template getData<PotentialField>(potential_id);
      WALBERLA_FOR_ALL_CELLS_XYZ(field, field->get(x, y, z) = 0.;)
    }
  }

  void add_charge_to_field(const BlockDataID &id, double valency,
                           bool is_double_precision) override {
    auto const factor = FloatType_c(valency) / FloatType_c(get_permittivity());
    // the FFT-solver re-uses the potential field for the charge
    const auto charge_id = get_potential_field_id();
    const auto &density_id = id;
    for (auto &block : *get_lattice()->get_blocks()) {
      auto charge_field = block.template getData<PotentialField>(charge_id);
      if (is_double_precision) {
        auto density_field =
            block.template getData<walberla::GhostLayerField<double, 1>>(
                density_id);
        WALBERLA_FOR_ALL_CELLS_XYZ(
            charge_field, charge_field->get(x, y, z) +=
                          factor * FloatType_c(density_field->get(x, y, z));)
      } else {
        auto density_field =
            block.template getData<walberla::GhostLayerField<float, 1>>(
                density_id);
        WALBERLA_FOR_ALL_CELLS_XYZ(
            charge_field, charge_field->get(x, y, z) +=
                          factor * FloatType_c(density_field->get(x, y, z));)
      }
    }
  }

  [[nodiscard]] BlockDataID get_potential_field_id() const noexcept override {
    return m_potential_field_id;
  }

  void solve() override {
    (*m_ft)();
    ghost_communication();
  }
  void ghost_communication() override { (*m_full_communication)(); }
};
} // namespace walberla

#endif // ESPRESSO_FFT_HPP
