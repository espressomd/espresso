#ifndef ESPRESSO_EKREACTIONS_HPP
#define ESPRESSO_EKREACTIONS_HPP

#include "EKinWalberlaBase.hpp"
#include "LatticeWalberla.hpp"

#include <memory>

namespace walberla {

void EKReactions(std::shared_ptr<LatticeWalberla> &lattice, double coeff,
                 std::size_t density_field_id, double order) {
  //    auto blocks = lattice->get_blocks();
  //
  //    const auto density_id = BlockDataID(density_field_id);

  // TODO: generate a kernel for this
  //    auto kernel = pystencils::DiffusiveFluxKernelWithElectrostatic(
  //            get_diffusion(), m_flux_field_flattened_id,
  //            BlockDataID(potential_id), m_density_field_flattened_id,
  //            ext_field[0], ext_field[1], ext_field[2], get_kT(),
  //            get_valency());
  //    for (auto &block : *blocks) {
  //      kernel.run(&block);
  //    }
}

template <typename FloatType> class EKReactant {
private:
  std::shared_ptr<EKinWalberlaBase<FloatType>> m_ekspecies;
  FloatType m_stoech_coeff;
  FloatType m_order;

public:
  EKReactant(std::shared_ptr<EKinWalberlaBase<FloatType>> ekspecies,
             FloatType stoech_coeff, FloatType order)
      : m_ekspecies(std::move(ekspecies)), m_stoech_coeff(stoech_coeff),
        m_order(order) {}

  void set_stoech_coefficient(FloatType stoech_coeff) {
    m_stoech_coeff = stoech_coeff;
  }
  [[nodiscard]] FloatType get_stoech_coeff() const { return m_stoech_coeff; }

  void set_order(FloatType order) { m_order = order; }
  [[nodiscard]] FloatType get_order() const { return m_order; }

  void set_species(std::shared_ptr<EKinWalberlaBase<FloatType>> ekspecies) {
    m_ekspecies = std::move(ekspecies);
  }
  [[nodiscard]] auto get_species() const { return m_ekspecies; }
};

template <typename FloatType> class EKReaction {
private:
  std::vector<EKReactant<FloatType>> m_reactants;
  FloatType m_coefficient;
};
} // namespace walberla
#endif // ESPRESSO_EKREACTIONS_HPP
