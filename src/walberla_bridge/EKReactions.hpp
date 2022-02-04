#ifndef ESPRESSO_EKREACTIONS_HPP
#define ESPRESSO_EKREACTIONS_HPP

#include "LatticeWalberla.hpp"

#include <memory>

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

#endif // ESPRESSO_EKREACTIONS_HPP
