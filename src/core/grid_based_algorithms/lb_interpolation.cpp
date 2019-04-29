#include <boost/mpi/collectives.hpp>

#include "communication.hpp"
#include "config.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lattice.hpp"
#include "grid_based_algorithms/lb_interpolation.hpp"
#include <utils/Vector.hpp>

#include "lb.hpp"
#include "lb_interface.hpp"
#include "lbgpu.hpp"

namespace {

InterpolationOrder interpolation_order = InterpolationOrder::linear;
}

void mpi_set_interpolation_order_slave(int, int) {
  boost::mpi::broadcast(comm_cart, interpolation_order, 0);
}

#if defined(LB) || defined(LB_GPU)

void lb_lbinterpolation_set_interpolation_order(
    InterpolationOrder const &order) {
  interpolation_order = order;
  mpi_call(mpi_set_interpolation_order_slave, 0, 0);
  boost::mpi::broadcast(comm_cart, interpolation_order, 0);
}

InterpolationOrder lb_lbinterpolation_get_interpolation_order() {
  return interpolation_order;
}

namespace {
template <typename Op>
void lattice_interpolation(Lattice const &lattice, Utils::Vector3d const &pos,
                           Op &&op) {
  Utils::Vector<std::size_t, 8> node_index{};
  Utils::Vector6d delta{};

  /* determine elementary lattice cell surrounding the particle
     and the relative position of the particle in this cell */
  lattice.map_position_to_lattice(pos, node_index, delta, my_left, local_box_l);
  for (int z = 0; z < 2; z++) {
    for (int y = 0; y < 2; y++) {
      for (int x = 0; x < 2; x++) {
        auto &index = node_index[(z * 2 + y) * 2 + x];
        auto const w = delta[3 * x + 0] * delta[3 * y + 1] * delta[3 * z + 2];

        op(index, w);
      }
    }
  }
}

Utils::Vector3d node_u(Lattice::index_t index) {
#ifdef LB_BOUNDARIES
  if (lbfields[index].boundary) {
    return lbfields[index].slip_velocity;
  }
#endif // LB_BOUNDARIES
  auto const modes = lb_calc_modes(index);
  auto const local_rho = lbpar.rho + modes[0];
  return Utils::Vector3d{modes[1], modes[2], modes[3]} / local_rho;
}

} // namespace

const Utils::Vector3d
lb_lbinterpolation_get_interpolated_velocity(const Utils::Vector3d &pos) {
  if (lattice_switch == ActiveLB::CPU) {
#ifdef LB
    Utils::Vector3d interpolated_u{};

    /* calculate fluid velocity at particle's position
       this is done by linear interpolation
       (Eq. (11) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
    lattice_interpolation(lblattice, pos,
                          [&interpolated_u](Lattice::index_t index, double w) {
                            interpolated_u += w * node_u(index);
                          });

    return interpolated_u;
  }
#endif
  return {};
}

const Utils::Vector3d lb_lbinterpolation_get_interpolated_velocity_global(
    const Utils::Vector3d &pos) {
  auto const folded_pos = folded_position(pos);
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    Utils::Vector3d interpolated_u{};
    switch (interpolation_order) {
    case (InterpolationOrder::linear):
      lb_get_interpolated_velocity_gpu<8>(folded_pos.data(),
                                          interpolated_u.data(), 1);
      break;
    case (InterpolationOrder::quadratic):
      lb_get_interpolated_velocity_gpu<27>(folded_pos.data(),
                                           interpolated_u.data(), 1);
      break;
    }
    return interpolated_u;
#endif
  }
  if (lattice_switch == ActiveLB::CPU) {
#ifdef LB
    switch (interpolation_order) {
    case (InterpolationOrder::quadratic):
      throw std::runtime_error("The non-linear interpolation scheme is not "
                               "implemented for the CPU LB.");
    case (InterpolationOrder::linear):
      auto const node = map_position_node_array(folded_pos);
      if (node == 0) {
        return lb_lbinterpolation_get_interpolated_velocity(folded_pos);
      }
      return mpi_recv_lb_interpolated_velocity(node, folded_pos);
    }
  }
#endif
  return {};
}

#ifdef LB
void lb_lbinterpolation_add_force_density(
    const Utils::Vector3d &pos, const Utils::Vector3d &force_density) {
  switch (interpolation_order) {
  case (InterpolationOrder::quadratic):
    throw std::runtime_error("The non-linear interpolation scheme is not "
                             "implemented for the CPU LB.");
  case (InterpolationOrder::linear):
    lattice_interpolation(lblattice, pos,
                          [&force_density](Lattice::index_t index, double w) {
                            auto &field = lbfields[index];
                            field.force_density += w * force_density;
                          });
    break;
  }
}
#endif

#endif
