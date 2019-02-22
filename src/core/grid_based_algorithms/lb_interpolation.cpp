#include "config.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "lattice.hpp"
#include "utils/Vector.hpp"

#include "lb.hpp"
#include "lb_interface.hpp"
#include "lbgpu.hpp"

#if defined(LB) || defined(LB_GPU)

namespace {

template <typename Op>
void lattice_interpolation(Lattice const &lattice, Vector3d const &pos,
                           Op &&op) {
  Vector<8, std::size_t> node_index{};
  Vector6d delta{};

  /* determine elementary lattice cell surrounding the particle
     and the relative position of the particle in this cell */
  lattice.map_position_to_lattice(pos, node_index, delta);
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

Vector3d node_u(Lattice::index_t index) {
#ifdef LB_BOUNDARIES
  if (lbfields[index].boundary) {
    return lbfields[index].slip_velocity;
  }
#endif // LB_BOUNDARIES
  auto const modes = lb_calc_modes(index);
  auto const local_rho = lbpar.rho + modes[0];
  return Vector3d{modes[1], modes[2], modes[3]} / local_rho * lb_lbfluid_get_agrid() / lb_lbfluid_get_tau();
}

} // namespace

const Vector3d lb_lbinterpolation_get_interpolated_velocity(const Vector3d &pos) {
  if (lattice_switch & LATTICE_LB) {
#ifdef LB
    Vector3d interpolated_u{};

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

const Vector3d lb_lbinterpolation_get_interpolated_velocity_global(const Vector3d &pos) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    Vector3d interpolated_u{};
    lb_get_interpolated_velocity_gpu(pos.data(), interpolated_u.data(), 1);
    return interpolated_u;
#endif
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    auto const node = map_position_node_array(pos);
    if (node == 0) {
      return lb_lbinterpolation_get_interpolated_velocity(pos);
    } else {
      return mpi_recv_lb_interpolated_velocity(node, pos);
    }
  }
#endif
  return {};
}

#ifdef LB
void lb_lbinterpolation_add_force_density(const Vector3d &pos,
                                  const Vector3d &force_density) {
  lattice_interpolation(lblattice, pos,
                        [&force_density](Lattice::index_t index, double w) {
                            auto &field = lbfields[index];
                            field.force_density += w * force_density;
                        });
}
#endif

#endif
