#include "config.hpp"
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
template <typename Op>
void new_lattice_interpolation(Vector3i const &ind, Vector6d const &delta,
                           Op &&op) {
  Vector3i neighbor_ind{};
  for (int z = 0; z < 2; z++) {
    for (int y = 0; y < 2; y++) {
      for (int x = 0; x < 2; x++) {
        auto const w = delta[3 * x + 0] * delta[3 * y + 1] * delta[3 * z + 2];
        neighbor_ind = ind + Vector3i{{x, y, z}};
        for (int i=0; i<3; ++i) {
          if (neighbor_ind[i] == box_l[i] / lb_lbfluid_get_agrid()) {
            neighbor_ind[i] = 0;
          }
        }
        op(neighbor_ind, w);
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

const Vector3d lb_lbfluid_get_interpolated_velocity(const Vector3d &pos) {
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

const Vector3d lb_lbfluid_get_interpolated_velocity_global(const Vector3d &pos) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    Vector3d interpolated_u{};
    lb_get_interpolated_velocity_gpu(pos.data(), interpolated_u.data(), 1);
    return interpolated_u;
#endif
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    /* calculate fluid velocity at particle's position
       this is done by linear interpolation
       (Eq. (11) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
    Vector3i ind{};
    Vector6d delta{};
    Lattice::map_position_to_lattice_global(pos, ind, delta, lb_lbfluid_get_agrid());
    Vector3d interpolated_u{};
    new_lattice_interpolation(ind, delta,
                          [&interpolated_u](Vector3i const& ind, double w) {
                            interpolated_u += w * lb_lbnode_get_u(ind);
                          });

    return interpolated_u;
  }
#endif
  return {};
}

#ifdef LB
void lb_lbfluid_add_force_density(const Vector3d &pos,
                                  const Vector3d &force_density) {
  lattice_interpolation(lblattice, pos,
                        [&force_density](Lattice::index_t index, double w) {
                          auto &node = lbfields[index];

                          node.force_density[0] += w * force_density[0];
                          node.force_density[1] += w * force_density[1];
                          node.force_density[2] += w * force_density[2];
                        });
}
#endif

#endif
