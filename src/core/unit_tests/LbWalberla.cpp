#define BOOST_TEST_MODULE Walberla node setters and getters test
#define BOOST_TEST_DYN_LINK
#include "config.hpp"

#ifdef LB_WALBERLA

#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <memory>
#include <vector>

#include "boost/mpi.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/LbWalberlaD3Q19MRT.hpp"
#include "grid_based_algorithms/LbWalberla_impl.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_walberla_instance.hpp"
#include "utils/Vector.hpp"
#include <iostream>

using Utils::Vector3d;
using Utils::Vector3i;
using walberla::LbWalberlaD3Q19MRT;

double viscosity = 3;
Vector3d box_dimensions = {10, 12, 14};
double agrid = 0.5;
Vector3i grid_dimensions{int(box_dimensions[0] / agrid),
                         int(box_dimensions[1] / agrid),
                         int(box_dimensions[2] / agrid)};
double tau = 0.34;
double density = 2.5;
Vector3i mpi_shape;


BOOST_AUTO_TEST_CASE(basic_params) {
  LbWalberlaD3Q19MRT lb = LbWalberlaD3Q19MRT(
      viscosity, density, agrid, tau, box_dimensions, mpi_shape, 2);
  BOOST_CHECK(lb.get_grid_dimensions() == grid_dimensions);
  BOOST_CHECK(lb.get_grid_spacing() == agrid);
  BOOST_CHECK(lb.get_tau() == tau);

  BOOST_CHECK_CLOSE(lb.get_viscosity(), viscosity, 1E-11);
  double new_viscosity = 2.0;
  lb.set_viscosity(new_viscosity);
  BOOST_CHECK_CLOSE(lb.get_viscosity(), new_viscosity, 1E-11);
}

BOOST_AUTO_TEST_CASE(boundary) {
  Vector3d vel = {0.2, 3.8, 4.2};
  LbWalberlaD3Q19MRT lb = LbWalberlaD3Q19MRT(
      viscosity, density, agrid, tau, box_dimensions, mpi_shape, 2);
  for (Vector3i node :

       std::vector<Vector3i>{{0, 0, 0}, {0, 1, 2}, {9, 9, 9}}) {
    if (lb.node_in_local_halo(node)) {
      BOOST_CHECK(lb.set_node_velocity_at_boundary(node, vel));
      auto vel_check = lb.get_node_velocity_at_boundary(node);
      // Do we have a value
      BOOST_CHECK(vel_check);
      // Check the value
      BOOST_CHECK_SMALL((*vel_check - vel).norm(), 1E-12);
      BOOST_CHECK(lb.remove_node_from_boundary(node));
      auto res = lb.get_node_is_boundary(node, true);
      // Did we get a value?
      BOOST_CHECK(res);
      // Should not be a boundary node
      BOOST_CHECK(*res == false);
    } else {
      // Not in the local halo.
      BOOST_CHECK(!lb.set_node_velocity_at_boundary(node, vel));
      BOOST_CHECK(!lb.get_node_velocity_at_boundary(node));
      BOOST_CHECK(!lb.remove_node_from_boundary(node));
      BOOST_CHECK(!lb.get_node_is_boundary(node));
    }
  }
}


BOOST_AUTO_TEST_CASE(boundary_flow_shear) {
  Vector3d vel = {0.2, -0.3, 0};
  LbWalberlaD3Q19MRT lb = LbWalberlaD3Q19MRT(
      viscosity, density, agrid, tau, box_dimensions, mpi_shape, 2);
  int lower = 1, upper = 12;
  for (int x = -1; x <= grid_dimensions[0]; x++) {
    for (int y = -1; y <= grid_dimensions[1]; y++) {
      Vector3i node{x, y, lower};
      if (lb.node_in_local_halo(node)) {
        BOOST_CHECK(lb.set_node_velocity_at_boundary(node, Vector3d{}));
        BOOST_CHECK(lb.get_node_is_boundary(node, true));
      }
      node[2] = upper;
      if (lb.node_in_local_halo(node)) {
        BOOST_CHECK(lb.set_node_velocity_at_boundary(node, vel));
      }
    }
  }

  for (int i = 0; i < 200; i++) {
    lb.integrate();
  }
  for (int j = lower + 1; j < upper; j++) {
    auto v = lb.get_node_velocity(Vector3i{0, 0, j});
    if (v) {
      double res = ((*v) - (j - lower) / double(upper - lower) * vel).norm();
      BOOST_CHECK_SMALL(res, 0.05 * vel.norm());
    }
  }
}

std::vector<Vector3i> all_nodes_incl_ghosts(int n_ghost_layers) {
  std::vector<Vector3i> res;
  for (int x = -n_ghost_layers; x < grid_dimensions[0] + n_ghost_layers; x++) {
    for (int y = -n_ghost_layers; y < grid_dimensions[1] + n_ghost_layers;
         y++) {
      for (int z = -n_ghost_layers; z < grid_dimensions[2] + n_ghost_layers;
           z++) {
        res.push_back(Vector3i{x, y, z});
      }
    }
  }
  return res;
}

BOOST_AUTO_TEST_CASE(domain_and_halo) {
  int n_ghost_layers = 2;
  LbWalberlaD3Q19MRT lb =
      LbWalberlaD3Q19MRT(viscosity, density, agrid, tau, box_dimensions,
                                 mpi_shape, n_ghost_layers);

  auto my_left = lb.get_local_domain().first;
  auto my_right = lb.get_local_domain().second;

  for (auto const &n : all_nodes_incl_ghosts(n_ghost_layers)) {
    const Vector3d pos =
        Vector3d{{double(n[0] + .5), double(n[1] + .5), double(n[2] + .5)}};
    int is_local = 0;
    // Nodes in local domain
    if ((n[0] >= my_left[0] and n[1] >= my_left[1] and n[2] >= my_left[2]) and
        (n[0] < my_right[0] and n[1] < my_right[1] and n[2] < my_right[2])) {
      BOOST_CHECK(lb.node_in_local_domain(n));
      BOOST_CHECK(lb.node_in_local_halo(n));

      BOOST_CHECK(lb.pos_in_local_domain(pos));
      BOOST_CHECK(lb.pos_in_local_halo(pos));
      is_local = 1;
    } else {
      // in local halo?
      if ((n[0] + n_ghost_layers >= my_left[0] and
           n[1] + n_ghost_layers >= my_left[1] and
           n[2] + n_ghost_layers >= my_left[2]) and
          (n[0] - n_ghost_layers < my_right[0] and
           n[1] - n_ghost_layers < my_right[1] and
           n[2] - n_ghost_layers < my_right[2])) {
        BOOST_CHECK(!lb.node_in_local_domain(n));
        BOOST_CHECK(lb.node_in_local_halo(n));

        BOOST_CHECK(!lb.pos_in_local_domain(pos));
        BOOST_CHECK(lb.pos_in_local_halo(pos));
      } else {
        // neither in domain nor in halo
        BOOST_CHECK(!lb.node_in_local_domain(n));
        BOOST_CHECK(!lb.node_in_local_halo(n));

        BOOST_CHECK(!lb.pos_in_local_domain(pos));
        BOOST_CHECK(!lb.pos_in_local_halo(pos));
      }
    }

    // If the cell is in the global physical domain
    // check that only one mpi rank said the node was local
    if (n[0] >= 0 and n[1] >= 0 and n[2] >= 0 and n[0] < grid_dimensions[0] and
        n[1] < grid_dimensions[1] and n[2] < grid_dimensions[2]) {
      MPI_Allreduce(MPI_IN_PLACE, &is_local, 1, MPI_INT, MPI_SUM,
                    MPI_COMM_WORLD);
      BOOST_CHECK(is_local == 1);
    }
  }
}

BOOST_AUTO_TEST_CASE(velocity) {
  int n_ghost_layers = 2;
  LbWalberlaD3Q19MRT lb =
      LbWalberlaD3Q19MRT(viscosity, density, agrid, tau, box_dimensions,
                                 mpi_shape, n_ghost_layers);

  // Values
  auto n_pos = [](Vector3i node) {
    return Vector3d{
        {double(node[0] + .5), double(node[1] + .5), double(node[2] + .5)}};
  };

  auto fold = [&](Vector3i n) {
    for (int i = 0; i < 3; i++) {
      if (n[i] < 0)
        n[i] += grid_dimensions[i];
      else if (n[i] >= grid_dimensions[i])
        n[i] -= grid_dimensions[i];
    }
    return n;
  };

  auto n_vel = [&fold](Vector3i node) {
    return fold(node) + Vector3d{{1, 2, -.5}};
  };

  // Assign velocities
  for (auto const &node : all_nodes_incl_ghosts(n_ghost_layers)) {
    if (lb.node_in_local_domain(node)) {
      BOOST_CHECK(lb.set_node_velocity(node, n_vel(node)));
    } else {
      // Check that access to node velocity is not possible
      BOOST_CHECK(!lb.set_node_velocity(node, Vector3d{}));
    }
  }

  lb.ghost_communication();

  // check velocities
  for (auto const &node : all_nodes_incl_ghosts(n_ghost_layers)) {
    double eps = 1E-8;

    if (lb.node_in_local_domain(node)) {
      auto res = lb.get_node_velocity(node);
      BOOST_CHECK(res);                               // value available
      BOOST_CHECK((*res - n_vel(node)).norm() < eps); // correct value?
    }
    if (lb.pos_in_local_halo(n_pos(node))) {
      // Check that the interpolated velocity at the node pos equals the node
      // vel
      auto res = lb.get_velocity_at_pos(n_pos(node), true);
      BOOST_CHECK(res);                                    // locally available
      BOOST_CHECK_SMALL((*res - n_vel(node)).norm(), eps); // value correct?
    } else {
      // Check that access to node velocity is not possible
      BOOST_CHECK(!lb.get_node_velocity(node));
      BOOST_CHECK(!lb.get_velocity_at_pos(n_pos(node), true));
    }
  }
}

BOOST_AUTO_TEST_CASE(total_momentum) {
  LbWalberlaD3Q19MRT lb = LbWalberlaD3Q19MRT(
      viscosity, density, agrid, tau, box_dimensions, mpi_shape, 2);
  auto v = Vector3d{1.5, 2.5, -2.2};
  lb.set_node_velocity(Vector3i{1, 1, 1}, v);
  lb.set_node_velocity(Vector3i{3, 5, 7}, v);
  auto mom = lb.get_momentum();
  auto mom_exp = 2 * density * v;
  MPI_Allreduce(MPI_IN_PLACE, mom.data(), 3, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  BOOST_CHECK_SMALL((mom - mom_exp).norm(), 1E-10);
}

BOOST_AUTO_TEST_CASE(integrate_with_volume_force) {
  LbWalberlaD3Q19MRT lb = LbWalberlaD3Q19MRT(
      viscosity, density, agrid, tau, box_dimensions, mpi_shape, 2);
  auto f = Vector3d{0.01, 0.2, -0.2};
  lb.set_external_force(f);
  BOOST_CHECK_SMALL(lb.get_momentum().norm(), 1E-10);

  int n_nodes = grid_dimensions[0] * grid_dimensions[1] * grid_dimensions[2];
  for (int i = 1; i < 30; i++) {
    lb.integrate();
    auto mom = lb.get_momentum() / double(n_nodes);
    auto mom_exp = (i + .5) * f;
    MPI_Allreduce(MPI_IN_PLACE, mom.data(), 3, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    BOOST_CHECK_SMALL((mom - mom_exp).norm(), 1E-7);

    auto v_exp = mom_exp/density;
    for (int i = 0; i < grid_dimensions[0]; i++)
      for (int j = 0; j < grid_dimensions[1]; j++)
        for (int k = 0; k < grid_dimensions[2]; k++) {
          auto v = lb.get_node_velocity(Vector3i{{i, j, k}});
          if (v) BOOST_CHECK_SMALL(((*v)- v_exp).norm(), 1E-10);
        }

  }
}

BOOST_AUTO_TEST_CASE(integrate_with_point_forces) {
  LbWalberlaD3Q19MRT lb = LbWalberlaD3Q19MRT(
      viscosity, density, agrid, tau, box_dimensions, mpi_shape, 2);
  // auto f = Vector3d{0.15, 0.25, -0.22};
  auto f = Vector3d{0.0006, -0.0013, 0.000528};
  auto f2 = Vector3d{0.095, 0.23, -0.52};
  lb.set_external_force(f);
  lb.add_force_at_pos(Utils::Vector3d{2, 2, 2}, f2);
  BOOST_CHECK_SMALL(lb.get_momentum().norm(), 1E-10);
  lb.integrate();
  auto mom = lb.get_momentum();
  auto mom_exp =
      1.5 * f * grid_dimensions[0] * grid_dimensions[1] * grid_dimensions[2] +
      1.5 * f2;
  MPI_Allreduce(MPI_IN_PLACE, mom.data(), 3, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  BOOST_CHECK_SMALL((mom - mom_exp).norm(), 4E-6);

  for (int i = 1; i < 30; i++) {
    lb.integrate();
    auto mom_exp = (i + 1.5) * (f * grid_dimensions[0] * grid_dimensions[1] *
                                grid_dimensions[2]) +
                   f2;
    auto mom = lb.get_momentum();
    MPI_Allreduce(MPI_IN_PLACE, mom.data(), 3, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    BOOST_CHECK_SMALL((mom - mom_exp).norm(), 2E-4);
  }
}

BOOST_AUTO_TEST_CASE(forces_initial_state) {
  LbWalberlaD3Q19MRT lb = LbWalberlaD3Q19MRT(
      viscosity, density, agrid, tau, box_dimensions, mpi_shape, 2);

  for (Vector3i n : all_nodes_incl_ghosts(1)) {
    if (lb.node_in_local_halo(n)) {
      auto res = lb.get_node_force_to_be_applied(n);
      BOOST_CHECK(res);
      BOOST_CHECK_SMALL((*res).norm(), 1E-10);
      res = lb.get_node_last_applied_force(n, true);
      BOOST_CHECK(res);
      BOOST_CHECK_SMALL((*res).norm(), 1E-10);
    }
  }
}

BOOST_AUTO_TEST_CASE(forces_interpolation) {
  LbWalberlaD3Q19MRT lb = LbWalberlaD3Q19MRT(
      viscosity, density, agrid, tau, box_dimensions, mpi_shape, 2);

  for (Vector3i n : all_nodes_incl_ghosts(1)) {
    if (lb.node_in_local_halo(n)) {
      Vector3d pos{double(n[0]), double(n[1]),
                   double(n[2])}; // Mid point between nodes
      Vector3d f = {1, 2, -3.5};
      lb.add_force_at_pos(pos, f);
      // Check neighboring nodes for force to be applied
      for (int x : {0, 1})
        for (int y : {0, 1})
          for (int z : {0, 1}) {
            Vector3i check_node{{n[0] - x, n[1] - y, n[2] - z}};
            if (lb.node_in_local_halo(check_node)) {
              auto res = lb.get_node_force_to_be_applied(check_node);
              BOOST_CHECK_SMALL(((*res) - f / 8.0).norm(), 1E-10);
            }
          }
      // Apply counter force to clear force field
      lb.add_force_at_pos(pos, -f);
    }
  }
}

BOOST_AUTO_TEST_CASE(force_book_keeping) {
  LbWalberlaD3Q19MRT lb = LbWalberlaD3Q19MRT(
      viscosity, density, agrid, tau, box_dimensions, mpi_shape, 2);

  // Forces added go to force_to_be_applied. After integration, they should be
  // in last_applied_force, where they are used for velocity calculation

  Vector3i origin{};
  Vector3i middle{
      {grid_dimensions[0] / 2, grid_dimensions[1] / 2, grid_dimensions[2] / 2}};
  Vector3i right = grid_dimensions - Vector3i{{1, 1, 1}};

  Vector3d f{{1, -2, 3.1}};

  for (auto n : {origin, middle, right}) {
    // Add force to node position
    if (lb.node_in_local_domain(n)) {
      lb.add_force_at_pos(Vector3d{{n[0] + .5, n[1] + .5, n[2] + .5}}, f);
      BOOST_CHECK_SMALL((*(lb.get_node_force_to_be_applied(n)) - f).norm(),
                        1E-10);
    }
    lb.integrate();
    // Check nodes incl some of the ghosts
    for (auto cn : {n, n + grid_dimensions, n - grid_dimensions,
                    n + Vector3i{{grid_dimensions[0], 0, 0}}}) {
      if (lb.node_in_local_halo(cn)) {
        BOOST_CHECK_SMALL(
            (*(lb.get_node_last_applied_force(cn, true)) - f).norm(), 1E-10);
        BOOST_CHECK_SMALL((*(lb.get_node_force_to_be_applied(cn))).norm(),
                          1E-10);
      }
    }
  }
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  walberla_mpi_init();
  int n_nodes;
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  MPI_Dims_create(n_nodes, 3, mpi_shape.data());

  auto res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

#else // ifdef LB_WALBERLA
int main(int argc, char **argv) {}
#endif
