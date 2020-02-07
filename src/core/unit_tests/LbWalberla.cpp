#define BOOST_TEST_MODULE Walberla Node setters / getters
#define BOOST_TEST_DYN_LINK
#include "config.hpp"

#ifdef LB_WALBERLA

#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <memory>
#include <vector>

#include "boost/mpi.hpp"
#include "grid_based_algorithms/LbWalberla.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_walberla_instance.hpp"
#include "grid.hpp"
#include "utils/Vector.hpp"
#include <iostream>

using Utils::Vector3d;
using Utils::Vector3i;

double viscosity = 3;
Vector3d box_dimensions = {10, 12, 14};
double agrid = 0.5;
Vector3i grid_dimensions{int(box_dimensions[0] / agrid),
                         int(box_dimensions[1] / agrid),
                         int(box_dimensions[2] / agrid)};
double tau = 0.34;
double skin = 0.01;
double density = 2.5;
Vector3i node_grid;

BOOST_AUTO_TEST_CASE(viscosity_test) {
  LbWalberla lb = LbWalberla(viscosity, density, agrid, tau, box_dimensions,
                             node_grid, skin);
  BOOST_CHECK(lb.get_grid_dimensions() == grid_dimensions);
  BOOST_CHECK(lb.get_grid_spacing() == agrid);
  BOOST_CHECK(lb.get_tau() == tau);

  BOOST_CHECK(fabs(lb.get_viscosity() - viscosity) <
              std::numeric_limits<double>::epsilon());
  double new_viscosity = 2.0;
  lb.set_viscosity(new_viscosity);
  BOOST_CHECK(fabs(lb.get_viscosity() - new_viscosity) <
              std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(boundary) {
  Vector3d vel = {0.2, 3.8, 4.2};
  LbWalberla lb = LbWalberla(viscosity, density, agrid, tau, box_dimensions,
                             node_grid, skin);
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
      auto res = lb.get_node_is_boundary(node);
      // Did we get a value?
      BOOST_CHECK(res);
      // Should not be a boundary node
      BOOST_CHECK(*res == false);
    } else {
      // Not on local domain. None of the following should succeed.
      BOOST_CHECK(!lb.set_node_velocity_at_boundary(node, vel));
      BOOST_CHECK(!lb.get_node_velocity_at_boundary(node));
      BOOST_CHECK(!lb.remove_node_from_boundary(node));
      BOOST_CHECK(!lb.get_node_is_boundary(node));
    }
  }
}

// BOOST_AUTO_TEST_CASE(boundary_flow_single_node) {
//  Vector3d vel = {0.2, 3.8, 0};
//  LbWalberla lb = LbWalberla(viscosity, density, agrid, tau, box_dimensions,
//                             node_grid, skin);
//    Vector3i node{5,5,5};
//    if (lb.node_in_local_domain(node)) {
//      BOOST_CHECK(lb.set_node_velocity_at_boundary(node, vel));
//    }
//    for( int i=0;i<10;i++) {
//      lb.integrate();
//    }
//    for (int j=0;j<9;j++) {
//        auto v = lb.get_node_velocity(Vector3i{j,node[1],node[2]});
//        if (v) {
//          if (j!=node[0]) {
//            BOOST_CHECK((*v)[0]>1E-4);
//            BOOST_CHECK((*v)[1]>1E-4);
//            BOOST_CHECK(fabs((*v)[2])<1E-12);
//            printf("%d,%g %g %g\n",j,(*v)[0],(*v)[1],(*v)[2]);
//          }
//        }
//  }
//}

BOOST_AUTO_TEST_CASE(boundary_flow_shear) {
  Vector3d vel = {0.2, -0.3, 0};
  LbWalberla lb = LbWalberla(viscosity, density, agrid, tau, box_dimensions,
                             node_grid, skin);
  for (int x = 1; x < grid_dimensions[0] - 1; x++) {
    for (int y = 1; y < grid_dimensions[1] - 1; y++) {
      Vector3i node{x, y, 1};
      if (lb.node_in_local_domain(node)) {
        BOOST_CHECK(lb.set_node_velocity_at_boundary(node, vel));
      }
      node[2] = grid_dimensions[2] - 4;
      if (lb.node_in_local_domain(node)) {
        BOOST_CHECK(lb.set_node_velocity_at_boundary(node, -vel));
      }
    }
  }

  for (int i = 0; i < 150; i++) {
    lb.integrate();
  }
  for (int j = 0; j < grid_dimensions[2]; j++) {
    auto v = lb.get_node_velocity(Vector3i{0, 0, j});
  }
}

BOOST_AUTO_TEST_CASE(velocity) {
  LbWalberla lb = LbWalberla(viscosity, density, agrid, tau, box_dimensions,
                             node_grid, skin);
  for (Vector3i node : std::vector<Vector3i>{
           {2, 2, 3}, {1, 0, 0}, {0, 1, 2}, {3, 2, 3}, {3, 2, 3}}) {
    const Vector3d pos = Vector3d{
        {double(node[0] + .5), double(node[1] + .5), double(node[2] + .5)}};
    if (lb.node_in_local_domain(node)) {
      auto res = lb.get_node_velocity(node);
      BOOST_CHECK(res);
      BOOST_CHECK_SMALL((*res - Vector3d{0, 0, 0}).norm(), 1E-10);
      const Vector3d v{{double(node[0]) + 1, -1, 2.5 - double(node[2])}};
      double eps = 1E-8;
      BOOST_CHECK(lb.set_node_velocity(node, v));
      res = lb.get_node_velocity(node);
      BOOST_CHECK(res);
      BOOST_CHECK((*res - v).norm() < eps);
      res = lb.get_velocity_at_pos(pos);
      BOOST_CHECK(res);
      BOOST_CHECK_SMALL((*res - v).norm(), 1E-10);
      BOOST_CHECK(lb.set_node_velocity(node, Vector3d{{0, 0, 0}}));
    } else {
      BOOST_CHECK(!lb.get_node_velocity(node));
      const Vector3d v{{double(node[0]) + 1, -1, 2.5 - double(node[2])}};
      BOOST_CHECK(!lb.set_node_velocity(node, v));
      BOOST_CHECK(!lb.get_node_velocity(node));
      BOOST_CHECK(!lb.get_velocity_at_pos(pos));
      BOOST_CHECK(!lb.set_node_velocity(node, Vector3d{{0, 0, 0}}));
    }
  }
  box_geo.set_length(box_dimensions);

  std::vector<Vector3d> corners;
  auto const local_domain =lb.get_local_domain();
  auto const my_left = local_domain.first;
  auto const my_right = local_domain.second;
  auto const local_domain_size = my_right - my_left;
  for (double x : {0,1}) 
    for (double y : {0,1}) 
      for (double z: {0,1}) {
        corners.push_back(my_left +hadamard_product(local_domain_size, Vector3d{x,y,z}));
      }
  
  // check corners
  for (Vector3d pos: corners) {
    // Left corners should be in local domain
    if (pos == my_left) {
      BOOST_CHECK(lb.pos_in_local_domain(pos));
    }
    // All corners should be in local halo
    BOOST_CHECK(lb.pos_in_local_halo(pos));
    // velocity should be accessible
    BOOST_CHECK(lb.get_velocity_at_pos(pos));
  }
  
  // Corners + offsets
  std::vector<Vector3d> offsets;
  for (double x : {-1, 1}) 
    for (double y : {-1, 1}) 
      for (double z: {-1, 1}) {
        offsets.push_back(0.99* Vector3d{x,y,z});
      }
  for (auto corner : corners)
    for (auto offset: offsets) {
      auto const pos = corner + offset;
      // Should be in local halo
      BOOST_CHECK(lb.pos_in_local_halo(pos));
      // Velocity should be accessible
      BOOST_CHECK(lb.get_velocity_at_pos(pos));
      // Positions >= my_left and < my_right
      bool in_local_domain = true;
      for (int i=0; i<3;i++) 
        if (pos[i] < my_left[i] or pos[i] >= my_right[i]) in_local_domain = false;
      
      if (in_local_domain) 
        BOOST_CHECK(lb.pos_in_local_domain(pos));
      else {
        BOOST_CHECK(! lb.pos_in_local_domain(pos));
        BOOST_CHECK(lb.pos_in_local_halo(pos));
      }
    }
}

BOOST_AUTO_TEST_CASE(total_momentum) {
  LbWalberla lb = LbWalberla(viscosity, density, agrid, tau, box_dimensions,
                             node_grid, skin);
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
  LbWalberla lb = LbWalberla(viscosity, density, agrid, tau, box_dimensions,
                             node_grid, skin);
  auto f = Vector3d{0.015, 0.25, -0.22};
  lb.set_external_force(f);
  BOOST_CHECK_SMALL(lb.get_momentum().norm(), 1E-10);

  for (int i = 1; i < 30; i++) {
    lb.integrate();
    auto mom = lb.get_momentum();
    auto mom_exp = (i + .5) * f * grid_dimensions[0] * grid_dimensions[1] *
                   grid_dimensions[2];
    MPI_Allreduce(MPI_IN_PLACE, mom.data(), 3, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    //    printf("%d, %g %g %g, %g %g
    //    %g\n",i,mom[0],mom[1],mom[2],mom_exp[0],mom_exp[1],mom_exp[2]);
    BOOST_CHECK_SMALL((mom - mom_exp).norm(), 1E-7);
  }
}

BOOST_AUTO_TEST_CASE(integrate_with_point_forces) {
  LbWalberla lb = LbWalberla(viscosity, density, agrid, tau, box_dimensions,
                             node_grid, skin);
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
    BOOST_CHECK_SMALL((mom - mom_exp).norm(), 8E-5);
  }
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  walberla_mpi_init();
  int n_nodes;
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  MPI_Dims_create(n_nodes, 3, node_grid.data());

  auto res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

#else // ifdef LB_WALBERLA
int main(int argc, char **argv) {}
#endif
