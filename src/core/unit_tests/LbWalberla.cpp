#define BOOST_TEST_MODULE Walberla Node setters / getters
#define BOOST_TEST_DYN_LINK
#include "config.hpp"

#ifdef LB_WALBERLA

#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <memory>
#include <vector>

#include "boost/mpi.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/LbWalberlaD3Q19TRT.hpp"
#include "grid_based_algorithms/LbWalberla_impl.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_walberla_instance.hpp"
#include "utils/Vector.hpp"
#include <iostream>

using Utils::Vector3d;
using Utils::Vector3i;
using walberla::LbWalberlaD3Q19TRT;

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
  LbWalberlaD3Q19TRT lb = LbWalberlaD3Q19TRT(viscosity, density, agrid, tau,
                                             box_dimensions, mpi_shape, 2);
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
  LbWalberlaD3Q19TRT lb = LbWalberlaD3Q19TRT(viscosity, density, agrid, tau,
                                             box_dimensions, mpi_shape, 2);
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
      // Not in the local halo.
      BOOST_CHECK(!lb.set_node_velocity_at_boundary(node, vel));
      BOOST_CHECK(!lb.get_node_velocity_at_boundary(node));
      BOOST_CHECK(!lb.remove_node_from_boundary(node));
      BOOST_CHECK(!lb.get_node_is_boundary(node));
    }
  }
}

BOOST_AUTO_TEST_CASE(boundary_flow_single_node) {
  Vector3d vel = {0.2, 3.8, 0};
  LbWalberlaD3Q19TRT lb = LbWalberlaD3Q19TRT(viscosity, density, agrid, tau,
  box_dimensions,
                             mpi_shape, 2);
    Vector3i node{5,5,5};
    if (lb.node_in_local_halo(node)) {
      BOOST_CHECK(lb.set_node_velocity_at_boundary(node, vel));
    }
    for( int i=0;i<10;i++) {
      lb.integrate();
    }
    for (int j=0;j<9;j++) {
        auto v = lb.get_node_velocity(Vector3i{j,node[1],node[2]});
        if (v) {
          if (j!=node[0]) {
            BOOST_CHECK((*v)[0]>1E-4);
            BOOST_CHECK((*v)[1]>1E-4);
            BOOST_CHECK(fabs((*v)[2])<1E-12);
          }
        }
  }
}

BOOST_AUTO_TEST_CASE(boundary_flow_shear) {
  Vector3d vel = {0.2, -0.3, 0};
  LbWalberlaD3Q19TRT lb = LbWalberlaD3Q19TRT(viscosity, density, agrid, tau,
                                             box_dimensions, mpi_shape, 2);
  int lower=1, upper=12;
  for (int x = -1; x <= grid_dimensions[0]; x++) {
    for (int y = -1; y <= grid_dimensions[1]; y++) {
      Vector3i node{x, y, lower};
      if (lb.node_in_local_halo(node)) {
        BOOST_CHECK(lb.set_node_velocity_at_boundary(node, Vector3d{}));
        BOOST_CHECK(lb.get_node_is_boundary(node));
      }
      node[2] =upper;
      if (lb.node_in_local_halo(node)) {
        BOOST_CHECK(lb.set_node_velocity_at_boundary(node, vel));
      }
    }
  }

  for (int i = 0; i < 200; i++) {
    lb.integrate();
  }
  for (int j = lower+1; j < upper; j++) {
    auto v = lb.get_node_velocity(Vector3i{0, 0, j});
    if (v) {
     double res = ((*v)-(j-lower)/double(upper-lower) *vel).norm();
     BOOST_CHECK_SMALL(res,0.05*vel.norm());
    }
       
  }
}



std::vector<Vector3i> all_nodes_incl_ghosts(int n_ghost_layers) {
  std::vector<Vector3i> res;
  for (int x = -n_ghost_layers; x<grid_dimensions[0]+n_ghost_layers; x++) {
    for (int y = -n_ghost_layers; y<grid_dimensions[1]+n_ghost_layers; y++) {
      for (int z = -n_ghost_layers; z<grid_dimensions[2]+n_ghost_layers; z++) {
        res.push_back(Vector3i{x,y,z});
      }
    }
  }
  return res;
}

BOOST_AUTO_TEST_CASE(domain_and_halo) {
  int n_ghost_layers =2;
  LbWalberlaD3Q19TRT lb = LbWalberlaD3Q19TRT(viscosity, density, agrid, tau,
                                             box_dimensions, mpi_shape, n_ghost_layers);
  
  auto my_left = lb.get_local_domain().first;
  auto my_right = lb.get_local_domain().second;

  for (auto const& n: all_nodes_incl_ghosts(n_ghost_layers)) {
    const Vector3d pos = Vector3d{
          {double(n[0] + .5), double(n[1] + .5), double(n[2] + .5)}};
    int is_local= 0;
    // Nodes in local domain
    if ((n[0]>=my_left[0] and n[1] >= my_left[1] and n[2] >= my_left[2]) and
        (n[0]<my_right[0] and n[1]<my_right[1] and n[2] < my_right[2])) {
          BOOST_CHECK(lb.node_in_local_domain(n));
          BOOST_CHECK(lb.node_in_local_halo(n));
          
          BOOST_CHECK(lb.pos_in_local_domain(pos));
          BOOST_CHECK(lb.pos_in_local_halo(pos));
          is_local = 1;
    } else {
      // in local halo?
      if ((n[0] + n_ghost_layers >=my_left[0] and n[1]+n_ghost_layers >= my_left[1] and  n[2] + n_ghost_layers >= my_left[2]) and
        (n[0]-n_ghost_layers <my_right[0] and n[1]-n_ghost_layers <my_right[1] and n[2]-n_ghost_layers < my_right[2])) {
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
    if (n[0]>=0 and n[1]>=0 and n[2]>=0 and
        n[0]<grid_dimensions[0] and n[1] < grid_dimensions[1] and n[2] < grid_dimensions[2]) {
      MPI_Allreduce(MPI_IN_PLACE, &is_local, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
      BOOST_CHECK(is_local==1);
    }
  }
}


BOOST_AUTO_TEST_CASE(velocity) {
  int n_ghost_layers = 2;
  LbWalberlaD3Q19TRT lb =
      LbWalberlaD3Q19TRT(viscosity, density, agrid, tau, box_dimensions,
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
      auto res = lb.get_velocity_at_pos(n_pos(node));
      BOOST_CHECK(res);                                    // locallly available
      auto v_exp = n_vel(node);
      printf("%d %d %d: %g %g %g | %g %g %g\n", node[0], node[1], node[2],
             (*res)[0], (*res)[1], (*res)[2], v_exp[0], v_exp[1], v_exp[2]);
      BOOST_CHECK_SMALL((*res - n_vel(node)).norm(), eps); // value correct?
    } else {
      // Check that access to node velocity is not possible
      BOOST_CHECK(!lb.get_node_velocity(node));
      BOOST_CHECK(!lb.get_velocity_at_pos(n_pos(node)));
    }
  }
};



BOOST_AUTO_TEST_CASE(total_momentum) {
  LbWalberlaD3Q19TRT lb = LbWalberlaD3Q19TRT(viscosity, density, agrid, tau,
                                             box_dimensions, mpi_shape, 2);
  auto v = Vector3d{1.5, 2.5, -2.2};
  lb.set_node_velocity(Vector3i{1, 1, 1}, v);
  lb.set_node_velocity(Vector3i{3, 5, 7}, v);
  auto mom = lb.get_momentum();
  auto mom_exp = 2 * density * v;
  MPI_Allreduce(MPI_IN_PLACE, mom.data(), 3, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  BOOST_CHECK_SMALL((mom - mom_exp).norm(), 1E-10);
};

BOOST_AUTO_TEST_CASE(integrate_with_volume_force) {
  LbWalberlaD3Q19TRT lb = LbWalberlaD3Q19TRT(viscosity, density, agrid, tau,
                                             box_dimensions, mpi_shape, 2);
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
  LbWalberlaD3Q19TRT lb = LbWalberlaD3Q19TRT(viscosity, density, agrid, tau,
                                             box_dimensions, mpi_shape, 2);
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

BOOST_AUTO_TEST_CASE(forces) {
  LbWalberlaD3Q19TRT lb = LbWalberlaD3Q19TRT(viscosity, density, agrid, tau,
                                             box_dimensions, mpi_shape, 2);
  Vector3d minus = grid_dimensions / 2 - Vector3d{0.1, 0.1, 0.1};
  Vector3d plus = grid_dimensions / 2 + Vector3d{0.1, 0.1, 0.1};
  for (Vector3d pos : std::vector<Vector3d>{{1.25, 1.25, 1.25},
                                            {1.0, 3.0, 1.0},
                                            {4.2, 1.7, 1.6},
                                            {3.2, 1.3, 3.1},
                                            minus,
                                            plus}) {
    if (lb.pos_in_local_halo(pos)) {
      auto res = lb.get_force_to_be_applied_at_pos(pos);
      BOOST_CHECK(res);
      BOOST_CHECK_SMALL((*res - Vector3d{0, 0, 0}).norm(), 1E-10);
      const Vector3d f{{pos[0] + 1, -1, 2.5 - pos[2]}};
      double eps = 1E-8;
      BOOST_CHECK(lb.add_force_at_pos(pos, f));
      *res = {0.0, 0.0, 0.0};
      for (int x = -1; x <= 1; x++)
        for (int y = -1; y <= 1; y++)
          for (int z = -1; z <= 1; z++) {
            auto res_temp = lb.get_force_to_be_applied_at_pos(
                {pos[0] + x, pos[1] + y, pos[2] + z});
            if (res_temp)
              *res += *res_temp;
          }
      if (lb.pos_in_local_domain(
              pos)) // If the force is applied on the halo region then the some
                    // of the nodes are outside of the ghost layer and no forces
                    // are added on these nodes, therefore the sum of all forces
                    // is not equal the injected force. TODO: maybe the forces
                    // distibuted on local nodes from the halo should to be
                    // checked further for correctness.
        BOOST_CHECK((*res - f).norm() < eps);
      BOOST_CHECK(lb.add_force_at_pos(pos, -1 * f));
      res = lb.get_force_to_be_applied_at_pos(pos);
      BOOST_CHECK(res);
      BOOST_CHECK_SMALL((*res - Vector3d{0, 0, 0}).norm(), 1E-10);
    } else {
      BOOST_CHECK(!lb.get_force_to_be_applied_at_pos(pos));
    }
  }
}

BOOST_AUTO_TEST_CASE(last_forces) {
  LbWalberlaD3Q19TRT lb = LbWalberlaD3Q19TRT(viscosity, density, agrid, tau,
                                             box_dimensions, mpi_shape, 2);
  auto positions = std::vector<Vector3d>{{1.25, 2.25, 1.25},
                                         {10.0, 3.0, 1.0},
                                         {16.2, 10.7, 7.6},
                                         {9.4, 2.3, 10.1},
                                         {10.6, 14.2, 10.2}};
  for (auto pos : positions) {
    if (lb.pos_in_local_halo(pos)) {
      auto res = lb.get_force_to_be_applied_at_pos(pos);
      BOOST_CHECK(res);
      BOOST_CHECK_SMALL((*res - Vector3d{0, 0, 0}).norm(), 1E-10);
      const Vector3d f{{pos[0] + 1, -1, 2.5 - pos[2]}};
      BOOST_CHECK(lb.add_force_at_pos(pos, f));
    }
  }

  lb.integrate();

  for (auto pos : positions) {
    if (lb.pos_in_local_halo(pos)) {
      const Vector3d f{{pos[0] + 1, -1, 2.5 - pos[2]}};
      double eps = 1E-8;
      Vector3d res = {0.0, 0.0, 0.0};
      for (int x = -1; x <= 1; x++)
        for (int y = -1; y <= 1; y++)
          for (int z = -1; z <= 1; z++) {
            auto res_temp = lb.get_force_last_applied_at_pos(
                {pos[0] + x, pos[1] + y, pos[2] + z});
            if (res_temp)
              res += *res_temp;
          }
      if (lb.pos_in_local_domain(pos))
        BOOST_CHECK((res - f).norm() < eps);
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
