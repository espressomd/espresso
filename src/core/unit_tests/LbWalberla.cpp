#define BOOST_TEST_MODULE Node setters/getters
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <vector>
#include <memory>
#include "config.hpp" 

#include "utils/Vector.hpp"
#include "grid_based_algorithms/LbWalberla.hpp"
#include <iostream>
#include "boost/mpi.hpp"

Vector3i grid_dimensions = {10,10,10};
double viscosity = 3;
Vector3d box_dimensions = {10,10,10};
double agrid=box_dimensions[0]/grid_dimensions[0];
double skin = 0.01;
Vector3i node_grid;


BOOST_AUTO_TEST_CASE(viscosity_test) {
LbWalberla lb = LbWalberla(viscosity,agrid,box_dimensions,node_grid,skin);
  BOOST_CHECK(abs(lb.get_viscosity()-viscosity) < std::numeric_limits<double>::epsilon());
  double new_viscosity = 2.0;
  lb.set_viscosity(new_viscosity);
  BOOST_CHECK(abs(lb.get_viscosity()-new_viscosity) < std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(boundary) {
Vector3d vel = {0.2,3.8,4.2};
LbWalberla lb = LbWalberla(viscosity,agrid,box_dimensions,node_grid,skin);
for (Vector3i node : std::vector<Vector3i>{{0,0,0},{0,1,2},{9,9,9}}) {
  if (lb.node_in_local_domain(node)) {
    BOOST_CHECK(lb.set_node_velocity_at_boundary(node,vel));
    auto vel_check = lb.get_node_velocity_at_boundary(node);
    // Do we have a value
    BOOST_CHECK(vel_check);
    // Check the value
    BOOST_CHECK_SMALL((*vel_check - vel).norm(),1E-12);
    BOOST_CHECK(lb.remove_node_from_boundary(node));
    auto res =lb.get_node_is_boundary(node);
    // Did we get a value?
    BOOST_CHECK(res);
    // Should not be a boundary node
    BOOST_CHECK(*res == false);
  } else {
    // Not on local domain. None of the following should succeed.
    BOOST_CHECK(!lb.set_node_velocity_at_boundary(node,vel));
    BOOST_CHECK(!lb.get_node_velocity_at_boundary(node));
    BOOST_CHECK(!lb.remove_node_from_boundary(node));
    BOOST_CHECK(!lb.get_node_is_boundary(node));
  }
}
}

BOOST_AUTO_TEST_CASE(velocity) {
LbWalberla lb = LbWalberla(viscosity,agrid,box_dimensions,node_grid,skin);
for (Vector3i node : std::vector<Vector3i>{{2,2,3},{1,0,0},{0,1,2},{3,2,3},{3,2,3}}) {
  const Vector3d pos =Vector3d{{double(node[0]+.5),double(node[1]+.5),double(node[2]+.5)}} ;
  if (lb.node_in_local_domain(node)) {
  auto res =lb.get_node_velocity(node);
  BOOST_CHECK(res);
  BOOST_CHECK_SMALL((*res-Vector3d{0,0,0}).norm(),1E-10);
  const Vector3d v{{double(node[0])+1,-1,2.5-double(node[2])}};
  double eps=1E-8;
  BOOST_CHECK(lb.set_node_velocity(node,v));
  res=lb.get_node_velocity(node);
  BOOST_CHECK(res);
  BOOST_CHECK((*res-v).norm() <  eps);
  res=lb.get_velocity_at_pos(pos);
  BOOST_CHECK(res);
  BOOST_CHECK_SMALL((*res-v).norm(),1E-10);
  BOOST_CHECK(lb.set_node_velocity(node,Vector3d{{0,0,0}}));
  } else {
  BOOST_CHECK(!lb.get_node_velocity(node));
  const Vector3d v{{double(node[0])+1,-1,2.5-double(node[2])}};
  BOOST_CHECK(!lb.set_node_velocity(node,v));
  BOOST_CHECK(!lb.get_node_velocity(node));
  BOOST_CHECK(!lb.get_velocity_at_pos(pos));
  BOOST_CHECK(!lb.set_node_velocity(node,Vector3d{{0,0,0}}));
}
}
}

int main(int argc, char **argv) {
  MPI_Init(&argc,&argv);
  walberla_mpi_init();
  int n_nodes;
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  MPI_Dims_create(n_nodes, 3, node_grid.data());


  auto res= boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}
