#define BOOST_TEST_MODULE Node setters/getters
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <vector>
#include <memory>
#include "config.hpp" 

#include "utils/Vector.hpp"
#include "grid_based_algorithms/LbWalberla.hpp"
#include <iostream>

Vector3i grid_dimensions = {10,10,10};
double viscosity = 3;
Vector3d box_dimensions = {10,10,10};
double agrid=box_dimensions[0]/grid_dimensions[0];
Vector3i node_grid = {1,1,1};
double skin = 0.01;


BOOST_AUTO_TEST_CASE(viscosity_test) {
LbWalberla lb = LbWalberla(viscosity,agrid,box_dimensions,node_grid,skin);
  BOOST_CHECK(abs(lb.get_viscosity()-viscosity) < std::numeric_limits<double>::epsilon());
  double new_viscosity = 2.0;
  lb.set_viscosity(new_viscosity);
  BOOST_CHECK(abs(lb.get_viscosity()-new_viscosity) < std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(boundary) {
Vector3d vel = {0.2,3.8,4.2};
Vector3d vel_check;
LbWalberla lb = LbWalberla(viscosity,agrid,box_dimensions,node_grid,skin);
for (Vector3i node : std::vector<Vector3i>{{0,0,0},{0,1,2},{9,9,9}}) {
  lb.set_node_velocity_at_boundary(node,vel);
  vel_check = lb.get_node_velocity_at_boundary(node);
  BOOST_CHECK(vel_check == vel);
  lb.remove_node_from_boundary(node);
//  BOOST_CHECK(lb.get_node_is_boundary(node)==0);
}
}

BOOST_AUTO_TEST_CASE(velocity) {
LbWalberla lb = LbWalberla(viscosity,agrid,box_dimensions,node_grid,skin);
for (Vector3i node : std::vector<Vector3i>{{2,2,3},{1,0,0},{0,1,2},{3,2,3},{3,2,3}}) {
  BOOST_CHECK_SMALL((lb.get_node_velocity(node)-Vector3d{0,0,0}).norm(),1E-10);
  const Vector3d v{{double(node[0])+1,-1,2.5-double(node[2])}};
  double eps=1E-8;
  lb.set_node_velocity(node,v);
  BOOST_CHECK((lb.get_node_velocity(node)-v).norm() <  eps);
  const Vector3d pos =Vector3d{{double(node[0]),double(node[1]),double(node[2])}} ;
  Vector3d vn=lb.get_node_velocity(node);
  printf("%d %d %d, %g %g %g, %g %g %g\n",node[0],node[1],node[2],pos[0],pos[1],pos[2],vn[0],vn[1],vn[2]);
  vn=lb.get_velocity_at_pos(pos);
  printf("%d %d %d, %g %g %g, %g %g %g\n",node[0],node[1],node[2],pos[0],pos[1],pos[2],vn[0],vn[1],vn[2]);
  BOOST_CHECK_SMALL((lb.get_velocity_at_pos(pos)-v).norm(),1E-10);
  lb.set_node_velocity(node,Vector3d{{0,0,0}});
}
}
