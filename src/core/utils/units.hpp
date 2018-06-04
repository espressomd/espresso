#ifndef UTILS_UNITS_HPP
#define UTILS_UNITS_HPP

#include <boost/units/base_unit.hpp>
#include <boost/units/derived_dimension.hpp>
#include <boost/units/dimensionless_type.hpp>
#include <boost/units/io.hpp>
#include <boost/units/make_system.hpp>
#include <boost/units/physical_dimensions/current.hpp>
#include <boost/units/physical_dimensions/electric_charge.hpp>
#include <boost/units/physical_dimensions/length.hpp>
#include <boost/units/physical_dimensions/mass.hpp>
#include <boost/units/physical_dimensions/time.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/static_constant.hpp>

namespace Utils {
namespace Units {
using namespace boost::units;
using boost::units::quantity;

struct sigma_unit : base_unit<sigma_unit, length_dimension, 1> {
  static std::string name() { return "lj_sigma"; }
  static std::string symbol() { return "lj_sigma"; }
};

struct mass_unit : base_unit<mass_unit, mass_dimension, 2> {
  static std::string name() { return "lj_mass"; }
  static std::string symbol() { return "lj_mass"; }
};

struct time_unit : base_unit<time_unit, boost::units::time_dimension, 3> {
  static std::string name() { return "lj_time"; }
  static std::string symbol() { return "lj_time"; }
};

struct electric_charge_unit
    : base_unit<electric_charge_unit, boost::units::electric_charge_dimension,
                4> {
  static std::string name() { return "lj_q"; }
  static std::string symbol() { return "lj_q"; }
};

typedef make_system<sigma_unit, mass_unit, time_unit,
                    electric_charge_unit>::type simulation_units;

typedef derived_dimension<length_base_dimension, -3>::type
    number_density_dimension;
typedef derived_dimension<mass_base_dimension, 1, length_base_dimension,
                          -3>::type mass_density_dimension;
typedef derived_dimension<mass_base_dimension, 1, length_base_dimension, 2,
                          time_base_dimension, -2>::type energy_dimension;
typedef derived_dimension<length_base_dimension, 1, time_base_dimension,
                          -2>::type acceleration_dimension;
typedef derived_dimension<mass_base_dimension, 1, length_base_dimension, 1,
                          time_base_dimension, -2>::type force_dimension;
typedef derived_dimension<length_base_dimension, 1, time_base_dimension,
                          -1>::type velocity_dimension;
typedef derived_dimension<mass_base_dimension, 1, length_base_dimension, -1,
                          time_base_dimension,
                          -1>::type dynamic_viscosity_dimension;
typedef derived_dimension<length_base_dimension, 2, time_base_dimension,
                          -1>::type kinematic_viscosity_dimension;
typedef derived_dimension<mass_base_dimension, 1, time_base_dimension, -1>::type
    friction_dimension;

typedef unit<dimensionless_type, simulation_units> dimensionless;
typedef unit<length_dimension, simulation_units> length;
typedef unit<mass_dimension, simulation_units> mass;
typedef unit<time_dimension, simulation_units> time;
typedef unit<electric_charge_dimension, simulation_units> charge;
typedef unit<number_density_dimension, simulation_units> number_density;
typedef unit<mass_density_dimension, simulation_units> mass_density;
typedef unit<energy_dimension, simulation_units> energy;
typedef unit<acceleration_dimension, simulation_units> acceleration;
typedef unit<force_dimension, simulation_units> force;
typedef unit<velocity_dimension, simulation_units> velocity;
typedef unit<current_dimension, simulation_units> electric_current;
typedef unit<dynamic_viscosity_dimension, simulation_units> dynamic_viscosity;
typedef unit<kinematic_viscosity_dimension, simulation_units>
    kinematic_viscosity;
typedef unit<friction_dimension, simulation_units> friction;

BOOST_UNITS_STATIC_CONSTANT(dimless, dimensionless);
BOOST_UNITS_STATIC_CONSTANT(sigma, length);
BOOST_UNITS_STATIC_CONSTANT(m, mass);
BOOST_UNITS_STATIC_CONSTANT(tau, time);
BOOST_UNITS_STATIC_CONSTANT(q, charge);
BOOST_UNITS_STATIC_CONSTANT(n_rho, number_density);
BOOST_UNITS_STATIC_CONSTANT(m_rho, mass_density);
BOOST_UNITS_STATIC_CONSTANT(epsilon, energy);
BOOST_UNITS_STATIC_CONSTANT(a, acceleration);
BOOST_UNITS_STATIC_CONSTANT(f, force);
BOOST_UNITS_STATIC_CONSTANT(v, velocity);
BOOST_UNITS_STATIC_CONSTANT(I, electric_current);
BOOST_UNITS_STATIC_CONSTANT(eta, dynamic_viscosity);
BOOST_UNITS_STATIC_CONSTANT(nu, kinematic_viscosity);
BOOST_UNITS_STATIC_CONSTANT(gamma, friction);
}
}

#endif
