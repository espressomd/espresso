#include "config.hpp"
#include "electrokinetics.hpp"

#include "particle_data.hpp"
#include "statistics.hpp"

#include "cuda_interface.hpp"
#include "lbgpu.hpp"

extern CUDA_particle_data *particle_data_host;
extern LB_parameters_gpu lbpar_gpu;

void ek_particles_add_momentum ( ekfloat velocity[3] ) {
  //std::cout << velocity[0] << velocity[1] << velocity[2] << std::endl;
  //std::cout << lbpar_gpu.number_of_particles << std::endl;
  //std::cout << particle_data_host[0].mass << std::endl;

  for(int i=0; i<lbpar_gpu.number_of_particles; i++){
    double mass;
#ifdef MASS
    mass = particle_data_host[i].mass;
#else
    mass = 1.;
#endif

    double new_velocity[3] = {
      particle_data_host[i].v[0] + velocity[0] * mass,
      particle_data_host[i].v[1] + velocity[1] * mass,
      particle_data_host[i].v[2] + velocity[2] * mass
    };
    //std::cout << new_velocity[0] << new_velocity[1] << new_velocity[2] << std::endl;
    set_particle_v( i, new_velocity );
  }
}

void ek_get_total_momentum ( ekfloat momentum[3] ) {
  std::vector<double> total_momentum;
  total_momentum = calc_linear_momentum(1,1);
  momentum[0] = total_momentum[0];
  momentum[1] = total_momentum[1];
  momentum[2] = total_momentum[2];
  return;
}
