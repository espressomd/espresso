#include "config.hpp"

#ifdef STOKESIAN_DYNAMICS

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <utils/Vector.hpp>

#include <boost/config.hpp>
#include <boost/serialization/split_free.hpp>

#include "communication.hpp"
#include "errorhandling.hpp"
#include "integrate.hpp"
#include "partCfg_global.hpp"
#include "particle_data.hpp"
#include "thermostat.hpp"

#if defined(BLAS) && defined(LAPACK)
#include "sd_cpu.hpp"
#endif

#ifdef CUDA
#include "sd_gpu.hpp"
#endif

#include "sd_interface.hpp"

namespace {
double sd_viscosity = -1.0;

enum { CPU, GPU, INVALID } device = INVALID;

std::unordered_map<int, double> radius_dict;

double sd_kT = 0.0;

std::size_t sd_seed = 0UL;

int sd_flags = 0;

/* type for particle data transfer between nodes */
struct SD_particle_data {
  int id = -1;
  
  /* particle radius */
  double r = -1;   // Where do I get this?  Particle type table?

  /* particle position */
  Utils::Vector3d pos = {0., 0., 0.};
  
  /* quaternion to define particle orientation */
  Utils::Vector4d quat = {1., 0., 0., 0.};  // TODO: make sure rotations work as they should
  
  /* external force */
  Utils::Vector3d ext_force = {0.0, 0.0, 0.0};  // double or float?
  
  /* external torque */
  Utils::Vector3d ext_torque = {0.0, 0.0, 0.0};
};

/* buffer that holds local particle data, and all particles on the master node 
   used for communication between nodes */
std::vector<SD_particle_data> parts_buffer{};

} // namespace



namespace boost {
namespace serialization {

template <typename Archive>
void load(Archive &ar, SD_particle_data &p,
          const unsigned int /* file_version */) {
  ar >> make_array(reinterpret_cast<char *>(&p), sizeof(SD_particle_data));
}

template <typename Archive>
void save(Archive &ar, SD_particle_data const &p,
          const unsigned int /* file_version */) {
  /* Cruel but effective */
  ar << make_array(reinterpret_cast<char const *>(&p),
                   sizeof(SD_particle_data));
}

template <class Archive>
void serialize(Archive &ar, SD_particle_data &p,
               const unsigned int file_version) {
  split_free(ar, p, file_version);
}

} // namespace serialization
} // namespace boost


/* packs selected properties of all local particles into a buffer */
void SD_pack_local_particle_data() {
  auto &parts = partCfg();
  std::size_t n_part = parts.size();
  parts_buffer.resize(n_part);
  std::size_t i = 0;
  
  for (auto const &p: parts) {
    parts_buffer[i].id = p.p.is_virtual ? -1 : p.p.identity;
    
    parts_buffer[i].pos[0] = p.r.p[0];
    parts_buffer[i].pos[1] = p.r.p[1];
    parts_buffer[i].pos[2] = p.r.p[2];
    
    // TODO: rotation
    parts_buffer[i].quat[0] = 1.0;
    parts_buffer[i].quat[1] = 0.0;
    parts_buffer[i].quat[2] = 0.0;
    parts_buffer[i].quat[3] = 0.0;
    
    parts_buffer[i].ext_force[0] = p.f.f[0];
    parts_buffer[i].ext_force[1] = p.f.f[1];
    parts_buffer[i].ext_force[2] = p.f.f[2];
    
    parts_buffer[i].ext_torque[0] = p.f.torque[0];
    parts_buffer[i].ext_torque[1] = p.f.torque[1];
    parts_buffer[i].ext_torque[2] = p.f.torque[2];
    
    //TODO: radius
    parts_buffer[i].r = radius_dict[p.p.type];
    
    i++;
  }
}

void SD_gather_particle_data() {
  /* Resizes the buffer always to local particle number. 
   * Not desirable on master node. TODO: make it more efficient.*/
  SD_pack_local_particle_data(); 
  if (this_node > 0) {
    Utils::Mpi::gather_buffer(parts_buffer, comm_cart);

    //Utils::Mpi::gather_buffer(buffer.data(), buffer.size(), comm_cart);
  } else {
    Utils::Mpi::gather_buffer(parts_buffer, comm_cart);

    //Utils::Mpi::gather_buffer(particle_data_host, n_part, comm_cart);
  }
}

void set_sd_viscosity(double eta) { sd_viscosity = eta; }

double get_sd_viscosity() { return sd_viscosity; }

void set_sd_device(std::string const &dev) {
  if (dev == "cpu") {
    device = CPU;
  } else if (dev == "gpu") {
    device = GPU;
  } else {
    device = INVALID;
  }
}

std::string get_sd_device() {
  switch (device) {
  case CPU:
    return "cpu";
  case GPU:
    return "gpu";
  default:
    return "invalid";
  }
}

void set_sd_radius_dict(std::unordered_map<int, double> const &x) {
  radius_dict = x;
}

std::unordered_map<int, double> get_sd_radius_dict() { return radius_dict; }

void set_sd_kT(double kT) { sd_kT = kT; }

double get_sd_kT() { return sd_kT; }

void set_sd_seed(std::size_t seed) { sd_seed = seed; }

std::size_t get_sd_seed() { return sd_seed; }

void set_sd_flags(int flg) { sd_flags = flg; }

int get_sd_flags() { return sd_flags; }

void propagate_vel_pos_sd() {
  if (thermo_switch & THERMO_SD) {
    if (BOOST_UNLIKELY(sd_viscosity < 0.0)) {
      runtimeErrorMsg() << "sd_viscosity has an invalid value: " +
                               std::to_string(sd_viscosity);
      return;
    }

    if (BOOST_UNLIKELY(sd_kT < 0.0)) {
      runtimeErrorMsg() << "sd_kT has an invalid value: " +
                               std::to_string(sd_viscosity);
      return;
    }

    // get particles data, does not work with parallel computing
    // parts is a C++ vector
    // should do something similar
    // auto &parts = partCfg();
    // std::size_t n_part = parts.size();
    
    SD_gather_particle_data();
    
    if (this_node == 0) {  // brauche ich das?

      std::vector<int> id(n_part);
      std::vector<double> x_host(6 * n_part);
      std::vector<double> f_host(6 * n_part);
      std::vector<double> a_host(n_part);

      std::size_t i = 0;
      for (auto const &p : parts_buffer) {
        id[i] = parts_buffer[i].id; // p.p.is_virtual ? -1 : p.p.identity;
        x_host[6 * i + 0] = p.pos[0]; // p.r.p[0];
        x_host[6 * i + 1] = p.pos[0]; // p.r.p[1];
        x_host[6 * i + 2] = p.pos[0]; // p.r.p[2];
        // TODO: Rotation
        x_host[6 * i + 3] = 0; 
        x_host[6 * i + 4] = 0;
        x_host[6 * i + 5] = 0;

        f_host[6 * i + 0] = p.ext_force[0]; // p.f.f[0];
        f_host[6 * i + 1] = p.ext_force[1]; // p.f.f[1];
        f_host[6 * i + 2] = p.ext_force[2]; // p.f.f[2];
        f_host[6 * i + 3] = p.ext_torque[0]; // p.f.torque[0];
        f_host[6 * i + 4] = p.ext_torque[1]; // p.f.torque[1];
        f_host[6 * i + 5] = p.ext_torque[2]; // p.f.torque[2];

        double radius = p.r; // radius_dict[p.p.type];
        
        // TODO: das hier fixen...
        /*
        if (BOOST_UNLIKELY(radius < 0.0)) {
          runtimeErrorMsg() << "particle of type " + std::to_string(p.p.type) +
                                   " has an invalid radius: " +
                                   std::to_string(radius);
          return;
        }
        */
        
        a_host[i] = radius;

        ++i;
      } // if (this_node==0)

      std::size_t offset = std::round(sim_time / time_step);
      std::vector<double> v_sd;
      switch (device) {

#if defined(BLAS) && defined(LAPACK)
      case CPU:
        v_sd = sd_cpu(x_host, f_host, a_host, n_part, sd_viscosity,
                      sd_kT / time_step, offset, sd_seed, sd_flags);
        break;
#endif

#ifdef CUDA
      case GPU:
        v_sd = sd_gpu(x_host, f_host, a_host, n_part, sd_viscosity,
                      sd_kT / time_step, offset, sd_seed, sd_flags);
        break;
#endif

      default:
        runtimeErrorMsg()
            << "Invalid device for Stokesian dynamics. Available devices:"
#if defined(BLAS) && defined(LAPACK)
               " cpu"
#endif
#ifdef CUDA
               " gpu"
#endif
            ;
        return;
      }

      for (int i = 0; i < n_part; ++i) {
        // skip virtual particles
        if (id[i] == -1) {
          continue;
        }

        // set new velocities and angular velocities
        double v[] = {v_sd[6 * i + 0], v_sd[6 * i + 1], v_sd[6 * i + 2]};
        set_particle_v(id[i], v);

        Utils::Vector3d omega = {v_sd[6 * i + 3], v_sd[6 * i + 4],
                                 v_sd[6 * i + 5]};
        set_particle_omega_body(id[i], omega);

        // integrate and set new positions
        double p[] = {x_host[6 * i + 0] + time_step * v[0],
                      x_host[6 * i + 1] + time_step * v[1],
                      x_host[6 * i + 2] + time_step * v[2]};
        place_particle(id[i], p);
      }
      
    } // if (this_node == 0)
  }
}

#endif // STOKESIAN_DYNAMICS
