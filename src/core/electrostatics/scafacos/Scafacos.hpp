#include <vector>
#include <string>
#include <list>
#include <fcs.h>
#include <mpi.h>

namespace Electrostatics {
namespace Scafacos {

/** \brief Abstraction of a method from the scafacos library */

struct Scafacos {
  Scafacos(const std::string &_method, MPI_Comm comm, const std::string &parameters);
  ~Scafacos();
  /** Parse parameter string */
  void parse_parameters(const std::string &s);
  /** Get the parameters from the library */
  std::string get_parameters();
  /** Set parameters common to all methods */
  void set_common_parameters(double *box_l, int *periodicity, int total_particles);
  /** Calulate short range pair force if supported by the method */
  inline double pair_force(double dist) const {     
    if(has_near) {
      fcs_float field;
      fcs_compute_near_field(handle, dist, &field);
      return field;
    }

    return 0.0;
  }
  /** Calculate the forces */
  void run(std::vector<double> &charges, std::vector<double> &positions,
           std::vector<double> &forces);
  /** Tune parameters */
  void tune(std::vector<double> &charges, std::vector<double> &positions);
  /** Get shortrange cutoff (0.0 if not supported) */
  double r_cut() const;
  /** Set cutoff */
  void set_r_cut(double r_cut); 

  /** Get a list of methods supported by the library */
  static std::list<std::string> available_methods();
  
  /** Handle from the library */
  FCS handle;
  /** Whether the method supports near field delegation */
  bool has_near;
  /** The scafacos method name of this instance */
  const std::string method;
  /** The last parameters set */
  std::string m_last_parameters;
};


}
}
