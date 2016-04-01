#include "Scafacos.hpp"
#include "errorhandling.hpp"

namespace Electrostatics {
namespace Scafacos {

#define handle_error(stmt) { const FCSResult res = stmt; if(res) runtimeError(fcs_result_get_message(res)); }

std::string Scafacos::get_parameters() {
  return method + " " + m_last_parameters;
}

std::list<std::string> Scafacos::available_methods() {
  std::list<std::string> methods;
  
#ifdef FCS_ENABLE_DIRECT
  methods.push_back("direct");
#endif
#ifdef FCS_ENABLE_EWALD
  methods.push_back("ewald");
#endif
#ifdef FCS_ENABLE_FMM
  methods.push_back("fmm");
#endif
#ifdef FCS_ENABLE_MEMD
  methods.push_back("memd");
#endif
#ifdef FCS_ENABLE_MMM1D
  methods.push_back("mmm1d");
#endif
#ifdef FCS_ENABLE_MMM2D
  methods.push_back("mmm2d");
#endif
#ifdef FCS_ENABLE_P2NFFT
  methods.push_back("p2nfft");
#endif
#ifdef FCS_ENABLE_P3M
  methods.push_back("p3m");
#endif
#ifdef FCS_ENABLE_PEPC
  methods.push_back("pepc");
#endif
#ifdef FCS_ENABLE_PP3MG
  methods.push_back("pp3mg");
#endif
#ifdef FCS_ENABLE_VMG
  methods.push_back("vmg");
#endif
#ifdef FCS_ENABLE_WOLF
  methods.push_back("wolf");
#endif

  return methods;
}

Scafacos::Scafacos(const std::string &_method, MPI_Comm comm, const std::string &parameters) : method(_method) {
  handle_error(fcs_init(&handle, method.c_str(), comm));
  
  int near_flag;
  fcs_get_near_field_delegation(handle, &near_flag);
  has_near = near_flag != 0;

  fcs_set_resort(handle, 0);
    
  parse_parameters(parameters);

}

Scafacos::~Scafacos() {
  fcs_destroy(handle);
}

void Scafacos::parse_parameters(const std::string &s) {
  m_last_parameters = s;
  handle_error(fcs_parser(handle, s.c_str(), 0));

}

double Scafacos::r_cut() const {
  if(has_near) {     
    fcs_float r_cut;
      
    fcs_get_r_cut(handle, &r_cut);      
    return r_cut;
  }
  else {
    return 0.0;
  }
}
  
void Scafacos::set_r_cut(double r_cut) {
  if(has_near) {
    fcs_set_r_cut(handle, r_cut);
  }
}

void Scafacos::run(std::vector<double> &charges, std::vector<double> &positions,
                   std::vector<double> &forces) {
  const int local_n_part = charges.size();

  forces.resize(3*local_n_part);

  handle_error(fcs_tune(handle, local_n_part, &(positions[0]), &(charges[0])));
  
  handle_error(fcs_run(handle, local_n_part, &(positions[0]), &(charges[0]), &(forces[0]), 0));
}

void Scafacos::tune(std::vector<double> &charges, std::vector<double> &positions) {
  handle_error(fcs_tune(handle, charges.size(), &(positions[0]), &(charges[0])));  
}

void Scafacos::set_common_parameters(double *box_l, int *periodicity, int total_particles) {
  double boxa[3] = { 0., 0., 0. }, boxb[3] = { 0., 0., 0. }, boxc[3] = { 0., 0., 0. }, off[3] = { 0., 0., 0. };
  boxa[0] = box_l[0];
  boxb[1] = box_l[1];
  boxc[2] = box_l[2];
  
  handle_error(fcs_set_common(handle, 0, boxa, boxb, boxc, off, periodicity, total_particles));
}

}
}
