#ifndef OBSERVABLES_PROFIELS_HPP
#define OBSERVABLES_PROFIELS_HPP
typedef struct { 
  IntList* id_list;
  double minx;
  double maxx;
  double miny;
  double maxy;
  double minz;
  double maxz;
  int xbins;
  int ybins;
  int zbins;
  void* container;
} profile_data;

typedef struct {
  IntList* id_list;
  double minr;
  double maxr;
  double minphi;
  double maxphi;
  double minz;
  double maxz;
  double center[3];
  double axis[3];
  int phibins;
  int rbins;
  int zbins;
  void* container;
} radial_profile_data;


#endif

