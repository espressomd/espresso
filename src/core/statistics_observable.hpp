/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef _STATISTICS_OBSERVABLE_H
#define _STATISTICS_OBSERVABLE_H

#include "config.hpp"
#include "utils.hpp"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <map>

#define CONST_UNITITIALIZED 1e-23

enum ObservableType { OBSERVABLE, AVERAGE, VARIANCE };

struct s_observable;

class Observable {
  public:
    Observable();
    int update();
    int calculate();
    virtual int actual_calculate(); 
    virtual int actual_update() ;


    /* IO functions for observables */
    int write(char *filename, bool binary);
    int read(char *filename, bool binary);
 

    // These should be made private with setters and getters.
    ObservableType type;
    char* obs_name;
    void* container;
    int n;
    double* last_value;
    double last_update;
    int autoupdate;
    double autoupdate_dt;
};


extern std::map<int,Observable*> observables;

// flag if autoupdates are necessary.
extern int observables_autoupdate;

void autoupdate_observables(); 


class ObservableParticleVelocities: public Observable {
  int actual_calculate();
};

class ObservableParticleBodyVelocities: public Observable {
  int actual_calculate();
};

class ObservableParticleAngularMomentum: public Observable {
  int actual_calculate();
};

class ObservableParticleBodyAngularMomentum: public Observable {
  int actual_calculate();
};

class ObservableComVelocity: public Observable { 
  int actual_calculate();
};

class ObservableBlockedComVelocity: public Observable { 
  int actual_calculate();
};

/** Obtain the particle positions.
 * TODO: Folded or unfolded?
 */ 
class ObservableParticlePositions: public Observable {
  int actual_calculate();
};

class ObservableParticleForces: public Observable {
  int actual_calculate();
};

class ObservableComForce: public Observable {
  int actual_calculate();
};


class ObservableBlockedComForce: public Observable {
  int actual_calculate();
};
class ObservableStressTensor: public Observable {
  int actual_calculate();
};
  
class ObservableStressTensorAcfObs: public Observable {
  int actual_calculate();
};

class ObservableComPosition: public Observable {
  int actual_calculate();
};

class ObservableBlockedComPosition: public Observable {
  int actual_calculate();
};

#ifdef ELECTROSTATICS

class ObservableParticleCurrents: public Observable {
  int actual_calculate();
};
class ObservableCurrents: public Observable {
  int actual_calculate();
};

class ObservableDipoleMoment: public Observable {
  int actual_calculate();
};

#endif

#ifdef DIPOLES
class ObservableComDipoleMoment: public Observable {
  int actual_calculate();
};
#endif

#ifdef LB
int mpi_observable_lb_radial_velocity_profile_parallel(void* pdata_, double* A, unsigned int n_A);
#endif

class ObservableAverage: public Observable {
  int actual_update();
  int reset();
};


typedef struct {
  Observable* reference_observable;
  unsigned int n_sweeps;
} observable_average_container;

/** Calculate structure factor from positions and scattering length */
class ObservableStructureFactor: public Observable {
  int actual_calculate();
};


/** Calculate structure factor from positions and scattering length */
class ObservableStructureFactorFast: public Observable {
  int actual_calculate();
};

typedef struct {
// FIXME finish the implementation of scattering length
  IntList* id_list;
  DoubleList *scattering_length; // Scattering lengths of particles
  int order;
  int dim_sf; // number of q vectors
  int *q_vals; // values of q vectors
  double *q_density; // number of q vectors per bin
  // entries for fast version
  //int num_k_vecs;
  int k_density;
} observable_sf_params;

/** See if particles from idList1 interact with any of the particles in idList2 
input parameters are passed via struct iw_params
*/
class ObservableInteractsWith: public Observable {
  int actual_calculate();
};


typedef struct {
  double cutoff;
  IntList *ids1;
  IntList *ids2;
} iw_params;


/** Do nothing */
class ObservableNothing: public Observable {
  int actual_calculate();
};

class ObservableFluxDensityProfile: public Observable {
  int actual_calculate();
};

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

class ObservableDensityProfile: public Observable {
  int actual_calculate();
};


class ObservableForceDensityProfile: public Observable {
  int actual_calculate();
};


class ObservableLbVelocityProfile: public Observable {
  int actual_calculate();
};


class ObservableRadialDensityProfile: public Observable {
  int actual_calculate();
};

class ObservableRadialFluxDensityProfile: public Observable {
  int actual_calculate();
};

class ObservableLbRadialVelocityProfile: public Observable {
  int actual_calculate();
};

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


void mpi_observable_lb_radial_velocity_profile_slave_implementation();

class ObservableRadialDensityDistribution: public Observable {
  int actual_calculate();
};


typedef struct { 
	IntList *id_list;
	int type;
	double minr;
	double maxr;
	int rbins;
	int start_point_id;
	int end_point_id;
	// id_flag == 0 : actual positions given, otherwise two particle ids for the start- and 
	// end-point are given
	int id_flag;
	double start_point[3];
	double end_point[3];
} radial_density_data;

class ObbservableSpatialPolymerProperties: public Observable {
  int actual_calculate();
};

typedef struct { 
	IntList *id_list;
	int npoly;
	int cut_off;
} spatial_polym_data;

class ObservablePersistenceLength: public Observable {
  int actual_calculate();
};


// uses the same data as spatial_polymer_properties

typedef struct {
	IntList *id_list;
	int poly_len;
	int npoly;
	int k;
	int n_bins;
	double r_min;
	double r_max;
} k_dist_data;
class ObservablePolymerKDistribution: public Observable {
  int actual_calculate();
};


class ObservableRdf: public Observable {
  int actual_calculate();
};

typedef struct {
  int *p1_types;
  int n_p1;
  int *p2_types;
  int n_p2;
  double r_min;
  double r_max;
  int r_bins;
} rdf_profile_data;


#endif
