/* vim: set ts=8 sts=2 sw=2 et: */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "electrokinetics_pdb_parse.hpp"
#include <iostream>
#include <string>
#include <sstream>

#ifdef EK_BOUNDARIES

/* Replacements for bool variables */
const int pdb_SUCCESS = 0;
const int pdb_ERROR = 1;

float* pdb_charge_lattice = NULL;
int* pdb_boundary_lattice = NULL;

typedef struct {
  int i; // index
  int m; // model index
  float x,y,z;
} pdb_ATOM;

typedef struct {
  int i;
  char type[2];
  float charge;
} itp_atoms;

typedef struct {
  char type[2];
  float sigma,epsilon;
} itp_atomtypes;

typedef struct {
  unsigned int pdb_n_particles;
  pdb_ATOM* pdb_array_ATOM;
  unsigned int itp_n_particles;
  itp_atoms* itp_array_atoms;
  unsigned int itp_n_parameters;
  itp_atomtypes* itp_array_atomtypes;
} particle_data;

typedef struct {
  float max_x;
  float max_y;
  float max_z;
  float min_x;
  float min_y;
  float min_z;
  float center[3];
} bounding_box;

/* BEGIN CODE */

void galloc(void** ptr, size_t size) {
  if (!*ptr) {
    if (size > 0) {
      *ptr = (void*) malloc(size);
    }
    else {
      printf("You cannot malloc to size 0\n");
    }
  }
  else {
    if (size > 0) {
      *ptr = (void*) realloc(*ptr, size);
    }
    else {
      free(*ptr);
    }
  }
}

unsigned int pdb_rhoindex_cartesian2linear(unsigned int x, unsigned int y, unsigned int z) {
  return z * ek_parameters.dim_y * ek_parameters.dim_x + y * ek_parameters.dim_x + x;
}

int print_charge_field(char* filename) {
  FILE* fp;
  if ((fp = fopen(filename,"w")) == NULL) return pdb_ERROR;
  
  if( fp == NULL ) {
    return 1;
  }
  
  fprintf( fp, "\
# vtk DataFile Version 2.0\n\
charge_density\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS charge_density float 1\n\
LOOKUP_TABLE default\n",
  ek_parameters.dim_x, ek_parameters.dim_y, ek_parameters.dim_z,
  ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f,
  ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
  ek_parameters.dim_x * ek_parameters.dim_y * ek_parameters.dim_z);
  
  for( unsigned int i = 0; i < (ek_parameters.dim_x * ek_parameters.dim_y * ek_parameters.dim_z); i++ ) {
    fprintf( fp, "%e ", pdb_charge_lattice[i] );
  }
  
  fclose( fp );
  return pdb_SUCCESS;
}

int print_boundary_lattice(char* filename) {
  FILE* fp;
  if ((fp = fopen(filename,"w")) == NULL) return pdb_ERROR;
  
  if( fp == NULL ) {
    return 1;
  }
  
  fprintf( fp, "\
# vtk DataFile Version 2.0\n\
boundary_flag\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS boundary_flag float 1\n\
LOOKUP_TABLE default\n",
  ek_parameters.dim_x, ek_parameters.dim_y, ek_parameters.dim_z,
  ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f, ek_parameters.agrid*0.5f,
  ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
  ek_parameters.dim_x * ek_parameters.dim_y * ek_parameters.dim_z);

  for( unsigned int i = 0; i < (ek_parameters.dim_x * ek_parameters.dim_y * ek_parameters.dim_z); i++ ) {
    fprintf( fp, "%d ", pdb_boundary_lattice[i] );
  }
  
  fclose( fp );
  return pdb_SUCCESS;
}

int pdb_parse_files(char* pdb_filename, char* itp_filename, particle_data* atom_data) {
  /*
   * This routine parses the pdb- and itp-file to extract
   * the relevant parameters. These are stored in arrays.
   */

  // Parse pdb-file
  int model = 0;
  char pdb_line[256];
  FILE* pdb_file;
  if ((pdb_file = fopen(pdb_filename,"r")) == NULL) return pdb_ERROR;
#ifdef DEBUG
  printf("### Reading pdb-file \"%s\" ###\n",pdb_filename);
#endif  
  while (fgets(pdb_line, sizeof(pdb_line), pdb_file)) {
    if (strncmp(pdb_line,"MODEL",5) == 0) {
      // read the MODEL identifier
      sscanf(pdb_line,"MODEL %d",&model);
#ifdef DEBUG
      printf("MODEL m=%d\n", model);
#endif
    }
    if ( strncmp(pdb_line,"ATOM",4) == 0) {
      // read all ATOMs
      galloc( (void**) &atom_data->pdb_array_ATOM , (atom_data->pdb_n_particles+1)*sizeof(pdb_ATOM) );
      pdb_ATOM* a = &atom_data->pdb_array_ATOM[atom_data->pdb_n_particles];
      // See http://deposit.rcsb.org/adit/docs/pdb_atom_format.html#ATOM for the meaning of the format string
      // sscanf(pdb_line,"ATOM %6d %*4s%*c%*4s%*c%*4d%*c %8f %8f %8f %*6f %*6f %*4s%*2s%*2s",&a->i,&a->x,&a->y,&a->z);
      std::istringstream str(pdb_line);

      std::string tmp;   
 
      str.ignore(246,' ');
      str >> a->i;
      str >> tmp >> tmp >> tmp >> tmp;    
      str >> a->x >> a->y >> a->z;

      a->x /= 10.0;
      a->y /= 10.0;
      a->z /= 10.0;
#ifdef DEBUG
      // Print all local variables
      printf("ATOM i=%d x=%f y=%f z=%f\n",a->i,a->x,a->y,a->z);
#endif
      atom_data->pdb_n_particles++;
    }
  }
  fclose(pdb_file);

  // Parse itp-file
  char itp_line[256];
  FILE* itp_file;
  if ((itp_file = fopen(itp_filename,"r")) == NULL) return pdb_ERROR;
#ifdef DEBUG
  printf("### Reading itp-file \"%s\" ###\n",itp_filename);
#endif  
  while (fgets(itp_line, sizeof(itp_line), itp_file)) {
    // get section
    char section[256];
    sscanf(itp_line,"[ %s ]",section);
    // only read non-comment, non-section and non-empty lines
    // TODO: Handling of lines consiting whitespace only (e.g. '\t' and ' ')
    if (itp_line[0] != '[' && itp_line[0] != ';' && itp_line[0] != '\r' && itp_line[0] != '\n') {
      if (strcmp(section,"atoms") == 0) {
        // section [ atoms ]
        galloc( (void**) &atom_data->itp_array_atoms , (atom_data->itp_n_particles+1)*sizeof(pdb_ATOM) );
        itp_atoms* a = &atom_data->itp_array_atoms[atom_data->itp_n_particles];
        // FIXME: no source :( Reverse engineered from the itp-file
        sscanf(itp_line," %d %2s %*d %*s %*s %*d %f %*f ; %*s %*f",&a->i,a->type,&a->charge);
#ifdef DEBUG
        // Print all local variables
        printf("[ atoms ] i=%d type=%s charge=%f\n",a->i,a->type,a->charge);
#endif
        atom_data->itp_n_particles++;
      }
      if (strcmp(section,"atomtypes") == 0) {
        // section [ atomtypes ]
        galloc( (void**) &atom_data->itp_array_atomtypes , (atom_data->itp_n_parameters+1)*sizeof(pdb_ATOM) );
        itp_atomtypes* a = &atom_data->itp_array_atomtypes[atom_data->itp_n_parameters];
        // FIXME: no source :( Reverse engineered from the itp-file
        sscanf(itp_line," %2s %*s %*f %*f %*c %f %f ; %*f %*f",a->type,&a->sigma,&a->epsilon);
#ifdef DEBUG
        // Print all local variables
        printf("[ atomtypes ] name=%s sigma=%f epsilon=%f\n",a->type,a->sigma,a->epsilon);
#endif
        atom_data->itp_n_parameters++;
      }
    }
  }
  fclose(itp_file);

  if (atom_data->pdb_n_particles != atom_data->itp_n_particles) return pdb_ERROR;
  return pdb_SUCCESS;
}

int calculate_bounding_box(bounding_box* bbox, particle_data* atom_data) {
  // prototype for joining the arrays
  if (atom_data->pdb_n_particles-1 == 0) return pdb_ERROR;
  pdb_ATOM* a = &atom_data->pdb_array_ATOM[0];
  bbox->max_x = a->x;
  bbox->max_y = a->y;
  bbox->max_z = a->z;
  bbox->min_x = a->x;
  bbox->min_y = a->y;
  bbox->min_z = a->z;

  for (unsigned int i = 1; i <= atom_data->pdb_n_particles-1; i++) {
    a = &atom_data->pdb_array_ATOM[i];
    if (bbox->max_x < a->x) bbox->max_x = a->x;
    if (bbox->max_y < a->y) bbox->max_y = a->y;
    if (bbox->max_z < a->z) bbox->max_z = a->z;
    if (bbox->min_x > a->x) bbox->min_x = a->x;
    if (bbox->min_y > a->y) bbox->min_y = a->y;
    if (bbox->min_z > a->z) bbox->min_z = a->z;
  }

  bbox->center[0] = ( bbox->max_x + bbox->min_x )/2;
  bbox->center[1] = ( bbox->max_y + bbox->min_y )/2;
  bbox->center[2] = ( bbox->max_z + bbox->min_z )/2;

  return pdb_SUCCESS;
}

int populate_lattice(particle_data* atom_data) {
  /*
   * This routine will populate the lattice using the
   * values read from the pdb and itp files.
   * WARNING: It contains much logic and interpolation stuff!
   */
#ifdef DEBUG
  printf("pdb_n_particles=%u, itp_n_particles=%u, itp_n_parameters=%u\n",atom_data->pdb_n_particles,atom_data->itp_n_particles,atom_data->itp_n_parameters);
#endif
  // TODO: Check if bounding box fits into simbox
  bounding_box bbox;
  calculate_bounding_box(&bbox, atom_data);

  // calculate the shift of the bounding box
  float shift[3];
  shift[0] = ek_parameters.agrid / 2.0 * ek_parameters.dim_x - bbox.center[0];
  shift[1] = ek_parameters.agrid / 2.0 * ek_parameters.dim_y - bbox.center[1];
  shift[2] = ek_parameters.agrid / 2.0 * ek_parameters.dim_z - bbox.center[2];

#ifdef DEBUG
  printf("bbox.max_x=%f, bbox.max_y=%f, bbox.max_z=%f, bbox.min_x=%f, bbox.min_y=%f, bbox.min_z=%f, bbox->center=[%f; %f; %f]\n", bbox.max_x, bbox.max_y, bbox.max_z, bbox.min_x, bbox.min_y, bbox.min_z, bbox.center[0], bbox.center[1], bbox.center[2]);
  printf("agrid=%f, dim_x=%d, dim_y=%d, dim_z=%d\n",ek_parameters.agrid, ek_parameters.dim_x, ek_parameters.dim_y, ek_parameters.dim_z);
  printf("shift=[%f; %f; %f]\n",shift[0], shift[1], shift[2]);
#endif

  // joining the array
  int lowernode[3];
  float cellpos[3];
  float gridpos;
  float a_x_shifted, a_y_shifted, a_z_shifted;

  for (unsigned int i = 0; i <= atom_data->pdb_n_particles-1; i++) {
    pdb_ATOM* a = &atom_data->pdb_array_ATOM[i];
    itp_atoms* b;
    itp_atomtypes* c;
    for (unsigned int j = 0; j <= atom_data->itp_n_particles-1; j++) {
      b = &atom_data->itp_array_atoms[j];
      if (a->i == b->i) {
        for (unsigned int k = 0; k <= atom_data->itp_n_parameters-1; k++) {
          c = &atom_data->itp_array_atomtypes[k];
          if (strcmp(b->type,c->type) == 0) {
#ifdef DEBUG
            printf("i=%d x=%f y=%f z=%f type=%s charge=%f sigma=%f epsilon=%f\n",a->i,a->x,a->y,a->z,b->type,b->charge,c->sigma,c->epsilon);
#endif

            // Interpolate the charge to the lattice
            gridpos      = (a->x + shift[0]) / ek_parameters.agrid - 0.5f;
            lowernode[0] = (int) floorf( gridpos );
            cellpos[0]   = gridpos - lowernode[0];
                                                
            gridpos      = (a->y + shift[1]) / ek_parameters.agrid - 0.5f;
            lowernode[1] = (int) floorf( gridpos );
            cellpos[1]   = gridpos - lowernode[1];
                                                
            gridpos      = (a->z + shift[2]) / ek_parameters.agrid - 0.5f;
            lowernode[2] = (int) floorf( gridpos );
            cellpos[2]   = gridpos - lowernode[2];
                                                
            lowernode[0] = (lowernode[0] + ek_parameters.dim_x) % ek_parameters.dim_x;
            lowernode[1] = (lowernode[1] + ek_parameters.dim_y) % ek_parameters.dim_y;
            lowernode[2] = (lowernode[2] + ek_parameters.dim_z) % ek_parameters.dim_z;

            pdb_charge_lattice[pdb_rhoindex_cartesian2linear( lowernode[0],lowernode[1],lowernode[2] )]
              = b->charge * ( 1 - cellpos[0] ) * ( 1 - cellpos[1] ) * ( 1 - cellpos[2] );

            pdb_charge_lattice[pdb_rhoindex_cartesian2linear( ( lowernode[0] + 1 ) % ek_parameters.dim_x,lowernode[1],lowernode[2] )]
              = b->charge * cellpos[0] * ( 1 - cellpos[1] ) * ( 1 - cellpos[2] );

            pdb_charge_lattice[pdb_rhoindex_cartesian2linear( lowernode[0],( lowernode[1] + 1 ) % ek_parameters.dim_y,lowernode[2] )]
              = b->charge * ( 1 - cellpos[0] ) * cellpos[1] * ( 1 - cellpos[2] );

            pdb_charge_lattice[pdb_rhoindex_cartesian2linear( lowernode[0],lowernode[1],( lowernode[2] + 1 ) % ek_parameters.dim_z )]
              = b->charge * ( 1 - cellpos[0] ) * ( 1 - cellpos[1] ) * cellpos[2];

            pdb_charge_lattice[pdb_rhoindex_cartesian2linear( ( lowernode[0] + 1 ) % ek_parameters.dim_x,( lowernode[1] + 1 ) % ek_parameters.dim_y,lowernode[2] )]
              = b->charge * cellpos[0] * cellpos[1] * ( 1 - cellpos[2] );

            pdb_charge_lattice[pdb_rhoindex_cartesian2linear( ( lowernode[0] + 1 ) % ek_parameters.dim_x,lowernode[1],( lowernode[2] + 1 ) % ek_parameters.dim_z )]
              = b->charge * cellpos[0] * ( 1 - cellpos[1] ) * cellpos[2];

            pdb_charge_lattice[pdb_rhoindex_cartesian2linear( lowernode[0],( lowernode[1] + 1 ) % ek_parameters.dim_y,( lowernode[2] + 1 ) % ek_parameters.dim_z )]
              = b->charge * ( 1 - cellpos[0] ) * cellpos[1] * cellpos[2];

            pdb_charge_lattice[pdb_rhoindex_cartesian2linear( ( lowernode[0] + 1 ) % ek_parameters.dim_x,( lowernode[1] + 1 ) % ek_parameters.dim_y,( lowernode[2] + 1 ) % ek_parameters.dim_z )]
              = b->charge * cellpos[0] * cellpos[1] * cellpos[2];
            // Interpolate lennard-jones parameters to boundary
            float r = pow(2,1./6.)*c->sigma;

            a_x_shifted = (a->x + shift[0]) / ek_parameters.agrid - 0.5f;
            a_y_shifted = (a->y + shift[1]) / ek_parameters.agrid - 0.5f;
            a_z_shifted = (a->z + shift[2]) / ek_parameters.agrid - 0.5f;

            for (float z = a->z - r; z <= a->z + r + ek_parameters.agrid; z += ek_parameters.agrid) {
              for (float y = a->y - r; y <= a->y + r + ek_parameters.agrid; y += ek_parameters.agrid) {
                for (float x = a->x - r; x <= a->x + r + ek_parameters.agrid; x += ek_parameters.agrid) {
                  gridpos      = (x + shift[0]) / ek_parameters.agrid - 0.5f;
                  lowernode[0] = (int) floorf( gridpos );

                  gridpos      = (y + shift[1]) / ek_parameters.agrid - 0.5f;
                  lowernode[1] = (int) floorf( gridpos );

                  gridpos      = (z + shift[2]) / ek_parameters.agrid - 0.5f;
                  lowernode[2] = (int) floorf( gridpos );

                  lowernode[0] = (lowernode[0] + ek_parameters.dim_x) % ek_parameters.dim_x;
                  lowernode[1] = (lowernode[1] + ek_parameters.dim_y) % ek_parameters.dim_y;
                  lowernode[2] = (lowernode[2] + ek_parameters.dim_z) % ek_parameters.dim_z;
#ifdef DEBUG
                  printf("shifted: %f %f %f\n", a_x_shifted, a_y_shifted, a_z_shifted);
                  printf("lowernode: %d %d %d\n", lowernode[0], lowernode[1], lowernode[2]);
                  printf("distance: %f %f %f\n", lowernode[0] - a_x_shifted, lowernode[1] - a_y_shifted, lowernode[2] - a_z_shifted);
                  printf("distance: %f <= %f\n\n", pow(lowernode[0] - a_x_shifted,2) + pow(lowernode[1] - a_y_shifted,2) + pow(lowernode[2] - a_z_shifted,2), pow(r/ek_parameters.agrid,2));
#endif
                  if ( pow(lowernode[0] - a_x_shifted,2) + pow(lowernode[1] - a_y_shifted,2) + pow(lowernode[2] - a_z_shifted,2) <= pow(r/ek_parameters.agrid,2) ) {
                    pdb_boundary_lattice[ek_parameters.dim_y*ek_parameters.dim_x*lowernode[2] + ek_parameters.dim_x*lowernode[1] + lowernode[0]] = 1;
                  }
                }
              }
            }

            break;
          }
        }
      }
    }
  }

  return pdb_SUCCESS;
}

int pdb_parse(char* pdb_filename, char* itp_filename) {
  /*
   * This is the main parsing routine, which is visible to the outside
   * through the header electrokinetics_pdb_parse.h. It doesn't contain any logic and just
   * deploys the input to the subroutines.
   */

  /* BEGIN DEPLOY */
  galloc( (void**) &pdb_charge_lattice, ek_parameters.dim_x * ek_parameters.dim_y * ek_parameters.dim_z * sizeof(float));
  galloc( (void**) &pdb_boundary_lattice, ek_parameters.dim_x * ek_parameters.dim_y * ek_parameters.dim_z * sizeof(int));
  for ( unsigned int i = 0; i < ek_parameters.dim_x * ek_parameters.dim_y * ek_parameters.dim_z; i++ ) {
    pdb_charge_lattice[i] = 0.0;
    pdb_boundary_lattice[i] = 0;
  }

  particle_data atom_data;
  atom_data.pdb_n_particles = 0;
  atom_data.pdb_array_ATOM = NULL;
  atom_data.itp_n_particles = 0;
  atom_data.itp_array_atoms = NULL;
  atom_data.itp_n_parameters = 0;
  atom_data.itp_array_atomtypes = NULL;

  pdb_parse_files(pdb_filename, itp_filename, &atom_data);

  populate_lattice(&atom_data);

  return pdb_SUCCESS;
}

#endif
