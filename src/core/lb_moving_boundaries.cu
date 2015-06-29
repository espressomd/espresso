#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <curand.h>
#include <cuda_runtime.h>
#include <cuda.h>

#define checkCudaErrors(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}



__constant__ int v[19][3] =        { {  0,  0,  0 },
                                     {  1,  0,  0 },
          			                 { -1,  0,  0 },
          		                     {  0,  1,  0 }, 
          			                 {  0, -1,  0 },
          			                 {  0,  0,  1 }, 
          			                 {  0,  0, -1 },
          			                 {  1,  1,  0 }, 
          			                 { -1, -1,  0 },
          			                 {  1, -1,  0 },
          			                 { -1,  1,  0 },
          			                 {  1,  0,  1 },
          			                 { -1,  0, -1 },
          			                 {  1,  0, -1 },
          			                 { -1,  0,  1 },
          			                 {  0,  1,  1 },
          			                 {  0, -1, -1 },
          			                 {  0,  1, -1 },
          			                 {  0, -1,  1 } } ;
__constant__ float dt = 1.f;
//----End Header----------------------------------------------------------------




struct nodes {
  int boundary_index;
  float n[19];
  int pos[3];
  float point_to_cent[3];
  float dist_from_cent;
};

struct params {
  int dimx;
  int dimy;
  int dimz;
  int node_count;
  float a_grid;
  int boundary_count;
  int* number_of_anchors;
};

struct boundaries {
  float velocity[3];
  float omega[3];
  float center[3];
  float orientation[3];
  float mass; //TOTAL mass (sum over anchor mass and central mass)
  float inertia_tensor[3][3];
  float radius; //for sphere
  float* anchors; //for concatenated spheres (rigid) they are stored as relative cartesians around center at zero orientation.
  float* anchor_r; //for radii.
  float* anchor_mass; //density may differ between anchor
  int number_of_anchors; // not needed for computation.
  int index; //which boundary this is (0 is fluid, then 1, 2, 3...)
};


//phi from -pi to pi





//host function to init boundary inertia

  
    
//used for the anchor functions 
//FIXME this is completely wrong!
__device__ void rotation_around_center(const float* xyz_in, const float* theta_phi, float* xyz_out){

  float tmpsint = sinf(theta_phi[0]); float tmpcost = cosf(theta_phi[0]);
  float tmpsinp = sinf(theta_phi[1]); float tmpcosp = cosf(theta_phi[0]);
  xyz_out[0] = xyz_in[0] * tmpsint * tmpcosp + xyz_in[1] * tmpcost * tmpcosp - xyz_in[2] * tmpsinp;
  xyz_out[1] = xyz_in[0] * tmpsint * tmpsinp + xyz_in[1] * tmpcost * tmpsinp + xyz_in[2] * tmpcosp;
  xyz_out[2] = xyz_in[0] * tmpcost - xyz_in[1] * tmpsint; 
}


//get inertia. host init function, doubles are necessary for big masses.
//TODO: get Hauptträgheitsachsen for diagonalised inertia tensor
void calc_inertia_tensor(struct boundaries* boundaries){
  
  double x, y, z;
  double delta_x, delta_y, delta_z, dist;
  double xmin, ymin, zmin;
  double xmax, ymax, zmax;
  double xy = 0, xz = 0, yz = 0, x_sq = 0, y_sq = 0, z_sq = 0;
  int hits, flag, ii, jj, kk;
  int precision = 100;
  double density;
  
  //determine xyz min max
  xmin = - boundaries->radius;
  ymin = - boundaries->radius;
  zmin = - boundaries->radius;
  xmax =  boundaries->radius;
  ymax =  boundaries->radius;
  zmax =  boundaries->radius;
#ifdef anchors  
  for(ii = 0; ii < boundaries->number_of_anchors; ii++){
    
    xmin = min(xmin, boundaries->anchors[ii*3 + 0] - boundaries->anchor_r[ii]);
    ymin = min(ymin, boundaries->anchors[ii*3 + 1] - boundaries->anchor_r[ii]);
    zmin = min(zmin, boundaries->anchors[ii*3 + 2] - boundaries->anchor_r[ii]);
    xmax = max(xmax, boundaries->anchors[ii*3 + 0] + boundaries->anchor_r[ii]);
    ymax = max(ymax, boundaries->anchors[ii*3 + 1] + boundaries->anchor_r[ii]);
    zmax = max(zmax, boundaries->anchors[ii*3 + 2] + boundaries->anchor_r[ii]);
  }
#endif
  xmax = (xmax - xmin) / precision; 
  ymax = (ymax - ymin) / precision;
  zmax = (zmax - zmin) / precision;
  
  for(ii = 0, hits = 0; ii < precision; ii++){
    x += xmax;
    y = ymin;
    for(jj = 0; jj < precision; jj++){
      y += ymax;
      z = zmin;
      for(kk = 0; kk < precision; kk++){
        z += zmax;
        flag = 0;
        //determine distance
        delta_x = x;
        delta_y = y;
        delta_z = z;
        
        dist = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
        if(dist < (boundaries->radius * boundaries->radius)){
          hits++;
          flag++;
          xy += x*y; xz += x * z; yz += y * z;
          x_sq += x*x; y_sq += y*y; z_sq += z*z;
        }
#ifdef anchors
        int ll;        
        for(ll = 0; ll < boundaries->number_of_anchors && !flag; ll++){
          //determine distance
          delta_x = x - boundaries->anchors[3 * ll + 0];
          delta_y = y - boundaries->anchors[3 * ll + 1];
          delta_z = z - boundaries->anchors[3 * ll + 2];
          
          dist = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
          if(dist < (boundaries->anchor_r[ll] * boundaries->anchor_r[ll])){
            hits++;
            flag++;
            xy += x*y; xz += x * z; yz += y * z;
            x_sq += x*x; y_sq += y*y; z_sq += z*z;
          }
        }
#endif
      }
    } 
  }
  density = boundaries->mass / hits;
  xy *= density; xz *= density; yz *= density; 
  x_sq *= density; y_sq *= density; z_sq *= density;  
  
  boundaries->inertia_tensor[0][0] = y_sq + z_sq;
  boundaries->inertia_tensor[0][1] = -xy;
  boundaries->inertia_tensor[0][2] = -xz;
  boundaries->inertia_tensor[1][0] = -xy;
  boundaries->inertia_tensor[1][1] = x_sq + z_sq;
  boundaries->inertia_tensor[1][2] = -yz;
  boundaries->inertia_tensor[2][0] = -xz;
  boundaries->inertia_tensor[2][1] = -yz;
  boundaries->inertia_tensor[2][2] = x_sq + y_sq;
     
        
}



__device__ void interpolate_n(struct nodes* l_node, const struct nodes* d_nodes, const int* s_dims){

  // add all neighbors
  short ii;
  short empty_neighbor_count;
#pragma unroll
  for(ii = 1, empty_neighbor_count = 0; ii < 19; ii++){
    short neighbor_index = (l_node->pos[0] + v[ii][0] + s_dims[0]) % s_dims[0]
                   + s_dims[0]* ((l_node->pos[1] + v[ii][1] + s_dims[1]) % s_dims[1])
                   + s_dims[0] * s_dims[1] * ((l_node->pos[2] + v[ii][2] + s_dims[2]) % s_dims[2]);
                   
    if(d_nodes[neighbor_index].boundary_index) empty_neighbor_count++;
    else {
#pragma unroll
      for (short jj = 0; jj < 19; jj++){
        l_node->n[jj] += d_nodes[neighbor_index].n[jj];
      }
    };
  }
  l_node->n[ii] /= 18 - empty_neighbor_count;
}


//this one requires atomics (within shared --> very small penalty) also: this is the smaller part of a thread divergence
//therefore the penalty might not even be noticeable, since the other path takes considerably longer.
__device__ void add_n_to_b(struct nodes* l_node, const struct boundaries* s_boundary, float* s_bound_delta_v){
  
  //create temp float to minimize atomics
  float v_temp[3] = {0, 0, 0};
  float n_per_bound_mass;
#pragma unroll
  for(short ii = 1; ii < 19; ii += 2) {
    n_per_bound_mass = (l_node->n[ii] - l_node->n[ii+1]) / s_boundary->mass;
#pragma unroll
    for(short jj = 0; jj < 3; jj++) {
      v_temp[jj] += n_per_bound_mass * v[ii][jj]; 
    }
  }
  atomicAdd(s_bound_delta_v, v_temp[0]); 
  atomicAdd(s_bound_delta_v + 1, v_temp[1]); 
  atomicAdd(s_bound_delta_v + 2, v_temp[2]);
}
  
  
  
  

//deprecated
__device__ void update_sphere_boundary(
struct nodes* l_node, const struct nodes* d_nodes, const int* s_dims, const struct boundaries* s_boundary, float* s_bound_delta_v){
  
  // calculate minimum image of periodic boundaries using delta as temporary var
  float delta = abs(s_boundary->center[0] - l_node->pos[0]);
  float delta_x_min = min(delta, s_dims[0] - delta);
  
  delta = s_boundary->center[1] - l_node->pos[1];
  float delta_y_min = min(delta, s_dims[1] - delta);
  
  delta = s_boundary->center[2] - l_node->pos[2];
  float delta_z_min = min(delta, s_dims[2] - delta);
  
  float dist = sqrtf(
  	    (delta_x_min)*(delta_x_min)
  	  + (delta_y_min)*(delta_y_min)
      + (delta_z_min)*(delta_z_min));

  //finishing all computation before if-branching to minimize time-penalty
  if(dist <= s_boundary->radius) {
    if(! l_node->boundary_index) add_n_to_b(l_node, s_boundary, s_bound_delta_v);
    l_node->boundary_index = 1;
#pragma unroll
    for(short ii = 0; ii < 19; ii++){
      l_node->n[ii] = 0.f;
    }
  } else {
      if(l_node->boundary_index) interpolate_n(l_node, d_nodes, s_dims);
      l_node->boundary_index = 0;
  }
  __syncthreads(); //debug
}


//deprecated
//__device__ void update_anchor_spheres(
//struct nodes* l_node, const struct nodes* d_nodes,const int* s_dims, const struct boundaries* s_boundary, float* s_bound_delta_v){
//  
//  float delta; float delta_x_min; float delta_y_min; float delta_z_min;
//  float dist;
//  short flag = 0;
//  short ii;
//  float xyz[3];
//  // calculate minimum image of periodic boundaries using delta as temporary var
//  //maybe bugged sphere coordinates.
//  //atomics for anchor??? if anchor stored as real angle and not relative angle
//  //as of now anchor is the relative r theta phi to center. may be subject to change.
//#pragma unroll
//  for(ii = 0; (ii < s_boundary->number_of_anchors) && !flag; ii++){
//  
//    //figure out real xyz of anchor center.
//    rotation_around_center(s_boundary->anchors + 3*ii, s_boundary->angle, xyz);
//    xyz[0] += s_boundary->center[0]; xyz[1] += s_boundary->center[1]; xyz[2] += s_boundary->center[2];
//    
//    delta = abs(xyz[0] - l_node->pos[0]);
//    delta_x_min = min(delta, s_dims[0] - delta);
//  
//    delta = abs(xyz[1] - l_node->pos[1]);
//    delta_y_min = min(delta, s_dims[1] - delta);
//    
//    delta = abs(s_boundary->anchors[3 * ii + 2] - l_node->pos[2]);
//    delta_z_min = min(delta, s_dims[2] - delta);
//  
//    dist = sqrtf(
//  	        (delta_x_min)*(delta_x_min)
//  	      + (delta_y_min)*(delta_y_min)
//          + (delta_z_min)*(delta_z_min));
//    if(dist <= s_boundary->anchor_r[ii]) flag++;
//  }
//  if(flag) {
//    if(! l_node->boundary_index) add_n_to_b(l_node, s_boundary, s_bound_delta_v);
//    l_node->boundary_index = 1;
//#pragma unroll
//    for(short ii = 0; ii < 19; ii++){
//      l_node->n[ii] = 0.f;
//    }
//  } else {
//      if(l_node->boundary_index) interpolate_n(l_node, d_nodes, s_dims);
//      l_node->boundary_index = 0;
//  }
//}

__device__ void update_boundaries(
struct nodes* l_node, const struct nodes* d_nodes, const int* s_dims, const struct boundaries* s_boundary, float* s_bound_delta_v){
  
  float delta; float delta_xyz_min[3]; //for center distance, stored seperately
  float delta_x_min, delta_y_min, delta_z_min; //for anchor distances
  float dist, cdist;
  short ii;
  float xyz[3];//for anchor positions
  // calculate minimum image in periodic boundaries using delta as temporary var
  delta = abs(s_boundary->center[0] - l_node->pos[0]);
  delta_xyz_min[0] = min(delta, s_dims[0] - delta);
  
  delta= abs(s_boundary->center[1] - l_node->pos[1]);
  delta_xyz_min[1] = min(delta, s_dims[1] - delta);
    
  delta = abs(s_boundary->center[2] - l_node->pos[2]);
  delta_xyz_min[2] = min(delta, s_dims[2] - delta);
  
  cdist = sqrtf(
  	      (delta_xyz_min[0])*(delta_xyz_min[0])
        + (delta_xyz_min[1])*(delta_xyz_min[1])
        + (delta_xyz_min[2])*(delta_xyz_min[2]));
        
  if(cdist <= s_boundary->radius){
    if(! l_node->boundary_index) add_n_to_b(l_node, s_boundary, s_bound_delta_v);
    
    l_node->dist_from_cent = cdist;
    l_node->boundary_index = s_boundary->index;
#pragma unroll
    for(ii = 0; ii < 19; ii++){
      l_node->n[ii] = 0.f;
    }
  } else {
      if(l_node->boundary_index) interpolate_n(l_node, d_nodes, s_dims);
      
      l_node->dist_from_cent = 0;
      l_node->boundary_index = 0;
  }
  //doing center sphere first before anchors, since center is generally the largest
#ifdef anchors
  short flag = 0;
#pragma unroll
  for(ii = 0; (ii < s_boundary->number_of_anchors) && !flag; ii++){
  
    //figure out relative xyz of anchor center in real space. (with bound center = 0)
    rotation_around_center(s_boundary->anchors + 3*ii, s_boundary->orientation, xyz);
    //add bound center to xyz to get absolute coordinates
    xyz[0] += s_boundary->center[0]; xyz[1] += s_boundary->center[1]; xyz[2] += s_boundary->center[2];+
    
    //minimum image
    delta = abs(xyz[0] - l_node->pos[0]);
    delta_x_min = min(delta, s_dims[0] - delta);
  
    delta = abs(xyz[1] - l_node->pos[1]);
    delta_y_min = min(delta, s_dims[1] - delta);
    
    delta = abs(xyz[2] - l_node->pos[2]);
    delta_z_min = min(delta, s_dims[2] - delta);
  
    dist = sqrtf(
  	        (delta_x_min)*(delta_x_min)
  	      + (delta_y_min)*(delta_y_min)
          + (delta_z_min)*(delta_z_min));
    if(dist <= s_boundary->anchor_r[ii]) flag++;
  }

  if(flag){
    if(! l_node->boundary_index) add_n_to_b(l_node, s_boundary, s_bound_delta_v);
    l_node->dist_from_cent = cdist;
    l_node->boundary_index = s_boundary->index;
#pragma unroll
    for(ii = 0; ii < 19; ii++){
      l_node->n[ii] = 0.f;
    }
  } else {
      if(l_node->boundary_index) interpolate_n(l_node, d_nodes, s_dims);
      l_node->boundary_index = 0;
  }
#endif //anchors
}      
//__device__ void init_bondary(int index,...)
__global__ void init_grid(
struct boundaries* d_boundaries, struct nodes* d_nodes, const struct params* d_params){

  //globaling and index for testing purposes, will change function to __device__ when finished
  int index = blockDim.x * blockIdx.x + threadIdx.x;
  
  // iterable
  int ii;
  //reading single node into local. TODO: check if more threads are launched than nodes available
  struct nodes l_node = d_nodes[index];
  
  //read boundary data to shared since it is used often
  //must be __shared__, not local, to minimize threadwise penalty for atomic operations.
  // IMPORTANT: SYNCTHREADS BEFORE REWRITING s_boundary!
  __shared__ struct boundaries s_boundary;
  //if (! threadIdx.x) s_boundary = d_boundaries[0];
  
  __shared__ float s_bound_delta_v[3]; //not used in init, but update_boundary requires it
  if(threadIdx.x == 1) {
    s_bound_delta_v[0] = 0; s_bound_delta_v[1] = 0; s_bound_delta_v[2] = 0;
  };
  
  //writing dimxyz to shared to avoid countless __global__ reads
  __shared__ int s_dims[3];//all three needed for periodic boundary condition
  if(threadIdx.x == 2) {
    s_dims[0] = d_params->dimx; 
    s_dims[1] = d_params->dimy; 
    s_dims[2] = d_params->dimz;
  }
  __syncthreads(); //sync after soon-to-be-read write to shared.
  
  //set x, y, z
  //not using modulo, since it may glitch with irregular lattice-sidelengths
  int z = index / (s_dims[0] * s_dims[1]);
  int y = (index - z * s_dims[1] * s_dims[0]) / s_dims[0];
  int x = index - z * s_dims[1] * s_dims[0] - y * s_dims[1]; // (...) / 1
  
  //write to local node
  l_node.pos[0] = x; l_node.pos[1] = y; l_node.pos[2] = z;
  
  for(ii = 0; ii < d_params->boundary_count; ii++){
  
    if (! threadIdx.x) s_boundary = d_boundaries[ii];
    __syncthreads();
    //divergence and atomicadding to delta_v starts now!
    //DO NOT initialize boundaries in grid with non-zero velocity density, n interpolation is active!
    update_sphere_boundary(&l_node, d_nodes, s_dims, &s_boundary, s_bound_delta_v);
    //atomicAdds stop here, it is theoretically possible to sync now
    //do the delta_v add
    __syncthreads();
  }
//set n to some non-zero values arbitrary
#pragma unroll
  for(ii = 0; ii < 19; ii++) {
    if(! l_node.boundary_index) {
      l_node.n[ii] = ii ;
    }
  }
  __syncthreads(); //divergence in update_boundary is fatal now
  
  //rewrite local global
  d_nodes[index] = l_node;

}



__global__ void propagate_boundary(
const struct boundaries* d_boundaries_in, struct boundaries* d_boundaries_out,
const struct nodes* d_nodes_in, struct nodes* d_nodes_out, 
const struct params* d_params){

  //globaling and index for testing purposes, will change function to __device__ when finished
  int index = blockDim.x * blockIdx.x + threadIdx.x;
  
  //iterable
  int ii;
  //reading single node into local
  struct nodes l_node = d_nodes_in[index];
  
  //read boundary data to shared since it is used often
  //must be __shared__, not local, to minimize threadwise penalty for atomic operations.
  // IMPORTANT: SYNCTHREADS BEFORE REWRITING s_boundary!
  __shared__ struct boundaries s_boundary;
  if (! threadIdx.x) s_boundary = *d_boundaries_in;
  
  __shared__ float s_bound_delta_v[3];

  
  //writing dimxyz to shared to avoid countless __global__ reads
  __shared__ int s_dims[3];
  if(threadIdx.x == 2) {
    s_dims[0] = d_params->dimx; 
    s_dims[1] = d_params->dimy; 
    s_dims[2] = d_params->dimz;
  }
  
  //iterate over different boundaries
  for(ii = 0; ii < d_params->boundary_count; ii++){


    if(! threadIdx.x){
      s_boundary = d_boundaries_in[ii];
      short jj;
#pragma unroll
      // the following operation must be done by threadIdx.x 0, since it is the 
      //only thread with guaranteed s_boundary before sync.
      //TODO: proper integration
      for(jj = 0; ii < 3; jj++){
        s_boundary.center[jj] += s_boundary.velocity[jj] * dt;
      }
    }
    if(threadIdx.x == 1) { //must be 1
      s_bound_delta_v[0] = 0; s_bound_delta_v[1] = 0; s_bound_delta_v[2] = 0;
    };
    
    __syncthreads(); //sync after write to shared.
    //divergence and atomicadding to delta_v starts now!
    //This function changes l_node, s_bound_delta_v
    update_sphere_boundary(&l_node, d_nodes_in, s_dims, &s_boundary, s_bound_delta_v);
    //atomicAdds stop here, it is theoretically possible to sync now, 
    //but not necessary until write to global, also hurts efficiency for the next if statement.
    __syncthreads();
    //rewrite shared to global out, sync before necessary.
    if(threadIdx.x == 1) { //must be 1
#pragma unroll
      for(ii = 0; ii < 3; ii++) {
        atomicAdd(d_boundaries_out->velocity + ii, s_bound_delta_v[ii]);
      }
    }
  }
  //stop iteration over all boundaries
  //set global node to new value
  d_nodes_out[index] = l_node;
}


 
int main(int argc, char** argv){

  //some testing stuff
  float center[3] = {10.121f, 10.33f, 10.623f};
  int center_grid[3] = {rint(center[0]), rint(center[1]), rint(center[2])};
  float boundary_displace[3] = {abs(center[0]-center_grid[0]),
  								abs(center[1]-center_grid[1]),
  								abs(center[2]-center_grid[2])};
  

  //iterator
  int ii;  
								
  //init params
  struct params* params;
  params = (struct params*) malloc(sizeof(struct params));
  
  //init number_of_anchors
  int* number_of_anchors;
  number_of_anchors = (int*) malloc(params->boundary_count * sizeof(int));
  
  //set params
  params->dimx = 128; params->dimy = 128; params->dimz = 128;
  params->a_grid = 1.f;
  params->node_count = params->dimx * params->dimy * params->dimz;
  params->boundary_count = 1;
  params->number_of_anchors = number_of_anchors;
  
  //set number_of_anchors (arbitrary atm)
  for(ii = 0; ii < params->boundary_count; ii++){
    number_of_anchors[ii] = ii;
  }
  
  //init boundary
  struct boundaries* boundaries;
  boundaries = (struct boundaries*) malloc(params->boundary_count * sizeof(struct boundaries));
  
  //init anchors
  float** anchorholder_pos;//only used for init, obsolete when boundaries are set
  float** anchorholder_r; //only used for init, obsolete when boundaries are set
  anchorholder_pos = (float**) malloc(params->boundary_count * sizeof(float*));
  anchorholder_r = (float**) malloc(params->boundary_count * sizeof(float*));
  
  for(ii = 0; ii < params->boundary_count; ii++){
    anchorholder_pos[ii] = (float*) malloc(params->number_of_anchors[ii] * 3 * sizeof(float));
    anchorholder_r[ii] = (float*) malloc(params->number_of_anchors[ii] * sizeof(float));
  }
  
  
  //set boundary (arbitrary atm)
  for(ii = 0; ii < params->boundary_count; ii++){
    
    boundaries[ii].center[0] = ii*center[0];
    boundaries[ii].center[1] = ii*center[1];
    boundaries[ii].center[2] = ii*center[2];
    boundaries[ii].mass = 3000.f;
    boundaries[ii].index = ii+1;
    boundaries[ii].radius = 6.1f;
    boundaries[ii].velocity[0] = 0.f;
    boundaries[ii].velocity[1] = 0.f;
    boundaries[ii].velocity[2] = 0.5f;
    boundaries[ii].anchors = anchorholder_pos[ii];
    boundaries[ii].anchor_r = anchorholder_r[ii];
    boundaries[ii].number_of_anchors = params->number_of_anchors[ii];
    //calc_inertia_tensor(boundaries + ii);
    //TODO: uncomment this ^
    //TODO: change r to r squared to save time. (no sqrt for distance calc)
  }
  
  free(anchorholder_r);
  free(anchorholder_pos);
  
  
  //TODO set anchors for testing purposes
  
  
  //copy
  //init nodes
  struct nodes* nodes;
  nodes = (struct nodes*) calloc(params->node_count, sizeof(struct nodes));
  
//debug
 int* zahl = (int*)malloc(sizeof(int));
 *zahl = 0;
 
  // ***begin init&set cuda memory***
  // create host copy of boundaries because I can't think of a better solution to cudacopy the anchors
  struct boundaries* boundaries_copy;
  boundaries_copy = (struct boundaries*) malloc(params->boundary_count * sizeof(struct boundaries));
  memcpy(boundaries_copy, boundaries, params->boundary_count * sizeof(struct boundaries));
  
  //set d_anchors and put them into the boundaries_copy
  float* d_anchors;
  float* d_anchor_r;
  float* d_anchor_mass;
  #pragma unroll
  for(ii = 0; ii < params->boundary_count; ii++){
  
    checkCudaErrors(cudaMalloc(&d_anchors, 3 * params->number_of_anchors[ii] * sizeof(float)));
    checkCudaErrors(cudaMemcpy(d_anchors, boundaries[ii].anchors, 3 * params->number_of_anchors[ii] * sizeof(float), cudaMemcpyHostToDevice));
    
    checkCudaErrors(cudaMalloc(&d_anchor_r, params->number_of_anchors[ii] * sizeof(float)));
    checkCudaErrors(cudaMemcpy(d_anchor_r, boundaries[ii].anchor_r, params->number_of_anchors[ii] * sizeof(float), cudaMemcpyHostToDevice));
    
    checkCudaErrors(cudaMalloc(&d_anchor_mass, params->number_of_anchors[ii] * sizeof(float)));
    checkCudaErrors(cudaMemcpy(d_anchor_mass, boundaries[ii].anchor_r, params->number_of_anchors[ii] * sizeof(float), cudaMemcpyHostToDevice));
    
    boundaries_copy[ii].anchors = d_anchors; 
    boundaries_copy[ii].anchor_r = d_anchor_r;
    boundaries_copy[ii].anchor_mass = d_anchors;
  }
  
  //cudacopy the modified boundaries
  struct boundaries* d_boundaries;
  checkCudaErrors(cudaMalloc(&d_boundaries, 2 * params->boundary_count * sizeof(struct boundaries)));
  checkCudaErrors(cudaMemcpy(d_boundaries, boundaries_copy, params->boundary_count * sizeof(struct boundaries), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_boundaries + params->boundary_count, boundaries_copy, params->boundary_count * sizeof(struct boundaries), cudaMemcpyHostToDevice));
  
  free(boundaries_copy);
  //double buffering of nodes
  struct nodes* d_nodes;
  checkCudaErrors(cudaMalloc(&d_nodes, 2*sizeof(struct nodes)*params->node_count));
  checkCudaErrors(cudaMemcpy(d_nodes, nodes, sizeof(struct nodes)*params->node_count, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_nodes + params->node_count, nodes, sizeof(struct nodes)*params->node_count, cudaMemcpyHostToDevice));
  
  //cudacopy params
  struct params* params_copy;
  params_copy = (struct params*) malloc(sizeof(struct params));
  memcpy(params_copy, params, sizeof(struct params));
  int* d_number_of_anchors;
  checkCudaErrors(cudaMalloc(&d_number_of_anchors, params->boundary_count * sizeof(int)));
  checkCudaErrors(cudaMemcpy(d_number_of_anchors, params->number_of_anchors, params->boundary_count * sizeof(int), cudaMemcpyHostToDevice));
  params_copy->number_of_anchors = d_number_of_anchors;
  
  struct params* d_params;
  checkCudaErrors(cudaMalloc(&d_params, sizeof(struct params)));
  checkCudaErrors(cudaMemcpy(d_params, params_copy, sizeof(struct params), cudaMemcpyHostToDevice));
  

  // **end init&set cuda memory
  
//  // Do the thing
//  init_grid<<<params->node_count / 64, 64>>>(d_boundaries, d_nodes, d_params);
//  cudaDeviceSynchronize(); //for GetLastError
//  checkCudaErrors(cudaGetLastError());


//  //update double buffer to transfer init data.
//  checkCudaErrors(cudaMemcpy(d_nodes + params->node_count, d_nodes, sizeof(struct nodes)*params->node_count, cudaMemcpyDeviceToDevice));
//  
//  //do the other thing. d_nodes input, d_nodes + params->node_count as output.
//  propagate_boundary<<<params->node_count / 64, 64>>>(d_boundaries, d_boundaries + params->boundary_count, d_nodes, d_nodes + params->node_count, d_params);
//  cudaDeviceSynchronize();//for getlasterror
//  checkCudaErrors(cudaGetLastError());
  
  //do the back and forth, now with d_nodes + params->node_count as input and d_nodes as output.

  
  
  
  

  
  
  
  
  
  
  
  
  
    //output
    

    

//  checkCudaErrors(cudaMemcpy(nodes, d_nodes + params->node_count, sizeof(struct nodes)*params->node_count, cudaMemcpyDeviceToHost));
//  int i;
//  printf("\n");
//  
//  for(i = 0; i < params->node_count/100; i++){
//  
//    printf("%d ", nodes[i* 100].boundary_index);
//    if(!(i % 30)) printf("\n");
//  }
//  

//  checkCudaErrors(cudaMemcpy(boundaries, d_boundaries, 2*sizeof(struct boundaries), cudaMemcpyDeviceToHost));
//  printf("\n");
//  printf("%f %f %f", boundaries[1].velocity[0], boundaries[1].velocity[1], boundaries[1].velocity[2]);
  
  return 0;
}  
