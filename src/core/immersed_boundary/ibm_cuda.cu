
// *******
// This is an internal file of the IMMERSED BOUNDARY implementation
// It should not be included by any main Espresso routines
// Functions to be exported for Espresso are in ibm_main.hpp

#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include "lbgpu.hpp"
#include "lb-boundaries.hpp"
#include "particle_data.hpp"
#include "cuda_utils.hpp"
#include "cuda_interface.hpp"
#include "immersed_boundary/ibm_main.hpp"
#include "immersed_boundary/ibm_cuda_interface.hpp"

// ****** Kernel functions for internal use ********
__global__ void ResetLBForces_Kernel(LB_node_force_gpu node_f, const LB_parameters_gpu *const paraP);
__global__ void ParticleVelocitiesFromLB_Kernel(LB_nodes_gpu n_curr, const IBM_CUDA_ParticleDataInput *const particles_input, IBM_CUDA_ParticleDataOutput *const particles_output, LB_node_force_gpu node_f, const float *const lb_boundary_velocity, const LB_parameters_gpu *const para);
__global__ void ForcesIntoFluid_Kernel(const IBM_CUDA_ParticleDataInput *const particle_input, LB_node_force_gpu node_f, const LB_parameters_gpu *const paraP);
__device__ inline void atomicadd( float* address,float value);

// ***** Other functions for internal use *****
void InitCUDA_IBM(const int numParticles);

// ***** Our own global variables ********
IBM_CUDA_ParticleDataInput *IBM_ParticleDataInput_device = NULL;
IBM_CUDA_ParticleDataOutput *IBM_ParticleDataOutput_device = NULL;

// ****** These variables are defined in lbgpu_cuda.cu, but we also want them here ****
extern LB_node_force_gpu node_f;
extern LB_nodes_gpu *current_nodes;

// ** These variables are static in lbgpu_cuda.cu, so we need to duplicate them here
// They are initialized in ForcesIntoFluid
// The pointers are on the host, but point into device memory
LB_parameters_gpu *paraIBM = NULL;
float* lb_boundary_velocity_IBM = NULL;

/****************
   IBM_ResetLBForces_GPU
Calls a kernel to reset the forces on the LB nodes to the external force
*****************/

void IBM_ResetLBForces_GPU()
{
  if ( this_node == 0 )
  {
    // Setup for kernel call
    int threads_per_block = 64;
    int blocks_per_grid_y = 4;
    int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
    dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

    KERNELCALL(ResetLBForces_Kernel, dim_grid, threads_per_block, (node_f, paraIBM));
  }
}

/******************
   IBM_ForcesIntoFluid_GPU
Called from integrate_vv to put the forces into the fluid
This must be the first CUDA-IBM function to be called because it also does some initialization
*******************/

void IBM_ForcesIntoFluid_GPU()
{
  // This function does
  // (1) Gather forces from all particles via MPI
  // (2) Copy forces to the GPU
  // (3) interpolate on the LBM grid and spread forces

  const int numParticles = gpu_get_global_particle_vars_pointer_host()->number_of_particles;

  // Storage only needed on master and allocated only once at the first time step
  if ( IBM_ParticleDataInput_host == NULL && this_node == 0 )
    InitCUDA_IBM(numParticles);

  // We gather particle positions and forces from all nodes
  IBM_cuda_mpi_get_particles();


  // ***** GPU stuff only on master *****
  if ( this_node == 0 && numParticles > 0)
  {

    // Copy data to device
    cuda_safe_mem(cudaMemcpy(IBM_ParticleDataInput_device, IBM_ParticleDataInput_host,  numParticles*sizeof(IBM_CUDA_ParticleDataInput), cudaMemcpyHostToDevice));

    // Kernel call for spreading the forces on the LB grid
    int threads_per_block_particles = 64;
    int blocks_per_grid_particles_y = 4;
    int blocks_per_grid_particles_x = (lbpar_gpu.number_of_particles + threads_per_block_particles * blocks_per_grid_particles_y - 1) / (threads_per_block_particles * blocks_per_grid_particles_y);
    dim3 dim_grid_particles = make_uint3(blocks_per_grid_particles_x, blocks_per_grid_particles_y, 1);

    KERNELCALL(ForcesIntoFluid_Kernel, dim_grid_particles, threads_per_block_particles, (IBM_ParticleDataInput_device, node_f, paraIBM));

  }
}

/***************
   InitCUDA_IBM
***************/

void InitCUDA_IBM(const int numParticles)
{

  if ( this_node == 0)    // GPU only on master
  {

    // Back and forth communication of positions and velocities
    IBM_ParticleDataInput_host = new IBM_CUDA_ParticleDataInput[numParticles];
    cuda_safe_mem(cudaMalloc((void**)&IBM_ParticleDataInput_device, numParticles*sizeof(IBM_CUDA_ParticleDataInput)));
    cuda_safe_mem(cudaMalloc((void**)&IBM_ParticleDataOutput_device, numParticles*sizeof(IBM_CUDA_ParticleDataOutput)));
    IBM_ParticleDataOutput_host = new IBM_CUDA_ParticleDataOutput[numParticles];


    // Copy parameters to the GPU
    LB_parameters_gpu *para_gpu;
    lb_get_para_pointer(&para_gpu);
    cuda_safe_mem(cudaMalloc((void**)&paraIBM, sizeof(LB_parameters_gpu)));
    cuda_safe_mem(cudaMemcpy(paraIBM, para_gpu, sizeof(LB_parameters_gpu), cudaMemcpyHostToDevice));

    // Copy boundary velocities to the GPU
    // First put them into correct format
#ifdef LB_BOUNDARIES_GPU
    float* host_lb_boundary_velocity = new float[3*(n_lb_boundaries+1)];

    for (int n=0; n<n_lb_boundaries; n++)
    {
      host_lb_boundary_velocity[3*n+0]=lb_boundaries[n].velocity[0];
      host_lb_boundary_velocity[3*n+1]=lb_boundaries[n].velocity[1];
      host_lb_boundary_velocity[3*n+2]=lb_boundaries[n].velocity[2];
    }

    host_lb_boundary_velocity[3*n_lb_boundaries+0] = 0.0f;
    host_lb_boundary_velocity[3*n_lb_boundaries+1] = 0.0f;
    host_lb_boundary_velocity[3*n_lb_boundaries+2] = 0.0f;
    cuda_safe_mem(cudaMalloc((void**)&lb_boundary_velocity_IBM, 3*n_lb_boundaries*sizeof(float)));
    cuda_safe_mem(cudaMemcpy(lb_boundary_velocity_IBM, host_lb_boundary_velocity, 3*n_lb_boundaries*sizeof(float), cudaMemcpyHostToDevice));
    delete[] host_lb_boundary_velocity;
#endif
  }
}

/**************
   Calc_m_from_n_IBM
This is our own version of the calc_m_from_n function in lbgpu_cuda.cu
It does exactly the same, but calculates only the first four modes
***************/

__device__ void Calc_m_from_n_IBM(const LB_nodes_gpu n_a, const unsigned int index, float *mode, const LB_parameters_gpu *const paraP)
{
  const LB_parameters_gpu &para = *paraP;
  const int ii = 0;
  // mass mode
  mode[0 + ii * LBQ] =   n_a.vd[( 0 + ii*LBQ ) * para.number_of_nodes + index]
  + n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index]
  + n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index]
  + n_a.vd[( 5 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 6 + ii*LBQ ) * para.number_of_nodes + index]
  + n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index]
  + n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index]
  + n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index]
  + n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index]
  + n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index]
  + n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index];

  // momentum modes

  mode[1 + ii * LBQ] =   (n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index])
  + (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
  + (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
  + (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
  + (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index]);

  mode[2 + ii * LBQ] =   (n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index])
  + (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
  - (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
  + (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
  + (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

  mode[3 + ii * LBQ] =   (n_a.vd[( 5 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 6 + ii*LBQ ) * para.number_of_nodes + index])
  + (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
  - (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index])
  + (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
  - (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);
}


/**************
   ParticleVelocitiesFromLB_GPU
Calls a kernel function to interpolate the velocity at each IBM particle's position
Store velocity in the particle data structure
**************/

void ParticleVelocitiesFromLB_GPU()
{
  // This function performs three steps:
  // (1) interpolate velocities on GPU
  // (2) transfer velocities back to CPU
  // (3) spread velocities to local cells via MPI


  const int numParticles = gpu_get_global_particle_vars_pointer_host()->number_of_particles;

  // **** GPU stuff only on master ****
  if (this_node == 0 && numParticles > 0)
  {
    // Kernel call
    int threads_per_block_particles = 64;
    int blocks_per_grid_particles_y = 4;
    int blocks_per_grid_particles_x = (lbpar_gpu.number_of_particles + threads_per_block_particles * blocks_per_grid_particles_y - 1) / (threads_per_block_particles * blocks_per_grid_particles_y);
    dim3 dim_grid_particles = make_uint3(blocks_per_grid_particles_x, blocks_per_grid_particles_y, 1);
    KERNELCALL(ParticleVelocitiesFromLB_Kernel, dim_grid_particles, threads_per_block_particles, (*current_nodes, IBM_ParticleDataInput_device, IBM_ParticleDataOutput_device, node_f, lb_boundary_velocity_IBM, paraIBM));


    // Copy velocities from device to host
    cuda_safe_mem(cudaMemcpy(IBM_ParticleDataOutput_host, IBM_ParticleDataOutput_device,  numParticles*sizeof(IBM_CUDA_ParticleDataOutput), cudaMemcpyDeviceToHost));

  }

  // ***** Back to all nodes ****
  // Spread using MPI
  IBM_cuda_mpi_send_velocities();

}

/***************
   ForcesIntoFluid_Kernel
****************/

__global__ void ForcesIntoFluid_Kernel(const IBM_CUDA_ParticleDataInput *const particle_input, LB_node_force_gpu node_f, const LB_parameters_gpu *const paraP)
{
  const unsigned int particleIndex = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  const LB_parameters_gpu &para = *paraP;

  if (particleIndex < para.number_of_particles && particle_input[particleIndex].isVirtual)
  {

    const float factor = powf( para.agrid,2)*para.tau*para.tau;
    const float particleForce[3] = { particle_input[particleIndex].f[0] * factor, particle_input[particleIndex].f[1] * factor, particle_input[particleIndex].f[2] * factor};
    const float pos[3] = { particle_input[particleIndex].pos[0], particle_input[particleIndex].pos[1], particle_input[particleIndex].pos[2] };

    // First part is the same as for interpolation --> merge into a single function
    float temp_delta[6];
    float delta[8];
    int my_left[3];
    int node_index[8];
    for(int i=0; i<3; ++i)
    {
      const float scaledpos = pos[i]/para.agrid - 0.5f;
      my_left[i] = (int)(floorf(scaledpos));
      temp_delta[3+i] = scaledpos - my_left[i];
      temp_delta[i] = 1.f - temp_delta[3+i];
    }

    delta[0] = temp_delta[0] * temp_delta[1] * temp_delta[2];
    delta[1] = temp_delta[3] * temp_delta[1] * temp_delta[2];
    delta[2] = temp_delta[0] * temp_delta[4] * temp_delta[2];
    delta[3] = temp_delta[3] * temp_delta[4] * temp_delta[2];
    delta[4] = temp_delta[0] * temp_delta[1] * temp_delta[5];
    delta[5] = temp_delta[3] * temp_delta[1] * temp_delta[5];
    delta[6] = temp_delta[0] * temp_delta[4] * temp_delta[5];
    delta[7] = temp_delta[3] * temp_delta[4] * temp_delta[5];

    // modulo for negative numbers is strange at best, shift to make sure we are positive
    const int x = my_left[0] + para.dim_x;
    const int y = my_left[1] + para.dim_y;
    const int z = my_left[2] + para.dim_z;

    node_index[0] = x%para.dim_x     + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*(z%para.dim_z);
    node_index[1] = (x+1)%para.dim_x + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*(z%para.dim_z);
    node_index[2] = x%para.dim_x     + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*(z%para.dim_z);
    node_index[3] = (x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*(z%para.dim_z);
    node_index[4] = x%para.dim_x     + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*((z+1)%para.dim_z);
    node_index[5] = (x+1)%para.dim_x + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*((z+1)%para.dim_z);
    node_index[6] = x%para.dim_x     + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z);
    node_index[7] = (x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z);

    for(int i=0; i<8; ++i)
    {

      // Atomic add is essential because this runs in parallel!
      atomicadd(&(node_f.force[0*para.number_of_nodes + node_index[i]]), (particleForce[0] * delta[i]));
      atomicadd(&(node_f.force[1*para.number_of_nodes + node_index[i]]), (particleForce[1] * delta[i]));
      atomicadd(&(node_f.force[2*para.number_of_nodes + node_index[i]]), (particleForce[2] * delta[i]));
    }
  }
}


/**************
   ParticleVelocitiesFromLB_Kernel
**************/

__global__ void  ParticleVelocitiesFromLB_Kernel(LB_nodes_gpu n_curr, const IBM_CUDA_ParticleDataInput *const particles_input, IBM_CUDA_ParticleDataOutput *const particles_output, LB_node_force_gpu node_f, const float *const lb_boundary_velocity, const LB_parameters_gpu *const paraP)
{

  const unsigned int particleIndex = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  const LB_parameters_gpu &para = *paraP;

  if (particleIndex < para.number_of_particles && particles_input[particleIndex].isVirtual)
  {

    // Get position
    float pos[3] = { particles_input[particleIndex].pos[0], particles_input[particleIndex].pos[1], particles_input[particleIndex].pos[2] };
    float v[3] = { 0 };

    // ***** This part is copied from get_interpolated_velocity
    // ***** + we add the force + we consider boundaries - we remove the Shan-Chen stuff

    float temp_delta[6];
    float delta[8];
    int my_left[3];
    int node_index[8];
    float mode[4];
    #pragma unroll
    for(int i=0; i<3; ++i)
    {
      const float scaledpos = pos[i]/para.agrid - 0.5f;
      my_left[i] = (int)(floorf(scaledpos));
      temp_delta[3+i] = scaledpos - my_left[i];
      temp_delta[i] = 1.f - temp_delta[3+i];
    }

    delta[0] = temp_delta[0] * temp_delta[1] * temp_delta[2];
    delta[1] = temp_delta[3] * temp_delta[1] * temp_delta[2];
    delta[2] = temp_delta[0] * temp_delta[4] * temp_delta[2];
    delta[3] = temp_delta[3] * temp_delta[4] * temp_delta[2];
    delta[4] = temp_delta[0] * temp_delta[1] * temp_delta[5];
    delta[5] = temp_delta[3] * temp_delta[1] * temp_delta[5];
    delta[6] = temp_delta[0] * temp_delta[4] * temp_delta[5];
    delta[7] = temp_delta[3] * temp_delta[4] * temp_delta[5];

    // modulo for negative numbers is strange at best, shift to make sure we are positive
    int x = my_left[0] + para.dim_x;
    int y = my_left[1] + para.dim_y;
    int z = my_left[2] + para.dim_z;

    node_index[0] = x%para.dim_x     + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*(z%para.dim_z);
    node_index[1] = (x+1)%para.dim_x + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*(z%para.dim_z);
    node_index[2] = x%para.dim_x     + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*(z%para.dim_z);
    node_index[3] = (x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*(z%para.dim_z);
    node_index[4] = x%para.dim_x     + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*((z+1)%para.dim_z);
    node_index[5] = (x+1)%para.dim_x + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*((z+1)%para.dim_z);
    node_index[6] = x%para.dim_x     + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z);
    node_index[7] = (x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z);

    for(int i=0; i<8; ++i)
    {
       double local_rho;
       double local_j[3];
#ifdef LB_BOUNDARIES_GPU
      if (n_curr.boundary[node_index[i]])
      {
        // Boundary node
        const int boundary_index= n_curr.boundary[node_index[i]];

        // lb_boundary_velocity is given in MD units --> convert to LB and reconvert back at the end of this function
        local_rho = para.rho[0]*para.agrid*para.agrid*para.agrid;
        local_j[0] = para.rho[0]*para.agrid*para.agrid*para.agrid*lb_boundary_velocity[3*(boundary_index-1)+0];
        local_j[1] = para.rho[0]*para.agrid*para.agrid*para.agrid*lb_boundary_velocity[3*(boundary_index-1)+1];
        local_j[2] = para.rho[0]*para.agrid*para.agrid*para.agrid*lb_boundary_velocity[3*(boundary_index-1)+2];

      }
      else
#endif
      {
        Calc_m_from_n_IBM(n_curr,node_index[i],mode, paraP);
        local_rho = para.rho[0]*para.agrid*para.agrid*para.agrid + mode[0];

        // Add the +f/2 contribution!!
        local_j[0] = mode[1] + node_f.force_buf[0*para.number_of_nodes + node_index[i]]/2.f;
        local_j[1] = mode[2] + node_f.force_buf[1*para.number_of_nodes + node_index[i]]/2.f;
        local_j[2] = mode[3] + node_f.force_buf[2*para.number_of_nodes + node_index[i]]/2.f;

      }


      // Interpolate velocity
      v[0] += delta[i]*local_j[0]/(local_rho);
      v[1] += delta[i]*local_j[1]/(local_rho);
      v[2] += delta[i]*local_j[2]/(local_rho);
    }


    // Rescale and store output
    particles_output[particleIndex].v[0] = v[0] * para.agrid / para.tau;
    particles_output[particleIndex].v[1] = v[1] * para.agrid / para.tau;
    particles_output[particleIndex].v[2] = v[2] * para.agrid / para.tau;

  }
}


/****************
   ResetLBForces_Kernel
*****************/

__global__ void ResetLBForces_Kernel(LB_node_force_gpu node_f, const LB_parameters_gpu *const paraP)
{

  const size_t index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  const LB_parameters_gpu &para = *paraP;

  if ( index < para.number_of_nodes )
  {
    const float force_factor=powf(para.agrid,2)*para.tau*para.tau;
    for(int ii=0;ii<LB_COMPONENTS;++ii)
    {
#ifdef EXTERNAL_FORCES
      if(para.external_force)
      {
        node_f.force[(0 + ii*3 ) * para.number_of_nodes + index] = para.ext_force[0 + ii*3 ]*force_factor;
        node_f.force[(1 + ii*3 ) * para.number_of_nodes + index] = para.ext_force[1 + ii*3 ]*force_factor;
        node_f.force[(2 + ii*3 ) * para.number_of_nodes + index] = para.ext_force[2 + ii*3 ]*force_factor;
      }
      else
#endif
      {
        node_f.force[(0 + ii*3 ) * para.number_of_nodes + index] = 0.0f;
        node_f.force[(1 + ii*3 ) * para.number_of_nodes + index] = 0.0f;
        node_f.force[(2 + ii*3 ) * para.number_of_nodes + index] = 0.0f;
      }
    }
  }
}

/************
   atomicadd
This is a copy of the atomicadd from lbgpu_cuda.cu
It seems that this function is re-implemented in all CUDA files: lbgpu_cuda, electrokinetics_cuda,...
Why is it not in a header?
*************/

__device__ inline void atomicadd( float* address, float value)
{

  #if !defined __CUDA_ARCH__ || __CUDA_ARCH__ >= 200 // for Fermi, atomicAdd supports floats
  atomicAdd(address, value);
  #elif __CUDA_ARCH__ >= 110

  #warning Using slower atomicAdd emulation

  //float-atomic-add from
  //[url="http://forums.nvidia.com/index.php?showtopic=158039&view=findpost&p=991561"]

  float old = value;
  while( ( old = atomicExch( address, atomicExch( address, 0.0f ) + old ) ) != 0.0f );

  #else
    #error CUDA compute capability 1.1 or higher required
  #endif
}


#endif
