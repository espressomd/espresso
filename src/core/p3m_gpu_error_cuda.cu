#include <cuda.h>
#include <thrust/reduce.h>
#include <thrust/device_vector.h>
#include <stdio.h>

#include "p3m_gpu_error.hpp"
#include "cuda_utils.hpp"

#define PI 3.14159265358979323846264338328

#define SQR(A) ((A)*(A))

/** @TODO: Extend to hight order. This comes from some 1/sin expansion in Hockney/Eastwood */

__device__ static double p3m_analytic_cotangent_sum(int n, double mesh_i, int cao)
{
  double c, res=0.0;
  c = SQR(cos(PI*mesh_i*(double)n));

  switch (cao) {
  case 1 : { 
    res = 1; 
    break; }
  case 2 : { 
    res = (1.0+c*2.0)/3.0; 
    break; }
  case 3 : { 
    res = (2.0+c*(11.0+c*2.0))/15.0; 
    break; }
  case 4 : { 
    res = (17.0+c*(180.0+c*(114.0+c*4.0)))/315.0; 
    break; }
  case 5 : { 
    res = (62.0+c*(1072.0+c*(1452.0+c*(247.0+c*2.0))))/2835.0; 
    break; }
  case 6 : { 
    res = (1382.0+c*(35396.0+c*(83021.0+c*(34096.0+c*(2026.0+c*4.0)))))/155925.0; 
    break; }
  case 7 : { 
    res = (21844.0+c*(776661.0+c*(2801040.0+c*(2123860.0+c*(349500.0+c*(8166.0+c*4.0))))))/6081075.0; 
    break; }
  }
  
  return res;
}


__device__ static inline double csinc(double d)
{
#define epsi 0.1

#define c2 -0.1666666666667e-0
#define c4  0.8333333333333e-2
#define c6 -0.1984126984127e-3
#define c8  0.2755731922399e-5

  double PId = PI*d, PId2;

  if (fabs(d)>epsi)
    return sin(PId)/PId;
  else {
    PId2 = SQR(PId);
    return 1.0 + PId2*(c2+PId2*(c4+PId2*(c6+PId2*c8)));
  }
}


__global__ void p3m_k_space_error_gpu_kernel_ik(int3 mesh, double3 meshi, int cao, double alpha_L, double * he_q)
{
  int nx = -mesh.x/2 + blockDim.x*blockIdx.x + threadIdx.x;
  int ny = -mesh.y/2 + blockDim.y*blockIdx.y + threadIdx.y;
  int nz = -mesh.z/2 + blockDim.z*blockIdx.z + threadIdx.z;

  if( (nx >= mesh.x/2) || (ny >= mesh.y/2) || (nz >= mesh.z/2))
    return;

  int lind = ( (nx + mesh.x/2) * mesh.y*mesh.z + (ny + mesh.y/2) * mesh.z + (nz + mesh.z/2));

  double alpha_L_i = 1./alpha_L;
  double n2, cs;
  double U2, ex, ex2;

  if((nx!=0) || (ny!=0) || (nz!=0)) {
    n2 = SQR(nx) + SQR(ny) + SQR(nz);

    cs = p3m_analytic_cotangent_sum(nz,meshi.z,cao)*p3m_analytic_cotangent_sum(nx,meshi.x,cao) * p3m_analytic_cotangent_sum(ny,meshi.y,cao);

    ex = exp(-(PI*alpha_L_i)*(PI*alpha_L_i)*n2);

    ex2 = SQR( ex );

    U2 = pow((double)csinc(meshi.x*nx)*csinc(meshi.y*ny)*csinc(meshi.z*nz), 2.0*cao);

    he_q[lind] = ex2/n2  -  SQR(U2*ex/cs) / n2;
    
  } else {
    he_q[lind] = 0;
  }
}

__global__ void p3m_k_space_error_gpu_kernel_ad(int3 mesh, double3 meshi, int cao, double alpha_L, double * he_q)
{


  int  nx = -mesh.x/2 + blockDim.x*blockIdx.x + threadIdx.x;
  int  ny = -mesh.y/2 + blockDim.y*blockIdx.y + threadIdx.y;
  int  nz = -mesh.z/2 + blockDim.z*blockIdx.z + threadIdx.z;

  if( (nx >= mesh.x/2) || (ny >= mesh.y/2) || (nz >= mesh.z/2))
    return;

  int lind = ( (nx + mesh.x/2) * mesh.y*mesh.z + (ny + mesh.y/2) * mesh.z + (nz + mesh.z/2));

  double alpha_L_i = 1./alpha_L;
  double n2;
  double U2, ex, ex2;
  int nmx, nmy, nmz;
  double alias1, alias2, alias3, alias4;

  alias1 = alias2 = alias3 = alias4 = 0;

  if((nx!=0) || (ny!=0) || (nz!=0)) {
    for(int mx = -1; mx <= 1; mx++) {
      nmx = nx + mx*mesh.x;
      for(int my = -1; my <= 1; my++) {
	nmy = ny + my*mesh.y;
	for(int mz = -1; mz <= 1; mz++) {
	  nmz = nz + mz*mesh.z;

	  n2 = SQR(nmx) + SQR(nmy) + SQR(nmz);

	  ex = exp(-(PI*alpha_L_i)*(PI*alpha_L_i)*n2);

	  ex2 = SQR( ex );

	  U2 = pow((double)csinc(meshi.x*nmx)*csinc(meshi.y*nmy)*csinc(meshi.z*nmz), 2.0*cao);

	  alias1 += ex2 / n2;
	  alias2 += U2 * ex;
	  alias3 += U2 * n2;
	  alias4 += U2;

	}
      }
    }

    if( (alias3 == 0.0) || (alias4 == 0.0) )
      he_q[lind] = 0;
    else
      he_q[lind] = alias1 - (alias2*alias2)/(alias3*alias4);
    
  } else {
    he_q[lind] = 0;
  }
}

__global__ void p3m_k_space_error_gpu_kernel_ik_i(int3 mesh, double3 meshi, int cao, double alpha_L, double * he_q)
{


  int  nx = -mesh.x/2 + blockDim.x*blockIdx.x + threadIdx.x;
  int  ny = -mesh.y/2 + blockDim.y*blockIdx.y + threadIdx.y;
  int  nz = -mesh.z/2 + blockDim.z*blockIdx.z + threadIdx.z;

  if( (nx >= mesh.x/2) || (ny >= mesh.y/2) || (nz >= mesh.z/2))
    return;

  int lind = ( (nx + mesh.x/2) * mesh.y*mesh.z + (ny + mesh.y/2) * mesh.z + (nz + mesh.z/2));

  double alpha_L_i = 1./alpha_L;
  double n2;
  double U2, ex, ex2;
  int nmx, nmy, nmz;
  double alias1, alias2, alias3, alias4;

  alias1 = alias2 = alias3 = alias4 = 0;

  if((nx!=0) || (ny!=0) || (nz!=0)) {
    for(int mx = -1; mx <= 1; mx++) {
      nmx = nx + mx*mesh.x;
      for(int my = -1; my <= 1; my++) {
	nmy = ny + my*mesh.y;
	for(int mz = -1; mz <= 1; mz++) {
	  nmz = nz + mz*mesh.z;

	  n2 = SQR(nmx) + SQR(nmy) + SQR(nmz);

	  ex = exp(-(PI*alpha_L_i)*(PI*alpha_L_i)*n2);

	  ex2 = SQR( ex );

	  U2 = pow((double)csinc(meshi.x*nmx)*csinc(meshi.y*nmy)*csinc(meshi.z*nmz), 2.0*cao);

	  alias1 += ex2 / n2;
	  alias2 += U2 * ex  * (nx*nmx + ny*nmy + nz*nmz) / n2;
	  alias3 += U2;

	  if (((mx+my+mz)%2)==0) {		//consider only even terms!
	    alias4 += U2;
	  } else {
	    alias4 -= U2;
	  }
	}
      }
    }

    he_q[lind] = alias1 - (alias2*alias2)/(0.5*(nx*nx + ny*ny + nz*nz)*(alias3*alias3 + alias4*alias4));
    
  } else {
    he_q[lind] = 0;
  }
}



__global__ void p3m_k_space_error_gpu_kernel_ad_i(int3 mesh, double3 meshi, int cao, double alpha_L, double * he_q)
{


  int  nx = -mesh.x/2 + blockDim.x*blockIdx.x + threadIdx.x;
  int  ny = -mesh.y/2 + blockDim.y*blockIdx.y + threadIdx.y;
  int  nz = -mesh.z/2 + blockDim.z*blockIdx.z + threadIdx.z;

  if( (nx >= mesh.x/2) || (ny >= mesh.y/2) || (nz >= mesh.z/2))
    return;

  int lind = ( (nx + mesh.x/2) * mesh.y*mesh.z + (ny + mesh.y/2) * mesh.z + (nz + mesh.z/2));

  double alpha_L_i = 1./alpha_L;
  double n2;
  double U2, ex, ex2;
  int nmx, nmy, nmz;
  double alias1, alias2, alias3, alias4, alias5, alias6;

  alias1 = alias2 = alias3 = alias4 = alias5 = alias6 = 0;

  if((nx!=0) && (ny!=0) && (nz!=0)) {
    for(int mx = -1; mx <= 1; mx++) {
      nmx = nx + mx*mesh.x;
      for(int my = -1; my <= 1; my++) {
	nmy = ny + my*mesh.y;
	for(int mz = -1; mz <= 1; mz++) {
	  nmz = nz + mz*mesh.z;

	  n2 = SQR(nmx) + SQR(nmy) + SQR(nmz);

	  ex = exp(-(PI*alpha_L_i)*(PI*alpha_L_i)*n2);

	  ex2 = SQR( ex );

	  U2 = pow((double)csinc(meshi.x*nmx)*csinc(meshi.y*nmy)*csinc(meshi.z*nmz), 2.0*cao);

	  alias1 += ex2 / n2;
	  alias2 += U2 * ex;
	  alias3 += U2 * n2;
	  alias4 += U2;

        if (((mx+my+mz)%2)==0) {					//even term
	   alias5 += U2 * n2;
	   alias6 += U2;
	 } else {						//odd term: minus sign!
	   alias5 -= U2 * n2;
	   alias6 -= U2;
	 }


	}
      }
    }

    he_q[lind] = (alias1  -  SQR(alias2) / (0.5*(alias3*alias4 + alias5*alias6)));
    
  } else {
    he_q[lind] = 0;
  }
}


static double *he_q = 0;
static size_t he_q_size = 0;

void p3m_gpu_error_reset() {
  if(he_q != 0)
    cudaFree(he_q);
  he_q_size = 0;
}

double p3m_k_space_error_gpu(double prefactor, int *mesh, int cao, int npart, double sum_q2, double alpha_L, double *box) {
  double he_q_final;
  size_t mesh_size = mesh[0]*mesh[1]*mesh[2];

  if(mesh_size >= he_q_size) {
    p3m_gpu_error_reset();
    cudaMalloc(&he_q, mesh_size*sizeof(double));
    he_q_size = mesh_size;
  }
  
  dim3 grid(max(1, mesh[0]/8 + 1),max(1, mesh[1]/8 + 1),max(1, mesh[2]/8 + 1));
  dim3 block(8,8,8);
  
  int3 mesh3;
  mesh3.x = mesh[0];
  mesh3.y = mesh[1];
  mesh3.z = mesh[2];
  double3 meshi;

  meshi.x = 1./mesh[0];
  meshi.y = 1./mesh[1];
  meshi.z = 1./mesh[2];
  
  KERNELCALL(p3m_k_space_error_gpu_kernel_ik,grid,block,(mesh3, meshi, cao, alpha_L, he_q));

  cudaError_t err;

  err = cudaGetLastError();
 
  if (err != cudaSuccess) { 
    printf("CUDA error: %s, %s line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
    return 0;
  }

  he_q_final = thrust::reduce(thrust::device_ptr<double>(he_q), thrust::device_ptr<double>(he_q + mesh_size),  0.0,  thrust::plus<double>());

  return 2.0*prefactor*sum_q2*sqrt( he_q_final / npart) / (box[1]*box[2]);
}
