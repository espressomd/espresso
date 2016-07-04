#include "utils.hpp"
#include "config.hpp"


typedef double real;

__device__ inline void get_mi_vector_real(real res[3], real a[3], real b[3],float box_l[3],int periodic[3])
{
  int i;

  for(i=0;i<3;i++) {
    res[i] = a[i] - b[i];
#ifdef PARTIAL_PERIODIC
    if (periodic[i])
#endif
      res[i] -= floor(res[i]/box_l[i]+0.5)*box_l[i];
  }
}


//__device__ real scalar(real a[3], real b[3])
//{
// real sum=0.;
// for (int i=0;i<3;i++)
//  sum+=a[i]*b[i];
// return sum;
//}

#define scalar(a,b) a[0]*b[0]+a[1]*b[1]+a[2]*b[2]



__device__ float dipole_ia(float pf, float* _r1, float *_r2, float* _dip1, float* _dip2, float* f1, float* torque1, int force_flag, float box_l[3], int periodic[3])
{
real dip1[3],dip2[3],r1[3],r2[3];
for (int i=0;i<3;i++)
{
 dip1[i]=_dip1[i];
 dip2[i]=_dip2[i];
 r1[i]=_r1[i];
 r2[i]=_r2[i];

}
  real u,r,pe1,pe2,pe3,pe4,r3,r5,r_sq,r7,a,b,cc,d,ab;
#ifdef ROTATION
  real bx,by,bz,ax,ay,az; 
#endif
  real dr[3];
 
	
  // Distance between particles
  get_mi_vector_real(dr,r1,r2,box_l,periodic);

  // Powers of distance
  r_sq=scalar(dr,dr);
  r=sqrt(r_sq);
  r3=r_sq*r;
  r5=r3*r_sq;
  r7=r5*r_sq;
 
  // Dot products
  pe1=scalar(dip1,dip2);
  pe2=scalar(dip1,dr);
  pe3=scalar(dip2,dr);
  pe4=3.0/r5;

  // Energy, if requested
  u= pf* ( pe1/r3 - pe4*pe2*pe3);

  // Force, if requested
  if(force_flag) { 
    a=pe4*pe1;
    b=-15.0*pe2*pe3/r7;
    ab =a+b;
    cc=pe4*pe3;
    d=pe4*pe2;
    
    //  Result
    f1[0]=(float)(pf*(ab*dr[0]+cc*dip1[0]+d*dip2[0]));
    f1[1]=(float)(pf*(ab*dr[1]+cc*dip1[1]+d*dip2[1]));
    f1[2]=(float)(pf*(ab*dr[2]+cc*dip1[2]+d*dip2[2]));
    
// Torques
#ifdef ROTATION
    ax=dip1[1]*dip2[2]-dip2[1]*dip1[2];
    ay=dip2[0]*dip1[2]-dip1[0]*dip2[2];
    az=dip1[0]*dip2[1]-dip2[0]*dip1[1];
    
    bx=dip1[1]*dr[2]-dr[1]*dip1[2];
    by=dr[0]*dip1[2]-dip1[0]*dr[2];
    bz=dip1[0]*dr[1]-dr[0]*dip1[1];
    
    torque1[0]=(float)(pf*(-ax/r3+bx*cc));
    torque1[1]=(float)(pf *(-ay/r3+by*cc));
    torque1[2]=(float)(pf *(-az/r3+bz*cc));
    
#endif
  }    
	
  // Return energy
  return u;
}
