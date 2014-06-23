#include "utils.hpp"
#include "grid.hpp"

double dipole_ia(float pf, float* r1, floar *r2, float* dip1, float* dip2, float* f1, float* torque1, int force_flag)
{
  double u,r,pe1,pe2,pe3,pe4,r3,r5,r2,r7,a,b,cc,d,ab;
#ifdef ROTATION
  double bx,by,bz,ax,ay,az; 
#endif
  double  ffx,ffy,ffz;
  double dr[3];
 
	
  // Distance between particles
  get_mi_vector(dr,r1,r2);

  // Powers of distance
  r2=dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
  r=sqrt(r2);
  r3=r2*r;
  r5=r3*r2;
  r7=r5*r2;
 
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
    f1[0]=pf*ab*dr[0]+cc*dip1[0]+d*dip2[0];
    f1[1]=pf*ab*dr[1]+cc*dip1[1]+d*dip2[1];
    f1[2]=pf*ab*dr[2]+cc*dip1[2]+d*dip2[2];
    
// Torques
#ifdef ROTATION
    ax=dip1[1]*dip2[2]-dip2[1]*dip1[2];
    ay=dip2[0]*dip1[2]-dip1[0]*dip2[2];
    az=dip1[0]*dip2[1]-dip2[0]*dip1[1];
    
    bx=dip1[1]*dr[2]-dr[1]*dip1[2];
    by=dr[0]*dip1[2]-dip1[0]*dr[2];
    bz=dip1[0]*dr[1]-dr[0]*dip1[1];
    
    p1->f.torque[0]+=pf*(-ax/r3+bx*cc);
    p1->f.torque[1]+=pf *(-ay/r3+by*cc);
    p1->f.torque[2]+=pf *(-az/r3+bz*cc);
    
#endif
  }    
	
  // Return energy
  return u;
}
