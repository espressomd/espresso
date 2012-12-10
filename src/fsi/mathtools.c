#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mathtools.h"

void double2double(int n, double* d, double* f) {
	int i;
	for (i=0;i<n;i++)
	  f[i] = (double) d[i];
}


int maximum(int a, int b) {
	if (a > b) return a;
	if (b > a) return b;
	return a;
}

void set0(int n, double* a) {
    int i;
    for (i=0;i<n;i++)
      a[i] = 0.0;
}

void set1(int n, double* a) {
    int i;
    for (i=0;i<n;i++)
      a[i] = 1.0;
}

double dot(int n, double* a, double* b) {
  double s;
  int i;
  s = 0.0;
  for (i=0;i<n;i++) {
	  s = s+a[i]*b[i];
  }
  return s;
}

double normr(int n, double* r) {
	//double x = sqrtf(dot(n,r,r));
	//if(x== 0.0 || x != x) printf("error normr\n");
	return sqrtf(dot(n,r,r));
}

void normv(int n, double* r, double* rdach) {
	double lengthr;
	lengthr = normr(n,r);
	bdivsx(n,lengthr,r,rdach);
}

void setvector3(double a, double b, double c, double* x) {
	x[0] = a; 
	x[1] = b;
	x[2] = c;
}

void bdivsa(int n, double s, double* a, double* b) {
	int i;
	for (i=0;i<n;i++) {
		b[i] = b[i]/(s*a[i]);
	}
}

void btimess(int n, double s, double* b) {
	int i;
	for (i=0;i<n;i++) {
		b[i] = b[i]*s;
	}
}

void copyv(int n, double* a, double* b) {
	int i;
	for (i=0;i<n;i++) b[i] = a[i];
}

void btimessx(int n, double s, double* b, double* x) {
	int i;
	for (i=0;i<n;i++) {
		//printf("(in btimessx 1) m%d= %e %e %e\n",i,x[3*i+0],x[3*i+1],x[3*i+2]);
		//printf("(in btimessx 1) b%d= %e %e %e\n",i,b[3*i+0],b[3*i+1],b[3*i+2]);
		//printf("(in btimessx 1) s= %e\n",s);
		x[i] = b[i]*s;
	}
}


void bdivs(int n, double s, double* b) {
	double si;
	si = 1.0/s;
	btimess(n,si,b);
}

void bdivsx(int n, double s, double* b, double* x) {
	double si;
	si = 1.0/s;
	btimessx(n,si,b,x);
}


void aminusb(int n, double* a, double* b, double* c) {
	int i;
	for (i=0;i<n;i++) {
		c[i] = a[i] - b[i];
	}
}

double angleab(double* a, double* b) {
	double d;
	double al,bl;
	double ac;
	d = dot(3,a,b);
	al = normr(3,a);
	bl = normr(3,b);
	// printf("(in angleab)d=%e al=%e bl=%e\n",d,al,bl);
	/* d = acos(alpha) * al * bl */
	//printf("(in angleab)acosf(d/(al*bl))=%e\n",acosf(d/(al*bl)));
	ac = d/(al*bl);
	if (ac>1.0) ac = 1.0;
	if (ac<-1.0) ac = -1.0;
	return ( acos(ac) );
}

void cisaxb(double* a, double* b, double* c) {
	c[0] = a[1]*b[2]-a[2]*b[1];
	c[1] = a[2]*b[0]-a[0]*b[2];
	c[2] = a[0]*b[1]-a[1]*b[0];
}

void cisaxaxb(double* a, double* b, double* c) {
	double axb[3];
	cisaxb(a,b,axb);
	cisaxb(a,axb,c);
}

void rotatecfromatob(double* m, double* h, double* c, double phi) {
	double aab[3];
	double x2[3];
	double x1[3];
	double cosphi,sinphi;
	int i;
	cisaxaxb(m,h,aab);
	normv(3,aab,x2);
	btimess(3,-1.0,x2); /* we need -mxmxh */
	normv(3,m,x1);
	cosphi = cos(phi);
	sinphi = sin(phi);
	//printf("(in rotatecfromatob) cosphi=%e, sinphi=%e, x1= %e %e %e, x2= %e %e %e\n",cosphi,sinphi,x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]);
	for (i=0;i<3;i++) c[i] = cosphi*x1[i] + sinphi*x2[i]; 
}

double distv(int n, double* a, double* b) {
   double d;
   int i;
   d = 0.0;
   for (i=0;i<n;i++) {
	   d = d + (a[i]-b[i])*(a[i]-b[i]);
   }
   d = sqrtf(d);
   return d;
}

double volume_sphere(double r) {
	return ((4*M_PI*r*r*r)/3.0);
}

double area_triangle_new(double *P1, double *P2) {
 // Computes the area of triangle between P1 and P2, by computing the crossproduct P1 x P2 and taking the half of its norm
 double area;
 double normal[3], n; //auxiliary variables
 cisaxb(P1,P2,normal); 
 n=normr(3,normal);
 area = 0.5*n;
 return area;
}

double area_triangle(double *P1, double *P2, double *P3) {
	// Computes the area of triangle P1,P2,P3 by computing the crossproduct P1P2 x P1P3 and taking the half of its norm
	double area;
	double u[3],v[3],normal[3], n; //auxiliary variables
	u[0] = P2[0] - P1[0]; // u = P1P2
	u[1] = P2[1] - P1[1]; 
	u[2] = P2[2] - P1[2]; 
	v[0] = P3[0] - P1[0]; // v = P1P3
	v[1] = P3[1] - P1[1]; 
	v[2] = P3[2] - P1[2]; 
	cisaxb(u,v,normal); 
	n=normr(3,normal);
	area = 0.5*n;
	return area;
}

//    This function returns the angle btw the triangle p1,p2,p3 and p2,p3,p4. 
//    Be careful, the angle depends on the orientation of the trianlges! 
//    You need to be sure that the orientation (direction of normal vector) 
//    of p1p2p3 is given by the cross product p2p1 x p2p3. 
//    The orientation of p2p3p4 must be given by p2p3 x p2p4. 

//    Example: p1 = (0,0,1), p2 = (0,0,0), p3=(1,0,0), p4=(0,1,0). 
//    The orientation of p1p2p3 should be in the direction (0,1,0) 
//    and indeed: p2p1 x p2p3 = (0,0,1)x(1,0,0) = (0,1,0)

//    This function is called in the beginning of the simulation when creating 
//    bonds depending on the angle btw the triangles, the bending_force.
//    Here, we determine the orientations by looping over the triangles 
//    and checking the correct orientation. So when defining the bonds by tcl command
//    "part p2 bond xxxx p1 p3 p4", we correctly input the particle id's.
//    So if you have the access to the order of particles, you are safe to call this
//    function with exactly this order. Otherwise you need to check the orientations.


double angle_btw_triangles_new(double *P1, double *P2, double *P3, double *P4) {
	double phi;
	double normal1[3],normal2[3]; //auxiliary variables

	cisaxb(P1,P2,normal1); 
 
	cisaxb(P3,P4,normal2); 

	double tmp11,tmp22,tmp33;
	// Now we compute the scalar product of n1 and n2 divided by the norms of n1 and n2
	tmp11 = dot(3,normal1,normal2);         // tmp11 = n1.n2
	tmp22 = normr(3,normal1);
	tmp33 = normr(3,normal2);
	tmp11 /= (tmp22*tmp33);  // tmp11 = n1.n2/(|n1||n2|)
	//if(tmp11<=-1. || tmp11>=1.)printf("x %e acos(x) %e dn1 %e dn2 %e ",tmp11,acos(tmp11),tmp22,tmp33);
	
	if(tmp11>=1.)tmp11=0.0;
	else if(tmp11<=-1.)tmp11=M_PI;
	phi = M_PI - acos(tmp11); 	// The angle between the faces (not considering the orientation, always less or equal to Pi) is
								// equal to Pi minus angle between the normals
	
	// Now we need to determine, if the angle btw two triangles is less than Pi or more than Pi. To do this we check, 
	// if the point P4 lies in the halfspace given by trianlge P1P2P3 and the normal to this triangle. If yes, we have 
	// angle less than Pi, if not, we have angle more than Pi.
	// General equation of the plane is n_x*x + n_y*y + n_z*z + d = 0 where (n_x,n_y,n_z) is the normal to the plane.
	// Point P1 lies in the plane, therefore d = -(n_x*P1_x + n_y*P1_y + n_z*P1_z)
	// Point P4 lies in the halfspace given by normal iff n_x*P4_x + n_y*P4_y + n_z*P4_z + d >= 0
	tmp11 = - (normal1[0]*P1[0] + normal1[1]*P1[1] + normal1[2]*P1[2]);
	if (normal1[0]*P4[0] + normal1[1]*P4[1] + normal1[2]*P4[2] + tmp11 < 0) phi = 2*M_PI - phi;
	//printf("%e\n",phi); 
	//printf("       %e %e %e\n",P2[0],P2[1],P2[2]);		
	//printf("       %e %e %e\n",P1[0],P1[1],P1[2]);		
	//printf("       %e %e %e\n",P3[0],P3[1],P3[2]);		
	//printf("       %e %e %e\n",P4[0],P4[1],P4[2]);		
	return phi;
}

double angle_btw_triangles(double *P1, double *P2, double *P3, double *P4) {
	double phi;
	double u[3],v[3],normal1[3],normal2[3]; //auxiliary variables
	u[0] = P1[0] - P2[0]; // u = P2P1
	u[1] = P1[1] - P2[1]; 
	u[2] = P1[2] - P2[2]; 
	v[0] = P3[0] - P2[0]; // v = P2P3
	v[1] = P3[1] - P2[1]; 
	v[2] = P3[2] - P2[2]; 
	cisaxb(u,v,normal1); 
	u[0] = P3[0] - P2[0]; // u = P2P3
	u[1] = P3[1] - P2[1]; 
	u[2] = P3[2] - P2[2]; 
	v[0] = P4[0] - P2[0]; // v = P2P4
	v[1] = P4[1] - P2[1]; 
	v[2] = P4[2] - P2[2]; 
	cisaxb(u,v,normal2); 

	double tmp11,tmp22,tmp33;
	// Now we compute the scalar product of n1 and n2 divided by the norms of n1 and n2
	tmp11 = dot(3,normal1,normal2);         // tmp11 = n1.n2
	tmp22 = normr(3,normal1);
	tmp33 = normr(3,normal2);
	tmp11 /= (tmp22*tmp33);  // tmp11 = n1.n2/(|n1||n2|)
	//if(tmp11<=-1. || tmp11>=1.)printf("x %e acos(x) %e dn1 %e dn2 %e ",tmp11,acos(tmp11),tmp22,tmp33);
	
	if(tmp11>=1.)tmp11=0.0;
	else if(tmp11<=-1.)tmp11=M_PI;
	phi = M_PI - acos(tmp11); 	// The angle between the faces (not considering the orientation, always less or equal to Pi) is
								// equal to Pi minus angle between the normals
	
	// Now we need to determine, if the angle btw two triangles is less than Pi or more than Pi. To do this we check, 
	// if the point P4 lies in the halfspace given by trianlge P1P2P3 and the normal to this triangle. If yes, we have 
	// angle less than Pi, if not, we have angle more than Pi.
	// General equation of the plane is n_x*x + n_y*y + n_z*z + d = 0 where (n_x,n_y,n_z) is the normal to the plane.
	// Point P1 lies in the plane, therefore d = -(n_x*P1_x + n_y*P1_y + n_z*P1_z)
	// Point P4 lies in the halfspace given by normal iff n_x*P4_x + n_y*P4_y + n_z*P4_z + d >= 0
	tmp11 = - (normal1[0]*P1[0] + normal1[1]*P1[1] + normal1[2]*P1[2]);
	if (normal1[0]*P4[0] + normal1[1]*P4[1] + normal1[2]*P4[2] + tmp11 < 0) phi = 2*M_PI - phi;
	//printf("%e\n",phi); 
	//printf("       %e %e %e\n",P2[0],P2[1],P2[2]);		
	//printf("       %e %e %e\n",P1[0],P1[1],P1[2]);		
	//printf("       %e %e %e\n",P3[0],P3[1],P3[2]);		
	//printf("       %e %e %e\n",P4[0],P4[1],P4[2]);		
	return phi;
}

void dipolefield(double* m, double* r, double *field) {
	/* this is for m as Ms times volume
	 * B(m,r) = (m0/(4pir^3))(3(m*r)r-m)+(2mu0/3)mdelta(r) 
	 * now for m as Js times volume 
	 * B(m,r) = (1/4pir^3))(3(m*r)r-m)+ (2/3)mdelta(r) 
	 *        = (3/4pir^3))((m*r)r-(1/3)m)+ (2/3)mdelta(r) */
	double mr;
	double v;
	double rdach[3];
	int i;
	normv(3,r,rdach);
	mr = dot(3,m,rdach);
	v = volume_sphere(normr(3,r));
	for (i=0;i<3;i++) {
		field[i] = (1.0/v)*(mr*rdach[i]-m[i]/3.0);
	}	
}

int meshsizesphere(double r, double dx) {
	int l;
	int i,j,k;
	int n;
	double y[3];
	l = (int)(r/dx);
	n = 0;
	for (i=-l;i<=l;i++) {
	for (j=-l;j<=l;j++) {
	for (k=-l;k<=l;k++) {
	   y[0] = i*dx; y[1] = j*dx; y[2] = k*dx;
       if (normr(3,y)<r) {  /* in sphere */
          n = n+1;
	   }
    }
    } 
    }
    return n;
}

int meshsphere(int n0, double r, double* center, double dx, double* x, int lx) {
	int i,j,k;
	int l;
	int n;
	double y[3];
	double volumecubes, volumesphere;
    double a,s;
    double mydx;
    if (n0<=0) {
      n = meshsizesphere(r,dx);
      if (3*n > lx) {
		printf("worskpace for meshing too small\n");
		return 0;
	  }
	}
	l = (int)(r/dx);
	n = 0;
	for (i=-l;i<=l;i++) {
	for (j=-l;j<=l;j++) {
	for (k=-l;k<=l;k++) {
	   y[0] = i*dx; y[1] = j*dx; y[2] = k*dx;
       if (normr(3,y)<r) {  /* in sphere */
          n = n+1;
	   }
    }
    } 
    }
    volumesphere = volume_sphere(r);
    volumecubes = n*dx*dx*dx;
    a = volumesphere/volumecubes;
    s = pow(a,0.33333333333333333333333333333);
    mydx = s*dx;
    n = 0;
	for (i=-l;i<=l;i++) {
	for (j=-l;j<=l;j++) {
	for (k=-l;k<=l;k++) {
	   y[0] = i*dx; y[1] = j*dx; y[2] = k*dx;
       if (normr(3,y)<r) {  /* in sphere */
          if (3*n+2 >= lx) {
		    printf("worskpace for meshing too small\n");
		    return 0;
	      } 
          x[3*n+0] = i*mydx + center[0]; 
          x[3*n+1] = j*mydx + center[1]; 
          x[3*n+2] = k*mydx + center[2];
          n = n+1;
	   }
    }
    } 
    }
    volumecubes = n*mydx*mydx*mydx;
    a = volumesphere/volumecubes;
    return n;
}

double interpolate(int n, double* x, double *y, double x0) {
	int i;
	double l;
	double lx,w;
	if (x0<=x[0]) return y[0];
	if (x0>=x[n-1]) return y[n-1];
	for (i=0;i<n-1;i++) {
		if ((x[i] <= x0) && (x0 <= x[i+1])) {
		   l = x[i+1]-x[i];
		   lx = x0 - x[i];
		   w = lx/l;
		   return (y[i]*(1.0-w)+y[i+1]*w);
		}
	}
    return y[n-1];
}

void get_n_triangle(double* p1, double* p2, double* p3, double* n){
	n[0]=(p2[1]-p1[1])*(p3[2]-p1[2])-(p2[2]-p1[2])*(p3[1]-p1[1]);
	n[1]=(p2[2]-p1[2])*(p3[0]-p1[0])-(p2[0]-p1[0])*(p3[2]-p1[2]);
	n[2]=(p2[0]-p1[0])*(p3[1]-p1[1])-(p2[1]-p1[1])*(p3[0]-p1[0]);
	
	//n1=(b2-a2)*(c3-a3)-(b3-a3)*(c2-a2);
	//n2=(b3-a3)*(c1-a1)-(b1-a1)*(c3-a3);
	//n3=(b1-a1)*(c2-a2)-(b2-a2)*(c1-a1);
}
