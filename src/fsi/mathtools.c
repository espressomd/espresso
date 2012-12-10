#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mathtools.h"


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


void aminusb(int n, double* a, double* b, double* c) {
	int i;
	for (i=0;i<n;i++) {
		c[i] = a[i] - b[i];
	}
}


void cisaxb(double* a, double* b, double* c) {
	c[0] = a[1]*b[2]-a[2]*b[1];
	c[1] = a[2]*b[0]-a[0]*b[2];
	c[2] = a[0]*b[1]-a[1]*b[0];
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







void get_n_triangle(double* p1, double* p2, double* p3, double* n){
	n[0]=(p2[1]-p1[1])*(p3[2]-p1[2])-(p2[2]-p1[2])*(p3[1]-p1[1]);
	n[1]=(p2[2]-p1[2])*(p3[0]-p1[0])-(p2[0]-p1[0])*(p3[2]-p1[2]);
	n[2]=(p2[0]-p1[0])*(p3[1]-p1[1])-(p2[1]-p1[1])*(p3[0]-p1[0]);
	
	//n1=(b2-a2)*(c3-a3)-(b3-a3)*(c2-a2);
	//n2=(b3-a3)*(c1-a1)-(b1-a1)*(c3-a3);
	//n3=(b1-a1)*(c2-a2)-(b2-a2)*(c1-a1);
}
