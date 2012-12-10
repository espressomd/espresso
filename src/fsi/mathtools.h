#ifndef included_mathtools_h
#define included_mathtools_h


/* one vectors */
double normr(int n, double* r);

/* two vectors */
double dot(int n, double* a, double* b);
void aminusb(int n, double* a, double* b, double* c);
void cisaxb(double* a, double* b, double* c);

/* geometry */
double area_triangle_new(double *P1, double *P2);
double area_triangle(double *P1, double *P2, double *P3);
double angle_btw_triangles(double *P1, double *P2, double *P3, double *P4);
void get_n_triangle(double* p1, double* p2, double* p3, double* n);  // calculate normal vector of triangle

#endif 
