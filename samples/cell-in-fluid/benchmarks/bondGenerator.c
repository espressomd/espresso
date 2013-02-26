#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Just ordinary cross product
int CrossProduct(double *u, double *v, double *res)
{
	res[0] = u[1]*v[2] - u[2]*v[1];
	res[1] = u[2]*v[0] - u[0]*v[2];
	res[2] = u[0]*v[1] - u[1]*v[0];
}

// calculates normal vector of triangle
double get_n_triangle(double *P1, double *P2, double *P3, double *normal){
	double u[3],v[3];
	u[0] = P1[0] - P2[0]; // u = P2P1
	u[1] = P1[1] - P2[1]; 
	u[2] = P1[2] - P2[2]; 
	v[0] = P3[0] - P2[0]; // v = P2P3
	v[1] = P3[1] - P2[1]; 
	v[2] = P3[2] - P2[2]; 
	CrossProduct(u,v,normal);
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
	double u[3],v[3];
	double normal1[3],normal2[3]; //auxiliary variables
	//u[0] = P1[0] - P2[0]; // u = P2P1
	//u[1] = P1[1] - P2[1]; 
	//u[2] = P1[2] - P2[2]; 
	//v[0] = P3[0] - P2[0]; // v = P2P3
	//v[1] = P3[1] - P2[1]; 
	//v[2] = P3[2] - P2[2]; 
	//CrossProduct(u,v,normal1); 
	get_n_triangle(P1,P2,P3,normal1);
	//u[0] = P3[0] - P2[0]; // u = P2P3
	//u[1] = P3[1] - P2[1]; 
	//u[2] = P3[2] - P2[2]; 
	//v[0] = P4[0] - P2[0]; // v = P2P4
	//v[1] = P4[1] - P2[1]; 
	//v[2] = P4[2] - P2[2]; 
	//CrossProduct(u,v,normal2); 
	get_n_triangle(P3,P2,P4,normal2);

	double tmp11;
	// Now we compute the scalar product of n1 and n2 divided by the norms of n1 and n2
	tmp11 = normal1[0]*normal2[0] + normal1[1]*normal2[1] + normal1[2]*normal2[2];         // tmp11 = n1.n2
	tmp11 /= sqrt(normal1[0]*normal1[0] + normal1[1]*normal1[1] + normal1[2]*normal1[2]);  // tmp11 = n1.n2/|n1|
	tmp11 /= sqrt(normal2[0]*normal2[0] + normal2[1]*normal2[1] + normal2[2]*normal2[2]);  // tmp11 = n1.n2/(|n1||n2|)
	
	if(tmp11>=1.)tmp11=0.0;
	else if(tmp11<=-1.)tmp11=M_PI;
	
	phi = M_PI - acos(tmp11); 	// The angle between the faces (not considering the orientation, always less or equal to Pi) is
								// equal to Pi minus angle between the normals

	// Now we need to determine, if the angle btw two triangles is less than Pi or more than Pi. To do this we check, if the point P4 lies in the halfspace given by trianlge P1P2P3 and the normal to this triangle. If yes, we have angle less than Pi, if not, we have angle more than Pi.
	// General equation of the plane is n_x*x + n_y*y + n_z*z + d = 0 where (n_x,n_y,n_z) is the normal to the plane.
	// Point P1 lies in the plane, therefore d = -(n_x*P1_x + n_y*P1_y + n_z*P1_z)
	// Point P4 lies in the halfspace given by normal iff n_x*P4_x + n_y*P4_y + n_z*P4_z + d >= 0
	tmp11 = - (normal1[0]*P1[0] + normal1[1]*P1[1] + normal1[2]*P1[2]);
	if (normal1[0]*P4[0] + normal1[1]*P4[1] + normal1[2]*P4[2] + tmp11 < 0) phi = 2*M_PI - phi;
	//printf("%d %d %d %d %e \n",p2->p.identity, p1->p.identity,p3->p.identity,p4->p.identity,phi);
	//printf("%e\n",phi); 
	//printf("       %e %e %e\n",P2[0],P2[1],P2[2]);		
	//printf("       %e %e %e\n",P1[0],P1[1],P1[2]);		
	//printf("       %e %e %e\n",P3[0],P3[1],P3[2]);		
	//printf("       %e %e %e\n",P4[0],P4[1],P4[2]);		

	return(phi);
}


double area_triangle(double *P1, double *P2, double *P3) {
	// Computes the area of triangle P1,P2,P3 by computing the crossproduct P1P2 x P1P3 and taking the half of its norm
	double area;
	//double u[3],v[3];
	double normal[3]; //auxiliary variables
	//u[0] = P2[0] - P1[0]; // u = P1P2
	//u[1] = P2[1] - P1[1]; 
	//u[2] = P2[2] - P1[2]; 
	//v[0] = P3[0] - P1[0]; // v = P1P3
	//v[1] = P3[1] - P1[1]; 
	//v[2] = P3[2] - P1[2]; 
	//CrossProduct(u,v,normal); 
	get_n_triangle(P1,P2,P3,normal);
	area = 0.5*sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
	return(area);
}

double discard_epsilon(double x)
{
    return ((-1e-10 < x) && (x < 1e-10)) ? 0.0 : x;
}

int main(int argc, char *argv[]){
	int i,j,k;
    
    printf("==== 1 ====\n");
	if(argc!=37){
		printf("arg needs 36 arguments (%d given): filenameNodes, filenameTriangles, nnodes, ntriangles,\
		rotate x y z, origin x y z, stretch x y z, \
		createPart, typePart, molPart, firstPartId, \
		bondS, bondB, bondAl, bondAg, bondV, bondVA, partS, partB, partAl, partAg, partV, partVA, ks, kb, kal, kag, kv, firstBondId, partmass\n", argc-1);
		for (i=0;i<argc;i++)printf("i:%d %s\n", i, argv[i]);
		return -1;
	}
	// input mesh
	char *filenamenodes = argv[1];
	char *filenametriangles = argv[2];
	int mesh_nnodes = (int)strtod(argv[3],NULL);
	int mesh_ntriangles = (int)strtod(argv[4],NULL);
	
	// pos, ori, size
	double rotate_X = strtod(argv[5], NULL);
	double rotate_Y = strtod(argv[6], NULL);
	double rotate_Z = strtod(argv[7], NULL);
	double origin_X = strtod (argv[8],NULL);
	double origin_Y = strtod (argv[9],NULL);
	double origin_Z = strtod (argv[10],NULL);
	double stretch_X = strtod (argv[11],NULL);
	double stretch_Y = strtod (argv[12],NULL);
	double stretch_Z = strtod (argv[13],NULL);
	
	// particle creation
	const char *createPart = argv[14];
	int typePart = (int)strtod(argv[15],NULL);
	int molPart = (int)strtod(argv[16],NULL);
	int firstPartId = (int)strtod(argv[17],NULL);
	
	// stretching, bending and area forces
	const char *bondsStretching = argv[18];
	const char *bondsBending = argv[19];
	const char *bondsAreaLocal = argv[20];
	const char *bondsAreaGlobal = argv[21];
	const char *bondsVolume = argv[22];
	const char *bondsVolumeAreaGlobal = argv[23];
	const char *partStretching = argv[24];
	const char *partBending = argv[25];
	const char *partAreaLocal = argv[26];
	const char *partAreaGlobal = argv[27];
	const char *partVolume = argv[28];
	const char *partVolumeAreaGlobal = argv[29];
	// Coefficients entering the definition of the corresponding bonds and forces, See Dupin 2007, Phys Rev E Stat Nonlin Soft Matter Phys. 2007 June ; 75(6 Pt 2): 066707
	double coeff_stretching = strtod (argv[30],NULL);
	double coeff_bending = strtod (argv[31],NULL);
	double coeff_area_local = strtod (argv[32],NULL);
	double coeff_area_global = strtod (argv[33],NULL);
	double coeff_volume = strtod (argv[34],NULL);
    int firstBondId = (int)strtod (argv[35],NULL); // Should be zero if no bonds were defined before, can be also nonzero if we have already defined other bond interactions
	double part_mass = strtod (argv[36],NULL);

	//const char *triIds = argv[31];
	
	//////////////////////////////////////////////
	//											//
	// End of user definitions				  	//
	//											//
	//////////////////////////////////////////////

	int firstID_StrBond = 0;		// will be defined later
	int firstID_BenBond = 0;		// will be defined later	
	int firstID_localAreaBond = 0;		// will be defined later
	int firstID_globalAreaBond = 0;		// will be defined later
	int firstID_VolumeBond = 0;		// will be defined later
	int n_StrBond = 0;		// will be defined later		
	int n_BenBond = 0;		// will be defined later
	int n_localAreaBond = 0;		// will be defined later
	int n_globalAreaBond = 0;		// will be defined later
	int n_VolumeBond = 0;		// will be defined later
	
	
	int MESH_MAX_NODES = 10000;         // Maximum number of nodes
	int MESH_MAX_EDGES = 10000;         // Maximum number of edges
	int MESH_MAX_TRIANGLES = 10000;     // Maximum number of triangles

	int mesh_nedges = 0;               // GID did not produce the list of edges. So we need to create one.

	FILE *f,*fpart;
	FILE *f2,*f2part;

	double mesh_nodes[MESH_MAX_NODES][3];              // Array of doubles storing the coordinates of the nodes
	int mesh_edges[MESH_MAX_EDGES][2];                 // Array of integers storing the couples of indices for nodes
	int mesh_triangles[MESH_MAX_TRIANGLES][3];         // Array of integers storing the triplets of indices for triangles

	                //calculate transformation matrix for rotation of the cell
	double ca = cos(rotate_X);
	double sa = sin(rotate_X);
	double cb = cos(rotate_Y);
	double sb = sin(rotate_Y);
	double cc = cos(rotate_Z);
	double sc = sin(rotate_Z);
	
	double rotation[3][3] = {
		{ cb * cc  ,   sa * sb * cc - ca * sc  ,   sc * sa + cc * sb * ca},
		{ cb * sc  ,   ca * cc - sa * sb * sc  ,   sc * sb * ca - cc * sa},
		{-sb       ,   cb * sa                 ,   ca * cb               }
	};
	

	printf("==== 2 ====\n");
	///////// Reading the nodes ////////////////
	f = fopen(filenamenodes,"rt");
	double xmin = 1000000.,xmax = -1000000.,ymin = 1000000.,ymax = -1000000.,zmin = 1000000.,zmax = -1000000.;
	for (i = 0; i < mesh_nnodes; i++) {
		fscanf(f,"%lf %lf %lf\n", &mesh_nodes[i][0], &mesh_nodes[i][1], &mesh_nodes[i][2]);
		// rotation of the object


		double xx = discard_epsilon(rotation[0][0] * mesh_nodes[i][0] + rotation[0][1] * mesh_nodes[i][1] + rotation[0][2] * mesh_nodes[i][2]);
		double yy = discard_epsilon(rotation[1][0] * mesh_nodes[i][0] + rotation[1][1] * mesh_nodes[i][1] + rotation[1][2] * mesh_nodes[i][2]);
		double zz = discard_epsilon(rotation[2][0] * mesh_nodes[i][0] + rotation[2][1] * mesh_nodes[i][1] + rotation[2][2] * mesh_nodes[i][2]);
		mesh_nodes[i][0] = xx;
		mesh_nodes[i][1] = yy;
		mesh_nodes[i][2] = zz;


		// stretch and translate to origin
		mesh_nodes[i][0] = mesh_nodes[i][0]*stretch_X + origin_X;
		mesh_nodes[i][1] = mesh_nodes[i][1]*stretch_Y + origin_Y;
		mesh_nodes[i][2] = mesh_nodes[i][2]*stretch_Z + origin_Z;
		
		
		if (mesh_nodes[i][0] < xmin) xmin = mesh_nodes[i][0];
		if (mesh_nodes[i][0] > xmax) xmax = mesh_nodes[i][0];
		if (mesh_nodes[i][1] < ymin) ymin = mesh_nodes[i][1];
		if (mesh_nodes[i][1] > ymax) ymax = mesh_nodes[i][1];
		if (mesh_nodes[i][2] < zmin) zmin = mesh_nodes[i][2];
		if (mesh_nodes[i][2] > zmax) zmax = mesh_nodes[i][2];
	}
	printf("Nodes have been read. Total %d nodes\n\n",mesh_nnodes);
	printf("xmin = %e, xmax = %e  nodes %d\n",xmin,xmax, mesh_nnodes);
	printf("ymin = %e, ymax = %e\n",ymin,ymax);
	printf("zmin = %e, zmax = %e\n",zmin,zmax);
	fclose(f);
	
	
	///////// Reading the triangles ////////////////
	f = fopen(filenametriangles,"rt");
	for (i = 0; i < mesh_ntriangles; i++) {
	fscanf(f,"%d %d %d\n", &mesh_triangles[i][0], &mesh_triangles[i][1], &mesh_triangles[i][2]);
	//mesh_triangles[i][0] -= 1;   //GID indexes the nodes from 1 to mesh_nnodes. We use the indexes from 0 to mesh_nnodes-1. Therefore the subtracting of 1
	//mesh_triangles[i][1] -= 1;
	//mesh_triangles[i][2] -= 1;

	// GID files corrected.

	//if (i >= mesh_ntriangles/2) {
	  //int tmp = mesh_triangles[i][1];
	  //mesh_triangles[i][1] = mesh_triangles[i][2];
	  //mesh_triangles[i][2] = tmp;   // Switching of the indexes because of the wrong orientation of the triangles obtained from GID - see notes.tex
	//}
	}
	printf("Triangles have been read. Total %d triangles\n\n",mesh_ntriangles);
	fclose(f);
	
	//////// Create corrected triangles files ///////////
  //f = fopen(triIds,"wt");

  //int n1,n2,n3;
  //for (int i = 0; i < mesh_ntriangles; i++) {
	  ////for(int j=0; j<mesh_ntriangles;j++){
		  ////if(j!=i){
			  ////for(int k=0;k<3;k++){
				////if(mesh_triangles[i][0]==mesh_triangles[j][(k+1)%3] && mesh_triangles[i][1]==mesh_triangles[j][k]) n1=j;  // AB i == BA,CB,AC j
				////if(mesh_triangles[i][1]==mesh_triangles[j][(k+1)%3] && mesh_triangles[i][2]==mesh_triangles[j][k]) n2=j;  // BC i == BA,CB,AC j
				////if(mesh_triangles[i][2]==mesh_triangles[j][(k+1)%3] && mesh_triangles[i][0]==mesh_triangles[j][k]) n3=j;  // CA i == BA,CB,AC j
			  ////}
		  ////}
	  ////}
    //fprintf(f,"%d %d %d\n", mesh_triangles[i][0],mesh_triangles[i][1],mesh_triangles[i][2]);  // A B C
  //}  
  //fclose(f); 

///////// Creating the list of edges
	for (i = 0; i < mesh_ntriangles; i++){
		int pa,pb,pc, is;
		pa = mesh_triangles[i][0];
		pb = mesh_triangles[i][1];
		pc = mesh_triangles[i][2];   // Take a triangle and copy the nodes of the triangle to pa,pb,pc (point A, point B, point C)
		is = 0;
		for (j = 0; j < mesh_nedges; j++){
		  if (mesh_edges[j][0] == pa && mesh_edges[j][1] == pb) is = 1;
		  if (mesh_edges[j][1] == pa && mesh_edges[j][0] == pb) is = 1;   // Chceck if the edge AB or BA is in the current list of edges
		}
		if (is == 0) {  // If AB nor BA is in the list then add the edge AB to the list
		  mesh_edges[mesh_nedges][0] = pa;
		  mesh_edges[mesh_nedges][1] = pb;
		  mesh_nedges++;
		}
		is = 0;
		for (j = 0; j < mesh_nedges; j++){
		  if (mesh_edges[j][0] == pb && mesh_edges[j][1] == pc) is = 1;
		  if (mesh_edges[j][1] == pb && mesh_edges[j][0] == pc) is = 1;  // Chceck if the edge BC or CB is in the current list of edges
		}
		if (is == 0) {  // If BC nor CB is in the list then add the edge BC to the list
		  mesh_edges[mesh_nedges][0] = pb;
		  mesh_edges[mesh_nedges][1] = pc;
		  mesh_nedges++;
		}
		is = 0;
		for (j = 0; j < mesh_nedges; j++){
		  if (mesh_edges[j][0] == pa && mesh_edges[j][1] == pc) is = 1;
		  if (mesh_edges[j][1] == pa && mesh_edges[j][0] == pc) is = 1;   // Chceck if the edge AC or CA is in the current list of edges
		}
		if (is == 0) {     // If AC nor CA is in the list then add the edge AC to the list
		  mesh_edges[mesh_nedges][0] = pa;
		  mesh_edges[mesh_nedges][1] = pc;
		  mesh_nedges++;
		}
	}    
	printf("Edges have been created. Total %d edges\n\n",mesh_nedges);

/////////// Generate createPart
  f = fopen(createPart,"wt");
  for (i = firstPartId; i < mesh_nnodes+firstPartId; i++){
    fprintf(f,"part %d pos %e %e %e type %d mol %d mass %e\n",i,mesh_nodes[i-firstPartId][0], mesh_nodes[i-firstPartId][1], mesh_nodes[i-firstPartId][2], typePart, molPart, part_mass);
  }  
  fclose(f); 
  
  
// Generation of stretching force bonds:
    printf("Generating stretching force bonds\n");
	firstID_StrBond = firstBondId;
	n_StrBond = mesh_nedges; // Stretching is coupled to the edges
    firstBondId += n_StrBond;

	double p1[3],p2[3],p3[3],p4[3];
	double dist = 0.;
	f = fopen(bondsStretching,"wt");
	fpart = fopen(partStretching,"wt");
	for (i = 0; i < n_StrBond; i++) {
		for(k = 0; k < 3; k++){		
			p1[k]=mesh_nodes[mesh_edges[i][0]][k];
			p2[k]=mesh_nodes[mesh_edges[i][1]][k];
		}
		dist = sqrt((p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]) + (p1[2]-p2[2])*(p1[2]-p2[2])); // We need just to compute the distance btw the vertexes
		fprintf(f,"inter %d stretching_force %e %e\n",firstID_StrBond + i, dist, coeff_stretching);
		fprintf(fpart,"part %d bond %d %d\n",mesh_edges[i][0]+firstPartId, firstID_StrBond + i, mesh_edges[i][1]+firstPartId);
	}
	fclose(f); 
	fclose(fpart); 
	printf("Done.\n\n");

// Generation of bending force bonds:
    printf("Generating bending force bonds\n");
	firstID_BenBond = firstBondId;
	n_BenBond = mesh_nedges; // Bending is coupled to the angles between triangles sharing the same edge -> 
    firstBondId += n_BenBond;

	//FILE *fpartpos;
	//fpartpos = fopen("/home/cimo/Dropbox/FH/espresso/testing/partpos.dat","wt");
	double phi = 0.;
	int p1id,p2id,p3id,p4id; // IDs of the vertexes
	double P1[3],P2[3],P3[3],P4[3]; // coordinates of the vertexes
	f = fopen(bondsBending,"wt");
	fpart = fopen(partBending,"wt");
	for (i = 0; i < n_BenBond; i++) { // Run over all edges
		p2id = mesh_edges[i][0]; // Put IDs of points to p2id,p3id
		p3id = mesh_edges[i][1]; 
		for(k = 0; k < 3; k++){	// Put coordinates of the edges's points
			p2[k]=mesh_nodes[p2id][k];
			p3[k]=mesh_nodes[p3id][k];
		}
		int detected = 0; // Number of detected triangles with current edge common
		// Algorithm is as follows: we run over all triangles and check is two vertices are those from current edge. If we find such triangle, we put the ID of the third vertex to p1id and moreover we check if the orientation p1id, p2id p3id is the same as was in the trianlge list (meaning, that we found one of the following three triples in the triangle list: p1id, p2id, p3id or p2id, p3id, p1id or p3id, p1id, p2id). If we have the same orientation, we set orient = 1, otherwise orient = -1.
		// Then we go further looking for the second triagle. The second triangle should have the opposite orientation.
		// The normal of the first triangle will be P1P2 x P1P3, of the second triangle will be P2P4 x P2P3
		
		int orient = 0; // We have two different orientations. I set orientation to 1 if the first triangle is oriented p1id

		for (k = 0; k < mesh_ntriangles; k++) { // Run over all triangles and determine the two triangles with the common current edge
			if ((mesh_triangles[k][0] == p2id) && (mesh_triangles[k][1] == p3id)) { // 
				 if (detected == 0) { // if no triangle was detected
					 p1id = mesh_triangles[k][2]; 
					 detected = 1; 
					 orient = 1;}
				else { // if already one triangle was detected - then also quit the k-loop
					p4id = mesh_triangles[k][2]; 
					k = mesh_ntriangles; } }
			if ((mesh_triangles[k][1] == p2id) && (mesh_triangles[k][2] == p3id)) {
				 if (detected == 0) {
					 p1id = mesh_triangles[k][0];
					 detected = 1;
					 orient = 1; }
				else {
					p4id = mesh_triangles[k][0];
					k = mesh_ntriangles; } }
			if ((mesh_triangles[k][2] == p2id) && (mesh_triangles[k][0] == p3id)) {
				 if (detected == 0) {
					 p1id = mesh_triangles[k][1];
					 detected = 1; 
					 orient = 1;}
				else {
					p4id = mesh_triangles[k][1];
					k = mesh_ntriangles; } }
			if ((mesh_triangles[k][1] == p2id) && (mesh_triangles[k][0] == p3id)) {
				 if (detected == 0) {
					 p1id = mesh_triangles[k][2];
					 detected = 1; 
					 orient = -1;}
				else {
					p4id = mesh_triangles[k][2];
					k = mesh_ntriangles; } }
			if ((mesh_triangles[k][2] == p2id) && (mesh_triangles[k][1] == p3id)) {
				 if (detected == 0) {
					 p1id = mesh_triangles[k][0];
					 detected = 1; 
					 orient = -1;}
				else {
					p4id = mesh_triangles[k][0];
					k = mesh_ntriangles; } }
			if ((mesh_triangles[k][0] == p2id) && (mesh_triangles[k][2] == p3id)) {
				 if (detected == 0) {
					 p1id = mesh_triangles[k][1];
					 detected = 1; 
					 orient = -1;}
				else {
					p4id = mesh_triangles[k][1];
					k = mesh_ntriangles; } }
		}
		if (orient == 1) { int tmp22 = p1id; p1id = p4id; p4id = tmp22;} //This is to have the correct orientation
		P1[0] = mesh_nodes[p1id][0]; P1[1] = mesh_nodes[p1id][1]; P1[2] = mesh_nodes[p1id][2];
		P2[0] = mesh_nodes[p2id][0]; P2[1] = mesh_nodes[p2id][1]; P2[2] = mesh_nodes[p2id][2];
		P3[0] = mesh_nodes[p3id][0]; P3[1] = mesh_nodes[p3id][1]; P3[2] = mesh_nodes[p3id][2];
		P4[0] = mesh_nodes[p4id][0]; P4[1] = mesh_nodes[p4id][1]; P4[2] = mesh_nodes[p4id][2];
		phi = angle_btw_triangles(P1,P2,P3,P4);
		fprintf(f,"inter %d bending_force %e %e\n",firstID_BenBond + i, phi, coeff_bending);
		fprintf(fpart,"part %d bond %d %d %d %d\n",p2id + firstPartId, firstID_BenBond + i, p1id+firstPartId,p3id+firstPartId,p4id+firstPartId);
		//fprintf(fpartpos,"%d %d %d %d %e\n",p2id + firstPartId, p1id+firstPartId, p3id+firstPartId, p4id+firstPartId, phi);
		//fprintf(fpartpos,"       %e %e %e\n",mesh_nodes[p2id][0],mesh_nodes[p2id][1],mesh_nodes[p2id][2]);		
		//fprintf(fpartpos,"       %e %e %e\n",mesh_nodes[p1id][0],mesh_nodes[p1id][1],mesh_nodes[p1id][2]);		
		//fprintf(fpartpos,"       %e %e %e\n",mesh_nodes[p3id][0],mesh_nodes[p3id][1],mesh_nodes[p3id][2]);		
		//fprintf(fpartpos,"       %e %e %e\n",mesh_nodes[p4id][0],mesh_nodes[p4id][1],mesh_nodes[p4id][2]);		
	}
	//fclose(fpartpos);

	fclose(f); 
	fclose(fpart); 
	printf("Done.\n\n");

// Generation of local area force bonds:
    printf("Generating local area force bonds\n");
	firstID_localAreaBond = firstBondId;
	n_localAreaBond = mesh_ntriangles; // Area is coupled to the triangles
    firstBondId += n_localAreaBond;

	double area = 0.;
	f = fopen(bondsAreaLocal,"wt");
	fpart = fopen(partAreaLocal,"wt");
	for (i = 0; i < n_localAreaBond; i++) {
		for(k = 0; k < 3; k++){		
			P1[k]=mesh_nodes[mesh_triangles[i][0]][k];
			P2[k]=mesh_nodes[mesh_triangles[i][1]][k];
			P3[k]=mesh_nodes[mesh_triangles[i][2]][k];
		}
		area = area_triangle(P1,P2,P3);
		fprintf(f,"inter %d area_force_local %e %e\n",firstID_localAreaBond + i, area, coeff_area_local);
		fprintf(fpart,"part %d bond %d %d %d\n",mesh_triangles[i][0]+firstPartId, firstID_localAreaBond + i, mesh_triangles[i][1]+firstPartId, mesh_triangles[i][2]+firstPartId);
	}
	fclose(f); 
	fclose(fpart); 
	printf("Done.\n\n");
	
// Generation of global area force bonds:
    printf("Generating global area force bonds\n");
	firstID_globalAreaBond = firstBondId;
	n_globalAreaBond = 1;
    firstBondId += n_globalAreaBond;

	area = 0.;
	double gl_area = 0.;
	f = fopen(bondsAreaGlobal,"wt");
	fpart = fopen(partAreaGlobal,"wt");

	for (i = 0; i < mesh_ntriangles; i++) {
		for(k = 0; k < 3; k++){		
			P1[k]=mesh_nodes[mesh_triangles[i][0]][k];
			P2[k]=mesh_nodes[mesh_triangles[i][1]][k];
			P3[k]=mesh_nodes[mesh_triangles[i][2]][k];
		}
		area = area_triangle(P1,P2,P3);
		gl_area += area;
		fprintf(fpart,"part %d bond %d %d %d\n",mesh_triangles[i][0]+firstPartId, firstID_globalAreaBond, mesh_triangles[i][1]+firstPartId, mesh_triangles[i][2]+firstPartId);
	}

	fprintf(f,"inter %d area_force_global %e %e\n",firstID_globalAreaBond, gl_area, coeff_area_global);

	fclose(f); 
	fclose(fpart); 
	printf("Done.\n\n");
	


// Generation of volume force bonds:
    printf("Generating volume force bonds\n");
	firstID_VolumeBond = firstBondId;
	n_VolumeBond = 1;
    firstBondId += n_VolumeBond;

	area = 0.;
	double volume = 0., hz = 0., norm[3], dn=0.;
	//double drmax=0., drtemp=0;
	f = fopen(bondsVolume,"wt");
	fpart = fopen(partVolume,"wt");
	f2 = fopen(bondsVolumeAreaGlobal,"wt");
	f2part = fopen(partVolumeAreaGlobal,"wt");
	for (i = 0; i < mesh_ntriangles; i++) {
		for(k = 0; k < 3; k++){		
			P1[k]=mesh_nodes[mesh_triangles[i][0]][k];
			P2[k]=mesh_nodes[mesh_triangles[i][1]][k];
			P3[k]=mesh_nodes[mesh_triangles[i][2]][k];
		}
		//drtemp=sqrt((P1[0]-P2[0])*(P1[0]-P2[0])+(P1[1]-P2[1])*(P1[1]-P2[1])+(P1[2]-P2[2])*(P1[2]-P2[2])); //distance P1P2
		//if(drmax < drtemp) drmax=drtemp;
		//drtemp=sqrt((P1[0]-P3[0])*(P1[0]-P3[0])+(P1[1]-P3[1])*(P1[1]-P3[1])+(P1[2]-P3[2])*(P1[2]-P3[2])); //distance P1P3
		//if(drmax < drtemp) drmax=drtemp;
		
		area = area_triangle(P1,P2,P3);
		get_n_triangle(P1,P2,P3,norm);
		dn=sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
		hz=1.0/3.0 *(P1[2]+P2[2]+P3[2]);
		volume += area * norm[2]/dn * hz;
		//printf("i %d bondGen: partVol=%e A=%e hz=%e dn=%e\n", i, volume, area, hz, dn);
		fprintf(fpart,"part %d bond %d %d %d\n",mesh_triangles[i][0]+firstPartId, firstID_VolumeBond, mesh_triangles[i][1]+firstPartId, mesh_triangles[i][2]+firstPartId);
		fprintf(f2part,"part %d bond %d %d %d\n",mesh_triangles[i][0]+firstPartId, firstID_VolumeBond, mesh_triangles[i][1]+firstPartId, mesh_triangles[i][2]+firstPartId);
	}
	
	//fprintf(f,"inter %d volume_force %e %e %e\n",firstID_VolumeBond, volume, coeff_volume, 2*drmax);
	fprintf(f,"inter %d volume_force %e %e\n",firstID_VolumeBond, volume, coeff_volume);
	fprintf(f2,"inter %d volume_areagl_force %e %e %e %e\n",firstID_VolumeBond, volume, coeff_volume, gl_area, coeff_area_global);
	fclose(f); 
	fclose(fpart); 
	fclose(f2); 
	fclose(f2part);
	printf("Done. Volume = %lf. \n END OF MOL ID %d.\n", volume, molPart);

  return 0;

}
