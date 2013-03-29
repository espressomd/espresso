#include "config.h"
#include "electrokinetics_tcl.h"
#include "electrokinetics.h"
#include "lb-boundaries.h"
#include "initialize.h"

#ifdef EK_BOUNDARIES

LB_Boundary* generate_lbboundary(float charge_density)
{
  n_lb_boundaries++;
  
  lb_boundaries = realloc(lb_boundaries,n_lb_boundaries*sizeof(LB_Boundary));
  lb_boundaries[n_lb_boundaries-1].type = LB_BOUNDARY_BOUNCE_BACK;
  
  lb_boundaries[n_lb_boundaries-1].velocity[0] =
  lb_boundaries[n_lb_boundaries-1].velocity[1] =
  lb_boundaries[n_lb_boundaries-1].velocity[2] = 0;
  
  lb_boundaries[n_lb_boundaries-1].force[0] =
  lb_boundaries[n_lb_boundaries-1].force[1] =
  lb_boundaries[n_lb_boundaries-1].force[2] = 0;
  
  lb_boundaries[n_lb_boundaries-1].charge_density = charge_density;
  
  return &lb_boundaries[n_lb_boundaries-1];
}

int tclcommand_lbboundary_wall(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv)
{
  int i;
  double norm;
  
  lbb->type = LB_BOUNDARY_WAL;
  
  /* invalid entries to start of */
  lbb->c.wal.n[0] = 
  lbb->c.wal.n[1] = 
  lbb->c.wal.n[2] = 0;
  
  lbb->c.wal.d = 0;
  
  while(argc > 0) {
    if(!strncmp(argv[0], "normal", strlen(argv[0]))) {
      if(argc < 4) {
				Tcl_AppendResult(interp, "lbboundary wall normal <nx> <ny> <nz> expected", (char *) NULL);
			  return (TCL_ERROR);
      }
    
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.wal.n[0])) == TCL_ERROR ||
	    	 Tcl_GetDouble(interp, argv[2], &(lbb->c.wal.n[1])) == TCL_ERROR ||
	    	 Tcl_GetDouble(interp, argv[3], &(lbb->c.wal.n[2])) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "dist", strlen(argv[0]))) {
      if (argc < 1) {
        Tcl_AppendResult(interp, "lbboundary wall dist <d> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if (Tcl_GetDouble(interp, argv[1], &(lbb->c.wal.d)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary wall type <t> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "velocity", strlen(argv[0]))) {
      if(argc < 4) {
	      Tcl_AppendResult(interp, "lbboundary wall velocity <vx> <vy> <vz> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(lbb->velocity[0])) == TCL_ERROR ||
      	 Tcl_GetDouble(interp, argv[2], &(lbb->velocity[1])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[3], &(lbb->velocity[2])) == TCL_ERROR)
	      return (TCL_ERROR);
      else {
#ifdef LB
        lbb->velocity[0]*=lbpar.tau/lbpar.agrid;
        lbb->velocity[1]*=lbpar.tau/lbpar.agrid;
        lbb->velocity[2]*=lbpar.tau/lbpar.agrid;
#endif
      }
      
      argc -= 4; argv += 4;
    }
    else
      break;
  }
  
  /* length of the normal vector */
  norm = SQR(lbb->c.wal.n[0])+SQR(lbb->c.wal.n[1])+SQR(lbb->c.wal.n[2]);
  
  if (norm < 1e-10) {
    Tcl_AppendResult(interp, "usage: lbboundary wall normal <nx> <ny> <nz> dist <d> type <t>", (char *) NULL);
    return (TCL_ERROR);    
  }
  
  /* normalize the normal vector */
  for (i=0;i<3;i++) lbb->c.wal.n[i] /= sqrt(norm);

  return (TCL_OK);
}

int tclcommand_lbboundary_sphere(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv)
{
  lbb->type = LB_BOUNDARY_SPH;

  /* invalid entries to start of */
  lbb->c.sph.pos[0] = 
  lbb->c.sph.pos[1] = 
  lbb->c.sph.pos[2] = 0;
  lbb->c.sph.rad = 0;
  lbb->c.sph.direction = -1;

  while (argc > 0) {
    if(!strncmp(argv[0], "center", strlen(argv[0]))) {
      if(argc < 4) {
	      Tcl_AppendResult(interp, "lbboundary sphere center <x> <y> <z> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.sph.pos[0])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[2], &(lbb->c.sph.pos[1])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[3], &(lbb->c.sph.pos[2])) == TCL_ERROR)
	      return (TCL_ERROR);
	      
        argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "radius", strlen(argv[0]))) {
      if(argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary sphere radius <r> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if (Tcl_GetDouble(interp, argv[1], &(lbb->c.sph.rad)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "direction", strlen(argv[0]))) {
      if (argc < 1) {
	      Tcl_AppendResult(interp, "-1/1 or inside/outside is expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(!strncmp(argv[1], "inside", strlen(argv[1])))
	      lbb->c.sph.direction = -1;
      else if(!strncmp(argv[1], "outside", strlen(argv[1])))
	      lbb->c.sph.direction = 1;
      else if(Tcl_GetDouble(interp, argv[1], &(lbb->c.sph.direction)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary sphere type <t> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  if(lbb->c.sph.rad < 0.) {
    Tcl_AppendResult(interp, "usage: lbboundary sphere center <x> <y> <z> radius <d> direction <direction> type <t>", (char *) NULL);
    return (TCL_ERROR);    
  }

  return (TCL_OK);
}

int tclcommand_lbboundary_cylinder(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv)
{
  double axis_len;
  int i;

  lbb->type = LB_BOUNDARY_CYL;
  
  /* invalid entries to start of */
  lbb->c.cyl.pos[0] = 
  lbb->c.cyl.pos[1] = 
  lbb->c.cyl.pos[2] = 0;
  
  lbb->c.cyl.axis[0] = 
  lbb->c.cyl.axis[1] = 
  lbb->c.cyl.axis[2] = 0;
  
  lbb->c.cyl.rad = 0;
  
  lbb->c.cyl.length = 0;
  
  lbb->c.cyl.direction = 0;
  
  while (argc > 0) {
    if(!strncmp(argv[0], "center", strlen(argv[0]))) {
      if(argc < 4) {
	      Tcl_AppendResult(interp, "lbboundary cylinder center <x> <y> <z> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.cyl.pos[0])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[2], &(lbb->c.cyl.pos[1])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[3], &(lbb->c.cyl.pos[2])) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "axis", strlen(argv[0]))) {
      if(argc < 4) {
	      Tcl_AppendResult(interp, "lbboundary cylinder axis <rx> <ry> <rz> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.cyl.axis[0])) == TCL_ERROR ||
  	     Tcl_GetDouble(interp, argv[2], &(lbb->c.cyl.axis[1])) == TCL_ERROR ||
    	   Tcl_GetDouble(interp, argv[3], &(lbb->c.cyl.axis[2])) == TCL_ERROR)
	      return (TCL_ERROR);

      argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "radius", strlen(argv[0]))) {
      if(argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary cylinder radius <rad> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.cyl.rad)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "length", strlen(argv[0]))) {
      if(argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary cylinder length <len> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.cyl.length)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "direction", strlen(argv[0]))) {
      if(argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary cylinder direction <dir> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if (!strncmp(argv[1], "inside", strlen(argv[1])))
	      lbb->c.cyl.direction = -1;
      else if (!strncmp(argv[1], "outside", strlen(argv[1])))
	      lbb->c.cyl.direction = 1;
      else if (Tcl_GetDouble(interp, argv[1], &(lbb->c.cyl.direction)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary cylinder type <t> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  axis_len=0.;
  
  for (i=0;i<3;i++)
    axis_len += SQR(lbb->c.cyl.axis[i]);

  if(lbb->c.cyl.rad < 0. || axis_len < 1e-30 ||
     lbb->c.cyl.direction == 0 || lbb->c.cyl.length <= 0) {
    Tcl_AppendResult(interp, "usage: lbboundary cylinder center <x> <y> <z> axis <rx> <ry> <rz> radius <rad> length <length> direction <direction> type <t>", (char *) NULL);
    return (TCL_ERROR);    
  }

  /*normalize the axis vector */
  axis_len = sqrt (axis_len);
  
  for (i=0;i<3;i++) {
    lbb->c.cyl.axis[i] /= axis_len;
  }
      
  return (TCL_OK);
}

int tclcommand_lbboundary_rhomboid(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv)
{
  double triple_product;
  double tmp[3];
  
  lbb->type = LB_BOUNDARY_RHOMBOID;
  
  lbb->c.rhomboid.pos[0] = 
  lbb->c.rhomboid.pos[1] = 
  lbb->c.rhomboid.pos[2] = 0;
  
  lbb->c.rhomboid.a[0] = 
  lbb->c.rhomboid.a[1] = 
  lbb->c.rhomboid.a[2] = 0;
  
  lbb->c.rhomboid.b[0] = 
  lbb->c.rhomboid.b[1] = 
  lbb->c.rhomboid.b[2] = 0;
  
  lbb->c.rhomboid.c[0] = 
  lbb->c.rhomboid.c[1] = 
  lbb->c.rhomboid.c[2] = 0;
  
  lbb->c.rhomboid.direction = 0;
  
  while(argc > 0) {
    if(!strncmp(argv[0], "a", strlen(argv[0]))) {
      if(argc < 4) {
				Tcl_AppendResult(interp, "lbboundary rhomboid a <ax> <ay> <az> expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.rhomboid.a[0])) == TCL_ERROR ||
	 			 Tcl_GetDouble(interp, argv[2], &(lbb->c.rhomboid.a[1])) == TCL_ERROR ||
	  		 Tcl_GetDouble(interp, argv[3], &(lbb->c.rhomboid.a[2])) == TCL_ERROR)
				return TCL_ERROR;
				
			argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "b", strlen(argv[0]))) {
      if(argc < 4) {
				Tcl_AppendResult(interp, "lbboundary rhomboid b <bx> <by> <bz> expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.rhomboid.b[0])) == TCL_ERROR ||
	 			 Tcl_GetDouble(interp, argv[2], &(lbb->c.rhomboid.b[1])) == TCL_ERROR ||
	  		 Tcl_GetDouble(interp, argv[3], &(lbb->c.rhomboid.b[2])) == TCL_ERROR)
				return TCL_ERROR;
				
			argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "c", strlen(argv[0]))) {
      if(argc < 4) {
				Tcl_AppendResult(interp, "lbboundary rhomboid c <cx> <cy> <cz> expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.rhomboid.c[0])) == TCL_ERROR ||
	 			 Tcl_GetDouble(interp, argv[2], &(lbb->c.rhomboid.c[1])) == TCL_ERROR ||
	  		 Tcl_GetDouble(interp, argv[3], &(lbb->c.rhomboid.c[2])) == TCL_ERROR)
				return TCL_ERROR;
				
			argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "corner", strlen(argv[0]))) { //this has to come after c
      if(argc < 4) {
				Tcl_AppendResult(interp, "lbboundary rhomboid corner <x> <y> <z> expected", (char *) NULL);
				return TCL_ERROR;
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.rhomboid.pos[0])) == TCL_ERROR ||
				 Tcl_GetDouble(interp, argv[2], &(lbb->c.rhomboid.pos[1])) == TCL_ERROR ||
	  		 Tcl_GetDouble(interp, argv[3], &(lbb->c.rhomboid.pos[2])) == TCL_ERROR)
				return TCL_ERROR;
				
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "direction", strlen(argv[0]))) {
      if (argc < 2) {
				Tcl_AppendResult(interp, "lbboundary rhomboid direction {inside|outside} expected", (char *) NULL);
				return (TCL_ERROR);
      }
      
      if(!strncmp(argv[1], "inside", strlen(argv[1])))
				lbb->c.rhomboid.direction = -1;
      else if(!strncmp(argv[1], "outside", strlen(argv[1])))
				lbb->c.rhomboid.direction = 1;
      else if(Tcl_GetDouble(interp, argv[1], &(lbb->c.rhomboid.direction)) == TCL_ERROR)
				return TCL_ERROR;
				
      argc -= 2; argv += 2;
    }
    else {
			Tcl_AppendResult(interp, "Error: Unknown parameter ", argv[0], " in lbboundary rhomboid", (char *) NULL);
			return TCL_ERROR;
    }
  }

  if( (lbb->c.rhomboid.a[0] == 0. && lbb->c.rhomboid.a[1] == 0. && lbb->c.rhomboid.a[2] == 0.) ||
  		(lbb->c.rhomboid.b[0] == 0. && lbb->c.rhomboid.b[1] == 0. && lbb->c.rhomboid.b[2] == 0.) ||
  		(lbb->c.rhomboid.c[0] == 0. && lbb->c.rhomboid.c[1] == 0. && lbb->c.rhomboid.c[2] == 0.) ||
  		lbb->c.rhomboid.direction == 0) {
    Tcl_AppendResult(interp, "usage: lbboundary rhomboid corner <x> <y> <z> a <ax> <ay> <az> b <bx> <by> <bz> c <cx> <cy> <cz> direction {inside|outside} [penetrable <0|1>] [reflecting <1|2>]", (char *) NULL);
    return TCL_ERROR;
  }
                     
  //If the trihedron a, b, c is left handed, then inside and outside will be exchanged since all normals will be reversed. This compensates  for that, so that the user doesn't have to take care of the handedness.
  triple_product = lbb->c.rhomboid.a[0]*( lbb->c.rhomboid.b[1]*lbb->c.rhomboid.c[2] - lbb->c.rhomboid.b[2]*lbb->c.rhomboid.c[1] ) +
                lbb->c.rhomboid.a[1]*( lbb->c.rhomboid.b[2]*lbb->c.rhomboid.c[0] - lbb->c.rhomboid.b[0]*lbb->c.rhomboid.c[2] ) + 
                lbb->c.rhomboid.a[2]*( lbb->c.rhomboid.b[0]*lbb->c.rhomboid.c[1] - lbb->c.rhomboid.b[1]*lbb->c.rhomboid.c[0] );
                
  if(triple_product < 0.)
  {    
    tmp[0] = lbb->c.rhomboid.a[0];
    tmp[1] = lbb->c.rhomboid.a[1];
    tmp[2] = lbb->c.rhomboid.a[2];
    
    lbb->c.rhomboid.a[0] = lbb->c.rhomboid.b[0];
    lbb->c.rhomboid.a[1] = lbb->c.rhomboid.b[1];
    lbb->c.rhomboid.a[2] = lbb->c.rhomboid.b[2];
    
    lbb->c.rhomboid.b[0] = tmp[0];
    lbb->c.rhomboid.b[1] = tmp[1];
    lbb->c.rhomboid.b[2] = tmp[2];
  }

  return TCL_OK;
}

int tclcommand_lbboundary_pore(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv)
{
  double axis_len;
  int i;

  lbb->type = LB_BOUNDARY_POR;
  
  /* invalid entries to start of */
  lbb->c.pore.pos[0] = 
  lbb->c.pore.pos[1] = 
  lbb->c.pore.pos[2] = 0;
  
  lbb->c.pore.axis[0] = 
  lbb->c.pore.axis[1] = 
  lbb->c.pore.axis[2] = 0;
  
  lbb->c.pore.rad_left = 0;
  
  lbb->c.pore.rad_right = 0;
  
  lbb->c.pore.length = 0;
  
  lbb->c.pore.reflecting = 0;
  
  lbb->c.pore.smoothing_radius = 1.;
  
  while (argc > 0) {
    if(!strncmp(argv[0], "center", strlen(argv[0]))) {
      if(argc < 4) {
	      Tcl_AppendResult(interp, "lbboundary pore center <x> <y> <z> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.pore.pos[0])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[2], &(lbb->c.pore.pos[1])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[3], &(lbb->c.pore.pos[2])) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "axis", strlen(argv[0]))) {
      if(argc < 4) {
	      Tcl_AppendResult(interp, "lbboundary pore axis <rx> <ry> <rz> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.pore.axis[0])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[2], &(lbb->c.pore.axis[1])) == TCL_ERROR ||
	       Tcl_GetDouble(interp, argv[3], &(lbb->c.pore.axis[2])) == TCL_ERROR)
	      return (TCL_ERROR);

      argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "radius", strlen(argv[0]))) {
      if(argc < 1) {
    	  Tcl_AppendResult(interp, "lbboundary pore radius <rad> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if(Tcl_GetDouble(interp, argv[1], &(lbb->c.pore.rad_left)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      lbb->c.pore.rad_right =  lbb->c.pore.rad_left; 
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "smoothing_radius", strlen(argv[0]))) {
      if (argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary pore smoothing_radius <smoothing_radius> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if (Tcl_GetDouble(interp, argv[1], &(lbb->c.pore.smoothing_radius)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "radii", strlen(argv[0]))) {
      if(argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary pore radii <rad_left> <rad_right> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if (Tcl_GetDouble(interp, argv[1], &(lbb->c.pore.rad_left)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      if (Tcl_GetDouble(interp, argv[2], &(lbb->c.pore.rad_right)) == TCL_ERROR)
	      return (TCL_ERROR);
	      
      argc -= 3; argv += 3;
    }
    else if(!strncmp(argv[0], "length", strlen(argv[0]))) {
      if (argc < 1) {
	      Tcl_AppendResult(interp, "lbboundary pore length <len/2> expected", (char *) NULL);
	      return (TCL_ERROR);
      }
      
      if (Tcl_GetDouble(interp, argv[1], &(lbb->c.pore.length)) == TCL_ERROR)
	      return (TCL_ERROR);

      argc -= 2; argv += 2;
    }
    else
      break;
  }

  axis_len=0.;
  
  for(i=0;i<3;i++)
    axis_len += SQR(lbb->c.pore.axis[i]);

  if(lbb->c.pore.rad_left < 0. || lbb->c.pore.rad_right < 0. || axis_len < 1e-30 || lbb->c.pore.length <= 0) {
    Tcl_AppendResult(interp, "usage: lbboundary pore center <x> <y> <z> axis <rx> <ry> <rz> radius <rad> length <length/2>", (char *) NULL);
    return (TCL_ERROR);
  }

  /*normalize the axis vector */
  axis_len = sqrt (axis_len);
  
  for (i=0;i<3;i++) {
    lbb->c.pore.axis[i] /= axis_len;
  }
  
  return (TCL_OK);
}

#endif /* EK_BOUNDARIES */

int tclcommand_electrokinetics(ClientData data, Tcl_Interp *interp, int argc, char **argv) {
#ifndef ELECTROKINETICS
  Tcl_AppendResult(interp, "Feature ELECTROKINETICS required", NULL);
  return TCL_ERROR;
#else

  argc--;
  argv++;

  int err = TCL_OK;
  int species;
  double floatarg;
  double vectarg[3];
  char int_buffer[TCL_INTEGER_SPACE];

  if(argc < 2) {
    Tcl_AppendResult(interp, "Usage of \"electrokinetics\":", (char *)NULL);
    Tcl_AppendResult(interp, "electrokinetics [agrid #float] [viscosity #float] [friction #float]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                [bulk_viscosity #float] [gamma_even #float] [gamma_odd #float]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                [print <density|velocity|potential|lbforce> vtk #string]\n", (char *)NULL);
    Tcl_AppendResult(interp, "electrokinetics boundary charge_density #float [wall ]\n", (char *)NULL); //TODO full description
    Tcl_AppendResult(interp, "                                               [sphere ]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                                               [cylinder ]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                                               [rhomboid ]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                                               [pore ]\n", (char *)NULL);
    Tcl_AppendResult(interp, "electrokinetics #int [density #float] [D #float] [T #float] [valency #float]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                     [ext_force #float #float #float]\n", (char *)NULL);
    Tcl_AppendResult(interp, "                     [print density vtk #string]\n", (char *)NULL);
    
    return TCL_ERROR;
  }
  else if(ARG0_IS_S("boundary")) {
#ifndef EK_BOUNDARIES
    Tcl_AppendResult(interp, "Feature EK_BOUNDARIES required", (char *) NULL);
    return (TCL_ERROR);
#else
    argc--;
    argv++;
    
    if(!ARG0_IS_S("charge_density") || !ARG1_IS_D(floatarg)) {
      Tcl_AppendResult(interp, "You need to specify the boundary charge density using\n", (char *) NULL);
      Tcl_AppendResult(interp, "electrokinetics boundary charge_density #float ...\n", (char *)NULL);
      return (TCL_ERROR);
    }
    
    argc -= 2;
    argv += 2;
    
    if(floatarg != 0.0) {
      species = -1;
      
      for(int i = 0; i < ek_parameters.number_of_species; i++)
        if(ek_parameters.valency[i] != 0.0) {
          species = i;
          break;
        }
      
      if(species == -1) {
        Tcl_AppendResult(interp, "You need to define at least one charged species in order to use charged walls\n", (char *) NULL);
        return (TCL_ERROR);
      }
    }
    
    int status = TCL_ERROR, c_num;
    
    if(ARG0_IS_S("wall"))
      status = tclcommand_lbboundary_wall(generate_lbboundary(floatarg), interp, argc - 1, argv + 1);
    else if(ARG0_IS_S("sphere"))
      status = tclcommand_lbboundary_sphere(generate_lbboundary(floatarg), interp, argc - 1, argv + 1);
    else if(ARG0_IS_S("cylinder"))
      status = tclcommand_lbboundary_cylinder(generate_lbboundary(floatarg), interp, argc - 1, argv + 1);
    else if(ARG0_IS_S("rhomboid"))
      status = tclcommand_lbboundary_rhomboid(generate_lbboundary(floatarg), interp, argc - 1, argv + 1);
    else if(ARG0_IS_S("pore"))
      status = tclcommand_lbboundary_pore(generate_lbboundary(floatarg), interp, argc - 1, argv + 1);
    else if(ARG0_IS_S("delete")) {
      if(argc < 3) {
        /* delete all */
        Tcl_AppendResult(interp, "Can only delete individual electrokinetics boundaries", (char *) NULL);
        status = TCL_ERROR;
      }
      else {
        if(Tcl_GetInt(interp, argv[2], &(c_num)) == TCL_ERROR)
          return (TCL_ERROR);
          
        if(c_num < 0 || c_num >= n_lb_boundaries) {
	        Tcl_AppendResult(interp, "Can not delete non existing electrokinetics boundary", (char *) NULL);
	        return (TCL_ERROR);
        }
      }
    }
    else {
      Tcl_AppendResult(interp, "possible electrokinetics boundary charge_density #float parameters: wall, sphere, cylinder, rhomboid, pore, delete {c} to delete a boundary", (char *) NULL);
      return (TCL_ERROR);
    }
        
    on_lbboundary_change();
#endif /* EK_BOUNDARIES */
  }
  else if(ARG0_IS_I(species)) {
    argc--;
    argv++;
    
    if(species < 0 || species > MAX_NUMBER_OF_SPECIES) {
      sprintf(int_buffer, "%d", MAX_NUMBER_OF_SPECIES);
      Tcl_AppendResult(interp, "electrokinetics #int requires a number between 0 and", int_buffer, "denoting the species\n", (char *)NULL);
      return TCL_ERROR;
    }
    
    while(argc > 0) {
      if(ARG0_IS_S("D")) {
        if(argc < 2 || !ARG1_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics #int D requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(floatarg < 0) {
          Tcl_AppendResult(interp, "electrokinetics #int D can not be negative\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_D(species, floatarg) == 0) {
            argc -= 2;
            argv += 2;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics #int D\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("density")) {
        if(argc < 2 || !ARG1_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics #int density requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(floatarg < 0) {
          Tcl_AppendResult(interp, "electrokinetics #int density must be positive\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_density(species, floatarg) == 0) {
            argc -= 2;
            argv += 2;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics #int density\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("ext_force")) {
  #ifndef EXTERNAL_FORCES
        Tcl_AppendResult(interp, "EXTERNAL_FORCES not compiled in\n", (char *)NULL);
         return TCL_ERROR;
  #else
        if(argc < 4 || !ARG_IS_D(1, vectarg[0]) || !ARG_IS_D(2, vectarg[1]) || !ARG_IS_D(3, vectarg[2])) {
          Tcl_AppendResult(interp, "electrokinetics #int ext_force requires three floating point numbers as arguments\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(ek_set_ext_force(species, vectarg[0], vectarg[1], vectarg[2]) == 0) {
          argc -= 4;
          argv += 4;
        }
        else {
          Tcl_AppendResult(interp, "Unknown error setting electrokinetics #int ext_force\n", (char *)NULL);
          return TCL_ERROR;
        }
  #endif
      }
      else if(ARG0_IS_S("valency")) {
        if(argc < 2 || !ARG1_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics #int valency requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_valency(species, floatarg) == 0) {
            argc -= 2;
            argv += 2;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics #int valency\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("print")) {
        argc--;
        argv++;
        
        if(argc != 3 || !ARG1_IS_S("vtk") || !ARG0_IS_S("density")) {
          Tcl_AppendResult(interp, "Wrong usage of electrokinetics #int print density vtk #string\n", (char *)NULL);
          return TCL_ERROR;
        }
        
        if(ARG0_IS_S("density")) {
          if(ek_print_vtk_density(species, argv[2]) == 0) {
            argc -= 3;
            argv += 3;

            if((err = gather_runtime_errors(interp, err)) != TCL_OK)
              return TCL_ERROR;
            else
              return TCL_OK;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error in electrokinetics #int print density vtk #string\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else {
    	  Tcl_AppendResult(interp, "unknown feature \"", argv[0],"\" of electrokinetics #int\n", (char *)NULL);
    	  return TCL_ERROR ;
      }
    }
  }
  else {
    Tcl_ResetResult(interp);
    
    while(argc > 0) {
      if(ARG0_IS_S("print")) {
        argc--;
        argv++;
        
        if(argc != 3 || !ARG1_IS_S("vtk") || (!ARG0_IS_S("velocity") && !ARG0_IS_S("density") && !ARG0_IS_S("boundary") && !ARG0_IS_S("potential") && !ARG0_IS_S("lbforce"))) {
          Tcl_AppendResult(interp, "Wrong usage of electrokinetics print <velocity|density|potential> vtk #string\n", (char *)NULL);
          return TCL_ERROR;
        }
        
        if(ARG0_IS_S("velocity")) {
          if(ek_lb_print_vtk_velocity(argv[2]) == 0) {
            argc -= 3;
            argv += 3;

            if((err = gather_runtime_errors(interp, err)) != TCL_OK)
              return TCL_ERROR;
            else
              return TCL_OK;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error in electrokinetics print velocity vtk #string\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
        else if(ARG0_IS_S("density")) {
          if(ek_lb_print_vtk_density(argv[2]) == 0) {
            argc -= 3;
            argv += 3;

            if((err = gather_runtime_errors(interp, err)) != TCL_OK)
              return TCL_ERROR;
            else
              return TCL_OK;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error in electrokinetics print density vtk #string\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
        else if(ARG0_IS_S("boundary")) {
#ifndef EK_BOUNDARIES
          Tcl_AppendResult(interp, "Feature EK_BOUNDARIES required", (char *) NULL);
          return (TCL_ERROR);
#else
          if(lb_lbfluid_print_vtk_boundary(argv[2]) == 0) {
            argc -= 3;
            argv += 3;

            if((err = gather_runtime_errors(interp, err)) != TCL_OK)
              return TCL_ERROR;
            else
              return TCL_OK;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error in electrokinetics print boundary vtk #string\n", (char *)NULL);
            return TCL_ERROR;
          }
#endif /* EK_BOUNDARIES */
        }
        else if(ARG0_IS_S("potential")) {
          if(ek_print_vtk_potential(argv[2]) == 0) {
            argc -= 3;
            argv += 3;

            if((err = gather_runtime_errors(interp, err)) != TCL_OK)
              return TCL_ERROR;
            else
              return TCL_OK;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error in electrokinetics print lbforce vtk #string\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
        else if(ARG0_IS_S("lbforce")) {
          if(ek_print_vtk_lbforce(argv[2]) == 0) {
            argc -= 3;
            argv += 3;

            if((err = gather_runtime_errors(interp, err)) != TCL_OK)
              return TCL_ERROR;
            else
              return TCL_OK;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error in electrokinetics print potential vtk #string\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
        else {
            Tcl_AppendResult(interp, "Unknown feature \"", argv[0], "\" in electrokinetics print\n", (char *)NULL);
            return TCL_ERROR;
        }
      }
      else if(ARG0_IS_S("T")) {
        argc--;
        argv++;
        
        if(argc < 1 || !ARG0_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics T requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(floatarg <= 0) {
          Tcl_AppendResult(interp, "electrokinetics T must be positive\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_T(floatarg) == 0) {
            argc --;
            argv ++;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics T\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("bjerrum_length")) {
        argc--;
        argv++;
      
        if(!ARG0_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics bjerrum_length requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(floatarg <= 0) {
          Tcl_AppendResult(interp, "electrokinetics bjerrum_length must be positive\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_bjerrumlength(floatarg) == 0) {
            argc--;
            argv++;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics bjerrum_length\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("agrid")) {
        argc--;
        argv++;
      
        if(!ARG0_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics agrid requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(floatarg <= 0) {
          Tcl_AppendResult(interp, "electrokinetics agrid must be positive\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_agrid(floatarg) == 0) {
            argc--;
            argv++;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics agrid\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("viscosity")) {
        if(argc < 2 || !ARG1_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics viscosity requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        } else if (floatarg <= 0) {
          Tcl_AppendResult(interp, "electrokinetics viscosity must be positive\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_viscosity(floatarg) == 0) {
            argc -= 2;
            argv += 2;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics viscosity\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("friction")) {
        if(argc < 2 || !ARG1_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics friction requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(floatarg <= 0) {
          Tcl_AppendResult(interp, "electrokinetics friction must be positive\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if (ek_set_friction(floatarg) == 0) {
            argc -= 2;
            argv += 2;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics friction\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("bulk_viscosity")) {
        if(argc < 2 || !ARG1_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics bulk_viscosity requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(floatarg <= 0) {
          Tcl_AppendResult(interp, "electrokinetics bulk_viscosity must be positive\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_bulk_viscosity(floatarg) == 0) {
            argc -= 2;
            argv += 2;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics bulk_viscosity\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("gamma_odd")) {
        if(argc < 2 || !ARG1_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics gamma_odd requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(fabs(floatarg) >= 1) {
          Tcl_AppendResult(interp, "electrokinetics gamma_odd must be smaller than 1\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_gamma_odd(floatarg) == 0) {
            argc -= 2;
            argv += 2;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics gamma_odd\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if(ARG0_IS_S("gamma_even")) {
        if(argc < 2 || !ARG1_IS_D(floatarg)) {
          Tcl_AppendResult(interp, "electrokinetics gamma_even requires one floating point number as argument\n", (char *)NULL);
          return TCL_ERROR;
        }
        else if(fabs(floatarg) >= 1) {
          Tcl_AppendResult(interp, "electrokinetics gamma_even must be smaller than 1\n", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          if(ek_set_gamma_even(floatarg) == 0) {
            argc -= 2;
            argv += 2;
          }
          else {
            Tcl_AppendResult(interp, "Unknown error setting electrokinetics gamma_even\n", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else {
    	  Tcl_AppendResult(interp, "unknown feature \"", argv[0],"\" of electrokinetics\n", (char *)NULL);
    	  return TCL_ERROR ;
      }
    }
  }

  if((err = gather_runtime_errors(interp, err)) != TCL_OK)
    return TCL_ERROR;
  
  err = ek_init();
  
  switch(err) {
    case 0:
      return TCL_OK;
      
    case 2:
      Tcl_AppendResult(interp, "electrokinetics agrid not compatible with box size\n", (char *)NULL);
      
    default:
      sprintf(int_buffer, "%d", err);
      Tcl_AppendResult(interp, "Error ", int_buffer, " during initialization of electrokinetics", (char *)NULL);
      break;
  }
  
  return TCL_ERROR;
  
#endif /* defined ELECTROKINETICS */
}
