/*
  Copyright (C) 2010,2011 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/** \file p3m-assignment.c  General functions related to the assignement functions for the P3M algorithms.
 *
*/



/************************************************************/
/*@{*/


#if ( defined(ELECTROSTATICS) || defined(MAGNETOSTATICS) ) && defined(ELP3M)

/* These functions are used for calculation of a charge assignment function */
double caf10(double x);

double caf20(double x);
double caf21(double x);

double caf30(double x);
double caf31(double x);
double caf32(double x);

double caf40(double x);
double caf41(double x);
double caf42(double x);
double caf43(double x);

double caf50(double x);
double caf51(double x);
double caf52(double x);
double caf53(double x);
double caf54(double x);

double caf60(double x);
double caf61(double x);
double caf62(double x);
double caf63(double x);
double caf64(double x);
double caf65(double x);

double caf70(double x);
double caf71(double x);
double caf72(double x);
double caf73(double x);
double caf74(double x);
double caf75(double x);
double caf76(double x);

/* the function to calculate caf without interpolation */
typedef double func(double);
func *int_caf_wi[7][7] = { {caf10},
                            {caf20, caf21},
			    {caf30, caf31, caf32},
			    {caf40, caf41, caf42, caf43},
			    {caf50, caf51, caf52, caf53, caf54},
			    {caf60, caf61, caf62, caf63, caf64, caf65},
			    {caf70, caf71, caf72, caf73, caf74, caf75, caf76}
                          };

/*@}*/

/************************************************************/

double caf10(double x)
{ double y;
   y = 1.0;
  return y;
}

double caf20(double x)
{ double y;
   y = 0.5-x;
  return y;
}

double caf21(double x)
{ double y;
   y = 0.5+x;
  return y;
}

double caf30(double x)
{ double y;
   y = 0.5*SQR(0.5 - x);
  return y;
}

double caf31(double x)
{ double y;
   y = 0.75 - SQR(x);
  return y;
}

double caf32(double x)
{ double y;
   y = 0.5*SQR(0.5 + x);
  return y;
}

double caf40(double x)
{ double y;
   y = ( 1.0+x*( -6.0+x*( 12.0-x* 8.0)))/48.0;
  return y;
}

double caf41(double x)
{ double y;
   y = (23.0+x*(-30.0+x*(-12.0+x*24.0)))/48.0;
  return y;
}

double caf42(double x)
{ double y;
   y = (23.0+x*( 30.0+x*(-12.0-x*24.0)))/48.0;
  return y;
}

double caf43(double x)
{ double y;
   y = ( 1.0+x*(  6.0+x*( 12.0+x* 8.0)))/48.0;
  return y;
}

double caf50(double x)
{ double y;
   y = (  1.0+x*( -8.0+x*(  24.0+x*(-32.0+x*16.0))))/384.0;
  return y;
}

double caf51(double x)
{ double y;
   y = ( 19.0+x*(-44.0+x*(  24.0+x*( 16.0-x*16.0))))/ 96.0;
  return y;
}

double caf52(double x)
{ double y;
   y = (115.0+x*       x*(-120.0+x*       x*48.0))  /192.0;
  return y;
}

double caf53(double x)
{ double y;
   y = ( 19.0+x*( 44.0+x*(  24.0+x*(-16.0-x*16.0))))/ 96.0;
  return y;
}

double caf54(double x)
{ double y;
   y = (  1.0+x*(  8.0+x*(  24.0+x*( 32.0+x*16.0))))/384.0;
  return y;
}

double caf60(double x)
{ double y;
   y = (  1.0+x*( -10.0+x*(  40.0+x*( -80.0+x*(  80.0-x* 32.0)))))/3840.0;
  return y;
}

double caf61(double x)
{ double y;
   y = (237.0+x*(-750.0+x*( 840.0+x*(-240.0+x*(-240.0+x*160.0)))))/3840.0;
  return y;
}

double caf62(double x)
{ double y;
   y = (841.0+x*(-770.0+x*(-440.0+x*( 560.0+x*(  80.0-x*160.0)))))/1920.0;
  return y;
}

double caf63(double x)
{ double y;
   y = (841.0+x*(+770.0+x*(-440.0+x*(-560.0+x*(  80.0+x*160.0)))))/1920.0;
  return y;
}

double caf64(double x)
{ double y;
   y = (237.0+x*( 750.0+x*( 840.0+x*( 240.0+x*(-240.0-x*160.0)))))/3840.0;
  return y;
}

double caf65(double x)
{ double y;
   y = (  1.0+x*(  10.0+x*(  40.0+x*(  80.0+x*(  80.0+x* 32.0)))))/3840.0;
  return y;
}

double caf70(double x)
{ double y;
   y = (    1.0+x*(   -12.0+x*(   60.0+x*( -160.0+x*(  240.0+x*(-192.0+x* 64.0))))))/46080.0;
  return y;
}

double caf71(double x)
{ double y;
   y = (  361.0+x*( -1416.0+x*( 2220.0+x*(-1600.0+x*(  240.0+x*( 384.0-x*192.0))))))/23040.0;
  return y;
}

double caf72(double x)
{ double y;
   y = (10543.0+x*(-17340.0+x*( 4740.0+x*( 6880.0+x*(-4080.0+x*(-960.0+x*960.0))))))/46080.0;
  return y;
}

double caf73(double x)
{ double y;
   y = ( 5887.0+x*          x*(-4620.0+x*         x*( 1680.0-x*        x*320.0)))   /11520.0;
  return y;
}

double caf74(double x)
{ double y;
   y = (10543.0+x*( 17340.0+x*( 4740.0+x*(-6880.0+x*(-4080.0+x*( 960.0+x*960.0))))))/46080.0;
  return y;
}

double caf75(double x)
{ double y;
   y = (  361.0+x*(  1416.0+x*( 2220.0+x*( 1600.0+x*(  240.0+x*(-384.0-x*192.0))))))/23040.0;
  return y;
}

double caf76(double x)
{ double y;
   y = (    1.0+x*(    12.0+x*(   60.0+x*(  160.0+x*(  240.0+x*( 192.0+x* 64.0))))))/46080.0;
  return y;
}


/************************************************************/



/** Computes the  assignment function of for the \a i'th degree
    at value \a x. */
double P3M_caf(int i, double x,int cao_value) {
  switch (cao_value) {
  case 1 : return 1.0;
  case 2 : {
    switch (i) {
    case 0: return 0.5-x;
    case 1: return 0.5+x;
    default:
      fprintf(stderr,"%d: Tried to access charge assignment function of degree %d in scheme of order %d.\n",this_node,i,cao_value);
      return 0.0;
    }
  } 
  case 3 : { 
    switch (i) {
    case 0: return 0.5*SQR(0.5 - x);
    case 1: return 0.75 - SQR(x);
    case 2: return 0.5*SQR(0.5 + x);
    default:
      fprintf(stderr,"%d: Tried to access charge assignment function of degree %d in scheme of order %d.\n",this_node,i,cao_value);
      return 0.0;
    }
  case 4 : { 
    switch (i) {
    case 0: return ( 1.0+x*( -6.0+x*( 12.0-x* 8.0)))/48.0;
    case 1: return (23.0+x*(-30.0+x*(-12.0+x*24.0)))/48.0;
    case 2: return (23.0+x*( 30.0+x*(-12.0-x*24.0)))/48.0;
    case 3: return ( 1.0+x*(  6.0+x*( 12.0+x* 8.0)))/48.0;
    default:
      fprintf(stderr,"%d: Tried to access charge assignment function of degree %d in scheme of order %d.\n",this_node,i,cao_value);
      return 0.0;
    }
  }
  case 5 : {
    switch (i) {
    case 0: return (  1.0+x*( -8.0+x*(  24.0+x*(-32.0+x*16.0))))/384.0;
    case 1: return ( 19.0+x*(-44.0+x*(  24.0+x*( 16.0-x*16.0))))/ 96.0;
    case 2: return (115.0+x*       x*(-120.0+x*       x*48.0))  /192.0;
    case 3: return ( 19.0+x*( 44.0+x*(  24.0+x*(-16.0-x*16.0))))/ 96.0;
    case 4: return (  1.0+x*(  8.0+x*(  24.0+x*( 32.0+x*16.0))))/384.0;
    default:
      fprintf(stderr,"%d: Tried to access charge assignment function of degree %d in scheme of order %d.\n",this_node,i,cao_value);
      return 0.0;
    }
  }
  case 6 : {
    switch (i) {
    case 0: return (  1.0+x*( -10.0+x*(  40.0+x*( -80.0+x*(  80.0-x* 32.0)))))/3840.0;
    case 1: return (237.0+x*(-750.0+x*( 840.0+x*(-240.0+x*(-240.0+x*160.0)))))/3840.0;
    case 2: return (841.0+x*(-770.0+x*(-440.0+x*( 560.0+x*(  80.0-x*160.0)))))/1920.0;
    case 3: return (841.0+x*(+770.0+x*(-440.0+x*(-560.0+x*(  80.0+x*160.0)))))/1920.0;
    case 4: return (237.0+x*( 750.0+x*( 840.0+x*( 240.0+x*(-240.0-x*160.0)))))/3840.0;
    case 5: return (  1.0+x*(  10.0+x*(  40.0+x*(  80.0+x*(  80.0+x* 32.0)))))/3840.0;
    default:
      fprintf(stderr,"%d: Tried to access charge assignment function of degree %d in scheme of order %d.\n",this_node,i,cao_value);
      return 0.0;
    }
  }
  case 7 : {
    switch (i) {
    case 0: return (    1.0+x*(   -12.0+x*(   60.0+x*( -160.0+x*(  240.0+x*(-192.0+x* 64.0))))))/46080.0;
    case 1: return (  361.0+x*( -1416.0+x*( 2220.0+x*(-1600.0+x*(  240.0+x*( 384.0-x*192.0))))))/23040.0;
    case 2: return (10543.0+x*(-17340.0+x*( 4740.0+x*( 6880.0+x*(-4080.0+x*(-960.0+x*960.0))))))/46080.0;
    case 3: return ( 5887.0+x*          x*(-4620.0+x*         x*( 1680.0-x*        x*320.0)))   /11520.0;
    case 4: return (10543.0+x*( 17340.0+x*( 4740.0+x*(-6880.0+x*(-4080.0+x*( 960.0+x*960.0))))))/46080.0;
    case 5: return (  361.0+x*(  1416.0+x*( 2220.0+x*( 1600.0+x*(  240.0+x*(-384.0-x*192.0))))))/23040.0;
    case 6: return (    1.0+x*(    12.0+x*(   60.0+x*(  160.0+x*(  240.0+x*( 192.0+x* 64.0))))))/46080.0;
    default:
      fprintf(stderr,"%d: Tried to access charge assignment function of degree %d in scheme of order %d.\n",this_node,i,cao_value);
      return 0.0;
    }
  }
  default :{
    fprintf(stderr,"%d: Charge assignment order %d unknown.\n",this_node,cao_value);
    return 0.0;
  }}}
}




#endif /* if ( defined(ELECTROSTATICS) || defined(MAGNETOSTATICS) ) && defined(ELP3M) */

