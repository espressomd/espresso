#ifndef LEES_EDWARDS_H
#define LEES_EDWARDS_H

#include "config.hpp"

/** \file lees_erdwards.hpp
*
*/

/** Enum for the different Lees Edwards Protocols: Off, Step. Steady Shear and Oscillatory shear  */
enum LeesEdwardsProtocolType {
    LEES_EDWARDS_PROTOCOL_OFF,
    LEES_EDWARDS_PROTOCOL_STEP,
    LEES_EDWARDS_PROTOCOL_STEADY_SHEAR,
    LEES_EDWARDS_PROTOCOL_OSC_SHEAR,
};

/** Struct holding all information concerning Lees Edwards  */
typedef struct {
/** Protocol type*/
  int type;
/** Time when Lees Edwards was started*/
  double time0;
/** Current offset*/
  double offset;
/** Current velocity*/
  double velocity;
/** Amplitude set via interface*/
  double amplitude;
/** Frequency set via interface*/
  double frequency;
/** Direction in which the velocity and position jump is applied*/
  int sheardir;
/** Direction in which get_mi_vector needs to be modified*/
  int shearplanenormal; 
} lees_edwards_protocol_struct;

extern lees_edwards_protocol_struct lees_edwards_protocol;
/** Function to determine the current offset / velocity depending on the protocol */
void setup_lees_edwards_protocol();
/** Calculation of current offset*/
double lees_edwards_get_offset(double time);
/** Calculation of current velocity*/
double lees_edwards_get_velocity(double time);

#endif
