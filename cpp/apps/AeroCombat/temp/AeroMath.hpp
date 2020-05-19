#pragma once

#include "wpch.hpp"
#include "ProkopMath.hpp"

#define G_CONST_ 9.8066f   // G_CONST is defined in vehicle.hpp ... I don't want make AeroMath.hpp dependent on vehicle.hpp

// Following functions a assuming Linear regime !!!; Drag independet on velocity; no supersonic, no Reynols number change )
// aerodynamic force nicely scales with V^2
inline float dragThrustToV2(float Drag, float Thrust)     { return Thrust / Drag;     };
inline float liftMassToTrunRadius(float Lift, float mass) { return mass / Lift;       };
inline float turnRadiusToTime(float R, float v)           { return 6.28318530718*R/v; };

// turn radius in flat turn with gravity
void evalThrustLimitedTurnGrav(float Lift, float Drag, float Thrust, float mass, float& v2, float& R, float& sinTheta);

void evalClimbRate(float Lift, float Drag, float Thrust, float mass, float& v, float& ca, float& sa);

inline float supersonicDrag(float m, float M0, float Mmin, float Mmax) 
{
  // Drag cureve should behave like this
  // http://adg.stanford.edu/aa241/drag/volumedrag.html
  // http://www.desktop.aero/appliedaero/compress3d/ssdragest.html
  // http://adg.stanford.edu/aa241/drag/SSDragCalc.html
  // http://aerorocket.com/HTV-3X.html
  // http://what-when-how.com/space-science-and-technology/rocket-propulsion-theory/
  //
  //  Mmin ... begining of supersonic transition
  //  Mmax ... end of      supersonic transition
  //  M0   ... position of wave drag divergence
  if (m < Mmin) 
    return 0.0;
  return cubicSmoothStep(m, Mmin, Mmax) / sqrt(m*m - M0 * M0);
}

namespace Atmosphere 
{
  // Reference data:
  // http://farside.ph.utexas.edu/teaching/sm1/lectures/node56.html
  // http://www.engineeringtoolbox.com/standard-atmosphere-d_604.html

  const int    _n = 47;
  const float  _dh = 2000.0; // [m]
  const float  _dh_inv = 1.0f / _dh;  // [1/m]
  const float  _density[_n] = {
    1.45,   // -2000.0 m
    1.1717, // 0.0 m
    0.97103,
    0.79314,
    0.64266,
    0.5214,
    0.41827,
    0.33109,
    0.25823,
    0.19214,
    0.13521,
    0.095807,
    0.06838,
    0.049144,
    0.035556,
    0.025891,
    0.018969,
    0.013981,
    0.010363,
    0.0077248,
    0.0057779,
    0.0043456,
    0.0032867,
    0.0025144,
    0.0019321,
    0.0014931,
    0.0011632,
    0.00091197,
    0.00071504,
    0.00056185,
    0.00044591,
    0.00035152,
    0.00027514,
    0.00021373,
    0.00016469,
    0.00012581,
    9.5234E-05,
    7.0578E-05,
    5.1557E-05,
    3.7497E-05,
    2.7148E-05,
    1.9564E-05,
    1.4031E-05,
    1.0013E-05,
    7.11E-06,
    5.0219E-06,
    3.5279E-06, // 90000 m
  };

  // ==== inline functions

  inline float GetDensity(float h) 
  {
    float u = h * _dh_inv;
    int   i = (int)u;
    u -= i;
    if ((i + 2) >= _n) 
      return 0;
    return (1 - u)*_density[i + 1] + u * _density[i + 2]; // TODO: cubic spline ?
    // return plante_rho0 * exp( planet_zrate * h );      // alternatively ?                                    
  };

};

