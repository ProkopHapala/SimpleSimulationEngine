#pragma once

#include "wpch.hpp"

/*
NOTE:
This is not directly AeroCraft related
I guess it should migrate somewhere else
Therefore I put it into single .h file and made it inline,
so that migration will be easier
*/


// complex multiplication  c = a*b 
inline void mulComplex(float ax, float ay, float bx, float by, float& cx, float& cy) { cx = ax * bx - ay * by;  cy = ax * by + ay * bx; }

// 2D rotation c = Rot(ang) * b 
inline void rot2D(float ang, float bx, float by, float& cx, float& cy)
{
  float ax = cos(ang); float ay = sin(ang);
  mulComplex(ax, ay, bx, by, cx, cy);
}

inline bool quadratic_roots(float a, float b, float c, float& x1, float& x2)
{
  float D = b * b - 4 * a*c;
  if (D < 0) 
    return false;
  float sqrtD = sqrt(D);
  float ia = -0.5 / a;
  if (a < 0) sqrtD *= -1; // make sure x1 < x2 
  x1 = (b + sqrtD)*ia;
  x2 = (b - sqrtD)*ia;
  //LogF( "  a,b,c, %f %f %f  x1,x2 %f %f \n", a,b,c, x1, x2 );
  return true;
}

// CubicSmoothStep function  https://en.wikipedia.org/wiki/Smoothstep
template <class TYPE>
inline TYPE cubicSmoothStep(TYPE x, TYPE x1, TYPE x2) {
  if (x < x1)
  {
    return 0.0;
  }
  else if (x > x2)
  {
    return 1.0;
  }
  else 
  { 
    TYPE a = (x - x1) / (x2 - x1); 
    return a * a*(3 - 2 * a); 
  }
}

// diagonalize 2x2 real symmetric matrix; rotation by complex number (ca,sa) =(cos(alfa),sin(alfa)) diagonalize the matrix 
void diagXY(float xx, float yy, float xy, float& ca, float& sa);  

// find maximum or minimum in array (sign=-1 minimum )
int findMax(int n, float * vals, float sign);                     

// =========================
//   RigidBody
//  ========================

class RigidBody 
{
public:
  // parameters
  float _invMass = 1.0;
  Matrix3	_invIbody = M3Identity;
  // State variables
  Vector3 _position = VZero;
  Vector3 _velocity = VZero;
  Vector3 _angMomentum = VZero;
  Matrix3 _orientation = M3Identity;
  Vector3 _angVelocity = VZero;

  //ClassIsGeneric(RigidBody); // https://docs.bisimulations.com/display/~Marek+Manukjan/Essence+library
  //ClassIsSimple(RigidBody);
  ClassIsBinary(RigidBody);

  // ==== function definitions
  void SetInertia(float mass, Vector3 diagI);
  void Move(double dt, Vector3 force, Vector3 torq);
  Vector3 GetControlTorq(float pitch, float roll, float yaw, float friction) const;

};

// =========================
//   RootFinder
//  ========================

// Iteratively search f(x) = 0
// 

class RootFinder
{
public:
  float   _x1, _x2;
  float   _y1;
  bool _bound = false;

  void  Init(float x1, float y1, float x2, float y2);
  float Bisect(float x, float y);  //  this is called when root is between bounds
  float Extent(float y);           //  this is called when root is not between bounds
  float Step(float x, float y);    //  returns next x 

};

// =========================
//   Optimizer1D_Golden
//  ========================

class Optimizer1D_Golden
{
public:
  static constexpr float invGoldenRatio = 0.61803398875;
  // see https://en.wikipedia.org/wiki/Golden-section_search
  float _yc, _yd;
  float _xc, _xd;
  float _xa, _xb;
  bool _lastC = false;
  bool _step0 = true;

  float init(float x1, float x2);
  float step(float y);

  inline float errX() const { return _xb - _xa; };

};