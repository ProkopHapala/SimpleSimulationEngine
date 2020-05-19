
#include "AeroMath.hpp"

void evalThrustLimitedTurnGrav(float Lift, float Drag, float Thrust, float mass, float& v2, float& R, float& sinTheta)
{
  v2 = dragThrustToV2(Drag, Thrust);
  float Fg = G_CONST_ * mass;
  float v2Lift = v2 * Lift;
  float cosTheta = 2;
  if (v2Lift != 0) cosTheta = Fg / (v2Lift);
  //LogF( "evalThrustLimitedTurnGrav: Lift %f Drag %f Thrust %f Fg %f v2 %f cosTheta %f ", Lift, Drag, Thrust, Fg, v2, cosTheta );
  if (cosTheta > 1) { sinTheta = NAN; R = NAN; return; }
  sinTheta = sqrt(1 - cosTheta * cosTheta);
  R = liftMassToTrunRadius(Lift*sinTheta, mass);
}

void evalClimbRate(float Lift, float Drag, float Thrust, float mass, float& v, float& ca, float& sa)
{
  // Derivation:
  //  Force balance (Pyctagorean triangle)
  //  (Thrust-FDrag)^2 + FLift^2 = Gravity^2  
  //  Aerodynamic force depend on velocity^2 (v2)
  //  (Thrust-Drag*v2)^2 + (Lift*v2)^2 = (mass*G_CONST)^2  
  //  (Thrust-Drag*v2)^2 + (Lift*v2)^2 = (mass*G_CONST)**2
  //  Thrust**2 - 2 * Drag*Thrust*v2 + (Drag*v2)**2 + (Lift*v2)**2 = (mass*G_CONST)**2
  // Quadratic equation to solve for v2 (velocity^2) 
  //  ((Drag + Lift)**2)*v2**2 - (2 * Drag*Ft)*v2 + (Thrust**2 - (mass*G_CONST)**2) = 0 
  float Fg = G_CONST_ * mass;
  float a = Drag * Drag + Lift * Lift;
  float b = -2 * Drag*Thrust;
  float c = Thrust * Thrust - Fg * Fg;
  float vv1, vv2;
  quadratic_roots(a, b, c, vv1, vv2);
  //LogF("evalClimbRate: Lift %f Drag %f Thrust %f Fg %f a,b,c(%f,%f,%f) vv2 %f ", Lift, Drag, Thrust, Fg,   a,b,c,  vv2 );
  float Fd = Drag * vv2;
  float Fl = Lift * vv2;
  sa = (Thrust - Fd) / Fg;
  ca = Fl / Fg;
  v = sqrt(vv2);
  //vx = ca * v;
  //vy = sa * v;
}