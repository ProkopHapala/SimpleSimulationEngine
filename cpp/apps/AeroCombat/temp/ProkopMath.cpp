
#include "ProkopMath.hpp"

int findMax(int n, float * vals, float sign)
{
  int imax = -1;
  float vmax = -INFINITY;
  for (int i = 0; i<n; i++)
  {
    float val = vals[i] * sign;
    if ((val > 0) || (val < 0))  // ignore NAN
      if (val > vmax)
      {
        vmax = val;
        imax = i;
      };
  }
  return imax;
}

void diagXY(float xx, float yy, float xy, float& ca, float& sa) 
{
  // see http://scipp.ucsc.edu/~haber/ph116A/diag2x2_11.pdf
  float u = xy + xy;
  float v = xx - yy;
  float c2a = v / sqrt(u*u + v * v);
  ca = sqrt((1 + c2a)*0.5);
  sa = sqrt((1 - c2a)*0.5);
  if (u > 0) { ca *= -1; }
  //LogF("diagXY xx,yy,xy %f,%f,%f u,v %f,%f c2a %f ca,sa %f,%f ", xx,yy,xy, u, v, c2a, ca, sa );
}

/*
void testMathSpeed(const Matrix3& rot, int n)
{
  Matrix3 rottot = M3Identity;
  Vector3 vec = VUp;
  long t1, t2;
  t1 = __rdtsc();
  for (int i = 0; i < n; i++)
    rottot = rot * rottot;
  t2 = __rdtsc();
  DIAG_MESSAGE_ID(100, 101, Format(" |M| %f, CPUticks/iter %f (%i) ", rottot.Determinant(), ((float)(t2 - t1)) / n, n));
  LogF(" |M| %f, CPUticks/iter %f  (%i) \n", rottot.Determinant(), ((float)(t2 - t1)) / n, n);

  t1 = __rdtsc();
  for (int i = 0; i < n; i++) { vec = rot * vec; }
  t2 = __rdtsc();
  DIAG_MESSAGE_ID(100, 102, Format(" |v| %f, CPUticks/iter %f (%i) ", vec.Size(), ((float)(t2 - t1)) / n, n));
  LogF(" |v| %f, CPUticks/iter %f  (%i) \n", vec.Size(), ((float)(t2 - t1)) / n, n);

  t1 = __rdtsc();
  float c = 0;
  for (int i = 0; i < n; i++) {
    c += rot.Direction().DotProduct(vec) + rot.DirectionUp().DotProduct(vec) + rot.DirectionAside().DotProduct(vec);
  }
  t2 = __rdtsc();
  DIAG_MESSAGE_ID(100, 103, Format(" c %f, CPUticks/iter %f (%i) ", c, ((float)(t2 - t1)) / n, n));
  LogF(" c %f, CPUticks/iter %f  (%i) \n", c, ((float)(t2 - t1)) / n, n);

  t1 = __rdtsc();
  for (int i = 0; i < n; i++) {
    vec += rot.Direction().CrossProduct(vec) + rot.DirectionUp().CrossProduct(vec) + rot.DirectionAside().CrossProduct(vec);
  }
  t2 = __rdtsc();
  DIAG_MESSAGE_ID(100, 104, Format(" v cross %f CPUticks/iter %f (%i) ", vec.Size(), ((float)(t2 - t1)) / n, n));
  LogF(" v cross %f, CPUticks/iter %f  (%i) \n", vec.Size(), ((float)(t2 - t1)) / n, n);
}
*/

// =========================
//   RigidBody
//  ========================

void RigidBody::SetInertia(float mass, Vector3 diagI)
{
  _invMass = 1 / mass;
  _invIbody = M3Identity;
  _invIbody.SetDirectionAside(Vector3(1 / diagI[0], 0.0, 0.0));
  _invIbody.SetDirectionUp(Vector3(0.0, 1 / diagI[1], 0.0));
  _invIbody.SetDirection(Vector3(0.0, 0.0, 1 / diagI[2]));
};

void RigidBody::Move(double dt, Vector3 force, Vector3 torq)
{
  _velocity += force * (dt * _invMass);
  _position += _velocity * dt;
  _angMomentum += torq * dt;
  Matrix3 invI = _orientation * _invIbody * _orientation.InverseRotation();
  _angVelocity = invI * _angMomentum;
  float r2Omega = _angVelocity.SquareSize();
  if (r2Omega>1e-8)
  {
    float rOmega = sqrt(r2Omega);
    Matrix3 drot = M3Identity;
    drot.SetRotationAxis(_angVelocity / rOmega, dt*rOmega);
    _orientation = drot * _orientation;
    _orientation.Orthogonalize(); // we probably do not need this each iteration
  }
  //LogF("force (%3.3f,%3.3f,%3.3f) vel (%3.3f,%3.3f,%3.3f) pos (%3.3f,%3.3f,%3.3f)\n", force.x,force.y,force.z, vel.x,vel.y,vel.z,  pos.x, pos.y, pos.z  );
  //LogF("L (%3.3f,%3.3f,%3.3f) omega (%3.3f,%3.3f,%3.3f) qrot (%3.3f,%3.3f,%3.3f,%3.3f)\n", L.x,L.y,L.z, omega.x,omega.y,omega.z,  qrot.x, qrot.y, qrot.z, qrot.w  );
};

Vector3 RigidBody::GetControlTorq(float pitch, float roll, float yaw, float friction) const
{
  Vector3 torq = VZero;
  float fa = _angVelocity.DotProduct(_orientation.DirectionAside());
  float fb = _angVelocity.DotProduct(_orientation.DirectionUp());
  float fc = _angVelocity.DotProduct(_orientation.Direction());
  torq += _orientation.Direction() * (roll - fc * friction);
  torq += _orientation.DirectionAside() * (pitch - fa * friction);
  torq += _orientation.DirectionUp() * (yaw - fb * friction);
  return torq;
};

// =========================
//   RootFinder
//  ========================

void RootFinder::Init(float x1, float y1, float x2, float y2)
{
  bool order;
  if ((y1*y2)<0) 
  {
    _bound = true;
    order = y1 > 0;
  }
  else
  {
    _bound = false;
    float dy = y2 - y1;
    order = (y1*dy) > 0;
  }
  if (order) {
    _x1 = x2;
    _x2 = x1;
    _y1 = y2;
  }
  else 
  {
    _x1 = x1;
    _x2 = x2;
    _y1 = y1;
  }
  //LogF("bound: %i x1: %f x2: %f \n", _bound, _x1, _x2);
}

float RootFinder::Bisect(float x, float y)
{
  if (y>0)
    _x2 = x;
  else
    _x1 = x;
  return 0.5*(_x1 + _x2);
}

float RootFinder::Extent(float y) {
  //float factor = 2.0;
  constexpr float factor = 1.61803398875;
  if (y*_y1>0) 
  {
    _x2 = _x1 + (_x2 - _x1)*factor;
    return _x2;
  }
  else 
  {
    _bound = true;
    if (_y1>0) {
      float tmp = _x1;
      _x1 = _x2;
      _x2 = tmp;
    };
    return 0.5*(_x1 + _x2);
  }
};

float RootFinder::Step(float x, float y) 
{
  if (_bound)
    return Bisect(x, y);
  else
    return Extent(y);
};

// =========================
//   Optimizer1D_Golden
//  ========================

float Optimizer1D_Golden::init(float x1, float x2)
{
  if (x1 > x2) { _xb = x1; _xa = x2; }
  else { _xa = x1; _xb = x2; }
  _xc = _xb - (_xb - _xa)*invGoldenRatio;
  return _xc;
};

float Optimizer1D_Golden::step(float y)
{
  if (_step0) { 
    _step0 = false; 
    _lastC = false; 
    _yc = y; 
    _xd = _xa + (_xb - _xa)*invGoldenRatio; 
    return _xd; 
  }
  if (_lastC) 
    _yc = y;
  else 
    _yd = y; ;
  if (_yc < _yd) 
  {
    _xb = _xd;
    _xd = _xc;
    _yd = _yc;
    _xc = _xb - (_xb - _xa)*invGoldenRatio;
    _lastC = true;
    return _xc;
  }
  else 
  {
    _xa = _xc;
    _xc = _xd;
    _yc = _yd;
    _xd = _xa + (_xb - _xa)*invGoldenRatio;
    _lastC = false;
    return _xd;
  }
};