
//#include "rString.hpp"

#include <math.h>
#include "AeroSurf.hpp"

#if _DEBUG_AEROSURF
void drawWing(Vector3 pos, Matrix3 rot, float a, float AspectRatio )
{
  Matrix4 pose;
  pose.SetPosition(Vector3D(pos.X(), pos.Y(), pos.Z()));
  pose.SetOrientation(rot);
  renderer.DrawBox(Vector3D(0, 0, 0), Vector3(AspectRatio*a*0.5, 0.025, a*0.5), PackedWhite, &pose);
}
void drawRigidBody(Vector3 pos, Matrix3 rot, float sz) 
{
  renderer.Add3dLine(pos, pos + rot.Direction()       * sz, PackedBlue  );
  renderer.Add3dLine(pos, pos + rot.DirectionUp()     * sz, PackedGreen );
  renderer.Add3dLine(pos, pos + rot.DirectionAside()  * sz, PackedRed   );

  Matrix4 pose;
  pose.SetPosition(Vector3D(pos.X(), pos.Y(), pos.Z()));
  pose.SetOrientation(rot);
  renderer.DrawBox(Vector3D(0, 0, 0), Vector3(0.5, 0.25, 1.0), PackedWhite, &pose);
};
#endif // _DEBUG_AEROSURF

//==============
//  PolarModel
//==============

inline void PolarModel::GetAeroCoefs(float ca, float sa, float& CD, float& CL) const {
  float abs_sa = (sa>0) ? sa : -sa;                                       // symmetric wing  
  float wS = cubicSmoothStep<float>(abs_sa, _sStall, _sStall + _wStall);  // stalled mixing factor
  float mS = 1 - wS;                                                      // not-stalled mixing factor
  CD = _CD0 + ( mS * _dCD * abs_sa + wS * _dCDS) * abs_sa;                // drag
  if (ca <0) {   // for Angle of Attack > 90 deg ( air comes from behind )
    ca = -ca;
    sa = -sa;
  };
  CL = ( mS* _dCL + wS * _dCLS*ca) * sa;   // lift
}

//=====================
//  PolarModelTorque
//=====================

inline void  PolarModelTorque::GetAeroCoefs(float ca, float sa, float flap, float& CD, float& CL, float& CM) const
{
  float abs_sa = (sa>0) ? sa : -sa;                                      // symmetric wing  
  float wS = cubicSmoothStep<float>(abs_sa, _sStall, _sStall + _wStall); // stalled mixing factor
  float mS = 1 - wS;                                                     // not-stalled mixing factor
  CD = _CD0 + (mS*_dCD*abs_sa + wS * _dCDS) * abs_sa;                    // drag
  if (ca <0) {     // for Angle of Attack > 90 deg ( air comes from behind )
    ca = -ca;
    sa = -sa;
  };
  CL = (mS*_dCL + wS * _dCLS*ca) * sa;    // lift
  float sa_, ca_;
  rot2D(flap, ca, sa, ca_, sa_);
  abs_sa = (sa_>0) ? sa_ : -sa_;
  //// We should use the same stall cutoff as main wing independent of flap (?)
  //wS = cubicSmoothStep<float>(abs_sa, _sStall, _sStall + _wStall); // stalled mixing factor
  //mS = 1 - wS;                                                     // not-stalled mixing factor
  CM = (mS*_dCM + wS * _dCMS*ca_) * sa_;                             // TODO: FIXME pitching moment 
}

// ======================================================
//             AeroSurface
// ======================================================

float AeroSurface::Control2tilt(float control) const 
{
  // if _maxTilt <> -_minTilt   this will produce piecewise linear function f(0)=0 with different slope for positive and negative value of control  
  if (control > 0) { return control * _maxTilt; }
  else { return control * -_minTilt; }
}

float AeroSurface::Controls2tilt(float elevator, float rudder, float aileron) const
{
  float control = _cElevator*elevator + _cRudder*rudder + _cAileron*aileron;
  if (control > 0) { return control * _maxTilt; }
  else { return control * -_minTilt; }
}

void AeroSurface::Tilt(float angle) 
{
  Matrix3 tiltMat = M3Identity;
  tiltMat.SetRotationX(angle);
  _lrot = _lrot * tiltMat;
};

void AeroSurface::ApplyForce(const Vector3& vair0, const Matrix3& craftRot, const Vector3& craftOmega, Vector3& force, Vector3& torq, float tilt, int i, bool debug, const Vector3& pos0 ) const
{
  // global coordinates
  Matrix3 grot; // rotation of panel in world coordinates
  if (_bControlable) { // Apply controls
    Matrix3 tiltMat = M3Identity;
    tiltMat.SetRotationX(tilt);
    grot = craftRot * _lrot * tiltMat;
  }
  else 
  {
    grot = craftRot * _lrot;
  }
  Vector3 gdpos = craftRot * _lpos;   // position of panel in world coordinates relative to COG

  Vector3 uair = vair0 + gdpos.CrossProduct(craftOmega); // velocity of panel due to rotation of aircraft
  //Vector3 uair = vair0;  // ignore rotation of aircraft

  float vrair2 = uair.SquareSize();
  if (vrair2 > lowSpeedCuoff) // for zero air-speed it would diverge
  {
    // unitary vector in direction of air flow
    float vrair = sqrt(vrair2);
    uair *= (1 / vrair);

    // decompose uair to panel coordinates
    float ca = grot.DirectionAside().DotProduct(uair);
    float cb = grot.DirectionUp().DotProduct(uair);
    float cc = grot.Direction().DotProduct(uair);

    // All aerodynamic forces are scalled by this factor
    float prefactor = vrair2 * _area;   // density not included

    force = VZero;

#if _DEBUG_AEROSURF
    Vector3 pos = pos0 + gdpos;       // position of panel in world coordinates; just for ploting
    if (debug)
      renderer.Add3dLine(pos, pos + uair * 5.0, PackedBlue);
#endif // _DEBUG_AEROSURF

      // get Lift and Drag coefs from polar (dimensionless)
      float CD, CL;
      _polar.GetAeroCoefs(-cc, cb, CD, CL);
      if (_bFiniteWingCorrection)
      {
        // NOTE: This is used if we use polar of infinite wing (of 2D airfoil )
        //       but finite wing effect can be already included in polar, than _bFiniteWingCorrection should be set to false
        // http://www.srmuniv.ac.in/sites/default/files/downloads/class4-2012.pdf
        // https://en.wikipedia.org/wiki/Lifting-line_theory#Useful_approximations
        CL *= _aspectRatio / ( 2 + _aspectRatio );
        CD += CL * CL / (3.1415926535 * _aspectRatio); // induced drag // https://en.wikipedia.org/wiki/Lift-induced_drag
      }
      // Scale Lift and Drag (no density yet)
      CL *= prefactor; 
      CD *= prefactor;
      if ((cb*cb) > 1e-8) // for zero Angle of Attack we cannot determine lift direction 
      {
        Vector3 airUp = grot.DirectionUp() + (uair * -cb); // component of grot.Up perpendicular to uair
        airUp.Normalize();
        force = airUp * CL + uair * CD;                    // combine Lift and Drag force

#if _DEBUG_AEROSURF
        if (debug) 
        {
          renderer.Add3dLine(pos, pos + airUp * 5.0, PackedGreen);
          char plotName[128];
          sprintf(plotName, "wing[%i]", i);
          //TO_DBG_MGRAPH(plotName, Lift, CL, 200);
          //TO_DBG_MGRAPH(plotName, Drag, CD, 200);
          /*
          if (debugPlot>0) {
            int N = 200;
            if (!liftPlot && GDebugGui) {
              liftPlot = GDebugGui->AddMultiGraph(plotName, "Lift", N);
              GEngine->ShowDiagGraph(liftPlot);
              iLift = liftPlot->AddGraph("Lift");
            }
            if (!dragPlot && GDebugGui) {
              dragPlot = GDebugGui->AddMultiGraph(plotName, "Drag", N);
              GEngine->ShowDiagGraph(dragPlot);
              iDrag = dragPlot->AddGraph("Drag");
            }
            if (debugPlot == 1) {
              if (-1 != iLift && liftPlot)  liftPlot->AddValue(iLift, CL);
              if (-1 != iDrag && dragPlot)  dragPlot->AddValue(iDrag, CD);
            }
            if (debugPlot == 2) {
              liftPlot->Clear();
              dragPlot->Clear();
              float dAoA = 6.28 / N;
              for (int i = 0; i < N; i++ ) {
                float AoA = dAoA * i - 3.14;
                float ca_ = cos(AoA);
                float sa_ = sin(AoA);
                polar->getAeroCoefs(ca_, sa_, CD, CL);
                liftPlot->AddValue(iLift, CL);
                dragPlot->AddValue(iDrag, CD);
              }
              float AoA = atan2( cb, -cc );
              liftPlot->HighlightB = AoA/6.28 + 0.5;
              DIAG_MESSAGE_ID(100, 30 + i, Format("AoA %f", AoA ) );
            }
          }
          */
          // DECL_DBG_2DSLIDER(coords, 1, 2, -5, 5, -5, 5); // Creates two floats, coord_X and coord_Y, with their initial values set to [1, 2] and the slider has a range of <-5,5> x <-5,5>
          //TO_DBG_MGRAPH(speed, z, FutureVisualState().RelativeSpeed().Z(), 200);
        }
#endif // _DEBUG_AEROSURF
      }
      else
      {
        force += uair * CD;   // for zero Angle of Attack only drag force
      }
    // apply anisotropic drag coef - Do we really need this here ?
    force += grot.DirectionAside() * (_cAnisoDrag[0] * ca * prefactor);
    force += grot.DirectionUp()    * (_cAnisoDrag[1] * cb * prefactor);
    force += grot.Direction()      * (_cAnisoDrag[2] * cc * prefactor);
    torq   = gdpos.CrossProduct(force);
#if _DEBUG_AEROSURF
    if (debug) {
      renderer.Add3dLine(pos, pos + grot.DirectionUp() * 5.0, PackedBlack);
      renderer.Add3dLine(pos, pos + grot.Direction()   * 5.0, PackedWhite);
      renderer.Add3dLine(pos, pos + force * 5.0, PackedRed);
    }
#endif // _DEBUG_AEROSURF
  }
};

void AeroSurface::FromString(const char * str) 
{
  float lrot_bx, lrot_by, lrot_bz,
        lrot_cx, lrot_cy, lrot_cz;
  sscanf(str, " %f %f %f    %f %f %f    %f %f %f    %f %f %f  %f %f",
    &_lpos[0], &_lpos[1], &_lpos[2],
    &lrot_bx, &lrot_by, &lrot_bz,
    &lrot_cx, &lrot_cy, &lrot_cz,
    &_cAnisoDrag[0], &_cAnisoDrag[1], &_cAnisoDrag[2],
    &_area, &_aspectRatio
  );
  _lrot.SetDirectionAndUp(Vector3(lrot_cx, lrot_cy, lrot_cz), Vector3(lrot_bx, lrot_by, lrot_bz));
};

void AeroSurface::FromStringPolarModel(const char * str) 
{
  float lrot_bx, lrot_by, lrot_bz,
        lrot_cx, lrot_cy, lrot_cz;
  sscanf(str, " %f %f %f    %f %f %f    %f %f %f    %f %f %f     %f %f   %f %f %f   %f %f   %f %f",
    &_lpos[0], &_lpos[1], &_lpos[2],
    &lrot_bx, &lrot_by, &lrot_bz,
    &lrot_cx, &lrot_cy, &lrot_cz,
    &_cAnisoDrag[0], &_cAnisoDrag[1], &_cAnisoDrag[2],
    &_area, &_aspectRatio,
    &(_polar._CD0), &(_polar._dCD), &(_polar._dCDS), &(_polar._dCL), &(_polar._dCLS), &(_polar._sStall), &(_polar._wStall)
  );
  _lrot.SetDirectionAndUp( Vector3(lrot_cx, lrot_cy, lrot_cz), Vector3(lrot_bx, lrot_by, lrot_bz) );
};

int AeroSurface::ToStringPolarModel(char * str) const {
  return sprintf(str, "%f %f %f    %f %f %f    %f %f %f    %f %f %f     %f %f   %f %f %f   %f %f   %f %f      %f %f   %i  \n",
    _lpos[0], _lpos[1], _lpos[2],
    _lrot.DirectionUp()[0], _lrot.DirectionUp()[1], _lrot.DirectionUp()[2],
    _lrot.Direction()[0], _lrot.Direction()[1], _lrot.Direction()[2],
    _cAnisoDrag[0], _cAnisoDrag[1], _cAnisoDrag[2],
    _area,_aspectRatio,
    _polar._CD0, _polar._dCD, _polar._dCDS, _polar._dCL, _polar._dCLS, _polar._sStall, _polar._wStall,
    _minTilt, _maxTilt, _bFiniteWingCorrection
  );
};

// =========================================
//             Propeler
// =========================================

float Propeler::GetBaypassThrust(float v0, float throtle) const
{
  // Thrust is limited by amount of air passing through propeler
  // Derived from those equations
  //   dm = S * v0
  //   F = dm * Dv
  //   Engine power accelerates both the aircraft and the air
  //   P = F * v0 + 0.5*dm*(Dv**2) = S*(v0**2)*Dv + 0.5*S*v0*(Dv**2)    // Quadratic equation
  //   v0      ... aircraft velocity
  //   vstatic ... we need non-zero flow through properler even when aircraft stand still on ground
  float dm = _area * (v0 + _vstatic); // mass of air-flow per second
  // coefs of quadratic equation
  float a = 0.5*dm;
  float b = dm * v0;
  float c = -_power * throtle;
  float Dv1, Dv2;
  quadratic_roots(a, b, c, Dv1, Dv2); // solve of square of air-speed
  return dm * Dv1 * _efficiency - dm * v0*_CD;
}

void Propeler::ApplyForce(const Vector3& vair0, float airDens, const Matrix3& craftRot, const Vector3& craftOmega, Vector3& force, Vector3& torq, float throtle, int i, bool debug, const Vector3& pos0) const
{
  // rotation of panel in world coordinates
  Vector3 gdir  = craftRot * _ldir;   // direction
  Vector3 gdpos = craftRot * _lpos;   // position
  Vector3 pos = pos0 + gdpos;  
  float v = 0.0;
  float Fout = _thrustRocket * throtle;                    // constant velocity independent thrust - like rocket
  if (_bSolveBypass) 
  {                                                        // thrust computed from velocity and engine power (ideal propeler)
    Vector3 vair = vair0 + gdpos.CrossProduct(craftOmega); // velocity of panel due to rotation of aircraft
    v = vair.Size() * gdir.DotProduct(vair);
    Fout += GetBaypassThrust(v,throtle)*airDens;
  }
  force = gdir * Fout;
  torq  = gdpos.CrossProduct(force);
};

float Propeler::AutoVStatic() const 
{
  return pow(4 * _power / _area, 0.333333); // this is just guess without good physical justification
};

void Propeler::FromString(const char * str)
{
  sscanf(str, " %f %f %f    %f %f %f    %f   %f %f %f %f  %f",
    &_lpos[0], &_lpos[1], &_lpos[2],
    &_ldir[0], &_ldir[1], &_ldir[2],
    &_thrustRocket, 
    &_area, &_power, &_efficiency, &_CD, &_vstatic
  );
  _ldir.Normalize();
  if (_vstatic < 0) _vstatic=AutoVStatic();
  printf("%lf %lf %lf %lf\n", _area, _power, _efficiency, _CD);
}

// ================================
//  ComposedAeroModel
// ================================

void ComposedAeroModel::ApplyForce(const Vector3& vair, const Matrix3& craftRot, const Vector3& craftOmega, Vector3& force, Vector3& torq, float elevator, float rudder, float aileron, float throtle, const Vector3& pos0) const 
{
  force = VZero;
  torq = VZero;
#if _DEBUG_AEROSURF
  bool debug = _DebugLevel >= 2;
  if(debug)
  {
    DIAG_MESSAGE_ID(100, 9, Format("vair = %f %f %f", vair.X(), vair.Y(), vair.Z()));
    renderer.Add3dLine(pos0, pos0 + vair * 10.0, PackedRed);
  }
#endif // _DEBUG_AEROSURF
  for (int iwing = 0; iwing < _wings.Size(); iwing++)
  {
    Vector3 torq_i, force_i;
    const AeroSurface& wing = _wings[iwing];

    // this is replaced by  AeroSurf::_cElevator,_cRudder,_cAileron
    //float tilt = 0.0;
    //if (iwing == _elevatorId)
    //  tilt = -pitch;
    //else if (iwing == _rudderId)
    //  tilt = -yaw;
    //else if (iwing == _leftAirelonId)
    //  tilt = roll;
    //else if (iwing == _rightAirelonId)
    //  tilt = -roll;
    //wing.ApplyForce(vair, craftRot, craftOmega, force_i, torq_i, wing.Control2tilt(tilt), iwing, debug, pos0);

    float tilt = 0.0;
    if(wing._bControlable)
      tilt = wing.Controls2tilt(elevator, rudder, aileron);
    wing.ApplyForce(vair, craftRot, craftOmega, force_i, torq_i, tilt, iwing, debug, pos0);

    force += force_i;
    torq += torq_i;
#if _DEBUG_AEROSURF
    if (debug)
    {
      Vector3 pos_i = pos0 + craftRot * wing._lpos;
      float a = sqrt(wing._area / wing._aspectRatio);
      drawWing(pos_i, craftRot*wing._lrot, a, wing._aspectRatio);
    }
#endif // _DEBUG_AEROSURF
  }

  if (_bSupersonic) 
  {
    // http://www.srmuniv.ac.in/sites/default/files/downloads/class4-2012.pdf
    float speedOfSound = 343.0f;  // TODO : speed of sound may change with altitude etc.  https://en.wikipedia.org/wiki/Speed_of_sound#/media/File:Comparison_US_standard_atmosphere_1962.svg
    float v2 = vair.SquareSize();
    if ( v2 > Square( speedOfSound*_MachMin ) ) 
    {
      float v = sqrt(v2);
      float wd = supersonicDrag( v/speedOfSound, _Mach0, _MachMin, _MachMax );
      float cdir = vair.DotProduct(force);
      force +=  vair * ( _cWaveDrag * wd * cdir / v2 ); // apply wave-drag-force along airflow direction;  1/v2 factor since (vair*cdir)~v2
    }
  }

  float fDens = 0.5*Atmosphere::GetDensity(pos0[1]); // rho/2;    hope that _position[1]==0 is sea level ?
#if _DEBUG_AEROSURF
  if (_DebugLevel >= 2) 
  {
    DIAG_MESSAGE_ID(15, 65, Format("Attitude %f AirDensity %f ", pos0[1], fDens*2.0));
  }
#endif // _DEBUG_AEROSURF
  force *= fDens;
  torq *= fDens;

  for (int i = 0; i < _propelers.Size(); i++) 
  {
    Vector3 torq_i, force_i;
    const Propeler& prop = _propelers[i];
    prop.ApplyForce(vair, fDens, craftRot, craftOmega, force_i, torq_i, throtle, i, debug, pos0);
    force += force_i;
    torq  += torq_i;
  }

};

void ComposedAeroModel::FromFile(FILE* pFile) 
{
    const int nbuf = 1024;
    char buf[nbuf];
    fgets(buf, nbuf, pFile);
    int nWings;
    sscanf(buf, "%i\n", &nWings);
    //_wings.Access(nWings);
    _wings.Resize(nWings);
    for (int i = 0; i < _wings.Size(); i++) 
    {
      fgets(buf, nbuf, pFile);
      _wings[i].FromStringPolarModel(buf);
    }
    //fgets(buf, nbuf, pFile);
    //sscanf(buf, "%i %i %i %i\n", &_leftAirelonId, &_rightAirelonId, &_elevatorId, &_rudderId);
    //_leftAirelonId--;
    //_rightAirelonId--;
    //_elevatorId--;
    //_rudderId--;
};

int ComposedAeroModel::LoadFile(const char * fname) 
{
  FILE * pFile;
  LogF(" AeroTestPlatform::InitWingsFile loading wings from: >>%s<<\n", fname);
  pFile = fopen(fname, "r");
  if (pFile == nullptr) 
  {
    LogF(" ComposedAeroModel::fromFile cannot open >>%s<< => call default InitWings() \n", fname);
    //exit(-1);
    return -1;
  }
  FromFile(pFile);
  LogF("AeroTestPlatform::InitWingsFile wings loaded from >>%s<<\n", fname);
  fclose(pFile);
  return 0;
};

char* ComposedAeroModel::ToString(char * str) const 
{
  str += sprintf(str, " ===== AeroTestPlatform::toString : _nWings %i \n", _wings.Size());
  for (int i = 0; i < _wings.Size(); i++) 
  {
    str += sprintf(str, " ===== wing[%i] ", i);
    str += _wings[i].ToStringPolarModel(str);
  };
  return str;
};

float ComposedAeroModel::WettedArea() const 
{
  float area = 0;
  for (int i = 0; i < _wings.Size(); i++) 
  {
    area += _wings[i]._area;
  }
  return area;
};

float ComposedAeroModel::MultArea(float f) 
{
  float area = 0;
  for (int i = 0; i < _wings.Size(); i++) 
  {
    float a = _wings[i]._area;
    a *= f;
    _wings[i]._area = a;
    area += a;
  }
  return area;
};

float ComposedAeroModel::GetMaxTorqScale() const 
{
  float tmax = 0;
  for (int i = 0; i < _wings.Size(); i++) 
  {
    float ti = _wings[i]._lpos.Size() * _wings[i]._area;
    tmax = fmax(tmax, ti);
  }
  return tmax;
};

// ======================================================
//              CompactAeroModel
// ======================================================

void CompactAeroModel::FromString(const char * str)
{
  sscanf(str, " %f   %f %f %f   %f %f %f   %f %f %f   %f %f %f    %f %f %f %f %f %f %f ",
    &_area,
    &_cAnisoDrag[0], &_cAnisoDrag[1], &_cAnisoDrag[2],
    &_cTorqDamp[0], &_cTorqDamp[1], &_cTorqDamp[2],
    &_cTorqControl[0], &_cTorqControl[1], &_cTorqControl[2],
    &_cTorqRelax[0], &_cTorqRelax[1], &_cTorqRelax[2],
    &(_polarFwUp._sStall), &(_polarFwUp._wStall), &(_polarFwUp._CD0), &(_polarFwUp._dCD), &(_polarFwUp._dCDS), &(_polarFwUp._dCL), &(_polarFwUp._dCLS)
  );
  LogF(" ===== CompactAeroModel::fromString \n");
  LogF(" area %f C %f %f %f D %f %f %f L %f %f %f S %f %f %f \n", _area,
    _cAnisoDrag[0], _cAnisoDrag[1], _cAnisoDrag[2],
    _cTorqDamp[0], _cTorqDamp[1], _cTorqDamp[2],
    _cTorqControl[0], _cTorqControl[1], _cTorqControl[2],
    _cTorqRelax[0], _cTorqRelax[1], _cTorqRelax[2]
  );
  LogF("polarFwUp   %f %f %f %f %f %f %f \n", _polarFwUp._sStall, _polarFwUp._wStall, _polarFwUp._CD0, _polarFwUp._dCD, _polarFwUp._dCDS, _polarFwUp._dCL, _polarFwUp._dCLS);
};

void CompactAeroModel::ApplyForce(const Vector3& vair0, const Matrix3& craftRot, const Vector3& craftOmega, Vector3& force, Vector3& torq, float elevator, float rudder, float aileron, const Vector3& pos0) const
{
  Vector3 uair = vair0;    // we shoud consider craftOmega
#if _DEBUG_AEROSURF
  renderer.Add3dLine(pos0, pos0 + uair * 5.0, PackedBlue);
  renderer.Add3dLine(pos0, pos0 + craftRot.DirectionAside()*5.0, PackedColor(0xff808080));
  renderer.Add3dLine(pos0, pos0 + craftRot.DirectionUp()*5.0, PackedBlack);
  renderer.Add3dLine(pos0, pos0 + craftRot.Direction()*5.0, PackedWhite);
#endif // _DEBUG_AEROSURF
  //Vector3 uair = vair0 + gdpos.CrossProduct(craftOmega);                                    
  float vrair2 = uair.SquareSize();
  torq = VZero;
  force = VZero;
  if (vrair2 > lowSpeedCuoff)   // for zero air-speed we would devide by 0 
  {
    float vrair = sqrt(vrair2);
    uair *= (1 / vrair);
    // decompose uair to panel coordinates
    float ca = uair.DotProduct(craftRot.DirectionAside());
    float cb = uair.DotProduct(craftRot.DirectionUp());
    float cc = uair.DotProduct(craftRot.Direction());
#if _DEBUG_AEROSURF
    DIAG_MESSAGE_ID(100, 11, Format("elevator %f rudder %f aileron %f ", elevator, rudder, aileron));
    DIAG_MESSAGE_ID(100, 12, Format("ca %f cb %f cc %f ", ca, cb, cc));
#endif // _DEBUG_AEROSURF
    float prefactor = vrair2 * _area; // scale all aerodynamic forces by this
    {
      float CD, CL;
      Vector3 airUp;
      _polarFwUp.GetAeroCoefs(-cc, cb, CD, CL);
      if ((cb*cb) > 1e-8)
      {
        airUp = craftRot.DirectionUp() + (uair * -cb);  // component of grot.Up perpendicular to uair
        airUp.Normalize();
      }
      else { airUp = craftRot.DirectionUp(); }
#if _DEBUG_AEROSURF
      renderer.Add3dLine(pos0, pos0 + airUp * 5.0, PackedGreen);
#endif // _DEBUG_AEROSURF
      force += (airUp * CL + uair * CD) * prefactor;
    }
#if _DEBUG_AEROSURF
    renderer.Add3dLine(pos0, pos0 + force * 5.0, PackedRed);
#endif // _DEBUG_AEROSURF

    float cc2 = cc * cc; // damp effect of control surfaces when aircraft is not alligned to airflow

    // Stabilization - aligns aircraft orientation toward airflow
    torq += craftRot.DirectionAside() * (_cTorqRelax[0] * cb*cc2* prefactor);
    torq += craftRot.DirectionUp()    * (_cTorqRelax[1] * -ca * cc2* prefactor);
    torq += uair.CrossProduct(craftRot.Direction())     * ( _cTorqRelax[2] * prefactor );

    // Control input
    torq += craftRot.DirectionAside() * (_cTorqControl[0] * elevator * prefactor * cc2);
    //torq += craftRot.DirectionUp()  * (_cTorqControl[1] * -rudder  * prefactor * cc2 );
    torq += craftRot.Direction()      * (_cTorqControl[2] * aileron  * prefactor * cc2);

    // damping of aircraft rotation due to drag of control surfaces
    torq += craftRot.DirectionAside() * ( craftOmega.DotProduct(craftRot.DirectionAside()) * -_cTorqDamp[0] * prefactor); // pitch drag
    torq += craftRot.DirectionUp()    * ( craftOmega.DotProduct(craftRot.DirectionUp())    * -_cTorqDamp[1] * prefactor); // yaw drag
    torq += craftRot.Direction()      * ( craftOmega.DotProduct(craftRot.Direction())      * -_cTorqDamp[2] * prefactor); // roll drag

    // apply anisotropic drag coef
    force += craftRot.DirectionAside() * (_cAnisoDrag[0] * ca*prefactor);
    force += craftRot.DirectionUp()    * (_cAnisoDrag[1] * cb*prefactor);
    force += craftRot.Direction()      * (_cAnisoDrag[2] * cc*prefactor);
#if _DEBUG_AEROSURF
    DIAG_MESSAGE_ID(100, 14, Format("force %f %f %f torq %f %f %f ", force.X(), force.Y(), force.Z(), torq.X(), torq.Y(), torq.Z()));
#endif // _DEBUG_AEROSURF
  }
  float fDens = 0.5*Atmosphere::GetDensity(pos0[1]); // hope that _position[1]==0 is ground level ?
  force *= fDens;
  torq *= fDens;
};

// ======================================================
//              CompactAeroMode2
// ======================================================

void CompactAeroModel2::FromString(const char * str)
{
  sscanf(str, " %f   %f %f %f   %f %f %f   %f %f %f   %f %f %f      %f %f %f %f %f %f %f %f %f      %f %f %f %f %f %f %f %f %f ",
    //&areaUp, &ARup, 
    //&areaSide, &ARside,
    &_area,
    &_cAnisoDrag[0], &_cAnisoDrag[1], &_cAnisoDrag[2],
    &_cTorqDamp[0], &_cTorqDamp[1], &_cTorqDamp[2],
    &_cTorqControl[0], &_cTorqControl[1], &_cTorqControl[2],
    &_cTorqRelax[0], &_cTorqRelax[1], &_cTorqRelax[2],
    &(_polarFwUp._sStall), &(_polarFwUp._wStall), &(_polarFwUp._CD0), &(_polarFwUp._dCD), &(_polarFwUp._dCDS), &(_polarFwUp._dCL), &(_polarFwUp._dCLS), &(_polarFwUp._dCM), &(_polarFwUp._dCMS),
    &(_polarFwSide._sStall), &(_polarFwSide._wStall), &(_polarFwSide._CD0), &(_polarFwSide._dCD), &(_polarFwSide._dCDS), &(_polarFwSide._dCL), &(_polarFwSide._dCLS), &(_polarFwSide._dCM), &(_polarFwSide._dCMS)
  );
  LogF(" ===== CompactAeroModel2::fromString \n");
  LogF(" area %f C %f %f %f D %f %f %f L %f %f %f S %f %f %f \n", _area,
    _cAnisoDrag[0], _cAnisoDrag[1], _cAnisoDrag[2],
    _cTorqDamp[0], _cTorqDamp[1], _cTorqDamp[2],
    _cTorqControl[0], _cTorqControl[1], _cTorqControl[2],
    _cTorqRelax[0], _cTorqRelax[1], _cTorqRelax[2]
  );
  LogF("polarFwUp   %f %f %f %f %f %f %f %f %f \n", _polarFwUp._sStall, _polarFwUp._wStall, _polarFwUp._CD0, _polarFwUp._dCD, _polarFwUp._dCDS, _polarFwUp._dCL, _polarFwUp._dCLS, _polarFwUp._dCM, _polarFwUp._dCM);
  LogF("polarFwSide %f %f %f %f %f %f %f %f %f \n", _polarFwSide._sStall, _polarFwSide._wStall, _polarFwSide._CD0, _polarFwSide._dCD, _polarFwSide._dCDS, _polarFwSide._dCL, _polarFwSide._dCLS, _polarFwSide._dCM, _polarFwSide._dCM);
};

void CompactAeroModel2::ApplyForce(const Vector3& vair0, const Matrix3& craftRot, const Vector3& craftOmega, Vector3& force, Vector3& torq, float elevator, float rudder, float aileron, const Vector3& pos0) const
{
  Vector3 uair = vair0;    // we shoud consider craftOmega
#if _DEBUG_AEROSURF
                           //LogF( "CompactAeroModel2::applyForce ID10: S %f %f %f L %f %f %f D %f %f %f ", S[0], S[1], S[2], L[0], L[1], L[2], D[0], D[1], D[2] );
                           //DIAG_MESSAGE_ID(100, 10, Format("ID10: S %f %f %f L %f %f %f D %f %f %f ", S[0], S[1], S[2], L[0], L[1], L[2], D[0], D[1], D[2]));
  renderer.Add3dLine(pos0, pos0 + uair * 5.0, PackedBlue);
  renderer.Add3dLine(pos0, pos0 + craftRot.DirectionAside() * 5.0, PackedColor(0xff808080));
  renderer.Add3dLine(pos0, pos0 + craftRot.DirectionUp()*5.0, PackedBlack);
  renderer.Add3dLine(pos0, pos0 + craftRot.Direction()*5.0, PackedWhite);
#endif // _DEBUG_AEROSURF
  //Vector3 uair = vair0 + gdpos.CrossProduct(craftOmega);                                    
  float vrair2 = uair.SquareSize();
  torq = VZero;
  force = VZero;
  if (vrair2 > lowSpeedCuoff)
  {
    float vrair = sqrt(vrair2);
    uair *= (1 / vrair);
    // decompose uair to panel coordinates
    float ca = uair.DotProduct(craftRot.DirectionAside());
    float cb = uair.DotProduct(craftRot.DirectionUp());
    float cc = uair.DotProduct(craftRot.Direction());
#if _DEBUG_AEROSURF
    DIAG_MESSAGE_ID(100, 11, Format("elevator %f rudder %f aileron %f ", elevator, rudder, aileron));
    DIAG_MESSAGE_ID(100, 12, Format("ca %f cb %f cc %f ", ca, cb, cc));
#endif // _DEBUG_AEROSURF

    float prefactor = vrair2 * _area;
    float CD, CL, CMup, CMside;
    Vector3 airUp;
    // ==== main wing polar (lift along aircraft up vector);     Pitch control
    { 
      _polarFwUp.GetAeroCoefs(-cc, cb, elevator*0.25, CD, CL, CMup);
      if ((cb*cb) > 1e-8)
      {
        airUp = craftRot.DirectionUp() + (uair * -cb); // component of grot.Up perpendicular to uair
        airUp.Normalize();
      }
      else
      {
        airUp = craftRot.DirectionUp();
      }
#if _DEBUG_AEROSURF
      renderer.Add3dLine(pos0, pos0 + airUp * 5.0, PackedGreen);
#endif // _DEBUG_AEROSURF
      force += (airUp * CL + uair * CD) * prefactor;
      torq += airUp.CrossProduct(uair) * -CMup * prefactor;
    }
    // ====  keel wing polar (lift along aircraft side vector);    Yaw control
    { 
      _polarFwSide.GetAeroCoefs(-cc, ca, 0, CD, CL, CMside);
      if ((ca*ca) > 1e-8)
      {
        airUp = craftRot.DirectionAside() + (uair * -ca); // component of grot.Up perpendicular to uair
        airUp.Normalize();
      }
      else
      {
        airUp = craftRot.DirectionAside();
      }
#if _DEBUG_AEROSURF
      renderer.Add3dLine(pos0, pos0 + airUp * 5.0, PackedColor(0xffff00ff));
#endif // _DEBUG_AEROSURF
      force += (airUp * CL + uair * CD) * prefactor * 0.2;
      torq += airUp.CrossProduct(uair) * -CMside * prefactor;
    }
#if _DEBUG_AEROSURF
    DIAG_MESSAGE_ID(100, 13, Format(" CMup %f CMside %f ", CMup, CMside));
    renderer.Add3dLine(pos0, pos0 + force * 5.0, PackedRed);
#endif // _DEBUG_AEROSURF

    float cc2 = cc * cc;  // factor used to damp torques for non-standard flight regimes ( large angle between vair and grot.Direction() )

    torq += craftRot.Direction() * (_cTorqControl[2] * aileron * prefactor * cc2);   // Roll control - we still don't do this by Polar

    // damping of aircraft rotation due to drag of control surfaces
    torq += craftRot.DirectionAside()* ( craftOmega.DotProduct(craftRot.DirectionAside()) * -_cTorqDamp[0] * prefactor); // pitch damp
    torq += craftRot.DirectionUp()*    ( craftOmega.DotProduct(craftRot.DirectionUp())    * -_cTorqDamp[1] * prefactor); // yaw   damp
    torq += craftRot.Direction()*      ( craftOmega.DotProduct(craftRot.Direction())      * -_cTorqDamp[2] * prefactor); // roll  damp

    // apply anisotropic drag coef
    force += craftRot.DirectionAside() * (_cAnisoDrag[0] * ca*prefactor);
    force += craftRot.DirectionUp()    * (_cAnisoDrag[1] * cb*prefactor);
    force += craftRot.Direction()      * (_cAnisoDrag[2] * cc*prefactor);

#if _DEBUG_AEROSURF
    DIAG_MESSAGE_ID(100, 14, Format("force %f %f %f torq %f %f %f ", force.X(), force.Y(), force.Z(), torq.X(), torq.Y(), torq.Z()));
#endif // _DEBUG_AEROSURF

  }
  float fDens = 0.5*Atmosphere::GetDensity(pos0[1]); // hope that _position[1]==0 is ground level ?
  force *= fDens;
  torq *= fDens;
};
