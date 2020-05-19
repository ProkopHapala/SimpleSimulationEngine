#pragma once

//#include "math3d.hpp"
#include "wpch.hpp"
#include "AeroMath.hpp"

#include "../../../VBS3_Git/lib/Debug/debugGui.hpp"

#define _DEBUG_AEROSURF (1 && _ENABLE_CHEATS)

float const lowSpeedCuoff = 1e-8;

#if _DEBUG_AEROSURF
void drawWing(Vector3 pos, Matrix3 rot, float a, float aspectRatio );
void drawRigidBody(Vector3 pos, Matrix3 rot, float sz);
#endif // _DEBUG_AEROSURF

//==============
//  PolarModel   - store List and Drag curve
//==============

class PolarModel 
{ 
public:
  
  float _CD0    = 0.02;           // residual drag for zero angle of attack [1]
  float _dCD    = 0.9;            // drag curve slope before stall  [1/Radian]
  float _dCDS   = 0.9;            // drag curve slope after stall   [1/Radian]
  float _dCL    = 6.28;           // lift curve slope before stall  [1/Radian]
  float _dCLS   = 2.82743338823;  // lift curve slope after stall   [1/Radian]
  float _sStall = 0.04;           // stall position                 [Radian]
  float _wStall = 0.4;            // width of stall transition      [Radian]
  
  //ClassIsGeneric(PolarModel); // https://docs.bisimulations.com/display/~Marek+Manukjan/Essence+library
  //ClassIsSimple(PolarModel); // https://docs.bisimulations.com/display/~Marek+Manukjan/Essence+library
  ClassIsBinary(PolarModel);

  // ==== inline functions

  inline void Set(float sStall, float wStall, float CD0, float dCD, float dCDS, float dCL, float dCLS) { _wStall = wStall; _sStall = sStall;_CD0  = CD0; _dCD  = dCD;_dCDS = dCDS;_dCL = dCL;_dCLS = dCLS;};
  inline void GetAeroCoefs(float ca, float sa, float& CD, float& CL) const;
  
};

//=====================
//  PolarModelTorque      - is polar including pitching moment CM, it is currently used just in CompatModel2 
//=====================

class PolarModelTorque : public PolarModel
{
public:
  float _dCM  = 5.0;     // pitching moment slope before stall
  float _dCMS = 2.0;     // with of stall transition     [1/Radian]

  //ClassIsGeneric(PolarModelTorque); // https://docs.bisimulations.com/display/~Marek+Manukjan/Essence+library
  //ClassIsSimple(PolarModelTorque);
  ClassIsBinary(PolarModelTorque);

  // ==== inline functions
  inline void Set(float sStall, float wStall, float CD0, float dCD, float dCDS, float dCL, float dCLS, float dCM, float dCMS) { 
    _wStall = wStall; _sStall = sStall; _CD0 = CD0; _dCD = dCD; _dCDS = dCDS; _dCL = dCL; _dCLS = dCLS; _dCM =dCM; _dCMS = dCMS;
  };

  inline void GetAeroCoefs(float ca, float sa, float flap, float& CD, float& CL, float& CM) const;
  
};

//=====================
//  AeroSurface          - a Wing exerting some aerodynamic force and torque depending on incident airflow
//=====================

//class AeroSurface : public KinematicBody {
class AeroSurface 
{ 
public:
  static constexpr float SAFETY_v = 1e-6;
  float _area = 1.0;
  float _aspectRatio = 6.0; 
  Vector3 _lpos = VZero;
  Matrix3 _lrot = M3Identity;
  Vector3 _cAnisoDrag = VZero; // additional anisotropic drag coefficient
  PolarModel _polar;

  bool _bFiniteWingCorrection = true;
  float _minTilt=-0.4;
  float _maxTilt=0.4;

  bool _bControlable=false;
  float _cRudder = 0.0, _cElevator = 0.0, _cAileron = 0.0; // this will allow us to bind any linear combination of controls

  // ======= DEBUG GUI
#if _DEBUG_AEROSURF
  int debugPlot = 0;
  Ref<DebugMultiGraph> liftPlot;
  Ref<DebugMultiGraph> dragPlot;
  int iLift=-1, iDrag=-1;
#endif // _DEBUG_AEROSURF
  
  //ClassIsGeneric(AeroSurface); // https://docs.bisimulations.com/display/~Marek+Manukjan/Essence+library
  //ClassIsSimple(AeroSurface);
  ClassIsBinary(AeroSurface);
  
  // ==== function declarations

  // ==== inline functions

  inline void GlobalPos(const Vector3& pos0, const Matrix3& rot0, Vector3& gpos) const { 
    gpos = rot0 * _lpos; 
    gpos += pos0;
  };

  // ==== function definitions

  // get aerodynamic force and torque (not yet scalled by air density)
  void ApplyForce(const Vector3& vair0, const Matrix3& craftRot, const Vector3& craftOmega, Vector3& force, Vector3& torq, float tilt, int i, bool debug, const Vector3& pos0 ) const;

  void FromString(const char * str);
  void FromStringPolarModel(const char * str);
  int  ToStringPolarModel(char * str) const;
  void Tilt(float angle);                       // rotate surface along its aside axis
  float Control2tilt(float control) const;      // tilt angle dependence of control input
  float Controls2tilt(float elevator, float rudder, float aileron) const;   // tilt angle dependence of control inputs
};

//=====================
//  Propeler            
//=====================

class Propeler
{
public:
  Vector3 _lpos = VZero;
  Vector3 _ldir = Vector3(0.0,0.0,1.0);
  float _thrustRocket=0.0;  // [N]

  const bool _bSolveBypass = false;
  float _area = 1.0;       // [m^2]
  float _power = 1.0;      // [W]
  float _efficiency = 1.0; // [1]
  float _CD = 0.05;        // [1]
  float _vstatic = 50.0;   // [m/s] velocity of stream for properler stationary with respect to air

  //ClassIsGeneric(Propeler); // https://docs.bisimulations.com/display/~Marek+Manukjan/Essence+library
  //ClassIsSimple(Propeler);
  ClassIsBinary(Propeler);

  // ==== function definitions

  // get aerodynamic force and torque (already scalled by air density)
  void  ApplyForce(const Vector3& vair0, float airDens, const Matrix3& craftRot, const Vector3& craftOmega, Vector3& force, Vector3& torq, float control, int i, bool debug, const Vector3& pos0) const;
  float GetBaypassThrust(float v0, float control) const;  // Solve propeler thrust for given velocity
  float AutoVStatic()const;                               // Get estimate of airflow speed through propeler at zero velocity 

  void FromString(const char * str);

};

//=====================
//  ComposedAeroModel
//=====================

class ComposedAeroModel
{
public:
  AutoArray<AeroSurface> _wings;
  AutoArray<Propeler> _propelers;

  // which AeroSurface is bind to which control ?
  //   this is replaced by  AeroSurf::_cElevator,_cRudder,_cAileron
  //int _leftAirelonId = -1;
  //int _rightAirelonId = -1; 
  //int _elevatorId = -1;
  //int _rudderId = -1;

  float _engineThrust = 80.64e+3;  // [N]

  // Supersonic correction - we do this on Aircraft level rather than wing level - i
  bool _bSupersonic = false;
  float _MachMin   = 0.8;
  float _MachMax   = 1.1;
  float _Mach0     = 0.7;
  float _cWaveDrag = 1.0;

  int _DebugLevel = 0;

  //ClassIsGeneric(ComposedAeroModel); // https://docs.bisimulations.com/display/~Marek+Manukjan/Essence+library
  //ClassIsSimple(ComposedAeroModel); // https://docs.bisimulations.com/display/~Marek+Manukjan/Essence+library
  ClassIsBinary(ComposedAeroModel);

  // ==== function definitions
  
  void ApplyForce(const Vector3& vair, const Matrix3& craftRot, const Vector3& craftOmega, Vector3& force, Vector3& torq, float elevator, float rudder, float aileron, float throtle, const Vector3& pos0) const;
  void FromFile(FILE* pFile);
  int LoadFile(const char * fname);
  char* ToString(char * str) const;

  float WettedArea() const;          // sum area of all wings
  float MultArea(float f);           // rescale all wing by the same factor
  float GetMaxTorqScale() const;     // maximim torque exerted by any wing

};

//=====================
//  CompactAeroModel
//=====================

// Fast Aerodynamic model #1 using only one wing
//  - control inputs exert same torq even in deep stall which is not realistic
//    => Aircraft is much more maneuverable in non-standard flight regimes than sould be in reality

class CompactAeroModel
{
public:
  float _area = 1.0;
  Vector3 _cAnisoDrag, _cTorqControl, _cTorqDamp, _cTorqRelax;
  PolarModel _polarFwUp;

  //ClassIsGeneric(CompactAeroModel); // https://docs.bisimulations.com/display/~Marek+Manukjan/Essence+library
  //ClassIsSimple(CompactAeroModel);
  ClassIsBinary(CompactAeroModel);

  // ==== function definitions

  // get aerodynamic force and torque (not yet scalled by air density)
  void ApplyForce(const Vector3& vair0, const Matrix3& craftRot, const Vector3& craftOmega, Vector3& force, Vector3& torq, float elevator, float rudder, float aileron, const Vector3& pos0) const;
  void FromString(const char * str);

};

//=====================
//  CompactAeroModel2
//=====================

// Fast Aerodynamic model #2 using only one wing
//  - Pitch and Yaw control (elevator and rudder) are solved using pitching moment coeffeciet (CM) in Polar 
//  - Therefore control is more realistic in stall and other non-standard flight-regimes (=> less controlable)

class CompactAeroModel2
{
public:
  float _area = 1.0;
  Vector3 _cAnisoDrag, _cTorqControl, _cTorqDamp, _cTorqRelax;
  PolarModelTorque _polarFwUp;
  PolarModelTorque _polarFwSide;

  //ClassIsGeneric(CompactAeroModel2); // https://docs.bisimulations.com/display/~Marek+Manukjan/Essence+library
  //ClassIsSimple(CompactAeroModel2);
  ClassIsBinary(CompactAeroModel2);

  // ==== function definitions

  // get aerodynamic force and torque (not yet scalled by air density)
  void ApplyForce(const Vector3& vair0, const Matrix3& craftRot, const Vector3& craftOmega, Vector3& force, Vector3& torq, float elevator, float rudder, float aileron, const Vector3& pos0) const;
  void FromString(const char * str);

};