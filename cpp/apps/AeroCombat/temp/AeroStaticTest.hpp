#pragma once

#include "wpch.hpp"
#include "../../../VBS3_Git/lib/Debug/debugGui.hpp"

#include "AeroSurf.hpp"

class ControlFeedback {
  public:
  float y0 = 0.0;
  float oy = 0.0;
  float Ky = 1.0;
  float Kdy = 1.0;

  inline void setup(float y0_, float Ky_, float Kdy_) { y0 = y0_; Ky = Ky_; Kdy = Kdy_; }

  inline float update(float y) 
  {
    float Dy = y - y0;
    float dy = y - oy;
    float dx = Dy * Ky + dy * Kdy;
    oy = y;
    return dx;
  }
  ControlFeedback() {};
  ControlFeedback(float y0, float Ky, float Kdy) { setup(y0, Ky, Kdy); };
};

class AeroPerformanceRecord 
{ 
  public:
  // I don't use auto arrays here intentionally, 
  // this way I can check if pointer == null, and more easily allocate/dealocate
  constexpr static int nLines = 13;
  int _n = 0;
  union 
  {
    float *_datalines[nLines];
    struct 
    {
      float *_elevs, *_AoAs, *_Ls, *_Ds, *_LDs, *_turnVs, *_turnTs, *_turnRs, *_turnSins, *_climbVs, *_climbVxs, *_climbVys, *_climbSins;
    };
  };

  //ClassIsGeneric(AeroPerformanceRecord); // https://docs.bisimulations.com/display/~Marek+Manukjan/Essence+library
  ClassIsBinary(AeroPerformanceRecord); // https://docs.bisimulations.com/display/~Marek+Manukjan/Essence+library

  AeroPerformanceRecord() { _n = 0; for (int i = 0; i < nLines; i++) _datalines[i] = 0; };
  AeroPerformanceRecord(int n) { Alloc(n); }
  ~AeroPerformanceRecord() { Dealloc(); };

  void Alloc(int n);
  void Dealloc();
  void ToFile(const char * fname) const;

};

class AeroTestPlatform : public RigidBody 
{ 
public:

  // values of control inputs
  float _c_elevator=0, _c_rudder=0, _c_aileron=0, _c_propel = 0, _c_break = 0, _c_throtle=1.0;

  // Deffierent aerodynamic modes - we use pointer so we can bind it to particular AirplaneType
  ComposedAeroModel* _compAeroModel = nullptr;
  CompactAeroModel*  _aeroModel     = nullptr;
  CompactAeroModel2* _aeroModel2    = nullptr;
  int _useAeroModelType = 1; // which aeroModel to use ?

  // Orientation
  int _nTrj=0,_iTrj=0;
  Vector3 * _trjPos = nullptr;
  Matrix3 * _trjRot = nullptr;

  //AutoArray<AeroPerformanceRecord> records; // this crashed for some reason
  AutoArray<AeroPerformanceRecord*> records;
  
  // Plots used to present aerodynamic preformace
  Ref<DebugGraphXY> _LDplot = 0;
  Ref<DebugGraphXY> _turnTimePlot = 0;
  Ref<DebugGraphXY> _turnRadiusPlot = 0;
  Ref<DebugGraphXY> _climbPlot = 0;
  
  // === function declaration
  void GetAeroForce(const Vector3& vwind, Vector3& force, Vector3& torq);

  int  InitWingsFile(const char * fname);                                                                     // load wings from text file

  // These functions does no ensure stable flight regime ( torque = 0;  )
  void SamplePolar(const char* fname, int nsamp, float angSpan, float speed, Vector3 v0, Vector3 axis);      // calculate for varying AngleOfAttack
  void FindLDMax(float dAoA0, float convAcc, float& AoA, float& L, float& D);                                // find AoA with maximum glide ratio (Lift/Drag) ( not necessarily stable torque != 0 )
  void SolveStraightFlight(int nmax, float convAcc, float Thrust, float Weight, float& speed, float& AoA);   // find AoA and Speed for which Lift==Weight     ( not necessarily stable torque != 0 )
  float ScanOmegaRoll(int n, float speed, float omega0, float dOmega, float aileron);                        // evaluate torque at given speed and angular velocity (assuming pre-set aileron)

  // Functions to relax pitch moment ( for pre-set elevator )
  float RelaxPitch(int nmax, float K, float maxAng, float convAcc, Vector3 uair, float speed);                // find elevator position for which pitching moment==0 ( torque==0 => stable flight regime )
  float RelaxPitchRF(float pitchMin, float pitchMax, float accConv, Vector3& force);                          //  ----,,-----  using root finder ( perhaps more robust (stable) algorithm )

  // These functions evaluet stable flight regimes ( torque = 0; for preset control )
  float RelaxRollRate(int nmax, float speed, float convAcc, float aileron);              // for given value of aileron find saturated roll rate at given speed (saturate = ignoring initial moment of intertia, just aerodynamics) 
  int   RunRelaxedTest(const char* fnameIn, const char* fnameOut);                       // for controls given in input file relax aircraft orientation to achieve stable flight and store results to output file

  // These function works assuming velocity independent drag and lift ( e.g. not Supersonic )
  float SamplePerformance(float elevMin, float elevMax, AeroPerformanceRecord* rec );                                        // Evaluate main performance characteristics of aircraft in stable fligt regimes ( climb rate, turn rate ... ) 
  void  FindMinTurnTime(float elevMin, float elevMax, float& AoA, float& Lift, float& Drag, float& TurnTime, float& TurnRadius, float& sinTheta, float& speed);  // find minimal turn time (flat turn with gravity) without need of precalculated polar
  void  FindMaxClimbRate(float elevMin, float elevMax, float& AoA, float& Lift, float& Drag, float& vx, float& ca, float& sa);                                   // find max climb rate (stright fight with gravity) without need of precalculated polar

  // These functions works with velocity dependent Drag and Lift coef ( e.g. Supersonic )
  void  SampleSpeed(int n, float vmin, float dv, float* Ls, float* Ds, float* vs = NULL);                     // Evaluate Lift and Drag for different speeds - this make sense only with non-linear corrections (e.g. Supersonic), otherwise it is just parabola ~v^2

  // Brute force rigid-body time-steping simulation  - should work for velocity dependnet Lift and Drag coef
  void NewTrj( int n );
  void Fly        ( int n, float dt );                          // fly n iterations with fixed (pre-set) values of controls
  void FlyStraight(int n, int nsub, float dt, Vector3 vGoal);   // fly n iterations with control feedback keeping it stright

#if _DEBUG_AEROSURF
  void DrawTrj(int i0, int n);
#endif //_DEBUG_AEROSURF

  ~AeroTestPlatform();
  
};