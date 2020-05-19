
#include <math.h>
#include "AeroStaticTest.hpp"

// =========================
//  AeroPerformanceRecord 
//  ========================

void AeroPerformanceRecord::Alloc(int n) 
{
  _n = n;
  for (int i = 0; i < nLines; i++) { _datalines[i] = new float[_n]; }
};

void AeroPerformanceRecord::Dealloc() 
{
  for (int i = 0; i < nLines; i++) { delete _datalines[i]; _datalines[i] = 0; }
}

void AeroPerformanceRecord::ToFile(const char * fname) const
{
  FILE * fout = fopen(fname, "w");
  for (int i = 0; i < _n; i++)
  {
    fprintf(fout, " %i  %f    %f %f %f    %f %f %f %f   %f %f \n", i, _elevs[i], _AoAs[i], _Ls[i], _Ds[i], _turnVs[i], _turnRs[i], _turnTs[i], _turnSins[i], _climbVs[i], _climbSins[i]);
  }
  fclose(fout);
}

// =========================
//   AeroTestPlatform
//  ========================


void AeroTestPlatform::GetAeroForce(const Vector3& vwind, Vector3& force, Vector3& torq) 
{
  switch (_useAeroModelType) 
  {
    case 1:  _compAeroModel->ApplyForce(vwind - _velocity, _orientation, _angVelocity, force, torq, _c_elevator, _c_rudder, _c_aileron, _c_throtle, _position); break;
    case 2:  _aeroModel->ApplyForce(vwind - _velocity, _orientation, _angVelocity, force, torq, _c_elevator, _c_rudder, _c_aileron, _position ); break;
    case 3:  _aeroModel2->ApplyForce(vwind - _velocity, _orientation, _angVelocity, force, torq, _c_elevator, _c_rudder, _c_aileron, _position ); break;
  }
};

int AeroTestPlatform::InitWingsFile(const char * fname) 
{
  const int nbuf = 1024;
  char buf[nbuf];
  FILE * pFile;
  LogF(" AeroTestPlatform::InitWingsFile loading wings from: >>%s<<\n", fname);
  pFile = fopen(fname, "r");
  if (pFile == nullptr) 
  {
    LogF(" AeroTestPlatform::InitWingsFile cannot open >>%s<< => call default InitWings() \n", fname);
    return -1;
  }
  _compAeroModel->FromFile(pFile);
  fgets(buf, nbuf, pFile); _aeroModel->FromString(buf); // CompactAeroModel  not used
  fgets(buf, nbuf, pFile); _aeroModel2->FromString(buf); // CompactAeroModel2 is  used
  LogF("AeroTestPlatform::InitWingsFile wings loaded from >>%s<<\n", fname);
  fclose(pFile);
  return 0;
  //exit(0);
}

void AeroTestPlatform::SamplePolar( const char * fname, int nsamp, float angSpan, float speed, Vector3 v0, Vector3 axis) 
{
  _velocity = VZero;
  Matrix3 drot;
  float dang = angSpan / nsamp;
  drot.SetRotationAxis( axis, dang );
  Vector3 uair = v0; uair.Normalize();
  _c_throtle = 0.0;
  
  LogF("samplePolar output file: %s \n", fname );
  FILE * pFile = fopen(fname, "w");
  float ang = 0.0;
  for (int i = 0; i < nsamp; i++) 
  {
    // get aerodynamioc force and torq
    Vector3 force, torq;
    GetAeroForce(uair*speed, force, torq); // Get force
    
    // store everything
    Vector3 airUp   = _orientation.DirectionUp() + (uair * -uair.DotProduct( _orientation.DirectionUp() ) );
    Vector3 airSide = uair.CrossProduct(airUp);
    float Drag  = force.DotProduct(uair);
    float Lift  = force.DotProduct(airUp);
    float Fside = force.DotProduct(airSide);
    float Tpitch = torq.DotProduct( _orientation.DirectionAside() );
    float Tyaw   = torq.DotProduct( _orientation.DirectionUp() );
    float Troll  = torq.DotProduct( _orientation.Direction() );
    LogF("samplePolar %i  %f  v:  %f %f %f  F:  %f %f %f  T:  %f %f %f  \n", i, ang, uair.X(), uair.Y(), uair.Z(), force.X(), force.Y(), force.Z(), torq.X(), torq.Y(), torq.Z() );
    fprintf( pFile, " %i  %f    %f %f %f    %f %f %f    %f %f %f  \n", i, ang, uair.X(), uair.Y(), uair.Z(), Drag, Lift, Fside, Tpitch, Tyaw, Troll );
    
    // update air direction
    ang += dang;
    uair = drot * uair; 
  };
  fclose(pFile);
};

float AeroTestPlatform::ScanOmegaRoll(int n, float speed, float omega0, float dOmega, float aileron)
{
  Vector3 uair = Vector3(0.0, 0.0, -1.0);
  Vector3 vair = uair * speed;
  _c_aileron = aileron;
  _c_rudder = 0.0;
  _c_elevator = 0.0;
  _velocity = VZero;
  _orientation = M3Identity;
  _angVelocity = Vector3(0, 0, omega0);
  LogF("AeroTestPlatform::scanOmegaRoll  airelon %f speed %f ", _c_aileron, speed);
  for (int i = 0; i < n; i++)
  {
    _angVelocity[2] += dOmega;
    Vector3 force, torq;
    GetAeroForce(vair, force, torq);
    float omega = _angVelocity.Size();
    float lTorq = torq.Size();
    LogF("%i omega |%f| %f %f %f torq |%f| %f %f %f ", i, omega, _angVelocity.X(), _angVelocity.Y(), _angVelocity.Z(), lTorq, torq.X(), torq.Y(), torq.Z());
    //if (lTorq < convAcc) {return omega;}
  }
  return NAN;
};

float AeroTestPlatform::RelaxPitch(int nmax, float K, float maxAng, float convAcc, Vector3 uair, float speed ) 
{
  for (int i = 0; i < nmax; i++) 
  {
    Vector3 force, torq;
    GetAeroForce(uair*speed, force, torq);
    Matrix3 rot;  
    float torqLength = torq.Size();
    if (torqLength > convAcc) 
    { // step of root finder - rotate in direction of torq, How much ?   torqLength*K
      Vector3 utorq = torq;
      utorq *= (1 / torqLength);
      float dang = torqLength*K;
      if (dang > maxAng) dang=maxAng;
      rot.SetRotationAxis(utorq, dang );
      //LogF(" |torq|=%f (%f,%f,%f) dang %f ", torqLength, torq.X(), torq.Y(), torq.Z(), dang );
      _orientation = rot * _orientation;
    }
    else 
    {
      //LogF(" |torq|(%f)< %f => relaxPitch DONE! ", torqLength, convAcc );
      return torqLength;
    }
  }
  LogF("AeroTestPlatform::relaxPitch not Converged in %i iterations ", nmax );
  exit(-1);
};

float AeroTestPlatform::RelaxPitchRF(float pitchMin, float pitchMax, float accConv, Vector3& force ) 
{
  Vector3 vair = Vector3(0.0, 0.0, -1.0);
  Vector3 torq;
  RootFinder rf;
  float x1 = pitchMin, x2 = pitchMax;
  float y1, y2;

  // eval min bound
  _orientation.SetRotationX(-x1); 
  GetAeroForce(vair, force, torq); 
  y1 = torq[0];

  // eval max bound
  _orientation.SetRotationX(-x2); 
  GetAeroForce(vair, force, torq);  
  y2 = torq[0];

  rf.Init(x1, y1, x2, y2);  // init root finder

  for (int j = 0; j<100; j++) 
  {
    // eval y=f(x)
    _orientation.SetRotationX(-x2); 
    GetAeroForce(vair, force, torq);  
    y2 = torq[0];

    // predict new 'x'
    x2 = rf.Step(x2, y2);  
    if (fabs(y2) < accConv) 
    {
      break;
    }
    //LogF("relaxPitchRF %i elev %f AoA %f torq.x %f fw.y %f ", j, _c_elevator, x2, y2, _orientation.Direction()[1]);
  }
  return x2;
};

float AeroTestPlatform::RelaxRollRate(int nmax, float speed, float convAcc, float aileron)
{
  convAcc *= speed * speed;
  Vector3 uair = Vector3(0.0, 0.0, -1.0);
  Vector3 vair = uair * speed;
  _c_aileron = aileron;
  _c_rudder = 0.0;
  _c_elevator = 0.0;
  _velocity = VZero;
  _orientation = M3Identity;
  _angVelocity = VZero;
  //LogF("AeroTestPlatform::relaxRollRate  airelon %f speed %f ", _c_aileron, speed );
  Vector3 force, torq;
  float x1 = -speed, x2 = speed;
  float y1, y2;

  // eval min bound
  _angVelocity[2] = x1;
  GetAeroForce(vair, force, torq);
  y1 = torq[2];

  // eval max bound
  _angVelocity[2] = x2;
  GetAeroForce(vair, force, torq);
  y2 = torq[2];

  RootFinder rf;
  rf.Init(x1, y1, x2, y2);  // init root finder
  for (int i = 0; i < nmax; i++)
  {
    // eval y=f(x)
    _angVelocity[2] = x2;
    GetAeroForce(vair, force, torq);
    y2 = torq[2];

    // predict new 'x'
    x2 = rf.Step(x2, y2);
    //LogF("%i omega %f torq %f maxErr(%f) ", i, x2, y2, convAcc );
    if (fabs(y2) < convAcc)
      return x2;
  }
  return NAN;
};

float AeroTestPlatform::SamplePerformance(float elevMin, float elevMax, AeroPerformanceRecord* rec)
{
  float torqSc = _compAeroModel->GetMaxTorqScale();
  float accConv = torqSc * 1e-6;
  //LogF("MaxTorqScale %f accConv %f", torqSc, accConv);
  Vector3 force;
  float pitchMin = -0.3, pitchMax = -pitchMin;
  _orientation = M3Identity;
  _angVelocity = Vector3(0.0, 0.0, 0.0);
  _velocity    = Vector3(0.0, 0.0, 0.0);
  _c_rudder    = 0.0;;
  _c_aileron   = 0.0;
  float delev = (elevMax - elevMin) / (rec->_n - 1); // elevator step
  float mass = 1 / _invMass;
  LogF("samplePerformanceRecord %i thrust %f mass %f torqSc %f ", rec->_n, _compAeroModel->_engineThrust, mass, torqSc );
  for (int i = 0; i < rec->_n; i++) 
  { 
    _c_elevator = elevMin + i * delev;     // set elevator
    float AoA = RelaxPitchRF( pitchMin, pitchMax, accConv, force);  // find AoA for given elevator 

    // Store everything
    float L= force[1], D= -force[2];
    rec->_elevs[i] = -_c_elevator;
    rec->_AoAs[i] = AoA;
    rec->_Ls[i] = L;
    rec->_Ds[i] = D;
    rec->_LDs[i] = L/D;
    float turnR = NAN, turnT = NAN, turnV = NAN, turnSin = NAN;
    float climbV= NAN, climbSin= NAN;
    if (rec->_turnVs) 
    {
      evalThrustLimitedTurnGrav(L, D, _compAeroModel->_engineThrust, mass, turnV, turnR, turnSin);
      turnV = sqrt(turnV);
      turnT = turnRadiusToTime(turnR, turnV);
      rec->_turnVs[i] = turnV;
      rec->_turnTs[i] = turnT;
      rec->_turnRs[i] = turnR;
      rec->_turnSins[i] = turnSin;
    }
    if (rec->_climbVs) 
    {
      float ca;
      evalClimbRate(L, D, _compAeroModel->_engineThrust, mass, climbV, ca, climbSin );
      rec->_climbVs[i]   = climbV;
      rec->_climbVxs[i]  = climbV * ca;
      rec->_climbVys[i]  = climbV * climbSin;
      rec->_climbSins[i] = climbSin;
    }

    LogF("samplePerformanceRecord>> %i elev %f (AoA,L,D) %f %f %f Turn(V,R,T,sa) %f %f %f %f climb(v,sa) %f %f ", i, rec->_elevs[i], rec->_AoAs[i], rec->_Ls[i], rec->_Ds[i], turnV, turnR, turnT, turnSin, climbV, climbSin);
  }
  return delev;
}

void AeroTestPlatform::FindMaxClimbRate(float elevMin, float elevMax, float& AoA, float& Lift, float& Drag, float& vx, float& vy, float& sinAlfa )
{
  const float pitchMin = -0.3, pitchMax = -pitchMin;
  const float mass = 1 / _invMass;
  const float Thrust = _compAeroModel->_engineThrust;  // FIXME: what if _compAeroModel==null and we use different AeroModel ?
  const float accConvTorq = _compAeroModel->GetMaxTorqScale() * 1e-5;
  const float accConvPitch = 1e-3;
  LogF("findMaxClimbRate Thrust %f mass %f accConvTorq %f ", Thrust, mass, accConvTorq);
  Vector3 force;
  Optimizer1D_Golden opt;
  _c_elevator = opt.init(elevMin, elevMax); // initialize optimized and get first estimate of optimum _c_elevator
  for (int i = 0; i<100; i++) 
  {
    AoA = RelaxPitchRF(pitchMin, pitchMax, accConvTorq, force); // for each value of elevator we must fully relax aircraft orientation to ensure stable flight regime (torq==0)
    
    // y=f(x)   climbRate=f(_c_elevator)
    float D = -force[2];
    float L =  force[1];
    evalClimbRate(L, D, Thrust, mass, vx, vy, sinAlfa );   

    // step estimating optimal _c_elevator
    _c_elevator = opt.step(-vy);                           
    float err = opt.errX();
    //LogF("findMaxClimbRate %i err %f elev %f AoAs %f TurnTime %f TurnRadius %f sinTheta %f speed %f ", i, err, _c_elevator, AoA, vx, vy, sinAlfa );
    if (err < accConvPitch) break;
  }
  LogF("findMaxClimbRate>> elev %f AoAs %f TurnTime %f TurnRadius %f sinTheta %f speed %f ", _c_elevator, AoA, vx, vy, sinAlfa );
};

void AeroTestPlatform::FindMinTurnTime(float elevMin, float elevMax, float& AoA, float& Lift, float& Drag, float& TurnTime, float& TurnRadius, float& sinTheta, float& speed) 
{
  const float mass = 1 / _invMass;
  const float Thrust = _compAeroModel->_engineThrust;  // FIXME: what if _compAeroModel==null and we use different AeroModel ?
  const float pitchMin = -0.3, pitchMax = -pitchMin;
  const float accConvTorq  = _compAeroModel->GetMaxTorqScale() * 1e-5;
  const float accConvPitch = 1e-3;
  LogF("findMinTurnTime Thrust %f mass %f accConvTorq %f ", Thrust, mass, accConvTorq);
  Vector3 force;
  Optimizer1D_Golden opt;
  _c_elevator = opt.init(elevMin, elevMax); // initialize optimized and get first estimate of optimum _c_elevator
  float v2;
  for (int i = 0; i < 100; i++) 
  {
    AoA = RelaxPitchRF(pitchMin, pitchMax, accConvTorq, force); // for each value of elevator we must fully relax aircraft orientation to ensure stable flight regime (torq==0)
    
     // y=f(x) : TurnTime=f(_c_elevator)
    float D = -force[2];
    float L =  force[1];
    evalThrustLimitedTurnGrav(L, D, Thrust, mass, v2, TurnRadius, sinTheta);
    speed    = sqrt(v2);
    TurnTime = turnRadiusToTime(TurnRadius, speed);   

    // step estimating optimal _c_elevator
    _c_elevator = opt.step(TurnTime);                
    float err = opt.errX();
    //LogF("findMinTurnTime %i err %f elev %f AoAs %f TurnTime %f TurnRadius %f sinTheta %f speed %f ", i, err, _c_elevator, AoA, TurnTime, TurnRadius, sinTheta, speed);
    if (err < accConvPitch) break;
  }
  LogF("findMinTurnTime>> elev %f AoAs %f TurnTime %f TurnRadius %f sinTheta %f speed %f ", _c_elevator, AoA, TurnTime, TurnRadius, sinTheta, speed );
};

void AeroTestPlatform::SampleSpeed(int n, float vmin, float dv, float* Ls, float* Ds, float* vs ) 
{
  Vector3 uair = Vector3(0.0, 0.0, -1.0);
  FILE* fout = fopen("SampleSpeed.log","w");
  for (int i = 0; i<n; i++) 
  {
    Vector3 force, torq;
    float speed = vmin + dv * i;
    GetAeroForce(uair*speed, force, torq);
    if (Ds) Ds[i] = -force[2];
    if (Ls) Ls[i] =  force[1];
    if (vs) vs[i] =  speed;
    LogF(         " SampleSpeed: %i |v| %f  L %f D %f  L/D %f  \n", i, speed, force[1], -force[2], force[1]/-force[2] );
    fprintf(fout,"%i  %f   %f %f  %f  \n", i, speed, force[1], -force[2], force[1] / -force[2] );
  }
  fclose(fout);
}

int AeroTestPlatform::RunRelaxedTest( const char* fnameIn, const char* fnameOut)
{
  FILE * fin  = fopen(fnameIn, "r");
  if (fin == nullptr) 
  {
    LogF("AeroTestPlatform::runRelaxedTest cannot open >>%s<< \n", fnameIn);
    //exit(-1);
    return -1;
  }
  FILE * fout = fopen(fnameOut, "w");
  int n;
  fscanf( fin, "%i\n", &n );
  float tmax = _compAeroModel->GetMaxTorqScale();
  LogF("MaxTorqScale %f", tmax);
  Vector3 uair = Vector3(0.0, 0.0, -1.0);
  float speed = 1.0;
  for (int i = 0; i<n; i++) 
  {
    fscanf( fin, "%f %f %f\n", &_c_elevator, &_c_rudder, &speed );
    float err = RelaxPitch( 100, 0.2 / tmax, 0.1, 1e-4*tmax, uair, 1.0 );
    Vector3 force, torq;
    GetAeroForce(uair*speed, force, torq);
    const Vector3& fw = _orientation.Direction();
    //LogF("runRelaxedTest %i %f %f %f err(%f) F: %f %f %f  fw: %f %f %f ", i, _c_elevator, _c_rudder, speed, err, force.X(), force.Y(), force.Z(), fw.X(), fw.Y(), fw.Z() );
    fprintf(fout,"%i  %f %f %f  %f   %f %f %f   %f %f %f \n", i, _c_elevator,_c_rudder,speed,   err,   force.X(),force.Y(),force.Z(), fw.X(), fw.Y(), fw.Z() );
  }
  fclose(fout);
  fclose(fin );
  return 0;
}

void AeroTestPlatform::FindLDMax(float dAoA0, float convAcc, float& AoA, float& L, float& D) 
{
  _orientation = M3Identity;
  _angVelocity = VZero;
  _velocity = Vector3(0.0, 0.0, 1.0);
  Vector3 force, torq;
  float LD,oLD=-10000.0;
  //LogF("AeroTestPlatform::findLDMax");
  for (int i = 0; i < 1000; i++) 
  {
    _orientation.SetRotationX(-AoA); 
    GetAeroForce(VZero, force, torq);  
    LD = force[1] / (-force[2]);
    //LogF("rought %i AoA %f dAoA(%f) | LD %f L %f D %f ", i, AoA, dAoA0, LD, force[1], force[2] );
    if (LD < oLD) 
      break;
    AoA += dAoA0;
    oLD = LD;
  };
  Optimizer1D_Golden opt;
  AoA = opt.init(AoA-dAoA0, AoA);
  for (int i = 0; i < 100; i++) 
  {
    _orientation.SetRotationX(-AoA); GetAeroForce(VZero, force, torq);   LD = force[1] / (-force[2]);
    AoA = opt.step( -LD );
    float err = opt.errX();
    //LogF("fine %i AoA %f err %f maxErr(%f) | LD %f L %f D %f ", i, AoA, err, convAcc, LD, force[1], force[2]);
    if ( err < convAcc) 
      break;
  }
  L =  force[1];
  D = -force[2];
};

void AeroTestPlatform::SolveStraightFlight(int nmax, float convAcc, float Thrust, float Weight, float& speed, float& AoA) 
{
  // for stationary stright flight hold { Lift == Weight; Drag = Thrust }   
  //    Lift = CL*rho*S*v^2; Drag = CD*rho*S*v^2
  //    Lift/Drag = CL/CD = Weight/Thrust
  float LDtarget = Weight / Thrust; // target lift to drag ratio
  RootFinder rf;
  _orientation = M3Identity;
  _angVelocity = VZero;
  _velocity    = Vector3(0.0, 0.0, 1.0);
  float x1 = -0.2, x2 = +0.2;
  float y1, y2;
  Vector3 force, torq;

  // eval min bound
  _orientation.SetRotationX(-x1);
  GetAeroForce(VZero, force, torq);   
  y1 = (force[1] / -force[2]) - LDtarget; 
  
  // eval max bound
  _orientation.SetRotationX(-x2);
  GetAeroForce(VZero, force, torq);   
  y2 = (force[1] / -force[2]) - LDtarget;
  
  // init root finder
  rf.Init(x1,y1, x2,y2);
  //LogF("init RootFinder AoA (%f,%f) err (%f,%f)|(%f,%f) LDtarget %f ", x1,x2, y1,y2, y1+LDtarget,y2+LDtarget,  LDtarget );
  for (int i = 0; i < nmax; i++) 
  {
    // y=f(x)
    _orientation.SetRotationX(-x2); 
    GetAeroForce(VZero, force, torq);   
    y2 = force[1]/(-force[2])-LDtarget;

    // predict new 'x'
    x2 = rf.Step(x2, y2);
    //LogF("%i AoA %f err %f maxErr(%f) | LD %f L %f D %f ", i, x2, y2, convAcc, force[1]/(-force[2]), force[1], force[2]);
    if (fabs(y2) < convAcc)
      break;
  }
  AoA   = x2;
  speed = sqrt( Weight / force[1] );
};

void AeroTestPlatform::NewTrj(int n) 
{
  _nTrj = n;
  if (_trjPos)
    delete[] _trjPos;  
  _trjPos = new Vector3[_nTrj];
  if (_trjRot)
    delete[] _trjRot;  
  _trjRot = new Matrix3[_nTrj];
} 

AeroTestPlatform::~AeroTestPlatform()
{
  if (_trjPos)
    delete[] _trjPos;
  if (_trjRot)
    delete[] _trjRot;
  for (AeroPerformanceRecord* rec: records )
  {
    delete rec;
  };
};

void AeroTestPlatform::Fly(int n, float dt)
{
  for (int i = 0; i < n; ++i)
  {
    Vector3 force, torq;
    torq = VZero; force = VZero;
    GetAeroForce(VZero, force, torq);
    force += VUp * -G_CONST_;
    RigidBody::Move(dt, force, torq);
    _trjPos[i] = _position;
    _trjRot[i] = _orientation;
    //LogF("AeroTestPlatform::fly %i pos %f %f %f ", _position[0], _position[1], _position[2] );
  }
};

void AeroTestPlatform::FlyStraight(int n, int nsub, float dt_, Vector3 vGoal ) 
{
  float speedGoal = vGoal.Size();
  Vector3 uGoal = vGoal / speedGoal;
  FILE * fout = fopen( "flyStraight.log", "w");
  if (fout == nullptr) { LogF("AeroTestPlatform::flyStraight cannot open >>%s<< \n", "flyStraight.log" ); exit(-1); }
  float t = 0;
  float dt = dt_/nsub;
  bool detailLog = true;
  Vector3 force, torq, wForce;
  int ii = 0;

  //ControlFeedback pitchFeedback(0.0,dt,1.0);
  ControlFeedback pitchFeedback(0.0, dt*5, 1.0*5);
  float BreakMax          = 160; // [m^2]
  float mass = 1 / _invMass;
  float GF = -G_CONST_ * mass ;
  //float ThrustMax = -Power2ThrustRatio*GF; // [N]
  float ThrustMax = _compAeroModel->_engineThrust;
  LogF("AeroTestPlatform::flyStraight GF %f TrustMax %f mass %f \n", GF, ThrustMax, mass );

  for (int i = 0; i < n; i++) 
  {
    for(int j = 0; j< nsub; j++) // dynamics substep allows us to use very short time step => very stable and precise simulation
    {
      float speed = _velocity.Size();
      float cvup = VUp.DotProduct(_velocity) / speed;
      float cpitch = _orientation.Direction()[1];
      float croll  = _orientation.DirectionAside()[1];
      torq = VZero; force = VZero; wForce = VZero;
      GetAeroForce(VZero, wForce, torq);
      
      force += wForce;
      force += _orientation.Direction() * ( ThrustMax * _c_propel);    // Propeler
      force += _velocity * (_velocity.Size() * -BreakMax * _c_break);  // Break
      force += VUp * GF; // Gravity
      float dpitch = pitchFeedback.update(cvup);
      _c_elevator += dpitch; saturate(_c_elevator, -1.0, 1.0);
      RigidBody::Move(dt, force, torq);
      if (detailLog || (j == 0)) 
      {
        //LogF("AeroTestPlatform::flyStraigh %i cpitch %f croll %f elevator %f speed %f ", i, cpitch, croll, _c_elevator, speed);
        fprintf(fout, "%i %f   %f %f %f %f %f   %f %f   %f %f %f   %f %f %f  %f %f %f \n", ii, t, _c_elevator, cpitch, cvup, dpitch, force[1], _c_aileron, croll, _c_propel, _c_break, speed, _velocity[0], _velocity[1], _velocity[2], _position[0], _position[1], _position[2]);
      }
      ii++;
      t += dt;
    }
    //_orientation = M3Identity;
    _trjPos[i] = _position;
    _trjRot[i] = _orientation;
    LogF("AeroTestPlatform::flyStraight wForce %f %f %f |%f| ", wForce[0], wForce[1], wForce[2], wForce.Size() );
  }
  fclose(fout);
  LogF("AeroTestPlatform::flyStraigh DONE! " );
};

#if _DEBUG_AEROSURF
void AeroTestPlatform::DrawTrj( int i0, int n ) 
{
  if (_trjPos == nullptr) return;
  if(n>(_nTrj - i0)) n=_nTrj - i0;
  for (int ii = 0; ii < n; ii++) 
  {
    int i = i0 + ii;
    if( i>0 )  renderer.Add3dLine( _trjPos[i-1],_trjPos[i], PackedBlack );
    if (_trjRot) 
    {
      Vector3& p = _trjPos[i];
      renderer.Add3dLine(p, p + _trjRot[i].DirectionAside(), PackedRed);
      renderer.Add3dLine(p, p + _trjRot[i].DirectionUp(), PackedGreen);
      renderer.Add3dLine(p, p + _trjRot[i].Direction(), PackedBlue);
    }
  }
}
#endif //_DEBUG_AEROSURF
