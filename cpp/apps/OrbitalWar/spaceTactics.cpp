
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "Draw3D.h"
#include "SDL_utils.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include "spline_hermite.h"
#include "SplineManager.h"

#include "SpaceBodies.h"
#include "SpaceWorld.h"
#include "SpaceDraw.h"
#include "raytrace.h"

#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "testUtils.h"

#include "appliedPhysics.h"
#include "spaceCombat.h"

#include "asteroidEngineering.h"



/*



ToDo:

Americium Weapons
    * Micro-Fission bombs
    * Micro-Fission engine / Rockets
Hafnium Weapons:      https://en.wikipedia.org/wiki/Hafnium_controversy
    * Hafnium rockets
    * Hafnium triggered fission weapons
    * Hafnium Lasers


#### Ideas:
 -

 - Optimal placement of units is given by
    - effective range of weapons and
    - narrow angle of shields

    *       *       *        *   //army 1
    ..
    . .
    .  .
    .   .
    .    .
    .     .
    .      .
    *       *       *        *   //army 1

*/

/*
void drawPolyLineRef( int n, Vec3d * ps, Vec3d * ps_ref ){   // closed=false
    //printf("%i %i\n", n, closed );
    glBegin(GL_LINE_STRIP);
    for(int i=0; i<n; i++){
        glVertex3d( ps[i].x-ps_ref[i].x, ps[i].y-ps_ref[i].y, ps[i].z-ps_ref[i].z );
    };
    glEnd();
};
*/

char strBuf[0x10000];

const int nThrusts = 6;
const double ThrustTimes[nThrusts  ] = {-10.0,35000.0,40000.0,60000.0,65000.0,100000.0};
const double Thrusts    [nThrusts*3] = {
    0.000, 0.0,0.0,
    0.001, 0.0,0.0,
    -0.200, 0.0,0.0,
    -0.100, 0.0,0.0,
    -0.001, 0.0,0.0,
    0.000, 0.0,0.0
};










void print_AblationRocket( double length, double pressure,  double caliber, double thickness, double dens, double molarMass, double temperature ){
    // (Force/area) / (mass/area)
    // s = 0.5 * a * t^2   ;   t = sqrt( 2*s/a )
    double area        = diskArea(0.5*caliber);
    double volume      = area*thickness;
    double mass        = volume*dens;

    double massPerArea = thickness*dens;
    double accel     = pressure / massPerArea;
    double t         = accelTime( accel, length );
    double vMuzzle   = accel * t;

    double vexh      = meanThermalVelocity( temperature, molarMass );
    //double areaMassFlow = pressure / vexh;

    double massRatio  = exp( vMuzzle/vexh );
    double mPropelant = (massRatio-1.) * mass;

    double Emuzzle    = 0.5*mass*vMuzzle*vMuzzle;
    double Ewaste     = 0.5*mPropelant*vexh*vexh;
    double Etot       = Emuzzle + Ewaste;
    double efficiency = Emuzzle/Etot;

    double power      = Etot/t;

    printf( "caliber %g [mm] thick %g [mm] dens %g [kg/l] \n", caliber*1000., thickness*1000.,  dens/1000. );
    printf( " -> area %g [m^2] volume %g [l] mass %g [kg] \n",  area, volume*1000., mass );
    printf( "areaMass %g [g/cm^2] pMax %g [GPa] length %g [m]\n", massPerArea/10., pressure*1e-9, length );
    printf( " -> accel %g [g] t %g [s] vMuzzle %g [km/s] \n", accel/9.81, t, vMuzzle*1e-3 );

    printf( " molarMass %g [g/mol] temperature %g [K] \n", molarMass, temperature );
    printf( " -> vexh %g [km/s] massRatio %g [1] EMuzzle %g [MJ] ETot %g [MJ] \n", vexh*1e-3, massRatio, Emuzzle*1e-6, Etot*1e-6 );
    printf( " -> efficiency %g [1] power %g [GW] \n", efficiency, power*1e-9 );

    const int nacc  = 5;
    const int ndist = 5;
    double dists [ndist]{1e+5, 1e+6, 1e+7, 1e+8, 1e+9}; // [m]
    double accels[nacc ]{0.01,0.1,1.0,10.0,100.0};      // [g]

    printf    ( "               ", dist );
    for(int j=0; j<nacc; j++){
        printf( "%3.3e ", accels[j] );
    }
    printf( "\n", dist );
    for(int i=0; i<ndist; i++){
        double dist = dists[i];
        double t_fly  = dist/vMuzzle + t;
        printf( "%3.3e[km] %3.3e[s] : ", dist/1000., t_fly );
        for(int j=0; j<nacc; j++){
            double spread = 0.5*accels[j]*sq(t_fly);
            //printf( "%3.3e|%3.3e ", spread, t_fly );
            //printf( "%3.3e ", spread );
            printf( "%3.3e[km] ", spread*1e-3 );
            //printf( "%3.3e ", t_fly );
        }
        printf( "\n", dist );
    }

}


void print_Laser( double wavelenght, double aperture, double distance, double power ){
    double r = difractionLimit_spot( wavelenght, aperture, distance );
    double intensity = power/(r*r);
    //double Efield = ;

    double time = distance/const_LightSpeed;
    printf( " wavelenght %g [um] aperture %g [m] distance %g [km] time %g [s] \n",  wavelenght*1e+6, aperture, distance*1e-3, time );
    printf( " r_spot %g [m] power %g [MW] -> Intensity %g [W/cm^2] \n", r, power*1e-6, intensity*1e-4  );
    printf(" ======================== \n ");
}




void print_WibleShield(double v0, double mass, double shieldMass, double standOff ){
    double tg        = kineticDispersion( v0, mass, shieldMass );
    double spread    = tg * standOff;
    double area      = M_PI * spread*spread;
    double E         = 0.5*mass*v0*v0;
    double intensity = E/area;
    printf( "v0 %g [km/s] mass %g [kg] shieldMass %g [kg] E %g [MJ] \n", v0*1e-3, mass, shieldMass, E*1e-6  );
    printf( "tg %g [1] spead %g [m] area %g [m^2] intensity %g [J/cm^2] \n", tg, spread, area, intensity*1e-4 );
}




class SpaceTactics : public AppSDL2OGL_3D { public:
    typedef AppSDL2OGL_3D Super;
    SplineManager splines;

    SpaceWorld world;

    bool dragging=false;
    int     nsamp=0;
    double tstart, tend, dt;


    std::vector<BodyInteraction> interactions;

    float scF = 1.0; float scSz = 0.2;
    //double view_scale = 1/1e+9;
    //Vec3d  view_shift = (Vec3d){0.0,0.0,0.0};


    //int iTrjMin=0,iTrjMax=0;
    //bool bRefTrj = false;
    //int trjStep = 10;
    SpaceBody* referenceBody = 0;

    int          fontTex;
    GUITextInput txtStatic;
    //GUITextInput txtPop;
    char curCaption[100];

    double timeCur = 0.0;
    double tmin    = 0;
    double tmax    = 100;

    //DEGUB
    //Vec3d *  dpos;
    //Vec3d * ddpos;

    int iShip=0;

    //void drawTrj    ( int n, Vec3d* ps );
    //void drawBodyTrj( SpaceBody& b );
    //void drawPlanet ( SpaceBody& b, int iTrj, double du);
    //void drawShip   ( SpaceBody& b, int iTrj, double du);
    //void drawInteraction( BodyInteraction& bi );
    //void drawInteractionTrj( BodyInteraction& bi );

	virtual void draw   ();
	virtual void drawHUD();
	virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	virtual void keyStateHandling( const Uint8 *keys );
    //virtual void mouseHandling( );

	SpaceTactics( int& id, int WIDTH_, int HEIGHT_ );

	void setTimeRange(double tmin_, double tmax_){ tmin=tmin_; tmax=tmax_; }

	inline double x2t(double x){ return x*(tmax-tmin)/WIDTH; };
	inline double t2x(double t){ return WIDTH*t/(tmax-tmin); };

};

SpaceTactics::SpaceTactics( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    //    "%i %lf %lf %lf", &n, &layerDens, &spacing, &critEdens
    whippleShieldType ws1; ws1.fromString("2 10.0 0.5 150000");
    //                                  &area, &HPs, &HPexponent
    ProjectedTarget tg1;   tg1.fromString( "10.0  1.5e+8 2.5"); tg1.wshield = &ws1;

    // sscanf( s,              "%lf %lf"  , &mass, &caliber );
    ProjectileType pt1;    pt1.fromString( "0.15 0.12" );
    //sscanf( s, "%lf %lf %lf %lf %lf"  , &length, &maxForce, &maxPower, &scatter, &fireRate  );
    SpaceGunType   sg1;    sg1.fromString( "800 60000 1e+9 2e-4 10" );
    SpaceGun g1( 1, &sg1, &pt1 );


    SpaceShipMobility mobil;
    //                  fuel_energy_density,  burnUp,  nozzle_efficiency,  pulse_fuel,  pulse_inert,  rate
    mobil.engine_hiIsp.fromString( "83140000e+6 0.15 0.8    0.150 0.5  0.5" );
    mobil.engine_loIsp.fromString( "83140000e+6 0.25 0.8    0.150 2.5  0.5" );
    mobil.mass_propletant = 1000.0*1e+3;
    mobil.mass_fuel       = 1000.0*1e+3;
    mobil.mass_empty      = 1000.0*1e+3;
    mobil.evalDeltaV( 1e+9, 0.0 );
    mobil.evalDeltaV( 1e+9, 1.0 );
    //sexit(0);



    DEBUG
    CombatAssembly battle;      DEBUG
    battle.addTarget( tg1 );    DEBUG
    battle.fireGun  ( g1, 1.0  );    DEBUG
    //      dist[m]  accel[m/s^2]
    battle.colide( 1e+5,    0.1 );
    //exit(0);



    //world.addPlanet( "Sun"  , 100000000000.0, 1.0, { 0.0, 0.0, 0.0}, {0.0,0.0,0.0} );
    //world.addPlanet( "Earth",   1000000000.0, 0.1, {10.0, 0.0, 0.0}, {0.0,1.0,0.0} );

    printf( "DEBUG 0 \n" );

    double vJup = 13.1e+3;
    double rJup = 7.784120E+011;

    world.addPlanet( "Sun"      ,   1.99E+030,  696.0000e+6,  { 0.0, 0.0, 0.0},          {0.0,0.0,0.0}    );
    world.addPlanet( "Jupiter"  ,   1.90E+27,     71.4920e+6, { 0.0,         rJup, 0.0}, {vJup,0.0  ,0.0} );
    world.addPlanet( "Io"       ,   8.93190E+23, 1.83000E+06, { 4.21700E+08, rJup, 0.0}, {vJup,17330,0.0} );
    world.addPlanet( "Europa"   ,   4.80000E+23, 1.56080E+06, { 6.71034E+08, rJup, 0.0}, {vJup,13740,0.0} );
    world.addPlanet( "Ganymede" ,   1.48190E+24, 2.63120E+06, { 1.07041E+09, rJup, 0.0}, {vJup,10880,0.0} );
    world.addPlanet( "Calisto"  ,   1.07590E+24, 2.41030E+06, { 1.88271E+09, rJup, 0.0}, {vJup,8204 ,0.0} );

    referenceBody = &world.planets[4];

    world.addShip( "ShipA1",   10.00E+3, 2.5E+03, (Vec3d){ -4e+6, 0, 0.0} +referenceBody->pos, (Vec3d){ 0.0, -5.6e+3,0.0}+referenceBody->vel );
    world.addShip( "ShipA2",   10.00E+3, 2.5E+03, (Vec3d){ -3e+6, 0, 0.0} +referenceBody->pos, (Vec3d){ 0.0, -5.6e+3,0.0}+referenceBody->vel );
    world.addShip( "ShipA3",   10.00E+3, 2.5E+03, (Vec3d){ -5e+6, 0, 0.0} +referenceBody->pos, (Vec3d){ 0.0, -5.6e+3,0.0}+referenceBody->vel );

    referenceBody = &world.planets[1];

    //world.addShip( "ShipB1",   10.00E+3, 2.5E+03, (Vec3d){ 3e+8, rJup, 0.0}, (Vec3d){ vJup, 26.6e+3,0.0} );
    //world.addShip( "ShipB2",   10.00E+3, 2.5E+03, (Vec3d){ 4e+8, rJup, 0.0}, (Vec3d){ vJup, 22.6e+3,0.0} );
    //world.addShip( "ShipB3",   10.00E+3, 2.5E+03, (Vec3d){ 5e+8, rJup, 0.0}, (Vec3d){ vJup, 20.6e+3,0.0} );

    world.addShip( "ShipB1",   10.00E+3, 2.5E+03, (Vec3d){ 0, rJup-4.5e+8, 0.0}, (Vec3d){ vJup+ 20.6e+3,0.0,0.0} );
    world.addShip( "ShipB2",   10.00E+3, 2.5E+03, (Vec3d){ 0, rJup-4.6e+8, 0.0}, (Vec3d){ vJup+ 20.6e+3,0.0,0.0} );
    world.addShip( "ShipB3",   10.00E+3, 2.5E+03, (Vec3d){ 0, rJup-4.7e+8, 0.0}, (Vec3d){ vJup+ 20.6e+3,0.0,0.0} );

    SpaceBody& jup = world.planets[1];
    world.intertialTansform( jup.pos*-1.0, jup.vel*-1.0 );
    SpaceDraw::bRefTrj = false;
    //referenceBody = &jup;
    world.allocateODE ();
    world.allocateTrjs( 500 );
    //world.allocateTrjs( 20 );
    interactions.push_back( {0,3, 0} );
    iShip = 3;

    double dt = 200.0;
    nonUni2spline( 0.0, dt, nThrusts, ThrustTimes, (Vec3d*)Thrusts, world.trj_n, world.ships[3].trjThrust );

    world.predictTrjs( 100000000, dt );

    setTimeRange(world.trj_t0, world.trj_t0+world.trj_dt*world.trj_n );
    //exit(0);


    SpaceDraw::trj_n = world.trj_n;

    //for( SpaceBody& b : world.ships ){
    //    for(int i=0; i<world.trj_n; i++){ printf("%s[%i] T= %f %f %f \n", b.name.c_str(), i, b.trjThrust[i].x, b.trjThrust[i].y, b.trjThrust[i].z ); }
    //}
    //exit(0);

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    //txtPop   .inputText = "txtPop";
    txtStatic.inputText = "txtStatic";

    /*
    double nsec=0.00015;
    for(int i=0; i<20; i++){
        printf( "nsec %e %e \n", nsec, Time::niceUnitAbove(nsec) );
        nsec*=10;
    }
    exit(0);
    */
    zoom = 5;





    double wavelenght = 1e-6;  // [m]
    double aperture   = 30;    // [m]
    double distance   = 1e+8;  // [m]
    double power      = 1e+11; // [W]
    print_Laser( wavelenght, aperture, distance, power );


    double v0         = 20e+3; // [m/s]
    double mass       = 1.0;   // [kg]
    double shieldMass = 0.2;   // [kg]
    double standOff   = 20.0;  // [m]
    print_WibleShield( v0, mass, shieldMass, standOff );

    //
    double length      = 100;   // [m]
    double maxPressure = 2e+9;  // [Pa]
    double caliber     = 0.2;   // [m]
    double thickness   = 0.01;  // [m]
    double dens        = 18000; // [kg/m^3]
    double molarMass   = 0.006; // [kg/mol]
    double temperature = 50000; // [K]
    print_AblationRocket( length, maxPressure, caliber, thickness, dens, molarMass, temperature );

    //exit(0);


    //                                               mass   caliber );
    world.projectileTypes.insert({"150g120mm",new ProjectileType("0.15   0.12")});
    //                                               length  maxForce maxPower scatter  fireRate  );
    world.gunTypes       .insert({"rail800m",new SpaceGunType  ("800     60000    1e+9     2e-4     10")});

    {
        //ProjectileType pt;
        //SpaceGunType   sg;
        SpaceCombatant  ship;
        SpaceGun        g;

        //               mass   caliber );
        //pt1.fromString( "0.15   0.12" );
        //               length  maxForce maxPower scatter  fireRate  );
        //sg1.fromString( "800     60000    1e+9     2e-4     10"       );
        //SpaceGun g1( 1, &sg1, &pt1 );
        //ship.guns.push_back( SpaceGun(1,world.gunTypes.get) );
        //ship.guns.push_back( SpaceGun(1,) );

        world.makeGun( g, 1, "rail800m", "150g120mm" );

        ship.faction = 1;
        for(int i=0; i<3; i++){
            ship.body = &world.ships[i];
            world.combatants.push_back( ship );
        }

    }

    SpaceCombatant target;
    target.faction = 2;
    target.targets.push_back( tg1 );
    double dmg = world.evalDamage( target, 10.0, 1, 0 );
    printf( "dmg %g \n", dmg );
    exit(0);

}

void SpaceTactics::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	//glEnable(GL_DEPTH_TEST);

    //bool bAnimate = false;

    /*
    bool bAnimate = true;
	if(bAnimate){
        int nPhases = 10;
        int nChunk = (world.trj_n/nPhases);
        iTrjMin = ( ((frameCount/10)%nPhases)*nChunk );
        //iTrjMax=_min(iTrjMin+nChunk,world.trj_n);
        iTrjMax=_min(iTrjMin+trjStep,world.trj_n);
	}else{
        iTrjMin=0; iTrjMax=world.trj_n;
	}
	*/

	world.trjTime2index(timeCur,SpaceDraw::iTrjMin); SpaceDraw::iTrjMax=SpaceDraw::iTrjMin+1;

	glColor3f(0.0f,0.0f,0.0f); for(SpaceBody& b : world.planets ){ SpaceDraw::bodyTrj( b ); SpaceDraw::planet( b, SpaceDraw::iTrjMin, 0, fontTex, zoom ); }
    glColor3f(0.0f,0.0f,1.0f); for(SpaceBody& b : world.ships   ){ SpaceDraw::bodyTrj( b ); SpaceDraw::ship  ( b, SpaceDraw::iTrjMin, 0 );  }

    glColor3f(0.0f,0.5f,0.0f); SpaceBody& b = world.planets[2];
    //printf("%s.trj \n", b.name.c_str() );
    //for(int i=0; i<world.trj_n; i++){ printf("%i %f %f %f \n ", i, b.trjPos[i].x, b.trjPos[i].y, b.trjPos[i].z ); };
    //exit(0);

    for( BodyInteraction& bi : interactions ){ SpaceDraw::interactionTrj( bi, world ); }

    for(SpaceBody& b : world.planets ){ b.getTrjPos(SpaceDraw::iTrjMin,0); };
    glColor3f(0.0f,0.0f,1.0f); for(SpaceBody& b : world.ships   ){ SpaceDraw::bodyTrj( b ); }

};

void SpaceTactics::drawHUD(){
    glDisable ( GL_LIGHTING );
    //glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );

    //int   nt   =(tmax-tmin)/Time::niceUnitBelow( (tmax-tmin)/40 );
    //float tstep=(tmax-tmin)/nt;
    float tstep=Time::niceUnitAbove( (tmax-tmin)/40.0 );
    int   nt = (int)((tmax-tmin)/tstep);
    //printf( "tstep %f nt %i   | %f %f \n", tstep, nt, tmax, tmin );

    //float xstep=WIDTH/(float)nt;
    //glPushMatrix();
    glColor3f(1.0f,1.0f,1.0f);
    glBegin(GL_LINES);
    float t,x;
    for( int i=0; i<nt; i++ ){
        t = i*tstep+tmin;
        x = t2x(t);
        glVertex3f(x,20,0); glVertex3f(x,0,0);
        for(int j=0; j<9; j++){ t+=0.1*tstep; x=t2x(t); glVertex3f(x,10,0); glVertex3f(x,0,0); }
    };
    glEnd();
    for( int i=0; i<nt; i++ ){
        t = i*tstep+tmin;
        x = t2x(t);
        Time::toStr( t, strBuf);
        Draw2D::drawText( strBuf, {x,10}, {80,15}, fontTex, 8 );
    }
    glColor3f(0.0,1.0,0.0);
    x=t2x(timeCur); Draw2D::drawLine({x,30},{x,0});
    Time::toStr( timeCur, strBuf);
    Draw2D::drawText( strBuf, {x,20}, {80,15}, fontTex, 8 );


    SpaceBody& thisShip = world.ships[iShip];
    if( thisShip.trjThrust ){
        //float scT = scF*100000;
        float scT = scF*100;
        //printf( "thisShip %i \n", thisShip.trjThrust );
        glBegin(GL_LINES);
        Vec3d oT;
        float ox;
        float T0=50.0;
        for(int i=0; i<world.trj_n; i++){
            double t = world.ind2time(i);
            double x = t2x(t);
            Vec3d T  = thisShip.getThrust(i,0);
            //printf( "T = (%f,%f,%f)\n", T.x,T.y,T.z);
            glColor3f(1.0,0.0,0.0); glVertex3f(ox,T0+oT.x*scT,0); glVertex3f(x,T0+T.x*scT,0);
            glColor3f(0.0,1.0,0.0); glVertex3f(ox,T0+oT.y*scT,0); glVertex3f(x,T0+T.y*scT,0);
            glColor3f(0.0,0.0,1.0); glVertex3f(ox,T0+oT.z*scT,0); glVertex3f(x,T0+T.z*scT,0);
            oT = T;
            ox=x;
        }
        glEnd();
    }


    //Draw::drawText( "abcdefghijklmnopqrstuvwxyz \n0123456789 \nABCDEFGHIJKLMNOPQRTSTUVWXYZ \nxvfgfgdfgdfgdfgdfgdfg", fontTex, 8, {10,5} );
    //glTranslatef( 10.0,HEIGHT-20,0.0  ); glColor3f(1.0,0.0,1.0); currentFormation->reportSetup (strBuf); Draw::drawText( strBuf, fontTex, 8, {80,15} );
    //glTranslatef( 300.0,      0.0,0.0 ); glColor3f(0.0,0.5,0.0); currentFormation->reportStatus(strBuf); Draw::drawText( strBuf, fontTex, 8, {80,15} );
    //glTranslatef( 300.0,      0.0,0.0 ); glColor3f(0.0,0.5,0.8); currentFormation->soldiers[0].type->toStrCaptioned(strBuf); Draw::drawText( strBuf, fontTex, 8, {80,15} );
    //glPopMatrix();


}

void SpaceTactics::mouseHandling( ){
    SDL_GetMouseState( &mouseX, &mouseY );   mouseY=HEIGHT-mouseY;
    mouse_begin_x = (2*mouseX-WIDTH )*zoom/HEIGHT;
    mouse_begin_y = (2*mouseY-HEIGHT)*zoom/HEIGHT;
    int mx,my; Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);
    //printf( " %i %i \n", mx,my );
    if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
        Quat4f q; q.fromTrackball( 0, 0, -mx*mouseRotSpeed, my*mouseRotSpeed );
        qCamera.qmul_T( q );
    }

    if ( buttons & SDL_BUTTON(SDL_BUTTON_LEFT)) {
        if( mouseY<20 ){ timeCur = x2t(mouseX); Time::toStr(timeCur,strBuf); printf( "%s\n", strBuf );  };
    }
    //Super::mouseHandling();
};

void SpaceTactics::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:  LMB = true;  break;
                case SDL_BUTTON_RIGHT: RMB = true;  break;
            };  break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:  LMB = false; break;
                case SDL_BUTTON_RIGHT: RMB = false; break;
            }; break;
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE:   quit(); break;
                //case SDLK_SPACE:    STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
                case SDLK_KP_MINUS: SpaceDraw::view_scale/=VIEW_ZOOM_STEP; break;
                case SDLK_KP_PLUS:  SpaceDraw::view_scale*=VIEW_ZOOM_STEP; break;
            } break;
        case SDL_QUIT: quit(); break;
    };
};

void SpaceTactics::keyStateHandling( const Uint8 *keys ){

    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -keyRotSpeed ); }

	double step = 0.05/SpaceDraw::view_scale;
    if( keys[ SDL_SCANCODE_W  ] ){ SpaceDraw::view_shift.y +=step; }
	if( keys[ SDL_SCANCODE_S  ] ){ SpaceDraw::view_shift.y -=step; }
	if( keys[ SDL_SCANCODE_A  ] ){ SpaceDraw::view_shift.x -=step; }
	if( keys[ SDL_SCANCODE_D  ] ){ SpaceDraw::view_shift.x +=step; }
	//printf(" view_shift %g %g %g \n", view_shift.x, view_shift.y, view_shift.z);
};

//void SpaceTactics::keyStateHandling( const Uint8 *keys ){ };

// ===================== MAIN

SpaceTactics * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
	thisApp = new SpaceTactics( junk , dm.w-150, dm.h-100 );
	//thisApp = new SpaceTactics( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















