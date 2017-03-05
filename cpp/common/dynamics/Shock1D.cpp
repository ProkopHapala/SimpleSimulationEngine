
#include "Shock1D.h" // THE HEADER

void ShockSystem1D::evalForce(){
    double volume0,area0,pressure0=0;
    boundProperties(0, volume0, area0);
    for(int i=0; i<nlayers; i++){
        double volume1,area1,pressure1;
        boundProperties( bounds[i], volume1, area1);
        //matData.evalEOS( 1/, E[i], pressure1 );
        pressure1 = cells[i].getPressure( (volume1-volume0) );
        bforce[i] = area0*(pressure1-pressure0); // force on inner boundary of this layer
        volume0=volume1; area0=area1;  pressure0=pressure1;
    }
    bforce[nlayers]= area0*(outer_pressure - pressure0);  // force on outer-most boundary
}

// problem is to move bonundary according to velocities of layer cog; this can be done acording to chain derivatives
// v_boundary = drb/dt;   v_cog = drc/dt; drc+ and drc- are cog of layer before and behind the boundary
//  drb/dt = (drb/drc+) * (drc+/dt)   +   (drb/drc-)*(drc-/dt) =  (drb/drc+)*vc+   +   (drb/drc-)*vc-   
void ShockSystem1D::move( double dt ){
    double bforce0=0,r0=0,dbound0=0;
    for(int i=0; i<nlayers; i++){
        double bforce1,r1;
        bforce1 = bforce[i];
        //double vel   = velocity[i] + (bforce0+bforce1)*imass[i]*dt;
        double vel     = velocity[i] + (bforce0+bforce1)*dt/cells[i].mass;
        velocity[i]    = vel;
        r1             = bounds[i];
        //double   pos   = getCOG( r0, r1 );
        //pos           += velocity*dt;
        double dr1,dr2;
        get_dR(r0,r1, dr1, dr2 );
        if(i>0) bounds[i-1] = r0 + dt*( vel*dr1 + dbound0 );
        dbound0     = vel * dr2;   // shift of outer bound due to cog movement 
        bforce0=bforce1;  r0=r1;
    }
    bounds[nlayers-1] += r0 + dbound0*dt;
}

void ShockSystem1D::update( double dt ){
    for( int i=0; i<nlayers; i++ ){
        cells[i].update(dt);
    }
    evalForce();
    move( dt );
}

void ShockSystem1D::allocate( int nlayers_ ){
    nlayers = nlayers_;
    if (bounds)   delete bounds;   bounds    = new double[nlayers+1];
    if (bforce)   delete bforce;   bforce    = new double[nlayers+1];
    //if (imass)    delete imass;    imass     = new double[nlayers];
    if (velocity) delete velocity; velocity  = new double[nlayers];
    if (cells)    delete cells;    cells     = new ShockVolume[nlayers]; 
}

void ShockSystem1D::init(){
    double volume0,area_;
    boundProperties(0, volume0, area_);
    for(int i=0; i<nlayers; i++){
        double volume1;
        boundProperties(bounds[1], volume1, area_);
        cells[i].init( volume1 - volume0 ); 
        //imass[i] = 1/cells[i].mass;
        volume0=volume1;
    }
}

