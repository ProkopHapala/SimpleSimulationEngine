
#ifndef MechPIC2D_Temperature_h
#define MechPIC2D_Temperature_h

#include "macroUtils.h"

#include "fastmath.h"
#include "Vec2.h"
//#include "Vec3.h"

#include "MechPIC2D.h"

class MechPIC2D_T{ public:

    // Particle In Cell Kind-off simultaion
    // - Material and Intertia (momentum) is stored per particle
    // - Pressure is stored per cell
    // - particle has point-size, but is interolated between boxes (4-boxes in 2D)
    // -

    // --- Cell Properties
    double step    = 1.0;
    double invStep = 1.0;
    Vec2i nc;
    int nctot=0;

    double* moles     =0;   // amount of material in cell
    double* moles_new =0;
    //double* Umol   =0;    // internal energy per mol of material due to pressure pV=nRT ( pressure normalized per mol off gas )
    double* temperature =0; // Temperature in cell

    // --- Particle Properties
    int np    = 0;
    int npMax = 0;
    int*     imats  = 0;   // materials
    //double* mass  =0;
    double*  pmoles  = 0;   // molar amount in cells
    double*  pVols   = 0;   // volume of particle (how much spac in cells it currently use)
    double   pTs     = 0;   // particle temperature
    //double*  pcmoles = 0; // molar distribution over individual cells
    Vec2d*   pos     = 0;   // position (global, or local within the box)
    Vec2d*   vel     = 0;   // velocity

    //inline void point2cell(const Vec2d& p, Vec2d& dp, Vec2i ip){
    //    dp.x=p.x*invStep;  ip.x=(int)dp.x; dp.x-=ip.x;
    //    dp.y=p.y*invStep;  ip.y=(int)dp.y; dp.y-=ip.y;
    //}

    CompressibleMaterial* materials = 0;

    // ========== Functions

    void realloc( int np_, Vec2i nc_ ){
        np    = np_;
        npMax = np*2;
        _realloc( imats,  npMax );
        _realloc( pmoles, npMax );
        _realloc( pos,    npMax );
        _realloc( vel,    npMax );
        //_realloc( pcmoles,npMax*4 );

        nc    = nc_;
        nctot = nc.x*nc.y;
        _realloc( temperature, nctot );
        _realloc( moles      , nctot );
        _realloc( moles_new  , nctot );
        //_realloc( Umol       , nctot );
    }

    inline void setStep(double step_){
        step=step_; invStep=1./step;
    }

    inline int point2cell(const Vec2d& p, Vec2d& d){
        int ix,iy;
        d.x=p.x*invStep;  ix=(int)d.x; d.x-=ix;
        d.y=p.y*invStep;  iy=(int)d.y; d.y-=iy;
        if(ix<0||ix>nc.x||iy<0||iy>nc.y) return -1;
        return nc.x*iy + ix;
    }

    void updateCellThermodynamics(){
        // Total Energy E = T + U   (Kinetic + potential ... we may ensure condensation using certain materials)
        double invV = invStep*invStep;
        double RT = const_Rgas * 1e+4; // [K]  - reference
        for(int i=0; i<nctot; i++){
            double T   = temperature[i];
            double N   = moles      [i];
            double N_  = moles_new [i];
            double T_  = T * pow( N_/N, 2./3. );
            // To Do : We must update temperature based on
            temperature[i] = T_;  // store updated temperature
            moles      [i] = N_;  // store updated molar amount
            // ToDo : consider also heat tranfer !!!!
            //Umol[i] =  RT * pow( rho, 5./3. ); // ToDo : check this equations of state ... they are probably wrong
        }
    }

    void clearCells(){ for(int i=0; i<nctot; i++){moles_new[i]=0;} }

    void particlesToCells(){
        // evaluated pressure in cells due to presence of particles (chunks of material)
        clearCells();
        const int io10=nc.x  ;
        const int io11=nc.x+1;
        for(int i=0; i<nctot; i++){ moles_new[i]=0; }
        for(int i=0; i<np; i++){
            //double nmol = ( mass[i] / materials[ imats[i] ].molarMass );
            double pN = pmoles[ i ];
            // ToDo : what if we want to implement different EOS for each material ?
            Vec2d d;
            int ic = point2cell( pos[i], d );
            if(ic<0) continue;
            double mx = 1-d.x;
            double my = 1-d.y;
            moles_new[ic      ] += pN*mx *my;
            moles_new[ic+1    ] += pN*d.x*my;
            moles_new[ic+io10 ] += pN*d.y*mx;
            moles_new[ic+io11 ] += pN*d.x*d.y;
            // To Do : We should consider also temperature transfer => we need to count hot material which left or entered the cell
        }
    }

    void moveMD( double dt ){
        // accelerate particles
        const int io10=nc.x  ;
        const int io11=nc.x+1;

        double damp = 1.0-100000.0*dt;
        if(damp<0)damp=0;
        damp=1; // damping off
        //printf("damp %g \n", damp);

        for(int i=0;i<np; i++){
            //double nmol = ( mass[i] / materials[ imats[i] ].molarMass ); // presure norm
            Vec2d d,f;
            int ic = point2cell( pos[i], d );
            if(ic<0) continue;
            double mx = 1-d.x;
            double my = 1-d.y;

            // prevent double counting
            double pN  = pmoles[ i ];
            double N00 = pN*mx * my;
            double N01 = pN*d.x* my;
            double N10 = pN*d.y* mx;
            double N11 = pN*d.x*d.y;
            int jc;
            jc=ic     ; double e00 = temperature[jc]*(moles[jc]-N00);
            jc=ic+   1; double e01 = temperature[jc]*(moles[jc]-N01);
            jc=ic+io10; double e10 = temperature[jc]*(moles[jc]-N10);
            jc=ic+io11; double e11 = temperature[jc]*(moles[jc]-N11);

            f.x = ( (e01-e00)*my + (e11-e10)*d.y )*-invStep;
            f.y = ( (e10-e00)*mx + (e11-e01)*d.x )*-invStep;
            double invM = materials[ imats[i] ].invMolarMass;

            // ToDo - velocity is damped when it goes along gradient of own molar mass
            Draw2D::drawVecInPos_d( f*1e-7, pos[i]*invStep );

            //printf(  "[%i] f(%g,%g) v(%g,%g) | e(%g,%g,%g,%g)\n", i, f.x,f.y, vel[i].x,vel[i].y, e00,e01,e10,e11 );
            vel[i].mul(damp);
            vel[i].add_mul( f,      dt*invM );

            pos[i].add_mul( vel[i], dt      );
        }
    }

    void update( double dt ){
        //clearPressure();
        particlesToCells();
        updateCellThermodynamics();
        moveMD( dt );
    }

    // ToDo : Merge & Split Particles
    // - During compression, many particles of the same meterial tend to be pressed into same box => we should merge them
    // - During expansion, single particle tend to expand over several cells => we should split them

};

#endif

