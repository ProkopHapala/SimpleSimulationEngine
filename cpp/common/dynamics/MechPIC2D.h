
#ifndef MechPIC2D_h
#define MechPIC2D_h

#include "macroUtils.h"

#include "fastmath.h"
#include "Vec2.h"
//#include "Vec3.h"

struct CompressibleMaterial{
    double molarMass    = 1.0;
    double invMolarMass = 1.0;

    void setMolarMass(double M){ molarMass=M; invMolarMass=1/M; }
};



const double const_Rgas    = 8.31446261815324; // J/(K*mol)
//const double const_Vsphere = M_PI*4./3.;

class MechPIC2D{ public:

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

    double* moles  =0; // amount of material in cell
    double* Umol   =0; // internal energy per mol of material due to pressure pV=nRT ( pressure normalized per mol off gas )

    // --- Particle Properties
    int np    = 0;
    int npMax = 0;
    int*     imats  = 0;   // materials
    //double* mass  =0;
    double*  pmoles = 0;   //
    Vec2d*   pos    = 0;   // position (global, or local within the box)
    Vec2d*   vel    = 0;   // velocity

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

        nc    = nc_;
        nctot = nc.x*nc.y;
        _realloc( moles, nctot  );
        _realloc( Umol,  nctot  );
    }

    inline void setStep(double step_){
        step=step_; invStep=1./step;
    }

    inline int point2cell(const Vec2d& p, Vec2d& d){
        int ix,iy;
        d.x=p.x*invStep;  ix=(int)d.x; d.x-=ix;
        d.y=p.y*invStep;  iy=(int)d.y; d.y-=iy;
        if(ix<0||ix>nc.x||iy<0||iy>nc.y) return 0;
        return nc.x*iy + ix;
    }

    void clearPressure(){
        for(int i=0; i<nctot; i++){
            moles[i]=0;
        }
    }

    void evalCellPressures(){
        // Total Energy E = T + U   (Kinetic + potential ... we may ensure condensation using certain materials)
        double invV = invStep*invStep;
        double RT = const_Rgas * 1e+4; // [K]  - reference
        for(int i=0; i<nctot; i++){
            double rho = moles[i] * invStep;
            Umol[i] =  RT * pow( rho, 5./3. ); // ToDo : check this equations of state ... they are probably wrong
        }
    }

    void particlesToCells(){
        // evaluated pressure in cells due to presence of particles (chunks of material)
        const int io10=nc.x  ;
        const int io11=nc.x+1;
        for(int i=0; i<np; i++){
            //double nmol = ( mass[i] / materials[ imats[i] ].molarMass );
            double nmol = pmoles[ i ];
            // ToDo : what if we want to implement different EOS for each material ?
            Vec2d d;
            int ic = point2cell( pos[i], d );
            if(ic<0) continue;
            double mx = 1-d.x;
            double my = 1-d.y;
            moles[ic      ] += nmol*mx *my;
            moles[ic+1    ] += nmol*d.x*my;
            moles[ic+io10 ] += nmol*d.y*mx;
            moles[ic+io11 ] += nmol*d.x*d.y;
        }
    }

    void moveMD( double dt ){
        // accelerate particles
        const int io10=nc.x  ;
        const int io11=nc.x+1;
        for(int i=0;i<np; i++){
            //double nmol = ( mass[i] / materials[ imats[i] ].molarMass ); // presure norm
            Vec2d d,f;
            int ic = point2cell( pos[i], d );
            if(ic<0) continue;
            double mx = 1-d.x;
            double my = 1-d.y;
            double e00 = Umol[ic     ];
            double e01 = Umol[ic+1   ];
            double e10 = Umol[ic+io10];
            double e11 = Umol[ic+io11];
            f.x = ( (e01-e00)*my + (e11-e10)*d.y )*-invStep;
            f.y = ( (e10-e00)*mx + (e11-e01)*d.x )*-invStep;
            double invM = materials[ imats[i] ].invMolarMass;

            Draw2D::drawVecInPos_d( f*1e-7, pos[i]*invStep );

            printf(  "[%i] f(%g,%g) v(%g,%g) | e(%g,%g,%g,%g)\n", i, f.x,f.y, vel[i].x,vel[i].y, e00,e01,e10,e11 );
            vel[i].add_mul( f,      dt*invM );

            pos[i].add_mul( vel[i], dt      );
        }
    }

    void update( double dt ){
        clearPressure();
        particlesToCells();
        evalCellPressures();
        moveMD( dt );
    }

    // ToDo : Merge & Split Particles
    // - During compression, many particles of the same meterial tend to be pressed into same box => we should merge them
    // - During expansion, single particle tend to expand over several cells => we should split them

};

#endif

