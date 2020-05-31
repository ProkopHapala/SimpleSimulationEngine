
#ifndef  IntersectionCurve_h
#define  IntersectionCurve_h

// Maybe rather call it "Implicit Curve Tracer" later
//
//  It traces 1D implicity function (a curve) emerging from intersection of several implicit surfaces defined as
//     F1(x,y,z) = 0    &    F2(x,y,z) = 0    &    F3(x,y,z) = 0    ...
//

//#include "Vec2.h"
#include "Vec3.h"


typedef  double (*FieldFunc)(const Vec3d& p, Vec3d& dFdp);

class Opt3d{ public:
// "Fast Inertial Realxation Engine" according to
// Bitzek, E., Koskinen, P., Gähler, F., Moseler, M. & Gumbsch, P. Structural relaxation made simple. Phys. Rev. Lett. 97, 170201 (2006).
// Eidel, B., Stukowski, A. & Schröder, J. Energy-Minimization in Atomic-to-Continuum Scale-Bridging Methods. Pamm 11, 509–510 (2011).
// http://users.jyu.fi/~pekkosk/resources/pdf/FIRE.pdf

    // parameters
    double fTinc = 1.1;   // factor by which time step is increased if going downhill
    double fTdec = 0.5;   // factor by which timestep is decreased if going uphill
    double fDamp = 0.99;  // rate of decrease of damping when going downhill
    //
    double dtmax   = 0.2;   // maximal timestep
    double dtmin   = 0.01;
    double dampMax = 0.1;   // default damping

    // variables
    double dt     = dtmax;    // time-step ( variable
    double damp   = dampMax;  // damping  ( variable

    inline void setup( double dtmax_, double dtmin_, double dampMax_){
        dtmax   = dtmax_;
        dtmin   = dtmin_;
        dampMax = dampMax_;
        dt      = dtmax;
        damp    = dampMax;
    }

    inline void  moveMD( const Vec3d& f, Vec3d& p, Vec3d& v ){
        v.mul( 1 - damp  );
        v.add_mul( f, dt );
        p.add_mul( v, dt );
    }

    // relaxation step using FIRE algorithm
    inline void move( const Vec3d& f, Vec3d& p, Vec3d& v ){
        double ff = f.norm2();
        double vv = v.norm2();
        double vf = f.dot(v);
        if( vf < 0 ){ // if velocity along direction of force
            v.set( 0.0d );
            dt   = fmax( dt * fTdec, dtmin );
            damp = dampMax;
        }else{       // if velocity against direction of force
            double cf  =     damp * sqrt(vv/ff);
            double cv  = 1 - damp;
            v.mul    ( cv    );
            v.add_mul( f, cf );	// v = cV * v  + cF * F
            dt    = fmin( dt * fTinc, dtmax );
            damp  = damp     * fDamp;
        }
        // normal leap-frog times step
        v.add_mul( f , dt );
        p.add_mul( v , dt );
    }

};

class IntersectionCurve{ public:
    constexpr static const int nFuncMax=4;
    int nfunc   =0;
    FieldFunc fields[nFuncMax]; // functions F(p) in formula (F(p) - C)^2 = 0
    double    Cs    [nFuncMax]; // constants C    in fomula  (F(p) - C)^2 = 0
    Opt3d opt;

    int nStepMax=0;
    int nStep=0;
    Vec3d* ps=0;

    //std::vector<Vec3d> ps;
    double fCurvPred=1.0;
    double Econv  = 1e-6*1000;
    double F2conv = 1e-6;
    int maxRelaxStep=100;
    int nevalTot=0;
    //Vec3d p;
    //Vec3d op;  // dimer direction
    //double l0=0.1; // distance between sampling points

    double relax(Vec3d& p, const Vec3d& op, double l0, bool bConstrain){
        // (F(p) - C)   =0
        // (F(p) - C)^2 =0
        // d_p (F(p) - C)^2 = 0
        //     (F(p) - C)*F(p) = 0
        Vec3d v =Vec3dZero;
        Vec3d hp; hp.set_sub( p, op ); hp.normalize();
        double E;
        double F2err;
        //glColor3b(1.0,0.0,0.0); Draw3D::drawLine( p, op ); //DEBUG
        //glColor3b(0.5,0.5,0.5);  //DEBUG
        //glBegin(GL_LINE_STRIP); //DEBUG
        for(int i=0; i<maxRelaxStep; i++){
            // -- constrain distance from previous point
            if(bConstrain){
                Vec3d dp;
                dp.set_sub( p, op );
                if( hp.dot(dp)<0 ){ dp.mul(-1.0); } // make sure we did not switched side accidently
                hp=dp;
                hp.normalize();
                p.set_add_mul(op,hp,l0);
            }
            // --- sum fields
            E=0;
            Vec3d f=Vec3dZero;
            for(int i=0;i<nfunc;i++){
                Vec3d  fi;
                double ei = fields[i]( p, fi ) - Cs[i];
                E += ei*ei;
                f.add_mul(fi,-2*ei);
            }
            nevalTot++;
            // --- termination
            F2err = f.norm();
            //printf( "relax[%i] p(%g,%g,%g) f(%g,%g,%g) E %g F2err %g \n", i, p.x,p.y,p.z,   f.x,f.y,f.z, E, F2err );
            //glVertex3f( p.x,p.y,p.z );
            if( (F2err<F2conv)&&(E<Econv) ){ break; }
            // --- constrain force and velocity
            if(bConstrain){
                f.makeOrthoU( hp );  // constrain (outproject force    along constrain)
                v.makeOrthoU( hp );  // constrain (outproject velocity along constrain)
            }
            // --- move
            opt.move( f, p, v );
        }
        //glEnd(); //DEBUG
        return E;
    }

    int trace( Vec3d p0, Vec3d dp, int nStepMax_ ){
        nevalTot=0;
        if(!ps){ nStepMax=nStepMax_; ps=new Vec3d[nStepMax]; };
        //ps.clear();
        Vec3d  op = p0;
        double l  = dp.norm();
        relax(op, op, l, false);
        p0=op;
        ps[0]=op;
        //return 0;
        Vec3d  odp= dp;
        Vec3d  p  = op+dp;
        double l2 = l*l;
        int i;
        for(i=1; i<nStepMax; i++){
            //printf( "DEBUG IntersectionCurve.trace[%i] \n", i );
            relax(p, op, l, true);
            //ps.push_back(p);
            //printf( "p[%i](%g,%g,%g) \n", i, p.x,p.y,p.z );
            ps[i]=p;
            Vec3d dtail; dtail.set_sub(p,p0);
            double l2tail = dtail.norm2();
            //printf( "l2tail %g \n", l2tail );
            if( (i>1) && ( l2tail<l2) ){ break; } // closed curve? .... catch on tail :-)
            Vec3d dp; dp.set_sub(p,op);
            op=p;
            //p.add(dp);                  // use previous step (without curvature predictor)
            //p.add_lincomb(2.,dp,-1.,odp); // with curvature predictor
            p.add_lincomb(1.+fCurvPred,dp,-fCurvPred,odp);
            odp=dp;
        }
        //DEBUG
        nStep=i;
        return nStep;
    }

};

#endif


