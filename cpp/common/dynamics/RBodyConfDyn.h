
#ifndef RBodyConfDyn_h
#define RBodyConfDyn_h

/*

Module for efficient exploration of configuration space of one R body with 6 DOFs (pos,rot)
  - the configuration udergo molecular dynamics relaxation under several force-fields:
    (1) froces form other sampled configurations ( in order to avoid visiting same configurations again )
    (2) derivatives of external objective function

*/

#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "DynamicOpt.h"

//typedef double (*FuncRBody2scalar)( const Vec3d& pos, const Quat4d& rot );
//typedef double (*FuncRBody2scalar)( const Vec3d& pos, const Quat4d& rot );

typedef double (*Func_array )(int n, double * arr1              );
typedef double (*Func_array2)(int n, double * arr1, double *arr2);

class RBodyPose{
    public:
    Vec3d   pos;
    Quat4d  rot;

    void print(){ printf("(%f %f %f)  (%f %f %f %f)\n", pos.x,pos.y,pos.z, rot.x,rot.y,rot.z,rot.w );}



};

class RBodyConfDyn{
	public:
	double RmutPos = 0.2;
	double RmutRot = 0.0;

	double RposMax    = 0.5; // maximum distance from which configurations feel each orther
	double RrotMax    = 0.2;
	//double RtotMax    =    ; // RtotMax**2 < Rrot**2 + Rpos**2
	//double RrotScale  = 1.0;
	double FposScale=30.0;
    double FrotScale=10.0;

    // axuliary
    double RposMax2;
	double RrotMax2;
	double RposeMax2;
    //double Fpos0;
	//double Frot0;

	double convF;
	int    nMaxSteps = 100;

	double E_last = 1e+300;

	//DynamicOpt * optimizer=NULL;
	DynamicOpt  optimizer;
	Vec3d  *pos=NULL,*vpos=NULL,*fpos=NULL;
	Quat4d *rot=NULL,*vrot=NULL,*frot=NULL;

    // configuration database - in future we may use more sophisticated detastructure with fast neighborhood search
    int nConfsMax=0;
    int nConfs=0;
    RBodyPose * confs = NULL;
    double    * confE = NULL;

    Func_array2 objectiveFuncDerivs = NULL;  // external objective function for the relaxation

    // ============= Functions

    void precompAux(){
        RposMax2  = RposMax*RposMax;
        RrotMax2  = RrotMax*RrotMax;
        RposeMax2 = RposMax2*0.2; // + RrotMax2

        FposScale=FposScale/RposMax2;
        FrotScale=FrotScale/RrotMax2;

        //Fpos0     = FposScale/RposMax2;
        //Frot0     = FrotScale/RrotMax2;
    }

    void mutate_near(){
        // there we do some initial mutation which should than relax
        // pos
        Vec3d  dpos;
        dpos.set(randf(-1,1),randf(-1,1),randf(-1,1));
        pos->add_mul( dpos, RmutPos/dpos.norm() );
        // rot
        dpos.set(randf(-1,1),randf(-1,1),randf(-1,1));
        rot->dRot_exact( RmutRot/dpos.norm(), dpos );
        //rot.dRot_taylor2 ( RmutRot/dpos.norm(), dpos ); // fast for small rotations
    };

    void mutate_far(){
        // there we do some initial mutation which should than relax
        // pos
        Vec3d  dpos;
        dpos.set(randf(-1,1),randf(-1,1),randf(-1,1));
        pos->add_mul( dpos, 1.0/dpos.norm() );
        // rot
        //dpos.set(randf(-1,1),randf(-1,1),randf(-1,1));
        //rot->dRot_exact( RmutRot/dpos.norm(), dpos );
        //rot.dRot_taylor2 ( RmutRot/dpos.norm(), dpos ); // fast for small rotations
    };

    void storeThisConf(){
        // this should be improved in future - if configuration database is changed

        for(int i=0; i<nConfs; i++){
            RBodyPose* conf_i = confs+i;
            double r2 = evalConfDist( *pos, *rot, conf_i->pos, conf_i->rot );
            if(r2 < RposeMax2){
                mutate_far();
                return;
            }
        }

        if(nConfs>=nConfsMax){ printf("error: unable to store configuration, database full !!! %i \n", nConfs ); return;}
        int i = nConfs;
        confs[i].pos=*pos;
        confs[i].rot=*rot;
        confE[i] = E_last;
        nConfs++;
    }

    double evalConfDist( const Vec3d& pos, const Quat4d rot, const Vec3d& pos0, const Quat4d rot0 ){
        Vec3d dpos; dpos.set_sub(pos,pos0);
        double r2 = dpos.norm2();
        //double rr2 += rot.ddist_cos( rot0, dRdq );
        return r2;
    }

    void forceFromConf( const Vec3d& pos, const Quat4d rot, const Vec3d& pos0, const Quat4d rot0, Vec3d& fpos, Quat4d& frot ){
        // NOTE :
        // Calculation of force between configurations is both (i) iefficient and (ii) useless
        // (i)  inefficient because hard and discontinuous constrain force will require fine time step
        // (ii) useless because we only need to know that a configuration hit an other to terminate the trijectory and generate new trial, repel it outward
        //      - generation of new trial can be done be extending in random direction by given length
        //  But dimer method is good to search for soft DOFs - this should be rather dynamics constrained on sphere
        constexpr double r2safety = 1e-8;
        // --- position force
        Vec3d dpos; dpos.set_sub(pos,pos0);
        double rp2 = dpos.norm2();
        if( rp2 > RposMax2 ) return;
        //double fpos_scale = FposScale/(rp2+r2safety)-Fpos0;   // we can put some 1D spline here in future
        double fpos_scale = FposScale*(RposMax2-rp2); // this makes problems because derivative discontinuity => unstable simulation
        //double fpos_scale = FposScale*sq(RposMax2-rp2); // to remove derivative discontinuity and make simulation smoother
        fpos.add_mul( dpos, fpos_scale );

        // --- rotation force
        /*
        Quat4d dRdq;
        double rr2 = rot.ddist_cos( rot0, dRdq );
        //printf("rp2 %f  rr2 %f \n", rp2, rr2 );
        //printf("dRdq (%g,%g,%g,%g)\n", dRdq.x, dRdq.y, dRdq.z, dRdq.w);
        if( rr2 > RrotMax2 ) return;
        double frot_scale = FrotScale/(rr2+r2safety)-Frot0;    // we can put some 1D spline here in future
        double frot_scale = FrotScale/(RrotMax2+rr2);
        frot.add_mul( dRdq, frot_scale );
        */
        //printf("FposScale %f  FrotScale %f \n", FposScale, FrotScale );
        //printf("fpos_scale %f  frot_scale %f \n", fpos_scale, frot_scale );
    }

    void assembleConfForces(){
        // in future we may use more sophisticated algorithm how to select just relevant nearest neighbor confs from database
        //Vec3d  fpos;
        //Quat4d frot;
        //Vec3d  *fpos =  (Vec3d* )(optimizer.pos  );
        //Quat4d *frot =  (Quat4d*)(optimizer.pos+3);
        for(int i=0; i<nConfs; i++){
            RBodyPose* conf_i = confs+i;
            forceFromConf( *pos, *rot, conf_i->pos, conf_i->rot, *fpos, *frot );
        }
        //*(Vec3d* )(optimizer.force  ) = fpos;
        //*(Quat4d*)(optimizer.force+3) = frot;
    }

    double optStep( ){
        optimizer.cleanForce( ); // set all forces to zero
        rot->normalize();        // keep quaternion normalized, otherwise unstable !!!
        if( objectiveFuncDerivs ){
            E_last = objectiveFuncDerivs( 7, optimizer.pos, optimizer.force );  // forces from derivative of external objective function
        }
        assembleConfForces();
        //printf("force"); ((RBodyPose*)optimizer.force)->print();              // forces form configurations stored in memory
        frot->sub_paralel_fast( *rot ); vrot->sub_paralel_fast( *rot );         // out-project component which would harm unitarity of quaternion
        //printf("force"); ((RBodyPose*)optimizer.force)->print();
        return optimizer.optStep();
        //return false;
    }

    void findNewConf(){
        mutate_near();
        for( int i=0; i<nMaxSteps; i++ ){
            double f = optStep( );
            if( f < convF ) break;
        }
        double storeThisConfiguration( );
    }


    void init( int nConfsMax_ ){
        optimizer.realloc( 7 );
        pos   = (Vec3d* )(optimizer.pos  );
        rot   = (Quat4d*)(optimizer.pos+3);
        vpos  = (Vec3d* )(optimizer.vel  );
        vrot  = (Quat4d*)(optimizer.vel+3);
        fpos  = (Vec3d* )(optimizer.force  );
        frot  = (Quat4d*)(optimizer.force+3);
        pos->set(0.0);
        rot->setOne();
        precompAux();
        optimizer.initOpt( 0.5, 0.1 );
        optimizer.setInvMass( 1.0 );
        nConfsMax=nConfsMax_;
        confs = new RBodyPose[nConfsMax];
        confE = new double[nConfsMax];
    }

};

#endif
