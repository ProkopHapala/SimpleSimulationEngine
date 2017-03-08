
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
};

class RBodyConfDyn{
	public:
	double RmutPos = 0.2;
	double RmutRot = 0.05;

	double RposMax    = 3.0; // maximum distance from which configurations feel each orther
	double RrotMax    = 0.1;
	//double RtotMax    =    ; // RtotMax**2 < Rrot**2 + Rpos**2
	//double RrotScale  = 1.0;

    double RposMax2;
	double RrotMax2;

	double FposScale;
	double Fpos0;
    double FrotScale;
	double Frot0;

	double convF;
	int    nMaxSteps = 100;

	double E_last = 1e+300;

	//DynamicOpt * optimizer=NULL;
	DynamicOpt  optimizer;
	Vec3d  *pos=NULL,*vpos=NULL,*fpos=NULL;
	Quat4d *rot=NULL,*vrot=NULL,*frot=NULL;

    // configuration database - in future we may use more sophisticated detastructure with fast neighborhood search
    int nConfsMax;
    int nConfs;
    RBodyPose * confs = NULL;
    double    * confE = NULL;

    Func_array2 objectiveFuncDerivs;  // external objective function for the relaxation

    // ============= Functions

    void precompAux(){
        RposMax2  = RposMax*RposMax;
        RrotMax2  = RrotMax*RrotMax;
        Fpos0     = FposScale/RposMax2;
        Frot0     = FrotScale/RrotMax2;
    }

    void storeThisConf(){
        // this should be improved in future - if configuration database is changed
        if(nConfs>=nConfsMax){ printf("error: unable to store configuration, database full !!! %i \n", nConfs ); return;}
        int i = nConfs;
        confs[i].pos=*pos;
        confs[i].rot=*rot;
        confE[i] = E_last;
    }

    void mutate(){
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

    void forceFromConf( const Vec3d& pos, const Quat4d rot, const Vec3d& pos0, const Quat4d rot0, Vec3d& fpos, Quat4d& frot ){
        // --- position force
        Vec3d dpos; dpos.set_sub(pos,pos0);
        double rp2 = dpos.norm2();
        if( rp2 > RposMax2 ) return;
        double fpos_scale = FposScale/rp2-Fpos0;   // we can put some 1D spline here in future
        fpos.add_mul( dpos, fpos_scale );

        // --- rotation force
        Quat4d dRdq;
        double rr2 = rot.ddist_cos( rot0, dRdq );
        if( rr2 > RrotMax2 ) return;
        double frot_scale = FrotScale/rr2-Frot0;    // we can put some 1D spline here in future
        frot.add_mul( dRdq, frot_scale );

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
        assembleConfForces();                                            // forces form configurations stored in memory
        E_last = objectiveFuncDerivs( 7, optimizer.pos, optimizer.force );         // forces from derivative of external objective function
        frot->sub_paralel_fast( *rot ); vrot->sub_paralel_fast( *rot );  // out-project component which would harm unitarity of quaternion
        return optimizer.optStep();
    }

    void findNewConf(){
        mutate();
        for( int i=0; i<nMaxSteps; i++ ){
            double f = optStep( );
            if( f < convF ) break;
        }
        double storeThisConfiguration( );
    }


    void init(){
        optimizer.allocate( 7 );
        pos   = (Vec3d* )(optimizer.pos  );
        rot   = (Quat4d*)(optimizer.pos+3);
        vpos  = (Vec3d* )(optimizer.vel  );
        vrot  = (Quat4d*)(optimizer.vel+3);
        fpos  = (Vec3d* )(optimizer.force  );
        frot  = (Quat4d*)(optimizer.force+3);
        pos->set(0.0);
        rot->setOne();
        precompAux();
        optimizer.initOpt( 0.01, 0.1 );
    }

};

#endif
