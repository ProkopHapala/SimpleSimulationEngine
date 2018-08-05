
#ifndef AeroControler1_h
#define AeroControler1_h

#include "AeroCraft.h"
#include "DynamicControl.h"
#include "AeroCraftWarrior.h"

class AeroControler1: public AnyControler{ public:

    AeroCraft* craft    =0;
    AeroCraft* craft_bak=0;

    DynamicControl roll;
    DynamicControl pitch;
    DynamicControl yaw;

    bool bActive=true, bUp=false, bDir = false;
    Vec3d goalUp  = (Vec3d){0.0,1.0,0.0};;
    Vec3d goalDir = (Vec3d){0.0,0.0,1.0};

    // ========= Functions

    void updateDirs(){
        goalDir = craft->rotMat.c;
        goalUp  = craft->rotMat.b;
    }

    void setup( bool bUp_, bool bDir_, AeroCraft* craft_, AeroCraft *craft_bak_ ){
        craft     = craft_;
        craft_bak = craft_bak_;
        goalDir   = craft->rotMat.c;
        bUp       = bUp_;
        bDir      = bDir_;

        updateDirs();
    }

    void controlUp( const Mat3d& rot, double dt ){
        //double roll = goalRoll.angleInPlane( rot.a*-1.0, rot.b );
        double droll = goalUp.angleInPlane( rot.b, rot.a );
        roll.x_O1( droll, dt );
        craft->leftAirelon ->lrot = craft_bak->leftAirelon ->lrot;
        craft->rightAirelon->lrot = craft_bak->rightAirelon->lrot;
        craft->leftAirelon ->lrot.rotate( roll.x,craft->leftAirelon ->lrot.a);
        craft->rightAirelon->lrot.rotate(-roll.x,craft->rightAirelon->lrot.a);
    }

    void controlDirTail( const Mat3d& rot, double dt ){

        double c_yaw   = -goalDir.dot( rot.a );
        double c_pitch = -goalDir.dot( rot.b );

        pitch.x_O1( c_pitch, dt );
        yaw  .x_O1( c_yaw  , dt );

        craft->rudder->lrot = craft_bak->rudder->lrot;
        craft->rudder->lrot.rotate( yaw.x, craft->rudder->lrot.a);

        craft->elevator->lrot = craft_bak->elevator->lrot;
        craft->elevator->lrot.rotate( pitch.x, craft->elevator->lrot.a);
    }

    /*
    virtual void control(void* obj, double dt ){
        AeroCraftWarrior* craft = (AeroCraftWarrior*)obj;
        Mat3d rot; rot.setT(craft->rotMat);
        if(bRoll)controlRoll( craft, rot, dt );
    };
    */

    virtual void update( double dt ){
        if(bActive){
            //Mat3d rot; rot.setT(craft->rotMat);
            Mat3d rot; rot.set(craft->rotMat);
            if(bUp )controlUp     ( rot, dt );
            if(bDir)controlDirTail( rot, dt );
        }
    }

};

#endif
