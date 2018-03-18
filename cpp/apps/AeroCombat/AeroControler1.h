
#ifndef AeroControler1_h
#define AeroControler1_h

#include "AeroCraft.h"
#include "DynamicControl.h"
#include "AeroCraftWarrior.h"

class AeroControler1: public AnyControler{ public:

    AeroCraft* craft_bak=0;

    bool bRoll = true;
    Vec3d goalRoll;
    DynamicControl rollControl;

    void controlRoll( AeroCraftWarrior* craft, const Mat3d& rot, double dt ){
        //double roll = goalRoll.angleInPlane( rot.a*-1.0, rot.b );
        double roll = goalRoll.angleInPlane( rot.b, rot.a );
        rollControl.x_O1( roll, dt );
        craft->leftAirelon ->lrot = craft_bak->leftAirelon ->lrot;
        craft->rightAirelon->lrot = craft_bak->rightAirelon->lrot;
        craft->leftAirelon ->lrot.rotate( rollControl.x,craft->leftAirelon ->lrot.a);
        craft->rightAirelon->lrot.rotate(-rollControl.x,craft->rightAirelon->lrot.a);
    };

    virtual void control(void* obj, double dt ){
        AeroCraftWarrior* craft = (AeroCraftWarrior*)obj;
        Mat3d rot; rot.setT(craft->rotMat);
        if(bRoll)controlRoll( craft, rot, dt );
    };

};

#endif
