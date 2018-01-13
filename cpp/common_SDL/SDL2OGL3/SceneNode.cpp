
#include "SceneOGL3.h" // THE HEADER

// ============== per frame

bool visibilityCheck( const Vec3d& camPos, const Mat3d& camRot, const Vec2d& tgFrustrum_ ){
    //double r2 = radius * radius;
    // FIXME TODO
    return true;
}

void SceneNode3D::render( const Vec3d& camPos_, const Mat3d& camRot_, const Vec2d& tgFrustrum, double scaleTot_ ){

    bool visible = visibilityCheck( camPos_, camRot_, tgFrustrum );

    if ( visible ) {
        Mat3d camRot;     camRot.set_mmul( rot, camRot_ );
        Vec3d camPos;     camPos.set_sub ( camPos_, pos );   camPos = rot.dot( camPos );
        double scaleTot = scaleTot_ * scale;

        for( GLObject* obj : objects ){
            obj->draw( );
        }
        for( SceneNode3D* scn : subnodes ){
            scn->render( camPos, camRot, tgFrustrum, scale );
        }
    }

};


