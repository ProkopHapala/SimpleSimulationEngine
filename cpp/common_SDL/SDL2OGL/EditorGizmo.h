
#ifndef  EditorGizmo_h
#define  EditorGizmo_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"
#include "Camera.h"

#include "Vec3.h"
#include "Vec2.h"
#include "raytrace.h"

//#include "Mesh.h"
// ToDo: ??? maybe use Mesh.h ???
// pickVertex     int pickVertex( const Vec3d &ray0, const Vec3d &hRay ){
// pickPolygon    double ray( const Vec3d &ray0, const Vec3d &hRay, Vec3d& normal, int& imin ){

//#include "Table.h"
//#include <string>

struct Pose3d{
    Vec3d pos;
    Mat3d rot;
};

class EditorGizmo{ public:

    //enum class MPick :char{ vertex='v', edge='e'               };
    //enum class MTrans:char{ move='m', rotate='r', scale='s';   };
    //enum class MOrig :char{ global='g', cluster='c', local='l' };

    int oglSphere=0;
    float Rhandle=0.1;  // size of handle
    float Rgizmo =1.0;  // size of gizmo
    //MPick  mPick;    // vertex, edge | polygon
    //MTrans mTransform;  // move, rotate, scale
    //MOrig  mOrigin;  // global, group
    char mPick ='v';   // vertex, edge | polygon
    char mTrans='m';   // move, rotate, scale
    char mOrig ='g';   // global, group
    bool bKeyboard   = true;
    bool bDragUpdate = true;
    bool bSelectAxis = true;
    bool dragged     = false;

    bool   axmask[3]; // ToDo: this could be bit-mask, but this is simpler
    Vec3d  axStart;   // start dragging on axes

    Pose3d pose = (Pose3d){Vec3dZero,Mat3dIdentity};
    Camera* cam;

    Vec3d rd,ro;
    Vec2f pix;

    int npoint;
    int nedge;
    Vec3d* points    =0;
    Vec3d* pointSizes=0;
    Vec2i* edges     =0; // edges between points
    //int*   groups    =0; // assign goup/cluster to each point
    // ToDo: Polygons ? editing polygons will take some time
    std::vector<Pose3d>         groupPose;
    std::unordered_map<int,int> selection; //   point_index -> group
    //std::unordered_set<int> selection;

    /*
    void applyTransform(){
        int nsel=selection.size();
        int sel     [nsel];
        int selgroup[nsel];
        int i=0; for( auto& s : selection ){ sel[i]=s.first; sel[i]=s.second; i++; }
        Vec3d shift;
        Vec3d center;
        Vec3d sc;
        Vec3d axis;
        double phi;
        switch(mTrans){
            case 'm':
                    Vec3d::move  ( nsel, sel, points, shift              );
                    // ToDo: consider different origin/orientation for different groups
                break;
            case 's': ;
                if(mOrig=='g'){
                    //Vec3d::scale ( n, selection, points, center, sc         );
                    pose.rot.scalePoints( nsel, sel, points, points, pose.pos, sc );
                }else{
                    for(int ii=0; ii<nsel; ii++){ int i=sel[ii]; int ig=selgroup[ii]; groupPose[ig].rot.scalePoint( points[i], points[i], groupPose[ig].pos, sc ); }
                }
                break;
            case 'r':
                Vec3d::rotate( nsel, sel, points, center, axis, phi );
                break;
        }
    }
    */

    void applyTranslation(const Vec3d& shift  ){
        int nsel=selection.size();
        int sel     [nsel];
        int selgroup[nsel];
        int i=0; for( auto& s : selection ){ sel[i]=s.first; sel[i]=s.second; i++; }
        //if(mOrig=='g'){
            Vec3d::move( nsel, sel, points, shift );
        //}else{
        //    for(int ii=0; ii<nsel; ii++){ int i=sel[ii]; int ig=selgroup[ii]; groupPose[ig].rot.scalePoint( points[i], points[i], groupPose[ig].pos, sc ); }
        //}
    };

    void applyScaling    (const Vec3d& sc     ){
        int nsel=selection.size();
        int sel     [nsel];
        int selgroup[nsel];
        int i=0; for( auto& s : selection ){ sel[i]=s.first; sel[i]=s.second; i++; }
        if(mOrig=='g'){
            //Vec3d::scale ( n, selection, points, center, sc         );
            pose.rot.scalePoints( nsel, sel, points, points, pose.pos, sc );
        }else{
            for(int ii=0; ii<nsel; ii++){ int i=sel[ii]; int ig=selgroup[ii]; groupPose[ig].rot.scalePoint( points[i], points[i], groupPose[ig].pos, sc ); }
        }
    };

    void applyRotation   (const Vec3d& uaxis, double phi ){
        int nsel=selection.size();
        int sel     [nsel];
        int selgroup[nsel];
        int i=0; for( auto& s : selection ){ sel[i]=s.first; sel[i]=s.second; i++; }
        // ToDo - not sure how to do rotation now
    };

    void clearAxMask(){ axmask[0]=0; axmask[1]=0; axmask[2]=0; }

    bool selectAxis(){
        double Rhandle2 = 0;
        double tmin;
        //ro, rd,
        for(int i=0; i<3; i++){
            Vec3d& ax     = pose.rot.vecs[i];
            Vec3d phandle = pose.pos+ax*Rgizmo;
            double t;
            double r2     = rayPointDistance2( ro, rd, phandle, t );
            if( (r2<Rhandle)&&(t<tmin)  ){
                clearAxMask();
                axmask[i]=true;
                tmin=t;
            }
        }
    }

    void projectToGizmo( Vec3d& cax ){
        Vec2f cxy= cam->ray2screen( (Vec3f)rd, (Vec3f)pose.pos ); // TodDo : shoudl we use double instead of float?
        // ToDo : we may perhaps simply use  pose.rot.dot( cxy )   and than just filter it
        const Mat3d camrot = (Mat3d)cam->rot;
        for(int i=0; i<3; i++){
            Vec3d& vax = pose.rot.vecs[i];
            if(axmask){ cax.array[i] = vax.dot(camrot.a)*cxy.a  +  vax.dot(camrot.b)*cxy.b; }
            else      { cax.array[i] = 0; }
        }
    }

    void mouseStart( ){
        if(bSelectAxis)selectAxis();
        Vec2f cxy = cam->ray2screen( (Vec3f)rd, (Vec3f)pose.pos ); // TodDo : shoudl we use double instead of float?
        projectToGizmo( axStart );
    }

    void mouseEnd( ){
        Vec3d axEnd;
        projectToGizmo( axEnd );
        switch(mTrans){
            case 'm':{
                Vec3d shift;
                pose.rot.dot_to_T( axEnd - axStart, shift );
                applyTranslation( shift );
                }break;
            case 's':{
                Vec3d sc; sc.set_div( axEnd, axStart );
                applyScaling( sc );
                }break;
            case 'r':{
                // -- ToDo: Not sure how to do rotations yet
                //void applyRotation(const Vec3d& uaxis, double phi );
                }break;
        }
    }

    void mousePick(){

    }

    //void onEvent( int mouseX, int mouseY, const SDL_Event& event ){
    //void onEvent( Vec3d mouseDir, const SDL_Event& event ){
    void onEvent( Vec2f pix_, const SDL_Event& event ){
        //printf("gizmo::onEvent \n");
        // ToDo : keyBinding !!!!!!
        //Vec3d mouseDir = mouseX + mouseY + cam.c;
        pix = pix_;
        Vec3f rd_,ro_;
        cam->pix2ray(pix,rd_,ro_);
        rd=(Vec3d)rd_; ro=(Vec3d)ro_;
        //printf(  "pix  (%g,%g) ro(%g,%g,%g) rd(%g,%g,%g) | zoom %g persc %i \n", pix.x,pix.y, ro.x,ro.y,ro.z,   rd.x,rd.y,rd.z, cam->zoom, cam->persp );

        switch( event.type ){
            case SDL_KEYDOWN:
                //printf("gizmo key(%i) | %i\n", event.key.keysym.sym, bKeyboard );
                if(bKeyboard){
                    switch( event.key.keysym.sym ){
                        // translation mode
                        case SDLK_t: mTrans='m';  break; // Translate
                        case SDLK_e: mTrans='s';  break; // Scale
                        case SDLK_r: mTrans='r';  break; // Rotate
                        // picking mode
                        case SDLK_i: mPick='e';  break; // edge
                        case SDLK_o: mPick='v';  break; // vertex
                        case SDLK_p: mPick='p';  break; // polygon
                        // Origin mode
                        case SDLK_g: mPick='g';  break; // global
                        case SDLK_h: mPick='c';  break; // cluster/group
                        case SDLK_j: mPick='l';  break; // individual/local
                    }
                }
                break;
            //case SDL_TEXTINPUT:
            //    break;
            //case SDL_MOUSEWHEEL:
            //    break;
            case SDL_MOUSEBUTTONDOWN:
                switch( event.button.button ){
                    case SDL_BUTTON_LEFT:
                        mouseStart();
                        dragged = 1;
                        break;
                    case SDL_BUTTON_RIGHT:
                        mousePick();
                        break;
                }; break;
            case SDL_MOUSEBUTTONUP:
                switch( event.button.button ){
                    case SDL_BUTTON_LEFT:
                        mouseEnd();
                        dragged = 0;
                        break;
                    case SDL_BUTTON_RIGHT:
                        break;
                }; break;
            case SDL_MOUSEMOTION:
                switch( event.button.button ){
                    case SDL_BUTTON_LEFT:
                        if(dragged && bDragUpdate){
                            mouseEnd();
                            //SDL_MouseMotionEvent* event_ = (SDL_MouseMotionEvent*)&event;
                            //dragged->moveBy( event_->xrel, -event_->yrel );
                        }
                        break;
                    case SDL_BUTTON_RIGHT:
                        break;
                }; break;
        }
        //return active;
    }

    void draw(){

        glColor3f(1.0,0.0,0.0); Draw3D::drawPointCross( ro   ,     0.5 );
        glColor3f(0.0,0.0,1.0); Draw3D::drawPointCross( ro+rd*0.0, 0.5 );

        Draw3D::drawMatInPos( pose.rot, pose.pos );
        //glCallList( oglSphere );
        //Draw3D::drawShape( pose.pos, pose.rot, shape );
        switch(mTrans){
            case 'm':
                glColor3f(1,0,0); Draw3D::drawConeFan( 4, Rhandle, pose.pos + pose.rot.a*Rgizmo, pose.pos+pose.rot.a*(Rgizmo+Rhandle) );
                glColor3f(0,1,0); Draw3D::drawConeFan( 4, Rhandle, pose.pos + pose.rot.b*Rgizmo, pose.pos+pose.rot.b*(Rgizmo+Rhandle) );
                glColor3f(0,0,1); Draw3D::drawConeFan( 4, Rhandle, pose.pos + pose.rot.c*Rgizmo, pose.pos+pose.rot.c*(Rgizmo+Rhandle) );
                break;
            case 's': ;
                glColor3f(1,0,0); Draw3D::drawConeFan( 4, Rhandle, pose.pos + pose.rot.a*Rgizmo, pose.pos+pose.rot.a*(Rgizmo-Rhandle) );
                glColor3f(0,1,0); Draw3D::drawConeFan( 4, Rhandle, pose.pos + pose.rot.b*Rgizmo, pose.pos+pose.rot.b*(Rgizmo-Rhandle) );
                glColor3f(0,0,1); Draw3D::drawConeFan( 4, Rhandle, pose.pos + pose.rot.c*Rgizmo, pose.pos+pose.rot.c*(Rgizmo-Rhandle) );
                break;
            case 'r': Draw3D::drawSphereOctLines ( 32, Rgizmo, pose.pos, pose.rot, true ); break;
        }
    }

};

#endif
