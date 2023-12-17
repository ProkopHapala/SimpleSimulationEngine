
#ifndef  SceneGraph_h
#define  SceneGraph_h

#include <string>
#include <vector>
#include <unordered_map>

#include <Draw3D.h>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"
#include "geom3D.h"
//#include "Camera.h"

class VarLink{
    // links two variables in scene
    void* from=0;
    void* to=0;
    int   nbytes=0;
    VarLink( void* from_, void* to_, int nbytes_ ):from(from_),to(to_),nbytes(nbytes_){};
    void apply(){
        for(int i=0;i<nbytes;i++){ ((char*)to)[i]=((char*)from)[i]; };
        // todo - we may make more sophisticated systems with makes arithmetic (add,sub,mul...) but for that we would need virtual functions
    }
};


namespace Scene{

    static int  nSphereRes =3;
    static bool bWire      =false;

// ============ SceneNodes

class Primitive{ public:
    virtual void render()=0;
	//virtual void (std::unordered_map<std::string,int>);
};

class Node: public Primitive { public:
    Vec3d pos   = Vec3dZero;
    Vec3d scale = Vec3dOne;
    Mat3d rot   = Mat3dIdentity;
    Primitive* obj=0; 
	// ==== functions
    Node( Primitive* obj_=0,Vec3d pos_=Vec3dZero, Mat3d rot_=Mat3dIdentity, Vec3d scale_=Vec3dOne):obj(obj_),pos(pos_),rot(rot_),scale(scale_){};
	virtual void render() override{ 
        //printf("DEBUG Scene::Node::render() \n");
        glPushMatrix();
        Draw3D::rigidTransform( pos, rot, scale );
        if(obj){ obj->render(); }else{ Draw3D::drawAxis(1.0); }
        glPopMatrix();
    };
};

class Group: public Node { public:
    std::vector<Primitive*> objs;
	// ==== functions
    Group(Vec3d pos_=Vec3dZero, Mat3d rot_=Mat3dIdentity, Vec3d scale_=Vec3dOne):Node(0,pos_,rot_,scale_){};
	virtual void render() override{
        //printf("DEBUG Scene::Group::render() \n");
        glPushMatrix();
        Draw3D::rigidTransform( pos, rot, scale );
        for( Primitive* o : objs ){
            o->render();
        }
        glPopMatrix();
    };
	//virtual void (std::unordered_map<std::string,int>);
};

// ============ Scene Primitives

class Line : public Primitive { public:
    Vec3d p0;
    Vec3d p1;
    virtual void render() override{
        //printf("DEBUG Scene::Line::render() \n");
        Draw3D::drawLine(p0,p1);
    }
    Line(Vec3d p0_,Vec3d p1_):p0(p0_),p1(p1_){};
};

class Sphere : public Primitive { public:
    double R;
    Vec3d  pos;
    Sphere(Vec3d pos_,double R_):pos(pos_),R(R_){};
    virtual void render() override{
        //printf("DEBUG Scene::Sphere::render() \n");
        Draw3D::drawSphere_oct( nSphereRes, R, pos, bWire );
    }
};

class Cone : public Primitive { public:
    int n;
    double R0,R1;
    Vec3d p0;
    Vec3d p1;
    Cone(int n_,Vec3d p0_,Vec3d p1_,double R0_,double R1_):n(n_),p0(p0_),p1(p1_),R0(R0_),R1(R1_){ }
    virtual void render() override{
        //printf("DEBUG Scene::Cone::render() \n");
        Draw3D::drawCone( n, 0, M_PI*2, R0, R1, p0, p1, false );
    }
};

class Capsula : public Primitive { public:
    int n;
    double R0,R1;
    double theta0,theta1;
    Vec3d p0;
    Vec3d p1;
    Capsula(int n_,Vec3d p0_,Vec3d p1_,double R0_,double R1_,double theta0_,double theta1_):n(n_),p0(p0_),p1(p1_),R0(R0_),R1(R1_),theta0(theta0_),theta1(theta1_){ }
    virtual void render() override{
        //printf("DEBUG Scene::Capsula::render() \n");
        Draw3D::drawCapsula( p0, p1, R0, R1, theta0, theta1, 0.2, n, true );
    }
};


// ============ Functions

void make_array( Group* g, Vec3i ns, Primitive* o, Vec3d pos0, Mat3d ds, Mat3d rot=Mat3dIdentity ){
    for(int ia=0;ia<ns.a; ia++){
        for(int ib=0;ib<ns.b; ib++){
            for(int ic=0;ic<ns.c; ic++){
                Vec3d p; ds.dot_to_T( {(double)ia,(double)ib,(double)ic}, p ); p.add(pos0);
                g->objs.push_back( new Node( o, p, rot ) );
            }
        }
    }
}

void make_circle( Group* g, int n, double phi, Primitive* o, Vec3d pos0, Vec3d shift, Vec3d ax, Mat3d rot=Mat3dIdentity){
    ax.normalize();
    double dphi = phi/n;
    double ca=cos(dphi);
    double sa=sin(dphi); 
    for(int i=0;i<n; i++){
        shift.rotate_csa( ca, sa, ax );
        rot  .rotate_csa( ca, sa, ax );
        g->objs.push_back( new Node( o, pos0+shift, rot ) );
    }
}

};

#endif
