
#ifndef AeroCraftDesign_h
#define AeroCraftDesign_h

#include <vector>;

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "ScreenSDL2OGL_3D.h"
#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "AeroCraft.h"

#include "quaternion.h"
#include "geom3D.h"

#include "PotentialFlow.h"

class PotentialFlowSystem{ public:
    std::vector<Vec2i>  wings;
    std::vector<Vec3d>  vorts;  // vortex line control points
    std::vector<Rayt3d> panels; // (pos,normal)
    double * derr      = 0;     // derivatives of panel strenght
    double * strenghts = 0;     // vortex strenghts
    double * areas     = 0;     // panel areas
    double vortexDamp2 = 0.01;  // damping of vortex with distance along traling filament - necessary to avoid singularities

    Vec3d evalVelocity( const Vec3d& R, const Vec3d& flightDir ){
        Vec3d B = Vec3dZero;
        int n   = vorts.size();
        for(int i=0; i<n; i++){
            horseshoeDecay( B, R, vorts[i], vorts[i+1], flightDir, strenghts[i], vortexDamp2 );
        }
        return B;
    }

    Vec3d evalErrGrad( Vec3d vair ){
        Vec3d vhat = vair; vhat.normalize();
        for(int i=0; i<panels.size(); i++){
            Vec3d v = evalVelocity( panels[i].p0, vhat );  // get air velocity at control point
            v.add(vair);
            derr[i] = v.dot( panels[i].hdir );                       // project along surface normal
        }
    }


};

// ===============
// ====  Wings
// ===============

class WingSection{ public:
    double x;
    double y;
    double z;
    double chord;
    double twist;
    //std::vector<double> section; // TODO LATER .. better as pointer to spline
    //std::vector<double> polar;
    int profile;

    //Vec2d  vtwist;    // twist by complex number
    //double thickness; // ? in profile - affects also bending strength
    //double camber;    // ? in profile

    void fromString(char * str ){ sscanf( str, "%lf %lf %lf %lf %lf %li\n", &x, &y, &z, &twist, &chord, &profile ); }

};

class GunSlot{ public:
    Vec3d  pos;
    Vec3d  dir;
    double length; // GunType can have variable length, - easy to calculate params for different length
    // group to cylinder ?
    int    type;
    Vec2i  nammo; // box nx*ny of ammo next to gun ... box computed from number of rounds

    double mass;  // gun + mount + ammo

    void set( const Vec3d& pos_, const Vec3d& dir_, int type_, Vec2i nammo_ ){
        pos=pos_; dir=dir_; length=dir.normalize(); type=type_; nammo=nammo_;
    };

};

class Radiator{ public:
    // NOTE - functionally radiator is very similar to Jet engine (heats air)
    //   also radial engine has integrated radiator
    //   therefore perhaps we should have common class for radiator and engine
    //   what if we use surface of wing as radiator ?
};

class WingDesign{ public:
    Vec3d pos    = Vec3dZero;
    double pitch = 0.0;
    double roll  = 0.0;
    //double length;
    bool   symmetric = false;
    std::vector<WingSection> sections;

    double area = 0.0;
    double eInducedDrag = 0.0; // shape factor of induced drag

    // ==== functions

    void addSection( double dx, double dy, double z, double chord, double twist, int profile ){
        //printf( "dx %g z %g chord %g twist %g profile %i \n", dx, z, chord, twist, profile );
        if  ( sections.size()>0 ){
            WingSection& ows = sections[sections.size()-1];
               sections.push_back( { ows.x+dx, ows.y+dy, z, chord, twist, profile } ); }
        else { sections.push_back( {        0,        0, z, chord, twist, profile } ); }
    };

    double getArea(){
        double ox = sections[0].x;
        double oc = sections[0].chord;
        double A = 0.0;
        for( int i=1; i<sections.size(); i++ ){
            double x = sections[i].x;
            double c = sections[i].chord;
            A += (x-ox)*(c+oc); // trapezoides
            ox=x; oc=c;
        }
        if( !symmetric ) A*=0.5;
        area = A;
        return A;
    };


    void fromString(char * str ){
        int n;
        sscanf( str, "%i %lf %lf %lf \n",   n,  &pos.x, &pos.y, &pos.z, &pitch, &roll );
        sections.resize(n);
    }

};


// ===============
// ====  Fuselage
// ===============

struct FuselageSection{
    Vec3d  pos; //= Vec3dZero;
    Vec2d  sz ; //= Vec2dOnes;
    double ax2; //= 0;
    double ay2; //= 0;
    double ay3; //= 0;
    //std::vector<double> Rs; // TODO LATER

    void fromString(char * str ){ sscanf( str, "%lf %lf   %lf %lf   %lf %lf %lf \n", &pos.x, &pos.y, &sz.x, &sz.y, &ax2, &ay2, &ay3 ); }

    inline double getX( const Vec2d& csa ){ return sz.x*(csa.x *( 1 + csa.x*csa.x*( ax2 + ay3*csa.y ) ) ) + pos.x; }
    inline double getY( const Vec2d& csa ){ return sz.y*(csa.y *( 1 + csa.y*csa.y*( ay2             ) ) ) + pos.y; }
    inline void getXY( const Vec2d& csa, double& x, double& y ){
        //x=sz.x*(csa.x *( 1 + csa.x*  ax2 ) )               + pos.x;
        //y=sz.y*(csa.y *( 1 + csa.y*( ay2 + ay3*csa.y ) ) ) + pos.y;
        x=getX( csa );
        y=getY( csa );
    }

};


class FuselageDesign{ public:
    Vec3d  pos = Vec3dZero;
    //double length;
    std::vector<FuselageSection> sections;

    void fromString(char * str ){
        int n;
        sscanf( str, "%i %lf %lf %lf \n",   n,  &pos.x, &pos.y, &pos.z  );
        sections.resize(n);
    }

    void addSection( Vec3d pos, Vec2d sz, double ax2, double ay2, double ax3 ){
        sections.push_back( { pos, sz, ax2, ay2, ax3 } );
    };


};


// ===============
// ====  Fuselage
// ===============

//class AeroCraftGUI : public ScreenSDL2OGL_3D {
class AeroCraftDesign{ public:
    std::vector<FuselageDesign> fuselages;
    std::vector<WingDesign>     wings;
    std::vector<GunSlot>        guns;

    double wettedArea = 0.0;
    double mass       = 0.0;
    Mat3d  Ibody      = Mat3dZero;

    void toFlowSystem( PotentialFlowSystem& pf ){

        for(WingDesign& wd : wings ){
            int i0 = pf.vorts.size();
            pf.wings.push_back( {i0,wd.sections.size()} );
            Vec3d op1,op2,nr;
            for(WingSection& ws: wd.sections ){
                double ca = cos(ws.twist + wd.pitch);
                double sa = sin(ws.twist + wd.pitch);
                double z2 = ws.z-ws.chord;
                Vec3d p1  = {ws.x, ws.y + ws.z*sa, ws.z*ca},
                      p2  = {ws.x, ws.y + z2  *sa,   z2*ca};
                Vec3d pv  =  p1*0.75 + p2*0.25;
                Vec3d pp  = (p1 + p2 + op1 + op2)*0.25;
                nr.set_cross( p1-op2, p2-op1 );
                double area = nr.normalize();
                //double strenght = area * ;
                pf.vorts.push_back( pv );
                //strenghts[]; // determined from size of panel ?
                // should we store panel area ?
                pf.panels.push_back( {pp,nr,area} );
                op1=p1; op2=p2;
            }
        }

        for(FuselageDesign& fd : fuselages ){
            for(FuselageSection& fs : fd.sections ){
                // TODO
            }
        }

    }

};


// ===============
// ====  draw
// ===============

void draw_( WingDesign& wd ){
    glBegin(GL_TRIANGLE_STRIP);
        //double ox = 0,oca=1.0,osa=0.0;
        for( const WingSection& ws : wd.sections ){
            double ca = cos(ws.twist + wd.pitch);
            double sa = sin(ws.twist + wd.pitch);
            //glVertex3f( ox  , ws.chord*ca, ws.chord*sa ); glVertex3f( ox  , ws.chord*ca, ws.chord*sa  );
            double z2 = ws.z-ws.chord;
            //printf( "section  z %g c %g p1(%g,%g,%g) p2(%g,%g,%g) \n",  ws.z, ws.chord,   ws.x, ws.z*ca, ws.z*sa, ws.x, z2*ca, z2*sa );
            glVertex3f( ws.x, ws.y + ws.z*sa, ws.z*ca ); glVertex3f( ws.x, ws.y + z2*sa, z2*ca );
            //ox = ws.x;
        }
    glEnd();
}

void draw_( FuselageDesign& fd ){
    int n = 16;
    double dphi=2*M_PI/n;
    Vec2d drot; drot.fromAngle( dphi );
    for( int i=1; i<fd.sections.size(); i++){
        Vec2d rot = (Vec2d){1.0,0.0};
        glBegin(GL_TRIANGLE_STRIP);
            for( int j=0; j<=n; j++ ){
                Vec2d p;
                //p = fd.sections[i-1].getXY( rot ); glVertex3f( p.x, p.y, fd.sections[i-1].pos.z );
                //p = fd.sections[i  ].getXY( rot ); glVertex3f( p.x, p.y, fd.sections[i  ].pos.z );
                fd.sections[i-1].getXY( rot, p.x, p.y ); glVertex3f( p.x, p.y, fd.sections[i-1].pos.z );
                fd.sections[i  ].getXY( rot, p.x, p.y ); glVertex3f( p.x, p.y, fd.sections[i  ].pos.z );
                //ox = ws.x;
                rot.mul_cmplx(drot);
            }
        glEnd();
    }
}

void drawFlat( WingDesign& wd ){
    glBegin(GL_QUADS);
        //double ox = 0,oca=1.0,osa=0.0;
        Vec3d op1,op2,nr;
        for( int i=0; i<wd.sections.size(); i++ ){
            const WingSection& ws = wd.sections[i];
            double ca = cos(ws.twist + wd.pitch);
            double sa = sin(ws.twist + wd.pitch);
            //glVertex3f( ox  , ws.chord*ca, ws.chord*sa ); glVertex3f( ox  , ws.chord*ca, ws.chord*sa  );
            double z2 = ws.z-ws.chord;
            //printf( "section  z %g c %g p1(%g,%g,%g) p2(%g,%g,%g) \n",  ws.z, ws.chord,   ws.x, ws.z*ca, ws.z*sa, ws.x, z2*ca, z2*sa );
            Vec3d p1 = {ws.x, ws.y + ws.z*sa, ws.z*ca},
                  p2 = {ws.x, ws.y + z2  *sa,   z2*ca};
            if(i>0){
                nr.set_cross( p1-op2, p2-op1 ); nr.normalize();
                glNormal3f( nr.x, nr.y, nr.z );
                glVertex3f( p1.x, p1.y, p1.z );
                glVertex3f( p2.x, p2.y, p2.z );
                glVertex3f( op2.x, op2.y, op2.z );
                glVertex3f( op1.x, op1.y, op1.z );
            }
            op1=p1; op2=p2;
            //ox = ws.x;
        }
    glEnd();
}

void drawWire( WingDesign& wd ){
    //double ox = 0,oca=1.0,osa=0.0;
    Vec3d op1,op2;
    for( int i=0; i<wd.sections.size(); i++ ){
        const WingSection& ws = wd.sections[i];
        double ca = cos(ws.twist + wd.pitch);
        double sa = sin(ws.twist + wd.pitch);
        //glVertex3f( ox  , ws.chord*ca, ws.chord*sa ); glVertex3f( ox  , ws.chord*ca, ws.chord*sa  );
        double z2 = ws.z-ws.chord;
        //printf( "section  z %g c %g p1(%g,%g,%g) p2(%g,%g,%g) \n",  ws.z, ws.chord,   ws.x, ws.z*ca, ws.z*sa, ws.x, z2*ca, z2*sa );
        Vec3d p1 = {ws.x, ws.y + ws.z*sa, ws.z*ca},
              p2 = {ws.x, ws.y + z2  *sa,   z2*ca};
        if(i>0){
            glBegin(GL_LINE_LOOP);
            glVertex3f( p1.x, p1.y, p1.z );
            glVertex3f( p2.x, p2.y, p2.z );
            glVertex3f( op2.x, op2.y, op2.z );
            glVertex3f( op1.x, op1.y, op1.z );
            glEnd();
        }
        op1=p1; op2=p2;
        //ox = ws.x;
    }
}


void drawFlat( FuselageDesign& fd ){
    int n = 16;
    double dphi=2*M_PI/n;
    Vec2d drot; drot.fromAngle( dphi );
    for( int i=1; i<fd.sections.size(); i++){
        Vec2d rot = (Vec2d){1.0,0.0};
        glBegin(GL_QUADS);
            Vec3d op1,op2;
            for( int j=0; j<=n; j++ ){
                Vec3d nr;
                Vec3d p1,p2;
                fd.sections[i-1].getXY( rot, p1.x, p1.y ); p1.z = fd.sections[i-1].pos.z;
                fd.sections[i  ].getXY( rot, p2.x, p2.y ); p2.z = fd.sections[i  ].pos.z;
                if(j>0){
                //p1 = (Vec3d){ fd.sections[i-1].getXY( rot ),  };
                //p2 = (Vec3d){ fd.sections[i  ].getXY( rot ), fd.sections[i  ].pos.z };
                nr.set_cross( p1-op2, p2-op1 ); nr.normalize();
                glNormal3f( nr.x, nr.y, nr.z );
                glVertex3f( p1.x, p1.y, p1.z );
                glVertex3f( p2.x, p2.y, p2.z );
                glVertex3f( op2.x, op2.y, op2.z );
                glVertex3f( op1.x, op1.y, op1.z );
                }
                op1=p1; op2=p2;
                //ox = ws.x;
                rot.mul_cmplx(drot);
            }
        glEnd();
    }
}

void draw( GunSlot& gun ){
    Draw3D::drawCylinderStrip_wire( 8, 0.1, 0.05, (Vec3f)gun.pos, (Vec3f)( gun.pos+gun.dir*gun.length ) );
    double caliber    = 0.02;
    double ammoLength = 0.2;
    Draw3D::drawBBox( gun.pos, gun.pos + ((Vec3d){caliber*gun.nammo.x,caliber*gun.nammo.y,ammoLength}) );
}

void draw_( AeroCraftDesign& ad ){
    for( WingDesign wd : ad.wings ){
        glPushMatrix();
            //double ca = cos(wd->roll);
            //double sa = sin(wd->roll);
            glTranslatef( wd.pos.x, wd.pos.y, wd.pos.z );
            glRotatef(wd.roll*(180.0/M_PI), 0.0, 0.0, 1.0 );
            //printf( "wing p (%g,%g,%g) roll %g \n", wd.pos.x, wd.pos.y, wd.pos.z, wd.roll  );
            //draw_( wd );

            glDisable( GL_LIGHTING ); glColor3f(0.0,0.0,0.0); drawWire( wd );
            glEnable ( GL_LIGHTING ); glColor3f(1.0,1.0,1.0); drawFlat( wd );
            if(wd.symmetric){
                glScalef(-1.0,1.0,1.0);
                //draw_   ( wd );

                glDisable( GL_LIGHTING ); glColor3f(0.0,0.0,0.0); drawWire( wd );
                glEnable ( GL_LIGHTING ); glColor3f(1.0,1.0,1.0); drawFlat( wd );
            }
        glPopMatrix();
    }
    //glEnable( GL_LIGHTING );
    for( FuselageDesign fd : ad.fuselages ){
        glPushMatrix();
        glTranslatef( fd.pos.x, fd.pos.y, fd.pos.z );
        //glRotatef(fd->pitch, 0, 0, 1 );
        //draw_( fd );
        drawFlat( fd );
        //drawWire( fd );
        glPopMatrix();
    }
    glDisable(GL_LIGHTING);
    glColor3f(0.5,0.0,0.0);
    for( GunSlot& gun : ad.guns ){
        draw( gun );
    }
    //exit(0);
}


#endif  // AeroCraftDesign_h

