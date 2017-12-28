
#ifndef AeroCraftDesign_h
#define AeroCraftDesign_h

#include <vector>;

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "ScreenSDL2OGL_3D.h"
#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "AeroCraft.h"

class WingSection{ public:
    double x;
    double twist;
    double chord;
    //std::vector<double> section; // TODO LATER .. better as pointer to spline
    //std::vector<double> polar;

    void fromString(char * str ){ sscanf( str, "%lf %lf %lf\n", &x, &twist, &chord); }

};

class WingDesign{ public:
    Vec3d pos;
    double pitch;
    double roll;
    //double length;
    bool   symmetric;
    std::vector<WingSection> sections;

    void fromString(char * str ){
        int n;
        sscanf( str, "%i %lf %lf %lf \n",   n,  &pos.x, &pos.y, &pos.z, &pitch, &roll );
        sections.resize(n);
    }

    void render(){
        double
        glBegin(GL_TRIANGLE_STRIP);
            //double ox = 0,oca=1.0,osa=0.0;
            for( WingSection ws : sections ){
                double ca = cos(ws.twist + pitch);
                double sa = sin(ws.twist + pitch);
                //glVertex3f( ox  , ws.chord*ca, ws.chord*sa ); glVertex3f( ox  , ws.chord*ca, ws.chord*sa  );
                glVertex3f( ws.x, ws.chord*ca, ws.chord*sa ); glVertex3f( ws.x, ws.chord*ca, ws.chord*sa );
                //ox = ws.x;
            }
        glEnd();
    }

};

class FuselageSection{ public:
    Vec2d  pos =(Vec2d){0.0,0.0};
    Vec2d  sz  =(Vec2d){1.0,1.0};
    double ax2=0;
    double ay2=0;
    double ay3=0;
    //std::vector<double> Rs; // TODO LATER

    void fromString(char * str ){ sscanf( str, "%lf %lf   %lf %lf   %lf %lf %lf \n", &pos.x, &pos.y, &sz.x, &sz.y, &ax2, &ay2, &ay3 ); }

    inline Vec2d getXY( const Vec2d& csa ){
        return {
            sz.x*(csa.x *( 1 + csa.x*  ax2 ) )              ,
            sz.y*(csa.y *( 1 + csa.y*( ay2 + ay3*csa.y ) ) ) + pos.y
        };
    }
};


class FuselageDesign{ public:
    Vec3d  pos;
    //double length;
    std::vector<FuselageSection> sections;

    void fromString(char * str ){
        int n;
        sscanf( str, "%i %lf %lf %lf \n",   n,  &pos.x, &pos.y, &pos.z  );
        sections.resize(n);
    }

    void render(){
        int n = 32;
        double dphi=(n-1);
        Vec2d drot; drot.fromAngle( dphi );
        for( int i=0; i<sections.size(); i++){
            Vec2d rot = (Vec2d){1.0,0.0};
            glBegin(GL_TRIANGLE_STRIP);
                for( int j=0; j<n; j++ ){
                    Vec2d p;
                    p = sections[i-1].getXY( rot ); glVertex3f( p.x, p.y, sections[i-1].pos.x );
                    p = sections[i  ].getXY( rot ); glVertex3f( p.x, p.y, sections[i  ].pos.x );
                    //ox = ws.x;
                }
            glEnd();
        }
    }
};


//class AeroCraftGUI : public ScreenSDL2OGL_3D {
class AeroCraftDesign{ public:
    std::vector<WingDesign*> wings;

    void render(){
        for( WingDesign* wd : wings ){
            glPushMatrix();
                //double ca = cos(wd->roll);
                //double sa = sin(wd->roll);
                glTranslatef( wd->pos.x, wd->pos.y, wd->pos.z );
                glRotatef(wd->roll, 0, 0, 1 );
                wd->render();
                if(wd->symmetric){
                    glScalef(-1.0,1.0,1.0);
                    wd->render();
                }
            glPopMatrix();
        }
    }

};

#endif  // AeroCraftDesign_h

