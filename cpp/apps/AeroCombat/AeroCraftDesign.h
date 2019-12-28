
#ifndef AeroCraftDesign_h
#define AeroCraftDesign_h

#include <vector>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "ScreenSDL2OGL_3D.h"
#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "AeroCraft.h"

#include "quaternion.h"
#include "geom3D.h"

#include "PotentialFlow.h"
#include "Lingebra.h"

class PotentialFlowSystem : public LinSolver{ public:

    // TODO: perhaps we need more control points (panels) per horse_shoe
    std::vector<Vec2i>    wings;
    std::vector<Vec3d>    vorts;   // vortex line control points
    std::vector<Rayt3d>   wpanels; // (pos,normal)
    std::vector<Sphere3d> sources;
    std::vector<Rayt3d>   panels; // (pos,normal)
    double *  errors    = 0;     // derivatives of panel strenght
    double *  strenghts = 0;     // vortex strenghts
    //double * dstrenghts = 0;     // vortex strenghts
    double * M = 0;
    //double * areas     = 0;   // panel areas now in panels[i].t;
    double vortexDamp2 = 0.01;  // damping of vortex with distance along traling filament - necessary to avoid singularities

    int nVorts  = 0;
    int nCoefs  = 0;
    double airSpeed;
    Vec3d  vair;
    Vec3d  hair = Vec3dZ*-1.0;
    int iter = 0;


    // ========== functions

    void setVair( const Vec3d& vair_ ){
        vair = vair_;
        hair = vair_*-1.0;
        airSpeed = hair.normalize();
    }

    void alloc(){
        nVorts = 0;
        for(Vec2i w : wings){ nVorts+=(w.y-1); }
        nCoefs = nVorts + sources.size()*3;

        printf( "PotentialFlowSystem::alloc: nPanels %i nVorts %i nCoefs %i \n", panels.size(), nVorts, nCoefs );
        _realloc( strenghts, nCoefs        );
        //_realloc(dstrenghts, nCoefs        );
        _realloc( errors,    panels.size() );
    }

    void allocCouplings(){ _realloc( M, panels.size()*nCoefs ); }

    void initialStrenghts(){
        int ip=0;
        //printf( "panels.size() %i \n", wings.size(), panels.size(), vorts.size() );
        for( const Vec2i& w : wings ){
            for(int ii=1; ii<w.y; ii++){
                //int iv=ii+w.x;
                //strenghts[ip] = -hair.dot( wpanels[ip].hdir ) * wpanels[ip].t;
                strenghts[ip] = 0;
                printf( "initial %i strenghts[%i] %g %g (%g,%g,%g)\n", ii, ip, strenghts[ip], wpanels[ip].t,  wpanels[ip].hdir.x, wpanels[ip].hdir.y, wpanels[ip].hdir.z );
                ip++;
            }
        }
        Vec3d* sp = (Vec3d*)(strenghts+ip);
        for(int i=0; i<sources.size(); i++){
            sp[i] = vair * sources[i].r * -0.1; // TODO: we don't realy know how to scale it
        }
    }

    void printVorts(){
        int ip=0;
        for(int iw=0; iw<wings.size(); iw++){
            const Vec2i& w = wings[iw];
            for(int ii=1; ii<w.y; ii++){
                int iv=ii+w.x;
                printf( "wing %i ii %i ip %i vr[%i] (%g,%g,%g) vr[%i] (%g,%g,%g) s %g\n",  iw, ii, ip, iv-1, iv, vorts[iv-1].x, vorts[iv-1].y, vorts[iv-1].z, vorts[iv].x, vorts[iv].y, vorts[iv].z, strenghts[ip] );
                ip++;
            }
        }
    }

    Vec3d evalVelocity( const Vec3d& R ) const {
        Vec3d B = Vec3dZero;
        int ip  = 0;
        for( const Vec2i& w : wings ){
            for(int ii=1; ii<w.y; ii++){
                int iv = w.x + ii;
                horseshoeDecay( B, R, vorts[iv-1], vorts[iv], hair, strenghts[ip], vortexDamp2 );
                //printf( "evalVelocity %i %i %g %g %g   %g\n", ii, i,  B.x, B.y, B.z, strenghts[i] );
                ip++;
            }
        }

        Vec3d* sp = (Vec3d*)(strenghts+ip);
        for(int i=0; i<sources.size(); i++){
            B.add( sourceDipol( R-sources[i].p, (Quat4d){sp[i].x,sp[i].y,sp[i].z, 0 } ) );
        }
        //exit(0);
        //printf( "%g %g %g \n", B.x, B.y, B.z );
        return B;
    }

    void evalErrors(){
        for(int i=0; i<panels.size(); i++){
            Vec3d v = evalVelocity( panels[i].p0 );  // get air velocity at control point
            v.add(vair);
            errors[i] = v.dot( panels[i].hdir );                       // project along surface normal
        }
    }

    void evalCouplings(){
        for(int j=0; j<panels.size(); j++){
            double* Mj = M + nCoefs*j;
            Vec3d p = panels[j].p0;
            for( const Vec2i& w : wings ){
                for(int ii=1; ii<w.y; ii++){
                    int iv = w.x + ii;
                    Vec3d B = Vec3dZero;
                    horseshoeDecay( B, p, vorts[iv-1], vorts[iv], hair, 1.0, vortexDamp2 );
                    *Mj = B.dot( panels[j].hdir );
                    /*
                    {
                        Vec3d v  = (vorts[iv-1]+vorts[iv])*0.5 - p;
                        Vec3d dp = v - p;
                        //printf( "Mji[%i,%i]: r=%g dp(%3.3e,%3.3e,%3.3e) v(%3.3e,%3.3e,%3.3e) p(%3.3e,%3.3e,%3.3e) \n", j, ip, dp.norm(), dp.x,dp.y,dp.z, v.x,v.y,v.z, p.x,p.y,p.z );
                        //printf( "Mji[%i,%i]: %g B(%5.5e,%5.5e,%5.5e) r %g \n", j, ip, *Mj, B.x,B.y,B.z, dp.norm() );
                        //printf( "Mji[%i,%i]: iv %i w.x %i \n", j, ip, iv, w.x );
                    }
                    */
                    Mj++;
                    //ip++;
                }
            }
            //Vec3d* sp = (Vec3d*)(Mj+ip);
            for(int i=0; i<sources.size(); i++){
                Vec3d dp,B;
                dp = p-sources[i].p;
                B = sourceDipol( dp, (Quat4d){1,0,0,0} ); *Mj = B.dot( panels[j].hdir ); Mj++;
                B = sourceDipol( dp, (Quat4d){0,1,0,0} ); *Mj = B.dot( panels[j].hdir ); Mj++;
                B = sourceDipol( dp, (Quat4d){0,0,1,0} ); *Mj = B.dot( panels[j].hdir ); Mj++;
            }
            //Vec3d v = evalVelocity( panels[i].p0 );  // get air velocity at control point
            //v.add(vair);
        }
    }

    double stepGD( double dt){
        double Etot=0,Ftot=0;
        //printf( "===== stepGD \n" );
        for(int j=0; j<panels.size(); j++){
            double* Mj = M+nCoefs*j;
            double err=0;
            for(int i=0; i<nCoefs; i++){
                //printf( "panel %i coef %i M %g \n", j, i, Mj[i] );
                //printf( "Mji[%i,%i]: %g s %g \n", j, i, Mj[i], strenghts[i]);
                err += Mj[i] * strenghts[i];
            }
            err += panels[j].hdir.dot( vair );

            err *= panels[j].t; // scale by panel area

            //printf( "err[%i] %g \n", j, err );
            Etot += err*err;
            errors[j] = err;
        }
        for(int i=0; i<nCoefs; i++){
            double* Mji  = M + i;
            double dEi = 0;
            for(int j=0; j<panels.size(); j++){
                dEi += (*Mji) * errors[j];
                Mji+=nCoefs; // FIXME : This is cache-unfriendly,   We should rather make Transpose(Mj)
            }
            Ftot += dEi*dEi;
            strenghts[i] += dEi * dt;
        }
        Ftot = sqrt(Ftot);
        printf( "iter %i E %g |F| %g \n", iter, Etot, Ftot );
        iter++;
        return Ftot;
       // exit(0);
    }

    virtual void dotFunc( int n, double * x, double * Ax ){
    };

};

void draw_( PotentialFlowSystem& fs, float psc, float vsc, bool bVCP ){
    int ip  = 0;

    for(Rayt3d& pr : fs.panels ){
        glColor3f(0.5,0.0,0.0); Draw3D::drawVecInPos( pr.hdir*pr.t*psc, pr.p0  );
        if( bVCP ){
            glColor3f(1.0,0.0,0.0);
            Vec3d v = fs.evalVelocity( pr.p0 );  // get air velocity at control point
            v.add( fs.vair );
            Draw3D::drawVecInPos( v*vsc, pr.p0 );
        }
    }
    for( Vec2i& iw : fs.wings){
        Vec3d op;
        int ip=0;
        for(int ii=0; ii<iw.y; ii++){
            int iv=ii+iw.x;
            Vec3d& p = fs.vorts[iv];
            glColor3f(0.0,0.5,0.0);
            Draw3D::drawVecInPos( fs.hair*-vsc, p );
            glColor3f(0.0,0.0,1.0);
            if(ii>0){
                Draw3D::drawLine( op, p );
                if( bVCP ){
                    Vec3d g; g.set_cross(fs.vair,p-op);
                    g.normalize();
                    Draw3D::drawVecInPos( g*fs.strenghts[ip], (p+op)*0.5 );
                }
                ip++;
            }
            op=p;
        }
    }
    glColor3f(0.5,0.5,0.0);
    for( Sphere3d& sph : fs.sources ){
        Draw3D::drawPointCross( sph.p, sph.r );
    }
}

void plotVecPlane( const PotentialFlowSystem& fs, Vec2i n, Vec3d p0, Vec3d a, Vec3d b, double sz, double dt ){
    glBegin(GL_LINES);
    for(int ia=0; ia<n.a; ia++ ){
        for(int ib=0; ib<n.b; ib++ ){
            Vec3d p = p0 + a*ia + b*ib;
            if( sz>0 ){
                glVertex3f( p.x-sz, p.y   , p.z    ); glVertex3f( p.x+sz, p.y,    p.z    );
                glVertex3f( p.x   , p.y-sz, p.z    ); glVertex3f( p.x,    p.y+sz, p.z    );
                glVertex3f( p.x   , p.y   , p.z-sz ); glVertex3f( p.x,    p.y,    p.z+sz );
            }
            glVertex3f( p.x, p.y, p.z );
            Vec3d v = fs.evalVelocity(p); //printf( "(%f,%f,%f) (%f,%f,%f)\n", p.x, p.y, p.z,  v.x, v.y, v.z );
            p.add_mul( v, dt);
            glVertex3f( p.x, p.y, p.z );
        }
    }
    glEnd();
}

void plotStreamLine( const PotentialFlowSystem& fs, int n, double dt, Vec3d p ){
    glBegin(GL_LINE_STRIP);
    for(int i=0; i<n; i++ ){
        Vec3d v = fs.evalVelocity( p );
        v.add(fs.vair);
        p.add_mul(v,dt);
        glVertex3f( p.x, p.y, p.z );
        //printf( "(%f,%f,%f) (%f,%f,%f)\n", p.x, p.y, p.z,  v.x, v.y, v.z );
    }
    glEnd();
    //exit(0);
}

void plotStreamLinePlane( const PotentialFlowSystem& fs, Vec2i n, int ns, Vec3d p0, Vec3d a, Vec3d b, double dt ){
    for(int ia=0; ia<n.a; ia++ ){
        for(int ib=0; ib<n.b; ib++ ){
            Vec3d p = p0 + a*ia + b*ib;
            //printf( "(%i,%i) (%f,%f,%f) \n", ia, ib, p.x, p.y, p.z );
            plotStreamLine( fs, ns, dt, p );
        }
    }
}

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


    void toFlowSystem( PotentialFlowSystem& pf, const WingSection& ws, const Vec3d& pos, const Mat3d& rot, Vec3d& op1, Vec3d& op2, bool& prev, const Vec2i& nsp ){
        double ca = cos(ws.twist);
        double sa = sin(ws.twist);
        double z2 = ws.z-ws.chord;
        Vec3d p1,p2;
        rot.dot_to_T( {ws.x, ws.y + ws.z*sa, ws.z*ca},p1);  p1.add(pos);
        rot.dot_to_T( {ws.x, ws.y + z2  *sa,   z2*ca},p2);  p2.add(pos);
        Vec3d pv  =  p1*0.75 + p2*0.25;
        //double strenght = area * ;
        //printf( "vorts.push[%i] (%g,%g,%g) \n", pf.vorts.size(), pv.x, pv.y, pv.z );
        pf.vorts.push_back( pv );
        //strenghts[]; // determined from size of panel ?
        // should we store panel area ?
        if( prev ){

            Vec3d nr; nr.set_cross( p1-op2, p2-op1 );
            double area = nr.normalize();
            //printf( "panels.push[%i] (%g,%g,%g) (%g,%g,%g) %g \n", pf.panels.size(), pp.x,pp.y,pp.z,   nr.x,nr.y,nr.z,  area  );
            pf.wpanels.push_back( { (p1+op2+ p2+op1)*0.5, nr, area } );
            double dx = 1.0/nsp.x;
            double dy = 1.0/nsp.y;
            for(int ix=0; ix<nsp.x; ix++ ){
                double fx = dx*(ix+0.5);
                double mx = 1-fx;
                for(int iy=0; iy<nsp.y; iy++ ){
                    double fy = dy*(iy+0.5);
                    //printf( " ixy (%i,%i) fxy (%g,%g) dxy (%g,%g) \n", ix,iy, fx,fy, dx, dy  );
                    Vec3d pp  = (op1*mx + p1*fx)*(1-fy) + (op2*mx + p2*fx)*fy;
                    // TODO : propper calculation of area ?


                    //Vec3d dx1 = ;
                    //Vec3d dy1 = ;
                    //Vec3d dx2 = ;
                    //Vec3d dy2 = ;


                    pf.panels.push_back( {pp,nr,area} );
                    //Rayt3d& pb=pf.panels.back(); printf( "panels.push[%i] (%g,%g,%g) (%g,%g,%g) %g \n", pf.panels.size()-1, pb.p0.x,pb.p0.y,pb.p0.z,   pb.hdir.x,pb.hdir.y,pb.hdir.z,  pb.t  );
                }
            }
        }else{
            prev=true;
        }
        op1=p1; op2=p2;
    }


    void toFlowSystem( PotentialFlowSystem& pf, Vec2i nsp ){
        for(WingDesign& wd : wings ){
            Mat3d rot = Mat3dIdentity;
            rot.rotate( wd.roll,  {0.0,0.0,1.0} );
            rot.rotate( wd.pitch, {1.0,0.0,0.0} );
            int ns = wd.sections.size();
            Vec3d op1,op2;
            if(wd.symmetric){
                rot.a.mul(-1);
                bool prev = false;
                int i0 = pf.vorts.size();
                for(int i=ns-1; i>=0; i-- ){ toFlowSystem( pf, wd.sections[i], wd.pos, rot, op1, op2, prev, nsp ); }
                pf.wings.push_back( {i0,pf.vorts.size()-i0} );
                rot.a.mul(-1);
            }
            bool prev = false;
            int i0 = pf.vorts.size();
            for(int i=0; i<ns; i++ )       { toFlowSystem( pf, wd.sections[i], wd.pos, rot, op1, op2, prev, nsp ); }
            pf.wings.push_back( {i0,pf.vorts.size()-i0} );
        }

        /*
        for(FuselageDesign& fd : fuselages ){
            for(int i=1; i<fd.sections.size(); i++ ){
                FuselageSection& fs1 = fd.sections[i-1];
                FuselageSection& fs2 = fd.sections[i  ];
                pf.sources.push_back( { fs1.pos*0.6666 + fs2.pos*0.3333, sqrt( (fs1.sz.x*0.6666+fs2.sz.x*0.3333)*(fs1.sz.y*0.6666+fs2.sz.x*0.3333)  ) } );
                pf.sources.push_back( { fs1.pos*0.3333 + fs2.pos*0.6666, sqrt( (fs1.sz.x*0.3333+fs2.sz.x*0.6666)*(fs1.sz.y*0.3333+fs2.sz.x*0.6666)  ) } );
            }
        }
        */

        printf( "toFlowSystem pf: nWings %i nPanels %i nVorts %i \n", pf.wings.size(), pf.panels.size(), pf.vorts.size() );

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


void drawWire( FuselageDesign& fd ){
    int n = 16;
    double dphi=2*M_PI/n;
    Vec2d drot; drot.fromAngle( dphi );
    for( int i=1; i<fd.sections.size(); i++){
        Vec2d rot = (Vec2d){1.0,0.0};
        Vec3d op1,op2;
        for( int j=0; j<=n; j++ ){
            Vec3d nr;
            Vec3d p1,p2;
            fd.sections[i-1].getXY( rot, p1.x, p1.y ); p1.z = fd.sections[i-1].pos.z;
            fd.sections[i  ].getXY( rot, p2.x, p2.y ); p2.z = fd.sections[i  ].pos.z;
            if(j>0){
                glBegin(GL_LINE_LOOP);
                glVertex3f( p1.x, p1.y, p1.z );
                glVertex3f( p2.x, p2.y, p2.z );
                glVertex3f( op2.x, op2.y, op2.z );
                glVertex3f( op1.x, op1.y, op1.z );
                glEnd();
            }
            op1=p1; op2=p2;
            //ox = ws.x;
            rot.mul_cmplx(drot);
        }
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
            //glEnable ( GL_LIGHTING ); glColor3f(1.0,1.0,1.0); drawFlat( wd );
            if(wd.symmetric){
                glScalef(-1.0,1.0,1.0);
                //draw_   ( wd );

                glDisable( GL_LIGHTING ); glColor3f(0.0,0.0,0.0); drawWire( wd );
                //glEnable ( GL_LIGHTING ); glColor3f(1.0,1.0,1.0); drawFlat( wd );
            }
        glPopMatrix();
    }
    //glEnable( GL_LIGHTING );

    for( FuselageDesign fd : ad.fuselages ){
        glPushMatrix();
        glTranslatef( fd.pos.x, fd.pos.y, fd.pos.z );
        //glRotatef(fd->pitch, 0, 0, 1 );
        //draw_( fd );
        //drawFlat( fd );
        drawWire( fd );
        glPopMatrix();
    }
    /*
    glDisable(GL_LIGHTING);
    glColor3f(0.5,0.0,0.0);
    for( GunSlot& gun : ad.guns ){
        draw( gun );
    }
    */
    //exit(0);
}


#endif  // AeroCraftDesign_h

