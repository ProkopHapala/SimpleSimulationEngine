
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"
#include "geom3D.h"
//#include "Body.h"

#include "Mesh.h"
//#include "CMesh.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

/*

ToDo:

For Optimal implementation of GPU / OpenCL we should use univorm datastructure => Triangle-surfaces
 * Each Triangle is a boundary between two material - Right-side / Left-side. Hiting the triangle means crossing this boundary
 * Triangles can be stored in acceleration Grid or Tree


*/


/*
enum class TRayObj{ sphere, cylinder, Tetrahedron, box, mesh };

struct RayedObject{
    int   kind;
    void* obj;
    int   imaterial;

    double ray( const Vec3d &ray0, const Vec3d &hRay ){
        // vtable vs switch :  https://stackoverflow.com/questions/4467485/vftable-performance-penalty-vs-switch-statement
        int   imin;
        Vec3d normal;
        switch(kind){
            case TRayObj::sphere :
                return ((Sphere3d*)obj)->ray(ray0,hRay);
                //inline double raySphere( const Vec3d& ray0, const Vec3d& hRay, double R, const Vec3d& center );
                break;
            case TRayObj::mesh :
                return ((Mesh*)obj)->ray( ray0, hRay, normal, imin );
                break;
            default:
        }
    }

}
*/

constexpr static const int nScatAng=16;

struct ScatterMatrial{
    //double scatterProb;     // [1/m  ]  scattering probability
    double mfp;             // [m]  mean free path
    double Sang[nScatAng];  // [1/rad]  angular scattering distribution function ( cumulative distribution, i.e. integral of density )

    inline double getScatterCos()const{
        double rnd = randf();
        double oSi=0;
        double t;
        for(int i=0; i<nScatAng; i++){ // ToDo : This can be improved by bisection-search
            double Si = Sang[i];
            //printf( "S[%i] %g <? %g \n", i, Si, rnd );
            if(rnd<Si){
                t= i + (1-(Si-rnd)/(Si-oSi));
                break;
            }
            oSi=Si;
        }
        double alpha = M_PI*(t/nScatAng);
        //printf( "getScatterCos t %g alpha %g cos(a) %g  rnd %g \n", t, alpha, cos( alpha ), rnd );
        return cos( alpha ); // ToDo : in future we can avoid using cos() by clever choice of intervals
    }

    void setLorenz(double width){
        double sum=0;
        double dx = 1./nScatAng;
        for(int i=0; i<nScatAng; i++){
            double val = 1/(1+sq(i*dx/width));
            //printf( "val[%i] %g \n", i , val );
            sum+=val;
            Sang[i]=sum;
        }
        double renorm=1/sum;
        for(int i=0; i<nScatAng; i++){  Sang[i]*=renorm; };
        double oSi=0;
        for(int i=0; i<nScatAng; i++){  printf("S[%i] %g %g \n", i, Sang[i], Sang[i]-oSi ); oSi=Sang[i]; };
    }

};

struct RayedTri{
    Vec3i points;
    Vec2i materials;   //
};

class RayWorld{ public:
    std::vector<Vec3d>           points;
    //std::vector<RayTri>          triangles;
    std::vector<Vec3i>           triangles;
    std::vector<Vec2i>           tri2mat;
    std::vector<ScatterMatrial>  materials;
    //std::vector<RayedObject>     objs;
    //Rayt3d

    int nScatMax=64;
    bool bDraw = false;

    void addMesh( CMesh& mesh, Vec2i mat){
        int ip0=points.size();
        points.resize( ip0 + mesh.nvert );
        for(int i=0; i<mesh.nvert; i++){
            points[ip0+i] = mesh.verts[i];
        }
        int it0 = triangles.size();
        triangles.resize( it0 + mesh.ntri );
        tri2mat  .resize( it0 + mesh.ntri );
        for(int i=0; i<mesh.ntri; i++){
            triangles[it0+i] = ( mesh.tris[i] + ip0 );
            tri2mat  [it0+i] =  mat ;
        }
    }

    double ray( Vec3d& ray0, Vec3d& hRay, int& imat)const{
        Vec3d normal;
        int imin;
        double t = rayTriangles( ray0, hRay, triangles.size(), &triangles[0], &points[0], normal, imin );
        double cdot = normal.dot( hRay );
        if( cdot>0 ){ imat = tri2mat[imin].b; }
        else        { imat = tri2mat[imin].a; }
        //printf( "cdot %g t %g imat %i \n", cdot, t, imat );
        return t;
    }

    double scatterInMat(int imat, double tmax, Vec3d& hRay )const{
        const ScatterMatrial& mat = materials[imat];
        double prob  = tmax/mat.mfp;  // scattering probability
        int nscat    = (int)(prob*4)+1; // number of sub-steps - higher number here makes it more precise
        if(nscat>nScatMax) nscat=nScatMax;
        double dt    = tmax/nscat;
        double dprob = prob/nscat;
        double t     = 0;
        //printf( "scatter nscat %i prob %g \n", nscat, prob );
        for(int i=0; i<nscat; i++){
            double rnd=randf();
            //printf( "scatter[%i] dprob %g rnd %g \n", i, dprob, rnd );
            if( rnd<dprob ){
                t+=(1-rnd/dprob)*dt;
                break;
            }
            t+=dt;
        }
        // generate scattering angle and vector  - ToDo : this is perhaps inefficient way
        double ca = mat.getScatterCos();
        double sa = sqrt(1-ca*ca);
        Vec3d  a,b; hRay.getSomeOrtho(a,b);
        double phi = randf()*2*M_PI;
        double cph = cos(phi);
        double sph = sin(phi);
        hRay.mul(ca);
        hRay.add_mul(a,sa*cph);
        hRay.add_mul(b,sa*sph);
        //printf( "ca %g \n", ca );
        return t;
    }

    void scatter( int n, Vec3d& ray0, Vec3d& hRay )const{
        if(bDraw){
            glBegin(GL_LINE_STRIP);
            Draw3D::vertex( ray0 );
        }
        double tmax;
        for(int i=0; i<n; i++){
            int imat;
            tmax = ray( ray0, hRay, imat);
            if(tmax>1e+100) break;
            //printf( "scatter[%i] imat %i \n", i, imat );
            Vec3d ohRay = hRay;
            double t;
            if( imat>=0 ){ t = scatterInMat( imat, tmax, hRay ); }
            else         { t = tmax;                             }
            ray0.add_mul(ohRay, t*(1+1e-8) );
            if(bDraw)Draw3D::vertex( ray0 );
        }
        if(bDraw){
            if(tmax>1e+100){
                Draw3D::vertex( ray0+hRay*2.0 );
            }
            glEnd();
        }
    }

};


// ============= Application

class TestAppRayScattererMMC : public AppSDL2OGL_3D { public:

    //int defaultObjectShape;
    RayWorld scene;
    int ogl;
    bool bIntegrate=false;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppRayScattererMMC( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppRayScattererMMC::TestAppRayScattererMMC( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    ScatterMatrial* mat;
    scene.materials.resize(2);
    //mat = &scene.materials[0]; mat->scatterProb=0.0; //mat->setLorenz(double width);
    //mat = &scene.materials[0]; mat->scatterProb=2.0; mat->setLorenz(0.25);
    mat = &scene.materials[0]; mat->mfp=0.5; mat->setLorenz(0.25);

    CMesh msh;
    ogl = Draw::list( 0 );
    msh = Solids::Tetrahedron;
    scene.addMesh( msh, {-1,0} );
    Draw3D::drawMesh( msh, 0, 0xFF00FF00 );
    glColor3f( 0.,0.,1.); Draw3D::drawTriangles( msh.ntri, (int*)msh.tris, msh.verts,  true );

    glColor3f( 1.,1.,1.);
    //msh = Solids::Tetrahedron;
    double tg=0.1;
    scene.bDraw = true;
    int nray=5;
    for(int i=0; i<nray; i++){
        printf( "============ ray[%i] \n", i );
        Vec3d ray0 = (Vec3d){10.0,0.0,0.0};
        Vec3d hRay = (Vec3d){-1.0,randf(-tg,tg),randf(-tg,tg)}; hRay.normalize();
        Draw3D::drawPointCross( ray0, 0.1  );
        Draw3D::drawVecInPos  ( hRay, ray0 );
        scene.scatter( 10, ray0, hRay );
    }

    glEndList();


}

void TestAppRayScattererMMC::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	//glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    //glEnable( GL_LIGHTING );
    //glEnable(GL_DEPTH_TEST);

    glDisable( GL_LIGHTING   );
    glDisable( GL_DEPTH_TEST );

    int    nray   = 50;
    int    perRay = 10;
    double tg   = 0.1;
    if(bIntegrate){
        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        glColor4f(1.0,1.0,1.0,0.05);
        //msh = Solids::Tetrahedron;
        scene.bDraw = true;
        for(int i=0; i<nray; i++){
            Vec3d ray0 = (Vec3d){10.0,0.0,0.0};
            Vec3d hRay = (Vec3d){-1.0,randf(-tg,tg),randf(-tg,tg)}; hRay.normalize();
            Draw3D::drawPointCross( ray0, 0.1  );
            Draw3D::drawVecInPos  ( hRay, ray0 );
            scene.scatter( perRay, ray0, hRay );
        }
    }else{
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
        glCallList(ogl);
    }

};


void TestAppRayScattererMMC::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_m:  bIntegrate=!bIntegrate; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}


void TestAppRayScattererMMC::drawHUD(){
    glDisable ( GL_LIGHTING );
    /*
    glColor3f( 0.0f, 1.0f, 0.0f );
    glBegin( GL_LINES );
        float whalf = WIDTH *0.5;
        float hhalf = HEIGHT*0.5;
        glVertex3f( whalf-10,hhalf, 0 ); glVertex3f( whalf+10,hhalf, 0 );
        glVertex3f( whalf,hhalf-10, 0 ); glVertex3f( whalf,hhalf+10, 0 );
    glEnd();
    */
}

// ===================== MAIN

TestAppRayScattererMMC * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppRayScattererMMC( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















