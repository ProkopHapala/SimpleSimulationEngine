
/// @file @brief This program simulates the growth of a fractal structure using the Diffusion-Limited Aggregation (DLA) algorithm. It starts with a seed particle, and new particles are added one by one, performing a random walk until they collide and stick to the growing cluster. The result is a natural, tree-like or coral-like structure. The simulation is accelerated using a `CubicRuler.h` spatial grid.
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <unordered_map>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"
#include "Solids.h"

#include "DynamicOpt.h"
#include "CubicRuler.h"

#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

#include "testUtils.h"

// ======================  TestApp

class TestAppSphereTree : public AppSDL2OGL_3D {
	public:

    bool running = false;
    int perFrame = 10;

	int nSpheres = 0;
	const static int nMaxSpheres = 1024*4; //*64;
	Vec3d sphere_pos[nMaxSpheres];
	//std::vector<Vec3d>

	std::unordered_multimap<int_fast64_t,int>  grid;

	CubicRuler ruler;

    int sphereShape, cursorShape;

    Vec3i iCurPos;

    double rSphere  = 0.25;
    double rSphere2 = rSphere*rSphere;
    //double collisionDist  = rSphere2;
    double collisionDist    = rSphere*1.5;
    //double collisionDist  = rSphere2*2.0;
    double collisionDist2   = collisionDist*collisionDist;

    const double rStartOff = 2.0;
    double rStart          = rStartOff;
    double driftSpeed      = 0.001;

    double maxStep         = 0.25;
    Vec3d maxPos,minPos;

    bool DLArunning = true;
    Vec3d DLA_pos;

    Vec3d hRay,pos0Ray;
    double tRay = 0;
    static const int nMaxRayShapes=256;
    int nRayShapes = 0;
    int rayShapes[nMaxRayShapes];

	// ---- function declarations
    void ray_next();
    int  raySpheresOnGrid( const Vec3d& hRay, const Vec3d& ray0, double tmax, double R, double& tmin );

    void printMapContent();
    void insertWarper( int ix, int iy, int iz, int val );
	bool insertSphere( const Vec3d& pos, int_fast64_t ind );
    bool insertSphere( const Vec3d& pos );
    bool insertSphere( const Vec3d& pos, double r );

	bool DLAstep ( Vec3d& pos );
	bool DLAstep_BruteForce ( Vec3d& pos );
	void DLAstart( );

	virtual void draw   ();
    //virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	TestAppSphereTree ( int& id, int WIDTH_, int HEIGHT_ );

};

void TestAppSphereTree::ray_next(){
    double dt = ruler.ray_next();

    printf( " iRay (%03i,%03i,%03i) mdRay (%3.3f,%3.3f,%3.3f) dt %3.3f tRay %3.3f \n",
           ruler.iRay.a, ruler.iRay.b, ruler.iRay.c,
           ruler.mdRay.a,ruler.mdRay.b,ruler.mdRay.c, dt, tRay );

    Vec3d lpos, p1,p2;
    p1.set_add_mul( pos0Ray, hRay ,  tRay    );
    p2.set_add_mul( pos0Ray, hRay ,  tRay+dt );
    ruler.index2pos( ruler.iRay, {0.0,0.0,0.0}, lpos );
    int rayShape = glGenLists(1);
    glNewList( rayShape , GL_COMPILE );
        glColor3f(0.0f,0.0f,0.0f); Draw3D::drawShape ( cursorShape, lpos, {1.0,0.0,0.0,  0.0,1.0,0.0, 0.0,0.0,1.0} );
        glColor3f(0.0f,0.0f,0.8f); Draw3D::drawLine      ( p1, p2 );
        glColor3f(0.8f,0.0f,0.0f); Draw3D::drawPointCross( p2, 0.2 );
    glEndList();
    rayShapes[nRayShapes] = rayShape;
    tRay+=dt;
    nRayShapes++;
}

int TestAppSphereTree::raySpheresOnGrid( const Vec3d& hRay, const Vec3d& ray0, double tmax, double R, double& tmin ){
    ruler.ray_start( hRay, ray0 );
    double t = 0;
    while( t < tmax ){
        double dt = ruler.ray_next();
        int_fast64_t ind =  xyz2id( ruler.iRay.a, ruler.iRay.b, ruler.iRay.c );
        auto range = grid.equal_range( ind );
        tmin=tmax;
        double imin=-1;
        for ( auto it = range.first; it != range.second; ++it ){
            int i = it->second;
            Vec3d posi = sphere_pos[i];
            //double r2 = rayPointDistance2( ray0, hRay, posi, thit );
            double t = raySphere( ray0, hRay, R, posi );
            if( t<tmin ){ tmin=t; imin=i; };
            i++;
        }
        if( imin>=0 ){
            return imin;
        }
    }
    return -1;
}

void TestAppSphereTree::printMapContent(){
    int ip = 0;
    for(auto p : grid ) {
        int_fast16_t ix,iy,iz;
        id2xyz( p.first, ix, iy, iz );
        if( p.first != 0 ) printf( " %i %i (%i,%i,%i) \n", p.first, p.second, ix, iy, iz );
        ip++;
    }
}

bool TestAppSphereTree::insertSphere( const Vec3d& pos, int_fast64_t ind ){
    if( nSpheres >= (nMaxSpheres-2) ) return false;
    sphere_pos[nSpheres] = pos;
    grid.insert({ind,nSpheres});
    nSpheres++;
    return true;
};

bool TestAppSphereTree::insertSphere( const  Vec3d& pos ){
    int_fast64_t ind = ruler.pos2index( pos );
    return insertSphere( pos, ind );
};


void TestAppSphereTree::insertWarper( int ix, int iy, int iz, int val ){
    int_fast64_t key = xyz2id( ix, iy, iz );
    grid.insert( {key,val} );
    //printf( " insert %i %i (%i,%i,%i) \n", key, val, ix,iy, iz );
}

bool TestAppSphereTree::insertSphere( const Vec3d& pos, double r ){
    // FIXME : This function will not work for non-cubic grid;
    // Does it mean we should work only on cubic grid ?
    if( nSpheres >= (nMaxSpheres-2) ) return false;
    sphere_pos[nSpheres] = pos;


    Vec3i iabc;
    Vec3d dabc;
    ruler.pos2index( pos, dabc, iabc );
    //printf( "== dabc (%3.3f,%3.3f,%3.3f) \n", dabc.x,dabc.y,dabc.z );
    int dix=0,diy=0,diz=0;
    int mask = 0;
    double mr = 1-r;
    insertWarper( iabc.x, iabc.y, iabc.z ,nSpheres );
    if     ( dabc.x<r  ){ insertWarper( iabc.x-1, iabc.y  , iabc.z   ,nSpheres ); dix=-1;                    }
    else if( dabc.x>mr ){ insertWarper( iabc.x+1, iabc.y  , iabc.z   ,nSpheres ); dix=+1; dabc.x = 1-dabc.x; }
    if     ( dabc.y<r  ){ insertWarper( iabc.x  , iabc.y-1, iabc.z   ,nSpheres ); diy=-1;                    }
    else if( dabc.y>mr ){ insertWarper( iabc.x  , iabc.y+1, iabc.z   ,nSpheres ); diy=+1; dabc.y = 1-dabc.y; }
    if     ( dabc.z<r  ){ insertWarper( iabc.x  , iabc.y  , iabc.z-1 ,nSpheres ); diz=-1;                    }
    else if( dabc.z>mr ){ insertWarper( iabc.x  , iabc.y  , iabc.z+1 ,nSpheres ); diz=+1; dabc.z = 1-dabc.z; }
    double r2  = r*r;
    double dx2 = dabc.x*dabc.x;
    double dy2 = dabc.y*dabc.y;
    double dz2 = dabc.z*dabc.z;
    if( r2>(dx2+dy2    ) ){ insertWarper( iabc.x+dix, iabc.y+diy, iabc.z     ,nSpheres ); };
    if( r2>(dx2+dz2    ) ){ insertWarper( iabc.x+dix, iabc.y    , iabc.z+diz ,nSpheres ); };
    if( r2>(dy2+dz2    ) ){ insertWarper( iabc.x    , iabc.y+diy, iabc.z+diz ,nSpheres ); };
    if( r2>(dx2+dy2+dz2) ){ insertWarper( iabc.x+dix, iabc.y+diy, iabc.z+diz ,nSpheres ); };

/*
    Vec3i iabc;
    Vec3d dabc;
    ruler.pos2index( pos, dabc, iabc );
    int dix=0,diy=0,diz=0;
    int mask = 0;
    double mr = 1-r;
    if     ( dabc.x<r  ){ grid.insert({xyz2id( iabc.x-1, iabc.y  , iabc.z   ),nSpheres}); dix=-1;                    printf( "\n" ); }
    else if( dabc.x>mr ){ grid.insert({xyz2id( iabc.x+1, iabc.y  , iabc.z   ),nSpheres}); dix=+1; dabc.x = 1-dabc.x; }
    if     ( dabc.y<r  ){ grid.insert({xyz2id( iabc.x  , iabc.y-1, iabc.z   ),nSpheres}); diy=-1;                    }
    else if( dabc.y>mr ){ grid.insert({xyz2id( iabc.x  , iabc.y+1, iabc.z   ),nSpheres}); diy=+1; dabc.y = 1-dabc.y; }
    if     ( dabc.z<r  ){ grid.insert({xyz2id( iabc.x  , iabc.y  , iabc.z-1 ),nSpheres}); diz=-1;                    }
    else if( dabc.z>mr ){ grid.insert({xyz2id( iabc.x  , iabc.y  , iabc.z+1 ),nSpheres}); diz=+1; dabc.z = 1-dabc.z; }
    double r2  = r*r;
    double dx2 = dabc.x*dabc.x;
    double dy2 = dabc.y*dabc.y;
    double dz2 = dabc.z*dabc.z;
    if( r2>(dx2+dy2    ) ){ grid.insert({xyz2id( iabc.x+dix, iabc.y+diy, iabc.z     ),nSpheres}); };
    if( r2>(dx2+dz2    ) ){ grid.insert({xyz2id( iabc.x+dix, iabc.y    , iabc.z+diz ),nSpheres}); };
    if( r2>(dy2+dz2    ) ){ grid.insert({xyz2id( iabc.x    , iabc.y+diy, iabc.z+diz ),nSpheres}); };
    if( r2>(dx2+dy2+dz2) ){ grid.insert({xyz2id( iabc.x-dix, iabc.y+diy, iabc.z+diz ),nSpheres}); };
*/


    /*
    Vec3d pmin,pmax, dmin,dmax;
    Vec3i imin,imax;
    pmin.set_add( pos, -r );   ruler.pos2index( pmin,dmin, imin );
    pmax.set_add( pos,  r );   ruler.pos2index( pmin,dmin, imin );
    for( int ix=imin.x; ix<=imax.x; ix++ ){
        for( int iy=imin.y; iy<=imax.y; iy++ ){
            for( int iz=imin.z; iz<=imax.z; iz++ ){

            }
        }
    }
    */
    nSpheres++;
    return true;
};

void TestAppSphereTree::DLAstart( ){
    double z   = randf( -1.0,1.0 );
    double xy  = rStart * sqrt (  1 - z*z );
    z         *= rStart;
    double phi = randf( 0.0,M_PI*2.0);
    DLA_pos.set( xy*cos(phi), xy*sin(phi), z );
}

bool TestAppSphereTree::DLAstep_BruteForce( Vec3d& pos ){

    Vec3d dpos;
    dpos.set( randf( -maxStep, maxStep ), randf( -maxStep, maxStep ), randf( -maxStep, maxStep ) );
    pos.add( dpos );
    double r = pos.norm();
    if( r > rStart ){ pos.mul( rStart/r ); }

    for( int i=0; i<nSpheres; i++ ){
        Vec3d posi = sphere_pos[i];
        double dr2 = pos.dist2( posi );
        if ( dr2 < collisionDist2 ){
            if( (r + rStartOff) > rStart ) rStart = r + rStartOff;
            //insertSphere( pos );
            insertSphere( pos, collisionDist );
            //printf( " insert (%3.3f,%3.3f,%3.3f) %i \n", pos.x,pos.y,pos.z, nSpheres );
            //printf( " insert %i (%3.3f,%3.3f,%3.3f) %f \n", nSpheres, pos.x,pos.y,pos.z, dr2  );
            return false;
        }
    }
    return true;
}

bool TestAppSphereTree::DLAstep( Vec3d& pos ){

    Vec3d dpos;
    dpos.set( randf( -maxStep, maxStep ), randf( -maxStep, maxStep ), randf( -maxStep, maxStep ) );
    pos.add( dpos );
    double r = pos.norm();
    if( r > rStart ){ pos.mul( rStart/r ); }

    // with HashMap acceleration
    Vec3d dabc;
    Vec3i iabc;
    ruler.pos2index( pos, dabc, iabc );
    int_fast64_t ind = xyz2id( iabc.x , iabc.y, iabc.z );
    //int_fast64_t ind = ruler.pos2index( pos );
    auto range = grid.equal_range( ind );
    //printf( " %i  (%i,%i,%i)  %i \n", ineigh, ix, iy, iz, ind );
    //printf( " %f  (%3.3f,%3.3f,%3.3f)  %i \n", r, pos.x, pos.y, pos.z, ind );
    for ( auto it = range.first; it != range.second; ++it ){
        Vec3d posi = sphere_pos[it->second];
        double dr2 = pos.dist2( posi );
        //printf( " > %f  (%3.3f,%3.3f,%3.3f) \n", dr2, posi.x, posi.y, posi.z );
        if ( dr2 < collisionDist2 ){
            if( (r + rStartOff) > rStart ) rStart = r + rStartOff;
            //insertSphere( pos, ind );
            insertSphere( pos, collisionDist );
            //printf( " insert %i (%3.3f,%3.3f,%3.3f) %f \n", nSpheres, pos.x,pos.y,pos.z, dr2  );
            //printf( " insert (%i,%i,%i) %i %i \n", iabc.x , iabc.y, iabc.z, ind, nSpheres );
            return false;
        }
    }

/*
    for( int i=0; i<nSpheres; i++ ){
        Vec3d posi = sphere_pos[i];
        double dr2 = pos.dist2( posi );
        if ( dr2 < collisionDist2 ){
            printf( " HashMap missed %i (%3.3f,%3.3f,%3.3f) %f \n", nSpheres, pos.x,pos.y,pos.z, dr2  );

            printf( " pos  missed %i (%3.3f,%3.3f,%3.3f) (%i,%i,%i) ind %i \n", nSpheres, pos.x,pos.y,pos.z, iabc.x,iabc.y,iabc.z,  ind  );
            ruler.pos2index( posi, dabc, iabc );
            ind = xyz2id( iabc.x , iabc.y, iabc.z );
            printf( " posi missed %i (%3.3f,%3.3f,%3.3f) (%i,%i,%i) ind %i \n", nSpheres, pos.x,pos.y,pos.z, iabc.x,iabc.y,iabc.z,  ind  );


            return false;
        }
    }
*/

    return true;
}


TestAppSphereTree::TestAppSphereTree( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    grid.reserve(1024);
    grid.max_load_factor ( 0.5 );

    maxPos.set(  100.0, 100.0, 100.0 );
    minPos.set( -100.0,-100.0,-100.0 );

    double gridStep = 1.0;

    ruler.setup( minPos*2.0, { gridStep,0.0,0.0,  0.0,gridStep,0.0,  0.0,0.0,gridStep } );

    double span = 10.0;



/*
    for(int i=0; i<32; i++){
        int_fast64_t id;
        //int          ix,iy,iz;
        int_fast64_t ix_,iy_,iz_;
        Vec3d pos,dabc;
        Vec3i iabc;
        pos.set( randf(-span,span), randf(-span,span), randf(-span,span) );
        //ix = rand()&0xFF;        //iy = rand()&0xFF;        //iz = rand()&0xFF;
        //id = xyz2id(     ix,  iy,  iz );
        ruler.pos2index( pos, dabc, iabc );
        id = xyz2id(     iabc.x,  iabc.y,  iabc.z );
        id2xyz     ( id, ix_, iy_, iz_ );
        //printf( " %i : (%i,%i,%i)  %i (%i,%i,%i) \n",i,  ix, iy, iz, id,  ix_, iy_, iz_ );
        printf( " %i : (%3.3f,%3.3f,%3.3f) (%i,%i,%i)  %i (%i,%i,%i) \n",i, pos.x,  pos.y,  pos.z,  iabc.x,  iabc.y,  iabc.z, id,  ix_, iy_, iz_ );
    }

*/
    //exit(0);

    insertSphere( {0.0,0.0,0.0}, rSphere  );
    //insertSphere( {1.0,-1.0,1.0}, rSphere );
    printMapContent();
    //exit(0);

    sphereShape = glGenLists(1);
    glNewList( sphereShape , GL_COMPILE );
        //glPushMatrix();
        //glDisable ( GL_LIGHTING );
        //Draw3D::drawAxis ( 3.0f );
        //glColor3f( 1.0f, 0.0f, 1.0f ); Draw3D::drawLines   ( Solids::Icosahedron_nedges, Solids::Icosahedron_edges, Solids::Icosahedron_verts                             );
        //glEnable( GL_LIGHTING );
        //glColor3f( 0.8f, 0.8f, 0.8f ); Draw3D::drawPolygons( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons, Solids::Icosahedron_faces,  Solids::Icosahedron_verts );
        glEnable( GL_LIGHTING );
        //glColor3f( 0.8f, 0.8f, 0.8f );
        Draw3D::drawSphere_oct( 2, (float)rSphere, {0.0,0.0,0.0} );
        //glPopMatrix();
    glEndList();

    cursorShape = glGenLists(1);
    glNewList( cursorShape , GL_COMPILE );
        glPushMatrix();
        glTranslatef(0.5f,0.5f,0.5f);
        glScalef( (float)gridStep*0.5, (float)gridStep*0.5, (float)gridStep*0.5 );
        glColor3f(0.01f,0.01f,0.01f); Draw3D::drawLines ( Solids::Cube_nedges, (int*)Solids::Cube_edges, Solids::Cube_verts );
        //Draw3D::drawAxis( 0.5f );
        glPopMatrix();
    glEndList();

    DLAstart( );

    Vec3d dabc;
    ruler.pos2index( {0.0,0.0,0.0}, dabc, iCurPos );


    srand(114545);

    pos0Ray.set( 0.5, 0.5, 0.5  );
    hRay   .set( 1.0, -0.5, 2.25 );
    hRay.normalize();
    ruler.ray_start( hRay, pos0Ray );
    tRay=0;

    printf( " hRay (%3.3f,%3.3f,%3.3f) \n", hRay.a, hRay.b, hRay.c );
    printf( " hRay (%3.3f,%3.3f,%3.3f) hRayInv (%3.3f,%3.3f,%3.3f) dt %3.3f \n",
       ruler.hRay.a, ruler.hRay.b, ruler.hRay.c,
       ruler.hRayInv.a, ruler.hRayInv.b, ruler.hRayInv.c );

    printf( " iRay (%03i,%03i,%03i) mdRay (%3.3f,%3.3f,%3.3f) dt %3.3f \n",
       ruler.iRay.a, ruler.iRay.b, ruler.iRay.c,
       ruler.mdRay.a,ruler.mdRay.b,ruler.mdRay.c );

};

void TestAppSphereTree::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glEnable(GL_DEPTH_TEST);

    float ambient  [] = { 0.8f, 0.8f, 0.8f, 1.0f };
	float diffuse  [] = { 0.2f, 0.2f,  0.2f,  1.0f };
	float specular [] = { 0.0f, 0.0f,  0.0f,  1.0f };
	//float shininess[] = { 80.0f                    };
	//glMaterialfv ( GL_FRONT_AND_BACK, GL_AMBIENT,   ambient);
	//glMaterialfv ( GL_FRONT_AND_BACK, GL_DIFFUSE,   diffuse);
	//glMaterialfv ( GL_FRONT_AND_BACK, GL_SPECULAR,  specular);
	glEnable     ( GL_COLOR_MATERIAL   );
	glLightfv    ( GL_LIGHT0, GL_AMBIENT,  ambient  );
	glLightfv    ( GL_LIGHT0, GL_DIFFUSE,  diffuse  );
	glLightfv    ( GL_LIGHT0, GL_SPECULAR, specular  );
	//glMaterialfv ( GL_FRONT_AND_BACK, GL_SHININESS, shininess);

    Vec3d lpos;
    Mat3d lrot; lrot.setOne();

    for(int i=0; i<nRayShapes; i++){
        glCallList( rayShapes[i] );
    }

    int perFrame = 400000;
    glColor3f(0.8f,0.2f,0.2f);
    long tcomp = getCPUticks();
    for( int i=0; i<perFrame; i++ ){
        if( nSpheres < (nMaxSpheres-2) ){
            if( DLArunning ){
                DLArunning = DLAstep( DLA_pos );
                //DLArunning = DLAstep_BruteForce( DLA_pos );
                //Draw3D::drawShape  ( DLA_pos, lrot, sphereShape );
            }else{
                DLAstart( );
                DLArunning = true;
            }
        }
    }
    tcomp = getCPUticks() - tcomp;

    //for( auto o : world.objects ) {
    //glColor3f(0.8f,0.8f,0.8f);
    long tview = getCPUticks();
    for( int i=0; i<nSpheres; i++ ){
        lpos = sphere_pos[i];
        double r = lpos.norm();
        float c  = (float)(r/rStart);
        glColor3f( c, 2*c*(1-c), 0 );
        Draw3D::drawShape( sphereShape, lpos, lrot );
    }
    tview = getCPUticks() - tview;

    double t;
    int isph = raySpheresOnGrid( (Vec3d)cam.rot.c, (Vec3d)cam.rot.c*(-10.0), 10.0, rSphere, t );
    if(isph>=0){
        lpos = sphere_pos[isph];
        glColor3f( 0.0, 1.0, 1.0 );
        Draw3D::drawShape ( sphereShape, lpos, {1.01,0.0,0.0,  0.0,1.01,0.0, 0.0,0.0,1.01} );
    }

    glColor3f( 1.0, 0.0, 1.0 );
    int_fast64_t ind = xyz2id(  iCurPos.x,  iCurPos.y,  iCurPos.z );
    auto range = grid.equal_range( ind );
    for ( auto it = range.first; it != range.second; ++it ){
        lpos = sphere_pos[it->second];
        Draw3D::drawShape ( sphereShape, lpos, {1.01,0.0,0.0,  0.0,1.01,0.0, 0.0,0.0,1.01} );
    }

//    printf( " frame %05i : alg %3.3f(%3.3f+e6) view %3.3f(%3.3f+e6) | nsph %05i \n", frameCount, (tcomp/(double)perFrame), tcomp*1e-6, (tview/(double)nSpheres), tview*1e-6, nSpheres );

    glColor3f( 1.0f, 1.0f, 1.0f );
	//glDisable(GL_DEPTH_TEST);
	ruler.index2pos( iCurPos, {0.0,0.0,0.0}, lpos );
	Draw3D::drawShape( cursorShape, lpos, lrot );
	//Draw3D::drawShape( lpos, lrot, sphereShape );

};

void TestAppSphereTree::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    int rot;
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p: printMapContent(); break;
                case SDLK_a:  iCurPos.x++; break;
                case SDLK_d:  iCurPos.x--; break;
                case SDLK_w:  iCurPos.y++; break;
                case SDLK_s:  iCurPos.y--; break;
                case SDLK_q:  iCurPos.z++; break;
                case SDLK_e:  iCurPos.z--; break;

                case SDLK_x:  qCamera.set( 0.0, 0.0, 0.0 ,1.0 ); break;
                case SDLK_y:  qCamera.set( 0.7071, 0.0, 0.0 ,0.7071 );; break;
                case SDLK_z:  qCamera.set( 0.0, 0.7071, 0.0 ,0.7071 ); break;

                case SDLK_n:  ray_next(); break;

                //case SDLK_a:  ix --; if( ix <  0            ) ix = builder.nMax-1;   break;
            }
            Vec3d pos;
            ruler.index2pos( iCurPos ,{0.0,0.0,0.0}, pos );
            int_fast64_t id = xyz2id(  iCurPos.x,  iCurPos.y,  iCurPos.z );
            //printf( " cursor (%i,%i,%i) %i (%3.3f,%3.3f,%3.3f) \n", iCurPos.x, iCurPos.y, iCurPos.z, id, pos.x, pos.y, pos.z );
            break;
        /*
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    pickParticle( world.picked );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    world.picked = NULL;
                    break;
            }
            break;
        */
    };


    AppSDL2OGL_3D::eventHandling( event );
}



// ===================== MAIN

TestAppSphereTree * testApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppSphereTree( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
