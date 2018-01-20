
// Copied tutorial from
//  http://www.opengl-tutorial.org/intermediate-tutorials/billboards-particles/particles-instancing/
//  https://github.com/opengl-tutorials/ogl/tree/master/tutorial18_billboards_and_particles

#include <stdlib.h>
#include <stdio.h>

#include <vector>
#include <algorithm>

#include <GL/glew.h>
#include <SDL2/SDL.h>

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Mat4.h"

#include "GL3Utils.h"
#include "GLObject.h"
#include "GLfunctions.h"
#include "GLobjects.h"
#include "GLInstances.h"
#include "IO_utils.h"
#include "Shader.h"


#include "DrawOGL3.h"

#include "Solids.h"
#include "CMesh.h"

#include "testUtils.h"

// =============== Global variables

Shader *shBranches,*shLeafs,*shTex,*shTexView;

//GLInstances instances;
//GLBillboards bilboards;
//GLuint texture_1;

GLMesh   *mshBranches,*mshLeafs,*msh1;

GLuint       texTest;
GLMesh       *glSprite,*glAtlas,*glAtlas_;
GLBillboards bilboards;
FrameBuffer  frameBuff1;


int frameCount = 0;
double lastTime = 0.0;

GLuint vao;     // vertex array object
int delay = 10;
int VSync = 0;
SDL_Window * window     = NULL;
SDL_GLContext   context = NULL;
int WIDTH  = 800;
int HEIGHT = 800;
float ASPECT_RATIO = HEIGHT/WIDTH;

int mouseX, mouseY;
Quat4f qCamera;
Mat3f  mouseMat;
Vec3f  camPos = (Vec3f){ 0.0f, 0.0f, -10.0f };

bool bTransparent  = false;

long lastCPUtick = 0;
double ticks_per_second=0;

// =============== Functions

GLMesh* makeQuad3D( Vec2f p0, Vec2f p1, Vec2f u0, Vec2f u1 ){
    //GLfloat verts = new GLfloat[3*2*3];
    //GLfloat vUVs  = new GLfloat[3*2*2];
    //delete [] verts;
    GLfloat verts[] = { p0.x,p0.y,0.0, p1.x,p1.y,0.0, p0.x,p1.y,0.0,   p0.x,p0.y,0.0, p1.x,p1.y,0.0, p1.x,p0.y,0.0 };
    Vec2f vUVs   [] = { u0.x,u0.y,     u1.x,u1.y,     u0.x,u1.y,       u0.x,u0.y,     u1.x,u1.y,     u1.x,u0.y     };
    GLMesh* glquad = new GLMesh();
    glquad->init( 6, 0, NULL, verts, NULL, NULL, vUVs );
    return glquad;
}

GLuint makeTestTextureRGBA( int W, int H ){
    double dx = 1.0d/W;
    double dy = 1.0d/H;
    uint32_t * c_img1 = new uint32_t[H*W];
    for( int iy=0; iy<H; iy++ ){
        for( int ix=0; ix<W; ix++ ){
            uint8_t r,g,b,a;
            r=ix; g=ix^iy; b=iy; a=255;
            c_img1[ iy*W + ix ] = (a<<24) | (b<<16) | (g<<8) | (r);
        }
    }
    GLuint texID;
    newTexture2D( texID, W, H, c_img1, GL_RGBA, GL_UNSIGNED_BYTE );
    return texID;
}

int pushCylinderTris( int n, float r1, float r2, Vec3f base, Vec3f tip, Vec3f up, std::vector<Vec3f>& verts, std::vector<Vec3f>* normals ){
	int nvert=0;

	Vec3f dir,left;
	dir.set_sub( tip, base );
	dir.normalize();
    left = dir.getOrtho(up);

    float alfa = 2*M_PI/n;
    Vec2f rot,drot;
    rot .set(1.0f,0.0f);
    drot.set( cos( alfa ), sin( alfa ) );

	Vec3f q; q.set(dir); q.add_mul( up, -(r1-r2) );
	float pnab =  dir.dot( q )/q.norm();
	float pnc  =  sqrt( 1 - pnab*pnab );

	Vec3f op,opn;
	for(int i=0; i<=n; i++ ){
		Vec3f p,pn;
		p .set( rot.x*up.x + rot.y*left.x, rot.x*up.y + rot.y*left.y, rot.x*up.z + rot.y*left.z );
		pn.set( pnab*p.x   + pnc*dir.x   , pnab*p.y   + pnc*dir.y   , pnab*p.z   + pnc*dir.z    );
		if( i>0 ){
            verts  .push_back(base+op*r1); verts  .push_back(tip+p*r2); verts  .push_back(tip+op*r2);
            verts  .push_back(base+op*r1); verts  .push_back(tip+p*r2); verts  .push_back(base+p*r1);
            if(normals){
                normals->push_back(opn); normals->push_back(pn); normals->push_back(opn);
                normals->push_back(opn); normals->push_back(pn); normals->push_back(pn);
            }
		}
		op=p; opn=pn;
        rot.mul_cmplx( drot );
	}
	return nvert;
};

void pushTris_Treestep( int level, Vec3f pos, Vec3f dir, Vec3f up, std::vector<Vec3f>& verts, std::vector<Vec3f>* normals ){
    float l = dir.norm();
    static const float drnd = 0.15;
    dir.add( randf(-drnd,drnd)*l, randf(-drnd,drnd)*l, randf(-drnd,drnd)*l );
    float f = randf(0.5,0.9);
    dir.mul( f );
    up=cross(dir,up); up.normalize();
    Vec3f pos_ = pos + dir;
    pushCylinderTris( 6, 0.1*l, 0.1*l*f, pos, pos_, (Vec3f){0.0f,1.0f,0.0f},verts,normals);
    //pushCylinderTris(verts,normals);
    if( level>0 ){
        level--;
        pushTris_Treestep( level, pos_, dir+up*0.3*l, up, verts, normals );
        pushTris_Treestep( level, pos_, dir-up*0.3*l, up, verts, normals );
    }
}

void tree_step( int level, Vec3f pos, Vec3f dir, std::vector<Vec3f>& branches, std::vector<Vec3f>& leafs ){
    static const float drnd = 0.6;
    //dir.x *= randf(1.0-drnd,1.0);
    //dir.y *= randf(1.0-drnd,1.0);
    //dir.z *= randf(1.0-drnd,1.0);
    float l = dir.norm();
    dir.add( randf(-drnd,drnd)*l, randf(-drnd,drnd)*l, randf(-drnd,drnd)*l );
    dir.mul( randf(0.5,0.9) );
    Vec3f pos_ = pos + dir;
    branches.push_back(pos );
    branches.push_back(pos_);
    if( level==0 ){
        leafs.push_back(pos_);
    }else{
        level--;
        tree_step( level, pos_, dir,  branches, leafs );
        tree_step( level, pos_, dir,  branches, leafs );
    }
}

void camera(float aspect, Mat4f& camMat){
    Mat4f mRot,mProj;
    qCamera.toMatrix(mouseMat);
    mRot.setOne(); mRot.setRot(mouseMat);
    // float fov = 3.0; mPersp.setPerspective( fov, fov*ASPECT_RATIO, 1.0, 1000.0 );
    mProj.setOrthographic( 8.0, 8.0*aspect, -1.0, -1000.0);
    //mProj.setOne();
    camMat.set_mmul_TN( mRot, mProj );
}


void renderPhases( Mat4f camMat ){

    //Mat4f camMat;  camera( 32, camMat);
    int nPhases = 32;
    Mat3f modelMat; modelMat.setOne(); modelMat.mul(0.3);
    //Vec3f modelPos; modelPos.set(-4.25f,-4.0f,10.0f);
    Vec3f modelPos; modelPos.set(-8.25f,-0.5f,10.0f);
    for(int iph=0; iph<nPhases; iph++){

        modelMat.rotate( 2*M_PI/nPhases, (Vec3f){0.0f,1.0f,0.0f});
        modelPos.add( 0.5f, 0.0f, 0.0f );

        Shader * sh;
        sh = shBranches;
        sh->use();
        sh->set_camPos  ( (GLfloat*)&camPos );
        sh->set_camMat  ( (GLfloat*)&camMat );
        sh->set_modelPos( (GLfloat*)&modelPos );
        sh->set_modelMat( (GLfloat*)&modelMat );

        glEnable(GL_DEPTH_TEST);
        mshBranches->draw(GL_TRIANGLES);
        //msh1->draw(GL_TRIANGLES);
    }
}


void renderPhases2D( Mat4f camMat, int nPhases ){
    //Mat4f camMat;  camera( 32, camMat);
    //int nPhases = 32;
    Mat3f modelMat;
    //Vec3f modelPos; modelPos.set(-4.25f,-4.0f,10.0f);
    //Vec3f modelPos; modelPos.set(-8.25f,-0.5f,10.0f);
    Vec3f modelPos;
    glEnable(GL_DEPTH_TEST);

    float dth =   0.5*M_PI/nPhases;
    for(int ith=-nPhases+1;ith<nPhases; ith++){
        int nph    = nPhases-ith;
        float dph  = 2*M_PI/nph;
        for(int iph=0; iph<nph; iph++){
            modelMat.setOne();
            modelMat.mul(0.25);
            /*
            modelMat.rotate( dph*iph, (Vec3f){0.0f,1.0f,0.0f});
            modelMat.rotate( dth*ith, (Vec3f){1.0f,0.0f,0.0f});
            modelPos=(Vec3f){ -7.5+1.0*iph, 7.5-1.0f*ith, 10.0f};
            */

            modelMat.rotate( dph*iph, (Vec3f){0.0f,1.0f,0.0f});
            modelMat.rotate( dth*ith, (Vec3f){1.0f,0.0f,0.0f});
            modelPos=(Vec3f){ -7.5+1.0*(iph+ith),-7.5-1.0f*(+ith-iph), 10.0f};

            Shader * sh;
            sh = shBranches;
            sh->use();
            sh->set_camPos  ( (GLfloat*)&camPos );
            sh->set_camMat  ( (GLfloat*)&camMat );
            sh->set_modelPos( (GLfloat*)&modelPos );
            sh->set_modelMat( (GLfloat*)&modelMat );

            mshBranches->draw(GL_TRIANGLES);
            msh1->draw(GL_TRIANGLES);
        }
    }
}

void renderPhasesOct( Mat4f camMat, int n ){
    Mat3f modelMat;
    Vec3f modelPos;
    glEnable(GL_DEPTH_TEST);

    Vec3f up=(Vec3f){0.0,1.0,0.0};
    Vec3f d;
    float step = 1.0/n;
    for(int ix=-n;ix<n; ix++){
        for(int iy=-n;iy<n; iy++){
            d.x=(ix+0.5);
            d.z=(iy+0.5);
            modelPos=(Vec3f){ 1.0*d.x,1.0f*d.z, 10.0f};

            //modelMat.setOne();

            d.mul(step);

            d.y = 1.0-fabs(d.x)-fabs(d.z);
            if( d.y<0 ){
                //continue;
                if(d.x<0.0){d.x+=1.0;}else{d.x-=1.0;}
                if(d.z<0.0){d.z+=1.0;}else{d.z-=1.0;}
            }
            //printf( " %i %i :  %f %f %f \n", ix,iy, d.x, d.y, d.z );

            d.normalize();
            modelMat.fromDirUp(d,up);
            modelMat.T();
            modelMat.mul(0.25);

            Shader * sh;
            sh = shBranches;
            sh->use();
            sh->set_camPos  ( (GLfloat*)&camPos );
            sh->set_camMat  ( (GLfloat*)&camMat );
            sh->set_modelPos( (GLfloat*)&modelPos );
            sh->set_modelMat( (GLfloat*)&modelMat );

            mshBranches->draw(GL_TRIANGLES);
            msh1->draw(GL_TRIANGLES);
        }
    }
}

Vec2f selectPhaseOct( Vec3f dir, int n ){
    //Vec3f up=(Vec3f){0.0,1.0,0.0};
    float step=0.5/n;

    float renorm = 1.0/(fabs(dir.x)+fabs(dir.y)+fabs(dir.z));
    //float x = (dir.x+1.0)+n;
    //float y = (dir.z+1.0)+n;
    Vec2f uv = (Vec2f){ (int)(dir.x*renorm*n+n)*step, (int)(dir.z*renorm*n+n)*step };
    if( dir.y<0 ){
        uv.x-=0.5;
        uv.y-=0.5;
    };
    printf( "  (%f,%f,%f)   (%f,%f) %f \n", dir.x, dir.y, dir.z,   uv.x, uv.y,  renorm );
    return uv;
    //glUniform2f(sh->getUloc("uv0s"), d.x*n, d.y*n );
}


int setup(){

    shBranches=new Shader();
    //shBranches->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
    shBranches->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/normal2color.glslf"   );
    shBranches->getDefaultUniformLocation();

    shLeafs=new Shader();
    shLeafs->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/pointSprite.glslf"   );
    shLeafs->getDefaultUniformLocation();

    shTexView=new Shader();
    shTexView->init( "common_resources/shaders/texture3D.glslv",   "common_resources/shaders/texture.glslf"   );
    shTexView->getDefaultUniformLocation();

    std::vector<Vec3f> branches;
    std::vector<Vec3f> leafs;
    tree_step( 5, (Vec3f){0.0f,0.0f,0.0f}, (Vec3f){0.0f,1.0f,0.0f}, branches, leafs );

    std::vector<Vec3f> verts;
    std::vector<Vec3f> normals;
    pushTris_Treestep( 5, (Vec3f){0.0f,0.0f,0.0f}, (Vec3f){0.0f,1.0f,0.0f}, (Vec3f){1.0f,0.0f,0.0f}, verts, &normals );

    //for(int i=0; i<branches.size(); i++){ printf("%i (%g,%g,%g)\n", i, branches[i].x, branches[i].y, branches[i].z ); }

    msh1 = hardTriangles2mesh( Solids::Cube );
    //msh1 = hardTriangles2mesh( Solids::Octahedron );


    mshBranches = new GLMesh();
    //mshBranches->init( branches.size(), 0, NULL, &branches[0],  NULL, NULL, NULL );
    mshBranches->init( verts.size(), 0, NULL, &verts[0],  &normals[0], NULL, NULL );

    mshLeafs = new GLMesh();
    mshLeafs->init( leafs.size(), 0, NULL, &leafs[0],  NULL, NULL, NULL );

    //ticks_per_second = calibrate_timer(100);
    //lastTime = glfwGetTime();
    lastTime = 0.0;
    qCamera.setOne();

    // ==== BEGIN : RENDER TO TEXTURE

    shTex=new Shader();
    //shTex->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
    //shTex->init( "common_resources/shaders/texture3D.glslv",   "common_resources/shaders/texture.glslf"   );
    //shTex->init( "common_resources/shaders/texture3D_anim.glslv",   "common_resources/shaders/texture.glslf"   );
    //shTex->init( "common_resources/shaders/texture3D_anim.glslv",   "common_resources/shaders/texture_anim.glslf"   );
    shTex->init( "common_resources/shaders/textureAtlas3D.glslv",   "common_resources/shaders/texture.glslf"   );
    shTex->getDefaultUniformLocation();

    glAtlas  = makeQuad3D( {0.0f,0.0f}, {32.0f,32.0f}, {0.0f,0.0f},  {1.f,1.0f} );
    glAtlas_ = makeQuad3D( {0.0f,0.0f}, {32.0f,32.0f}, {-0.5f,-0.5f},{0.5f,0.5f} );

    //glquad = makeQuad3D( {0.0f,0.0f}, {0.5f,1.0f}, {0.0f,0.0f}, {0.125f,1.0f} );
    glSprite = makeQuad3D( {0.0f,0.0f}, {8.0f,8.0f}, {0.0f,0.0f}, {1.0f,1.0f} );

    //glquad =new GLMesh();
    //glquad->init( 6, 0, NULL, DEFAULT_Bilboard_verts, NULL, NULL, DEFAULT_Bilboard_UVs );
    //glquad->init( 6, 0, NULL, DEFAULT_Bilboard_verts, DEFAULT_Bilboard_verts, NULL, NULL );

    //texTest = makeTestTextureRGBA( 256, 256);

    //frameBuff1.init( 2048, 64 );
    //frameBuff1.init( 2048, 128 );
    frameBuff1.init( 2048, 2048 );

    // ===== Prepare texture by rendering

    //glBindFramebuffer(GL_FRAMEBUFFER, frameBuff1.buff );
    frameBuff1.bind();

    //glClearColor(0.0, 0.0, 0.8, 1.0);
    glClearColor(0.9, 0.9, 0.9, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    Mat4f camMat;
    //camera( 1/16.0, camMat);
    camera( 1.0, camMat);
    //renderPhases( camMat);
    //renderPhases2D( camMat, 16 );
    renderPhasesOct( camMat, 8 );

    // ==== END   : RENDER TO TEXTURE

    return 0;
};

void draw( ){
    glBindFramebuffer(GL_FRAMEBUFFER, 0); glViewport(0,0,WIDTH,HEIGHT);
    glClearColor(0.8, 0.8, 0.8, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );
	// Simulate all particles

	Mat4f camMat;  camera( 1.0, camMat);
	//renderPhases(camMat);

    Mat3f modelMat; modelMat.setOne(); modelMat.mul(0.5f);
    //Vec3f modelPos;

    //glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, texTest );
    glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, frameBuff1.texRGB );

    // ======= Bare Texture Atlas
    Shader * sh;


    //sh = shTex;
    sh = shTexView;
    sh->use();
    sh->set_modelMat( (GLfloat*)&modelMat );
    sh->set_camPos  ( (GLfloat*)&camPos   );
    sh->set_camMat  ( (GLfloat*)&camMat   );
    glUniform1i(sh->getUloc("texture_1"), 0);

    sh->set_modelPos( (const GLfloat[]){-4.0f,-16.0f,5.0f} );
    glAtlas->draw(GL_TRIANGLES); glAtlas->drawPoints(10.0);

    sh->set_modelPos( (const GLfloat[]){-4.0f-8.0f,-24.0f,5.0f} );
    glAtlas_->draw(GL_TRIANGLES); glAtlas_->drawPoints(10.0);

    // ======= Animated Texture
    //sh = shTex;


    sh = shTex;
    sh->use();
    sh->set_modelMat( (GLfloat*)&modelMat );
    sh->set_camPos  ( (GLfloat*)&camPos   );
    sh->set_camMat  ( (GLfloat*)&camMat   );

    glUniform1i(sh->getUloc("texture_1"), 0);
    //glUniform1f(sh->getUloc("nPhases"), 32.0f);

    //glUniform2f(sh->getUloc("uv0s"), 0.0   ,0.0);

    //float the = +1.5;
    float the = -0.2;
    float phi = frameCount*0.01;
    Vec2f cst; cst.fromAngle(the);
    Vec2f uv = selectPhaseOct( (Vec3f){cst.a*cos(phi),cst.b,cst.a*sin(phi)}, 8 );
    glUniform2f(sh->getUloc("uv0"), uv.x, uv.y );
    glUniform2f(sh->getUloc("du" ), 1/16.0,0.0 );
    glUniform2f(sh->getUloc("dv" ), 0.0   ,1/16.0);

    sh->set_modelPos( (const GLfloat[]){0.0f,0.0,10.0f} );
    glSprite->draw(GL_TRIANGLES);
    glSprite->drawPoints(10.0);


    sh = shBranches;
    sh->use();
    sh->set_camPos  ( (GLfloat*)&camPos );
    sh->set_camMat  ( (GLfloat*)&camMat );
    sh->set_modelMat( (GLfloat*)&modelMat );

    //camMat.m
    Vec3f mpos = (Vec3f){uv.x*16-4+0.5,uv.y*16-16+0.5,10.0};
    sh->set_modelPos( (GLfloat*)&mpos );
    //mshBranches->draw(GL_TRIANGLES);
    msh1->draw(GL_TRIANGLES);



}

void init();
void quit();
void die( char const *msg );
void inputHanding();

int main(int argc, char *argv[]){
    init();
	setup();
    for ( frameCount=1; frameCount<1000000; frameCount++)    {

        draw(); SDL_GL_SwapWindow(window);
 		//if( !STOP ) draw();
		inputHanding();
        SDL_Delay(delay);
    }
    quit();
    return 0;
}

// FUNCTION ======	inputHanding
void inputHanding(){

    //float posstep = 0.1f; if(RayTerrain){ posstep = 2.0f; }
    float step          = 0.1f;
    float keyRotSpeed   = 0.01f;

    const Uint8 *keys = SDL_GetKeyboardState(NULL);

    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.roll2  (  (float)keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.roll2  ( -(float)keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch2(  (float)keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch2( -(float)keyRotSpeed ); }

    if( keys[ SDL_SCANCODE_Q  ] ){ camPos.add_mul( mouseMat.c, +step ); }
	if( keys[ SDL_SCANCODE_E  ] ){ camPos.add_mul( mouseMat.c, -step ); }
    if( keys[ SDL_SCANCODE_W  ] ){ camPos.add_mul( mouseMat.b, +step ); }
	if( keys[ SDL_SCANCODE_S  ] ){ camPos.add_mul( mouseMat.b, -step ); }
	if( keys[ SDL_SCANCODE_A  ] ){ camPos.add_mul( mouseMat.a, -step ); }
	if( keys[ SDL_SCANCODE_D  ] ){ camPos.add_mul( mouseMat.a, +step ); }

	//if( keys[SDL_SCANCODE_W]||keys[SDL_SCANCODE_S]||keys[SDL_SCANCODE_A]||keys[ SDL_SCANCODE_D  ]  )printf( "camPos (%g,%g,%g)\n", camPos.x, camPos.y, camPos.z );

	SDL_Event event;
	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE: quit(); break;
                case SDLK_t: bTransparent   =!bTransparent;    break;
                //case SDLK_KP_PLUS:  terrain_size[0] *=1.1; terrain_size[2] *=1.1; break;
                //case SDLK_KP_MINUS: terrain_size[0] /=1.1; terrain_size[2] /=1.1; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            printf( "" );
		}
		if( event.type == SDL_QUIT){ quit();  };
	}

	int dmx,dmy;
	SDL_GetMouseState( &mouseX, &mouseY );
    Uint32 buttons = SDL_GetRelativeMouseState( &dmx, &dmy);
    //printf( " %i %i \n", mx,my );
    float mouseRotSpeed = 0.002;
    if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
        Quat4f q; q.fromTrackball( 0, 0, -dmx*mouseRotSpeed, dmy*mouseRotSpeed ); qCamera.qmul_T( q );
        //qCamera.dyaw2(-dmx*mouseRotSpeed); qCamera.dpitch2(-dmy*mouseRotSpeed);
        //qCamera.dpitch2(-dmy*mouseRotSpeed); qCamera.dyaw2(-dmx*mouseRotSpeed);
        //qCamera.normalize();

        //pitch +=  dmy*mouseRotSpeed;
        //yaw   +=  dmx*mouseRotSpeed;
        qCamera.toMatrix(mouseMat);
        printf("mouseMat:\n");
        mouseMat.print();
    }
}

void init(){
    if (SDL_Init(SDL_INIT_VIDEO) < 0) die( "Unable to initialize SDL" );
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3); // Opengl 3.2
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3); // Opengl 3.2
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    //SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 32);
    window = SDL_CreateWindow("Tutorial2", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
    if ( !window ) die("Unable to create window");
    context = SDL_GL_CreateContext( window );
    //SDL_GL_SetSwapInterval(1); // VSync On
    SDL_GL_SetSwapInterval(VSync);

    glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		getchar();
		quit();
		//return -1;
	}

	// vertex array object
	glGenVertexArrays(1, &vao);  				// Allocate and assign a Vertex Array Object to our handle
	glBindVertexArray(vao); 					// Bind our Vertex Array Object as the current used object
}

void quit(){
	glDeleteVertexArrays(1, &vao);
    if( context != NULL ) SDL_GL_DeleteContext( context );
    if( window  != NULL ) SDL_DestroyWindow   ( window  );
    SDL_Quit();
	exit(0);
};

void die( char const *msg ){
    printf("%s: %s\n", msg, SDL_GetError());
    SDL_Quit();
    exit(1);
}
