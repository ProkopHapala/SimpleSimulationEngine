#ifndef  MusicRendererOGL3_h
#define  MusicRendererOGL3_h

#include "containers.h"

#include <GL/glew.h>

//#include "Mesh.h"
//#include "Solids.h"
#include "GLfunctions.h"
#include "GLobjects.h"
#include "GLObject.h"
#include "Shader.h"
#include "DrawOGL3.h"

#include "Noise.h"

#include "SDL_utils.h"
//#include "IO_utils.h"


/*

TODO : Share Camera between many shaders  (and other uniforms as well)

Sharing Uniforms (like Camera) between many shaders
Uniform Buffer Objects (OpenGL)
https://gamedev.stackexchange.com/questions/153896/camera-and-multiple-shaders

*/


void textureFillRandomRGB(int W, int H, uint tex ){
    //int W = layers.buffers[0]->W;
    //int H = layers.buffers[0]->H;
    const int nbuf = W*H;
    uint8_t buff[nbuf*3];
    //if(frameCount==1){
    for(int i=0;i<nbuf*3;i++) buff[i]=rand()&0xFF;
    //glActiveTexture(GL_TEXTURE0);
    glBindTexture  (GL_TEXTURE_2D, tex );
    glTexSubImage2D(GL_TEXTURE_2D,0,0,0,W,H,GL_RGB,GL_UNSIGNED_BYTE,buff);
    //glFlush();
    //}
}


void renderTexture( Shader* sh, GLMesh* mesh, Camera& cam, uint tx, int iTxUnit=0){
    sh->use();
    glActiveTexture(GL_TEXTURE0 + iTxUnit);
    glBindTexture  (GL_TEXTURE_2D, tx   );
    setCamera(*sh, cam );
    sh->setModelPoseT( (Vec3d){-4.,-4.,0.0}, Mat3dIdentity*8.0 );
    mesh->draw();
}

void renderTexture( GLMesh* mesh, uint tx, int iTxUnit=0){
    glActiveTexture(GL_TEXTURE0 + iTxUnit);
    glBindTexture  (GL_TEXTURE_2D, tx   );
    mesh->draw();
}


void plotBuff( GLMesh& mesh, int n, double* buff, float dx, float dy ){
    Vec3f ps[3*n];
    for(int i=0; i<n; i++){ ps[i].set( i*dx, buff[i]*dy, 0 ); }
    glBindBuffer(GL_ARRAY_BUFFER, mesh.vpos );
    glBufferSubData(GL_ARRAY_BUFFER, 0, 3*n*sizeof(float), ps   );
    mesh.draw();
}

void plotBuffStereo( GLMesh& mesh, Shader& sh, int n, double* buff, float dx, float dy ){
    Vec3f ps[3*n];

    glBindBuffer(GL_ARRAY_BUFFER, mesh.vpos );

    GLuint ucolor = sh.getUloc("baseColor");

    //printf( "plotBuffStereo n %i \n", n );
    for(int i=0; i<n; i++){ ps[i].set( i*dx, buff[i*2]*dy, 0 ); }
    glBufferSubData(GL_ARRAY_BUFFER, 0, 3*n*sizeof(float), ps   );
    //sh.setUniformVec4f( "baseColor", (Quat4f){1.f,0.f,0.f,1.f} );
    glUniform4f( ucolor, 0.0f, 0.0f, 1.0f, 1.0f );
    mesh.draw();

    for(int i=0; i<n; i++){ ps[i].set( i*dx, buff[i*2+1]*dy, 0 ); }
    glBufferSubData(GL_ARRAY_BUFFER, 0, 3*n*sizeof(float), ps   );
    //sh.setUniformVec4f( "baseColor", (Quat4f){0.f,0.f,1.f,1.f} );
    glUniform4f( ucolor, 1.0f, 0.0f, 0.0f, 1.0f );
    mesh.draw();
}


struct RenderPassRecord{
    int shader;    // index of program
    int output;    // output textures (FrameBuffers)
    int inputs[4]; // input  textures
    int nin=0;
    // Transforms ?
};

struct RenderScript{
    int npass=0;
    RenderPassRecord* passes=0;
};

class RenderStack{ public:
    // This class let you render use multiple FrameBuffers as input-output textures with multople shaders like multiple buffers on ShaderToy
    // This is usefull for grid-based physical simulations on GPU (like fluid simulations) or like texture dynamic wraping/transforms like movement module in WinAmp AVS

    std::vector<FrameBuffer*> buffers;
    std::vector<Shader*>      shaders;

    GLMesh* screenQuad=0;

    //bool bDrawRaw = true;
    bool bDrawRaw  = false;
    bool bFlushing = false;

    void makeBuffers(int n, int width, int height ){
        for(int i=0; i<n; i++){
            buffers.push_back( new FrameBuffer(width, height, true) );
        }
    }

    void bindOutput(int ibuf){
        if( (ibuf<0)||(ibuf>buffers.size()) ) { glBindFramebuffer(GL_FRAMEBUFFER, 0                );  }
        else                                  { glBindFramebuffer(GL_FRAMEBUFFER, buffers[ibuf]->buff );  }
    };
    void unbindOutput(){ glBindFramebuffer(GL_FRAMEBUFFER, 0); }

    void bindInput(int ibuf, int islot=0){
        glActiveTexture(GL_TEXTURE0+islot);
        if(ibuf<0){ glBindTexture  (GL_TEXTURE_2D, buffers[-ibuf]->texZ   ); } // for negative buffer index we take Z-buffer
        else      { glBindTexture  (GL_TEXTURE_2D, buffers[ ibuf]->texRGB ); } // for positive buffer index we take RGB-buffer
    }

    Shader* render( int ish, int iout, int nin, const int* ins, const int* texUnits=0 ){
        if(bFlushing)glFlush();
        //glBindFramebuffer(GL_FRAMEBUFFER, iout );
        bindOutput(iout);
        for(int i=0; i<nin; i++){
            //printf( "i %i nin %i \n", i, nin );
            int islot = i;
            if(texUnits) islot = texUnits[i];
            //glActiveTexture(GL_TEXTURE0+itex);
            //int ibuf = ins[i];
            //if(ibuf<0){ glBindTexture  (GL_TEXTURE_2D, buffers[-ibuf]->texZ   ); } // for negative buffer index we take Z-buffer
            //else      { glBindTexture  (GL_TEXTURE_2D, buffers[ ibuf]->texRGB ); } // for positive buffer index we take RGB-buffer
            bindInput( ins[i], islot=islot);
        }
        Shader* sh = shaders[ish];
        sh->use();
        // ToDo - would be nice to set shader parameters here
        if(bDrawRaw){ screenQuad->drawRaw(); }else{ screenQuad->draw(); }
        return sh;
    }

    Shader* render( const  RenderPassRecord& rc  ){ return render( rc.shader, rc.output, rc.nin, rc.inputs ); }
    void render   ( int n, RenderPassRecord* rcs ){ for(int i=0;i<n;i++){ render(rcs[i]); }; };

    void fillRandomRGB(int ibuf){ textureFillRandomRGB( buffers[ibuf]->W, buffers[ibuf]->H, buffers[ibuf]->texRGB ); }
};


class RenderStackManager{ public:
    RenderStack* layers=0;
    Dictionary   shaderDict;
    Dictionary   scriptDict;
    std::vector<RenderScript> scripts;
    std::vector<std::string> shaderNames;

    int maxTextures=0;

    void loadShader( const char* vertName, const char* fragName, bool bFile = true ){
        //Shader* sh = new Shader( "common_resources/shaders/texture3D.glslv" , "common_resources/shaders/texture.glslf" , true );
        Shader* sh = new Shader(vertName , fragName , bFile );
        layers->shaders.push_back(sh);
    }

    Shader* getShader(std::string name){
        auto got  =  shaderDict.find (name);
        if  ( got == shaderDict.end() ){ return 0; }
        else                           { return layers->shaders[got->second]; }
    }

    /*
    void parseScriptLine(char* s, char* name, int& output, int* inputs ){
        int in0,in1;
        int i=0;
        while(s[i]==' ')i++; in0-i;
        while(s[i]!=' ')i++; in1=i;
        while(s[i]!='=')i++; in1=i;
    }
    */

    void loadShaders(){
        const char prefix[]{"common_resources/shaders/"};
        char vertPath[1024];
        for(int i=0; i<shaderNames.size(); i++){
            sprintf(vertPath,"%s%s",prefix, shaderNames[i].c_str() );
            printf( "Loading %s \n", vertPath );
            loadShader( "common_resources/shaders/texture3D.glslv", vertPath );
        }
    }


    void prepare( int width, int height ){
        int nbuf = maxTextures - layers->buffers.size();
        layers->makeBuffers( nbuf, width, height );
        loadShaders();
    };

    void addScriptLine(char* s){
        // parse script line
        RenderPassRecord rc;
        char cname[64];
        int nread = sscanf( s, "%s %i = %i %i %i %i ", cname, rc.output, rc.inputs[0], rc.inputs[1], rc.inputs[2], rc.inputs[3] );
        if(nread<2){ printf( "ERORR : render pass needs to specify shader name & output at least ! \n %s", s ); exit(0); }
        rc.nin = nread - 2;
        // update required number of textures
        _setmax( maxTextures, rc.output );
        for(int i=0; i<rc.nin; i++){ _setmax( maxTextures, rc.inputs[i] ); }
        // register new shader name
        std::string name(cname);
        //int ish;
        auto got  =  shaderDict.find (name);
        if  ( got == shaderDict.end() ){ rc.shader = shaderDict.size(); shaderDict.insert( {name,rc.shader} ); shaderNames.push_back(name); }
        else                           { rc.shader = got->second; }
    }

    /*
    Script looks like:
    FluidPDE   2 = 1    // render FluidPDE from texture[1] to texture[2]
    FluidPDE   1 = 2    // render FluidPDE from texture[2] to texture[1] // Flip-Flopping
    FluidDrift 4 = 3 1  // render FluidDrift to texture[4] using texture[3] and texture[1] as input
    FluidPDE   2 = 1    //
    FluidPDE   1 = 2    //
    FluidDrift 3 = 4  1 //
    tx         0 = 3    // output to screen
    */
    void loadScriptString( char* str, const char* sep="\n"){
       char* token;
       token = strtok(str, sep);
       while( token != NULL ) {
          addScriptLine( token );
          token = strtok(NULL, sep);
       }
    }

    void loadScriptFile( const char* fname, char sep='\n'){
        /*
        FILE * pFile;
        const int nbuff = 4096;
        char str[nbuff];
        pFile = fopen ( fname , "r");
        if (pFile == NULL){ printf("ERROR : RenderStackManager: script file not found: %s \n", fname ); return(-1); };
        int n=0;
        while ( fgets( str , nbuff, pFile) != NULL ){
            if (str[0]=='#') continue;
            addScriptLine( str );
            n++;
        }
        fclose(pFile);
        return n;
        */
        processFileLines( fname, [this](char* line){ addScriptLine(line); } );
    }


};


class ParticleFlow{ public:
    int np=0;
    Vec2d* ps=0; // ToDo : we may do in also in 3D ( using procedural noise )
    Vec2d* vs=0;
    double bmix = 0.2;
    double resetProb = 0.01;
    Vec2d pmin=Vec2dZero;
    Vec2d pmax=Vec2dOnes;

    void init(Vec2d pmin_, Vec2d pmax_){
        pmin=pmin_;
        pmax=pmax_;
        for(int i=0; i<np; i++){
            ps[i].set( randf(pmin.x,pmax.x), randf(pmin.y,pmax.y) );
            vs[i] = Vec2dZero;
        }
    }

    void realloc(int np_){
        np=np_;
        _realloc( ps, np );
        _realloc( vs, np );
    }

    inline void flowField(const Vec2d& p, Vec2d& v){
        Noise::simplexNoise2D( p*25.0, v );
        //Noise::warpNoise3R( const Vec2d& pos0, const Vec2d& rot, double fdown, double strenght, int n, Vec2d& dpos_ );
    }

    void update(double dt){
        Vec2d vel;
        for(int i=0; i<np; i++){
            if(resetProb>1e-12){
                if(randf()<resetProb){
                    ps[i].set( randf(pmin.x,pmax.x), randf(pmin.y,pmax.y) );
                    vs[i] = Vec2dZero;
                }
            }
            Vec2d pos = ps[i];
            flowField(pos,vel);
            vel = vel*bmix + vs[i]*(1-bmix);
            pos.add_mul(vel,dt);
            ps[i] = pos;
            vs[i] = vel;
        }
    }

};


class ParticleFlowRender{ public:
    int np=0;
    Vec2d*  ps=0;
    Vec3f*  ops=0;
    GLMesh* mesh=0;
    Vec2f p0=Vec2fZero;
    Vec2f sc=Vec2fOnes;
    double max_dist = 0.1;

    void init(int np_, Vec2d* ps_, GLMesh* mesh_ = 0){
        np=np_;
        ps=ps_;
        _realloc(ops, 2*np);
        mesh=mesh_;
        if(mesh==0){
            if(mesh==0) mesh=new GLMesh();
            mesh->init( np*2, 0, 0, ops, 0, NULL, NULL, GL_STREAM_DRAW );
            mesh->draw_mode = GL_LINES;
        }
    }

    void fillBuff(){
        for(int i=0; i<np; i++){
            int i2=i*2;
            Vec3f op = ops[i2+1];
            Vec3f p  = (Vec3f){ ps[i].x*sc.x+p0.x, ps[i].y*sc.y+p0.y, 0. };
            if( p.dist2(op)>(max_dist*max_dist) ){ op=p; };
            ops[i2  ] = op;
            ops[i2+1] = p;
        }
    };

    void plotBuff(){

        fillBuff();

        glBindBuffer(GL_ARRAY_BUFFER, mesh->vpos );
        glBufferSubData(GL_ARRAY_BUFFER, 0, 6*np*sizeof(float), ops );
        mesh->draw();
    }

};

#endif
