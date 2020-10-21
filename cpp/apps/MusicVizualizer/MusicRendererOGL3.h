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




#endif
