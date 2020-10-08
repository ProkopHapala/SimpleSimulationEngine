#ifndef  MusicRendererOGL3_h
#define  MusicRendererOGL3_h

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



class ShaderStack{ public:
    // This class let you render use multiple FrameBuffers as input-output textures with multople shaders like multiple buffers on ShaderToy
    // This is usefull for grid-based physical simulations on GPU (like fluid simulations) or like texture dynamic wraping/transforms like movement module in WinAmp AVS

    std::vector<FrameBuffer*> buffers;
    std::vector<Shader*>      shaders;
    GLMesh* screenQuad=0;

    //bool bDrawRaw = true;
    bool bDrawRaw = false;

    void makeBuffers(int n, int width, int height ){
        for(int i=0; i<n; i++){
            buffers.push_back( new FrameBuffer(width, height) );
        }
    }

    void bindOutput(int i){
        if( (i<0)||(i>buffers.size()) ) { glBindFramebuffer(GL_FRAMEBUFFER, 0               );  }
        else                            { glBindFramebuffer(GL_FRAMEBUFFER, buffers[i]->buff );  }
    };
    void unbindOutput(){ glBindFramebuffer(GL_FRAMEBUFFER, 0); }

    Shader* render( int ish, int iout, int nin, const int* ins, const int* texUnits=0 ){
        //glBindFramebuffer(GL_FRAMEBUFFER, iout );
        bindOutput(iout);
        for(int i=0; i<nin; i++){
            printf( "i %i nin %i \n", i, nin );
            int itex = i;
            if(texUnits) itex = texUnits[i];
            glActiveTexture(GL_TEXTURE0+i);
            int ibuf = ins[i];
            if(ibuf<0){ glBindTexture  (GL_TEXTURE_2D, buffers[-ibuf]->texZ   ); } // for negative buffer index we take Z-buffer
            else      { glBindTexture  (GL_TEXTURE_2D, buffers[ ibuf]->texRGB ); } // for positive buffer index we take RGB-buffer
        }
        Shader* sh = shaders[ish];
        sh->use();
        // ToDo - would be nice to set shader parameters here
        if(bDrawRaw){ screenQuad->drawRaw(); }else{ screenQuad->draw(); }
        return sh;
    }

};



#endif
