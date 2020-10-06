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
    for(int i=0; i<n/2; i++){ ps[i].set( i*dx, buff[i*2]*dy, 0 ); }

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





#endif
