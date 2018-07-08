
#ifndef  GLUtils_h
#define  GLUtils_h

#include <GL/gl.h>
#include "Vec3.h"

inline float * toFloat(int n, double * arr_){
	float * arr = new GLfloat[n];
	for(int i=0; i<n; i++){ arr[i] = arr_[i];  }
	return arr;
}

inline int countVerts( int n, int * ngons ){
    int nVert = 0;
    for(int i=0; i<n; i++){ nVert += 3*(ngons[i]-2); }
    return nVert;
}

inline void hardFace( int n, int * ngons, int * faces, Vec3d * points, GLfloat* verts, GLfloat* normals ){
    int * facei = faces;
    GLfloat * verti = verts;
    GLfloat * normi = normals;
    for(int i=0; i<n; i++){
        int ni = ngons[i];
        Vec3d a = points[facei[0]];
        Vec3d b = points[facei[1]];
        //printf( " i %i ni %i \n", i, ni );
        for(int j=2; j<ni; j++){
            Vec3d c = points[facei[j]];
            Vec3d n; n.set_cross(a-c,b-a);
            double l = n.normalize();
            //printf( "%i %i %g \n", i, j, l );
            verti[0]=(GLfloat)a.x; verti[1]=(GLfloat)a.y; verti[2]=(GLfloat)a.z;
            verti[3]=(GLfloat)b.x; verti[4]=(GLfloat)b.y; verti[5]=(GLfloat)b.z;
            verti[6]=(GLfloat)c.x; verti[7]=(GLfloat)c.y; verti[8]=(GLfloat)c.z;
            normi[0]=(GLfloat)n.x; normi[1]=(GLfloat)n.y; normi[2]=(GLfloat)n.z;
            normi[3]=(GLfloat)n.x; normi[4]=(GLfloat)n.y; normi[5]=(GLfloat)n.z;
            normi[6]=(GLfloat)n.x; normi[7]=(GLfloat)n.y; normi[8]=(GLfloat)n.z;
            b = c;
            verti+=9;
            normi+=9;
        }
        facei += ni;
    }
}

static constexpr GLuint  QUAD_nVert = 4;
GLfloat QUAD_vertexes[4][2] = {
	{  -0.9f,  -0.9f  },
	{  -0.9f,   0.9f  },
	{   0.9f,  -0.9f  },
	{   0.9f,   0.9f  } };

#endif
