#include "MeshRenderOGL3.h"

#include "IO_utils.h"

void MeshRenderOGL3::initDefaultShaders(){
    const char* namesShade[] = { "shade3D.vert", "shade3D.frag" };
    char** srcShade = fileGetSections( "common_resources/shaders/Basic.glslf", 2, (char**)namesShade, (char*)"//>>" );
    sh_solid.init_str( srcShade[0], srcShade[1], nullptr );
    sh_solid.getDefaultUniformLocation();

    sh_const.init_default();
    sh_const.getDefaultUniformLocation();

    modelPos.set(0.0f,0.0f,0.0f);
    modelMat.setOne();
}

void MeshRenderOGL3::setModelPos(const Vec3f& pos){
    modelPos = pos;
}

void MeshRenderOGL3::setModelMat(const Mat3f& mat){
    modelMat = mat;
}

void MeshRenderOGL3::uploadMesh_d( int nVerts, int nTris, const int* tris, const double* verts, const double* nors ){
    if(!mesh_tri) mesh_tri = new GLMesh();
    mesh_tri->init_d( nVerts, nTris*3, (int*)tris, (double*)verts, (double*)nors, nullptr, nullptr );
    mesh_tri->draw_mode = GL_TRIANGLES;
}

void MeshRenderOGL3::uploadLines_d( int nVerts, int nEdges, const int* edges, const double* verts ){
    if(!mesh_lines) mesh_lines = new GLMesh();
    mesh_lines->init_d( nVerts, nEdges*2, (int*)edges, (double*)verts, nullptr, nullptr, nullptr );
    mesh_lines->draw_mode = GL_LINES;
}

void MeshRenderOGL3::draw(const Camera& cam){
    // Build camera matrix from Camera parameters (similar to setCameraPersp in SceneOGL3.h)
    Mat4f camMat, mRot, mPersp;
    if(cam.persp){
        mPersp.setPerspective(cam.aspect*cam.zoom, cam.zoom, cam.zmin, cam.zmax);
        mRot.setOne();
        mRot.setRot(cam.rot);
        camMat.set_mmul_TN(mRot, mPersp);
    }else{
        camMat.setOrthographic(cam.zoom, cam.zoom*cam.aspect, cam.zmin, cam.zmax);
    }

    // Debug: print camera matrix
    printf("MeshRenderOGL3::draw() camMat: \n"); camMat.print();

    if(mesh_tri){
        sh_solid.use();
        sh_solid.set_modelPos( (GLfloat*)&modelPos );
        sh_solid.set_modelMat( (GLfloat*)&modelMat );
        sh_solid.set_camPos  ( (GLfloat*)&cam.pos );
        sh_solid.set_camMat  ( (GLfloat*)&camMat );
        mesh_tri->draw();
    }
    if(mesh_lines){
        sh_const.use();
        sh_const.set_modelPos( (GLfloat*)&modelPos );
        sh_const.set_modelMat( (GLfloat*)&modelMat );
        sh_const.set_camPos  ( (GLfloat*)&cam.pos );
        sh_const.set_camMat  ( (GLfloat*)&camMat );
        const float colRed[4]   = {1.0f,0.0f,0.0f,1.0f};
        const float colGreen[4] = {0.0f,1.0f,0.0f,1.0f};
        sh_const.set_baseColor( colGreen );
        glLineWidth(3.0f);
        mesh_lines->draw();
    }
}
