#ifndef _GLMesh_H_
#define _GLMesh_H_

#include <utility>
#include "GLattribs.h"
#include "GLuniform.h"
#include "Shader.h"
#include "Camera.h"
#include "quaternion.h"
#include "GLvbo.h"

template<attrib ... attribs>
class GLMeshBase{
private:
static_assert(GLattrib::check_attribs<attribs...>(), "ERROR: attribute list cannot contain duplicate names.");
public:
    using vertex = GLvertex<typename decltype(attribs)::type ...>;
    using attrIdxSeq = std::make_index_sequence<sizeof...(attribs)>;
    using attr_monostate = attribs_monostate<attribs...>;

private:
    GLuint vao = 0;

    inline void ensure_vao(){
        if (vao) return;
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);
        
        verts->setup_vertex_attribs(
            [this](GLattrib::Name name){return shader->attrName2Loc(name);}
        );

        glBindVertexArray(0);
    }

public:
    GLenum drawMode;
    GLvbo<attribs...>* verts;
    Shader* shader;
    GLuniformSet uniforms;

    GLMeshBase(GLenum drawMode=GL_TRIANGLES, GLenum usage=GL_STATIC_DRAW, Shader* shader=defaultcolorShader<attribs...>)
        : drawMode(drawMode), shader(shader), verts(new GLvbo<attribs...>(usage)) {}
    GLMeshBase(GLvbo<attribs...>* vbo, GLenum drawMode=GL_TRIANGLES, Shader* shader=defaultcolorShader<attribs...>)
        : drawMode(drawMode), shader(shader), verts(vbo) {}

    ~GLMeshBase(){
        if (vao) glDeleteVertexArrays(1, &vao);
        vao = 0;
    }

    void clear(){ verts->clear(); }

    void addVertex(typename decltype(attribs)::type ... args){ verts->push_back(args...); }
    void addVertex_strip(typename decltype(attribs)::type ... args){
        vertex v = vertex(args...);

        if (vertexCount() >= 3){
            vertex v1 = verts[verts->size()-1];
            vertex v2 = verts[verts->size()-2];

            verts->push_back(v2);
            verts->push_back(v1);
        }

        verts->push_back(v);
    }
    inline int vertexCount() const { return verts->size(); }

    inline void setUniformName(std::string name, GLuniform u){ uniforms.set(name, u); }
    void setUniform1f(std::string name, GLfloat v)        {setUniformName(name, {.type=GLuniform::f1, .data={.f1=v}}); }
    void setUniform2f(std::string name, Vec2T<GLfloat> v) {setUniformName(name, {.type=GLuniform::f2, .data={.f2=v}}); }
    void setUniform3f(std::string name, Vec3T<GLfloat> v) {setUniformName(name, {.type=GLuniform::f3, .data={.f3=v}}); }
    void setUniform4f(std::string name, Vec4T<GLfloat> v) {setUniformName(name, {.type=GLuniform::f4, .data={.f4=v}}); }
    void setUniform1i(std::string name, GLint v)          {setUniformName(name, {.type=GLuniform::i1, .data={.i1=v}}); }
    void setUniform2i(std::string name, Vec2T<GLint> v)   {setUniformName(name, {.type=GLuniform::i2, .data={.i2=v}}); }
    void setUniform3i(std::string name, Vec3T<GLint> v)   {setUniformName(name, {.type=GLuniform::i3, .data={.i3=v}}); }
    void setUniform4i(std::string name, Vec4T<GLint> v)   {setUniformName(name, {.type=GLuniform::i4, .data={.i4=v}}); }
    void setUniformMatrix3f(std::string name, Mat3T<GLfloat> v) {setUniformName(name, {.type=GLuniform::m3, .data={.m3=v}}); }
    void setUniformMatrix4f(std::string name, Mat4T<GLfloat> v) {setUniformName(name, {.type=GLuniform::m4, .data={.m4=v}}); }  
    void setUniformTex(std::string name, GLTexture* tex) {setUniformName(name, {.type=GLuniform::tex, .data={.tex=tex}});} 

    inline void draw(GLenum drawMode=0){
        if (drawMode == 0) drawMode = this->drawMode;
        ensure_vao();

        glBindVertexArray(vao);
        verts->sync();

        shader->setUniforms(uniforms);
        shader->use();

        glDrawArrays(drawMode, 0, vertexCount());
        glBindVertexArray(0);
    }
};


template<attrib...attribs>
class GLMesh : public GLMeshBase<attribs...>{
public:
    using base = GLMeshBase<attribs...>;
    GLMesh(GLenum drawMode=GL_TRIANGLES, GLenum usage=GL_STATIC_DRAW, Shader* shader=defaultcolorShader<attribs...>)
        : base(drawMode, usage, shader) {}
    GLMesh(GLvbo<attribs...>* vbo, GLenum drawMode=GL_TRIANGLES, Shader* shader=defaultcolorShader<attribs...>)
        : base(vbo, drawMode, shader) {}
    GLMesh(Shader* shader)
        : base(GL_TRIANGLES, GL_STATIC_DRAW, shader) {}

    void drawMVP(Mat4f mvp, GLenum drawMode=0){
        base::uniforms.set4m("uMVPMatrix", mvp);
        base::draw(drawMode);
    }

    inline void draw(GLenum drawMode){draw(Vec3fZero, Vec3fOne, drawMode);}
    inline void draw(Vec3f position, float scale, GLenum drawMode=0){draw(position, (Vec3f){scale, scale, scale}, drawMode);}
    void draw(Vec3f position=Vec3fZero, Vec3f scale={1, 1, 1}, GLenum drawMode=0){
        if (GLES::active_camera == nullptr){
            printf("Warning: No active camera - skipping rendering.\n");
            return;
        }
    
        // scale + translation
        Mat4f mvpMatrix;
        mvpMatrix.setDiag(scale.x, scale.y, scale.z, 1);
        mvpMatrix.setPos(position);

        // view + projection
        mvpMatrix.mmulL(GLES::active_camera->viewProjectionMatrix());
    
        drawMVP(mvpMatrix, drawMode);
    }

    void draw(Vec3f position, Quat4f rotation, Vec3f scale={1, 1, 1}, GLenum drawMode=0){
        if (GLES::active_camera == nullptr){
            printf("Warning: No active camera - skipping rendering.\n");
            return;
        }

        // TODO: optimize this
    
        // scale + translation
        Mat4f ScaleTranslationMatrix;
        ScaleTranslationMatrix.setDiag(scale.x, scale.y, scale.z, 1);
        ScaleTranslationMatrix.setPos(position);
        
        // rotation
        Mat3f modelRotMatrix3;
        rotation.toMatrix(modelRotMatrix3);
        Mat4f modelRotMatrix4 = Mat4fIdentity;
        modelRotMatrix4.setRot(modelRotMatrix3);
    
        // model matrix
        Mat4f modelMatrix = modelRotMatrix4;
        modelMatrix.mmulL(ScaleTranslationMatrix);

        // MVP matrix
        Mat4f mvpMatrix = modelMatrix;
        mvpMatrix.mmulL(GLES::active_camera->viewMatrix());
        mvpMatrix.mmulL(GLES::active_camera->projectionMatrix());
    
        drawMVP(mvpMatrix, drawMode);
    }

    inline void draw2D(Vec3f pos, float scale, GLenum drawMode=0){draw2D(pos, (Vec2f){scale, scale}, drawMode);}
    void draw2D(Vec3f pos=Vec3fZero, Vec2f scale={1, 1}, GLenum drawMode=0){
        // convert from screen space ((0, 0)  to (WIDHT, HEIGHT)) to NDC ((-1, -1) to (1, 1))

        const int WIDTH = GLES::screen_size.x;
        const int HEIGHT = GLES::screen_size.y;

        pos.x = pos.x*2/WIDTH - 1;
        pos.y = pos.y*2/HEIGHT - 1;
    
        scale.x = scale.x*2/WIDTH;
        scale.y = scale.y*2/HEIGHT;

        draw2D_NDC(pos, scale, drawMode);
    }

    void draw2D_NDC(Vec3f pos=Vec3fZero, Vec2f scale={1, 1}, GLenum drawMode=0){
        Mat4f mvp;
        mvp.array[ 0] = scale.x; mvp.array[ 4] = 0;       mvp.array[ 8] = 0; mvp.array[12] = pos.x;
        mvp.array[ 1] = 0;       mvp.array[ 5] = scale.y; mvp.array[ 9] = 0; mvp.array[13] = pos.y;
        mvp.array[ 2] = 0;       mvp.array[ 6] = 0;       mvp.array[10] = 1; mvp.array[14] = pos.z;
        mvp.array[ 3] = 0;       mvp.array[ 7] = 0;       mvp.array[11] = 0; mvp.array[15] = 1;

        drawMVP(mvp, drawMode);
    }
};


#endif // _GLMesh_H_

