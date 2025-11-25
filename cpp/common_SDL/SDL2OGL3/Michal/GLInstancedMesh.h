#ifndef GLINSTANCEDMESH_H
#define GLINSTANCEDMESH_H

#include "GLES.h"
#include "GLattribs.h"
#include "GLvbo.h"
#include "Shader.h"
#include "Camera.h"

template <class T, template<attrib...> class Template>
struct is_specialization : std::false_type {};

template <template<attrib...> class Template, attrib... Args>
struct is_specialization<Template<Args...>, Template> : std::true_type {};

template <class T, template<attrib...> class Template>
constexpr bool is_specialization_v = is_specialization<T, Template>::value;








template<class VertVboType, attrib...instAttribs>
class GLInstancedMeshBase {
public:
    using vertex = GLvertex<typename decltype(instAttribs)::type ...>;
    using attrIdxSeq = std::make_index_sequence<sizeof...(instAttribs)>;
    using attr_monostate = attribs_monostate<instAttribs...>;

private:
    static_assert(is_specialization_v<VertVboType, GLvbo>, "VertVboType must be specialization of GLvbo");
    static_assert(GLattrib::check_attribs<instAttribs...>(), "Attribute list cannot contain duplicate names.");
    static_assert(GLattrib::check_attribs(typename VertVboType::attr_monostate{}, attr_monostate{}), "vertex attribs and instAttribs must have disjunct attribute names.");


    template<attrib...vertAttribs>
    static constexpr Shader* get_default_shader(attribs_monostate<vertAttribs...>){
        return defaultShader<vertAttribs..., instAttribs...>;
    }

    GLuint vao = 0;

    inline void ensure_vao(){
        if (vao) return;
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);
        
        verts->setup_vertex_attribs(
            [this](GLattrib::Name name){return shader->attrName2Loc(name);}
        );
        instances->setup_vertex_attribs(
            [this](GLattrib::Name name){return shader->attrName2Loc(name);},
            1
        );

        glBindVertexArray(0);
    }

public:
    GLenum drawMode = GL_TRIANGLES;
    VertVboType* verts;
    GLvbo<instAttribs...>* instances;
    Shader* shader;
    GLuniformSet uniforms;

    GLInstancedMeshBase(VertVboType* vbo, GLenum draw_mode=GL_TRIANGLES, Shader* shader=get_default_shader(typename VertVboType::attr_monostate{}))
        : shader(shader), verts(vbo), drawMode(draw_mode), instances(new GLvbo<instAttribs...>()) {}
    GLInstancedMeshBase(GLenum draw_mode=GL_TRIANGLES, Shader* shader=get_default_shader(typename VertVboType::attr_monostate{}))
        : shader(shader), verts(new VertVboType()), drawMode(draw_mode), instances(new GLvbo<instAttribs...>()) {}

    void addInstance(typename decltype(instAttribs)::type ... args){
        instances->push_back(vertex(args...));
    }

    inline void draw(GLenum drawMode=0){
        if (drawMode == 0) drawMode = this->drawMode;
        if (GLES::active_camera == nullptr){
            printf("Warning: No active camera - skipping rendering.\n");
            return;
        }

        ensure_vao();
        glBindVertexArray(vao);
        
        verts->sync();
        instances->sync();

        shader->setUniforms(uniforms);
        shader->setUniform4m("uMVPMatrix", GLES::active_camera->viewProjectionMatrix());
        shader->use();
        glDrawArraysInstanced(drawMode, 0, verts->size(), instances->size());

        glBindVertexArray(0);
    }
};



#endif // GLINSTANCEDMESH_H
