#ifndef _GLvbo_H_
#define _GLvbo_H_

#include "GLES.h"
#include "GLattribs.h"
#include <cstddef>
#include <functional>
#include <utility>

template<attrib ... attribs>
class GLvbo{
    static_assert(GLattrib::check_attribs<attribs...>(), "ERROR: attribute list cannot contain duplicate names.");
    static_assert(sizeof...(attribs) <= 16, "ERROR: OpenGL ES 3.0 supports a maximum of 16 vertex attributes.");
public:
    using vertex = GLvertex<typename decltype(attribs)::type ...>;
    using attrIdxSeq = std::make_index_sequence<sizeof...(attribs)>;
    using attr_monostate = attribs_monostate<attribs...>;
private:
    std::vector<vertex> elements;
    mutable GLuint id = 0;
    mutable bool _sync = false;
    GLenum usage = GL_STATIC_DRAW;

    template<size_t ... attrIdx>
    inline void _setup_vertex_attribs_impl(std::index_sequence<attrIdx...> seq, const std::function<GLint(GLattrib::Name)> attrName2locFunc, GLuint divisor) const {
        ensure_id();

        if (GLES::currentGL_ARRAY_BUFFER == id) return;

        glBindBuffer(GL_ARRAY_BUFFER, id);
        (glEnableVertexAttribArray(attrName2locFunc(attribs.name))
        ,...);
        (glVertexAttribPointer(
            attrName2locFunc(attribs.name),
            GLattrib::type2count <typename decltype(attribs)::type>(),
            GLattrib::type2GLenum<typename decltype(attribs)::type>(),
            GL_FALSE,
            sizeof(vertex),
            (void*)vertex::template get_offset<attrIdx>())
        ,...);
        (glVertexAttribDivisor(
            attrName2locFunc(attribs.name),
            divisor)
        ,...);
    }
    template<GLattrib::Name target_name, size_t target_idx, attrib...As>
    void set_vertex_attr_from_list(vertex& v, const typename decltype(As)::type ... list){
        ((
            (target_name == As.name) ? v.template set<target_idx>(list),0 : 0
        ),...);
    }

    template<attrib...As,size_t ... attrIdx>
    void _push_back_t_impl(std::index_sequence<attrIdx...> seq, const typename decltype(As)::type ... args){
        vertex v;
        (set_vertex_attr_from_list<attribs.name, attrIdx, As...>(v, args...)
        ,...);
        push_back(v);
    }

    inline void sync_raw() const { // assumes that the vbo is bound to GL_ARRAY_BUFFER
        glBufferData(GL_ARRAY_BUFFER, elements.size() * sizeof(vertex), elements.data(), usage);
        _sync = true;
    }

public:
    GLvbo(GLenum usage=GL_STATIC_DRAW) : usage(usage){} 

    inline void ensure_id() const {
        if (id) return;
        glGenBuffers(1, &id);
    }

    inline GLuint get_id() const {return id;}
    inline size_t size()   const {return elements.size();}
    inline const vertex& operator[] (size_t i) const {return elements[i];}

    void push_back(const vertex v){
        elements.push_back(v);
        _sync = false;
    }
    inline void push_back(const typename decltype(attribs)::type ... args){ push_back(vertex(args...)); }

    template<attrib...As>
    void push_back_t(const typename decltype(As)::type ... args){
        static_assert(GLattrib::is_subset(attr_monostate{}, attribs_monostate<As...>{}), "Error: call to 'push_back_t' does not supply all attributes used by this GLvbo.");
        _push_back_t_impl<As...>(attrIdxSeq{}, args...);
    }

    void clear() {elements.clear();}

    void addVertex_strip(typename decltype(attribs)::type ... args){
        vertex v = vertex(args...);

        if (elements.size() >= 3){
            vertex v1 = elements[elements.size()-1];
            vertex v2 = elements[elements.size()-2];

            elements.push_back(v2);
            elements.push_back(v1);
        }

        elements.push_back(v);
    }

    inline void setup_vertex_attribs(const std::function<GLint(GLattrib::Name)> attrName2locFunc, GLuint divisor=0)
        { _setup_vertex_attribs_impl(attrIdxSeq{}, attrName2locFunc, divisor); }

    inline void bind() const{
        ensure_id();
        glBindBuffer(GL_ARRAY_BUFFER, id);
    }
    inline void bind_sync() const{
        bind();
        if (!_sync) sync_raw();
    }
    inline void sync() const {
        if (_sync) return;
        bind();
        sync_raw();
    }
    inline void full_bind(const std::function<GLint(GLattrib::Name)> attrName2locFunc, GLuint divisor=0){ // binds, syncs, and sets up vertex attributes
        setup_vertex_attribs(attrName2locFunc, divisor);
        if (!_sync) sync_raw();
    }


    ~GLvbo(){
        if (id) glDeleteBuffers(1, &id);
        id = 0;
    }
};

#endif // _GLvbo_H_
