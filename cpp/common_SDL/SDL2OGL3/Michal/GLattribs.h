#ifndef _GLATTRIBS_H_
#define _GLATTRIBS_H_

#include "GLES.h"
#include "quaternion.h"
#include <stdexcept>
#include <utility>
#include <cstddef>
#include <variant>

namespace GLattrib{
    enum Name {
        Position  = 0,
        Normal    = 1,
        Color     = 2,
        UV        = 3,
        PosOffset = 4,
        Radius    = 5,
    ATTRIB_NAME_MAX}; // ATTRIB_NAME_MAX must be last
}


template<typename T>
struct attrib{
    using type = T;
    GLattrib::Name name;

    constexpr attrib(GLattrib::Name name) : name(name) {}
};

template<attrib...attribs>
struct attribs_monostate{};


template <typename T, typename ... Ts>
struct GLvertex{
    T first;
    GLvertex<Ts...> next;

    GLvertex(T first, Ts... next) : first(first), next(next...) {}
    GLvertex() {}

    template<size_t i>
    using get_type = std::conditional_t<i==0, T, typename decltype(next)::template get_type<i-1>>;

    template<size_t i> auto& get(){
        if constexpr(i==0) return first;
        else return next.template get<i-1>();
    }

    template <size_t i> void set(get_type<i> value){
        if constexpr(i==0) first = value;
        else next.template set<i-1>(value);
    }

    template<size_t i> static consteval size_t get_offset() {
        if constexpr(i==0) return __builtin_offsetof(GLvertex<T, Ts...>, first); // TODO: using just offsetof() throws error for some reason, so we use __builtin_offsetof() instead - but it isn't portable
        else return __builtin_offsetof(GLvertex<T, Ts...>, next) + decltype(next)::template get_offset<i-1>();
    }
};
template<typename T>
struct GLvertex<T>{
    T first;

    GLvertex(T first) : first(first) {}
    GLvertex() {}

    template<size_t i>
    using get_type = std::conditional_t<i==0, T, std::monostate>; // TODO: replace std::monostate with some sort of error

    template <size_t i> auto& get(){
        static_assert(i==0, "Error: index out of bounds");
        return first;
    }

    template <size_t i> void set(T value){
        static_assert(i==0, "Error: index out of bounds");
        first = value;
    }

    template<size_t i> static consteval size_t get_offset() {
        static_assert(i==0, "Error: index out of bounds");
        if constexpr(i==0) return offsetof(GLvertex<T>, first);
    }
};

namespace GLattrib{
    template<typename T> consteval GLenum type2GLenum();
    template<typename T> consteval GLint type2count();
    template<Name> consteval const char* name2str();

    template<attrib...attribs, size_t...i>
    static consteval bool _check_attribs_impl(std::index_sequence<i...>){
        Name name;
        size_t idx;
        bool invalid = false;
        ((
            name = attribs.name,
            idx = i,
            invalid = invalid || ((name == attribs.name && idx != i) || ...)
        ),...);
        return !invalid;
    }
    template<attrib...attribs>
    static consteval bool check_attribs(){
        return _check_attribs_impl<attribs...>(std::make_index_sequence<sizeof...(attribs)>{});
    }
    template<attrib...As>
    static consteval bool check_attribs(attribs_monostate<As...>){
        return check_attribs<As...>();
    }
    template<attrib...A1s, attrib...A2s>
    static consteval bool check_attribs(attribs_monostate<A1s...>, attribs_monostate<A2s...>){
        return check_attribs<A1s...,A2s...>();
    }

    template<attrib...BaseAttribs, attrib...As>
    static consteval bool is_subset(attribs_monostate<BaseAttribs...>, attribs_monostate<As...>){
        GLattrib::Name needs;
        bool foundAll = ((needs = BaseAttribs.name, ((As.name == needs)||...))&&...);
        return foundAll;
    }


    template<> consteval GLenum type2GLenum<float>() { return GL_FLOAT; }
    template<> consteval GLenum type2GLenum<Vec2f>() { return GL_FLOAT; }
    template<> consteval GLenum type2GLenum<Vec3f>() { return GL_FLOAT; }
    template<> consteval GLenum type2GLenum<Quat4f>(){ return GL_FLOAT; }
    template<> consteval GLenum type2GLenum<int>()   { return GL_INT;   }
    template<> consteval GLenum type2GLenum<Vec2i>() { return GL_INT;   }
    template<> consteval GLenum type2GLenum<Vec3i>() { return GL_INT;   }
    template<> consteval GLenum type2GLenum<unsigned int>(){ return GL_UNSIGNED_INT; }

    template<> consteval GLint type2count<float>() { return 1; }
    template<> consteval GLint type2count<Vec2f>() { return 2; }
    template<> consteval GLint type2count<Vec3f>() { return 3; }
    template<> consteval GLint type2count<Quat4f>(){ return 4; }
    template<> consteval GLint type2count<int>()   { return 1; }
    template<> consteval GLint type2count<Vec2i>() { return 2; }
    template<> consteval GLint type2count<Vec3i>() { return 3; }
    template<> consteval GLint type2count<unsigned int>(){ return 1; }

    template<> consteval const char* name2str<Position>() { return "vPosition" ; }
    template<> consteval const char* name2str<Normal>()   { return "vNormal"   ; }
    template<> consteval const char* name2str<Color>()    { return "vColor"    ; }
    template<> consteval const char* name2str<UV>()       { return "vUV"       ; }
    template<> consteval const char* name2str<PosOffset>(){ return "vPosOffset"; }
    template<> consteval const char* name2str<Radius>()   { return "vRadius"   ; }

    static inline const char* name2str(Name name){
        switch(name){
            case Position:  return "vPosition";
            case Normal:    return "vNormal";
            case Color:     return "vColor";
            case UV:        return "vUV";
            case PosOffset: return "vPosOffset";
            case Radius:    return "vRadius";
            default:       throw std::runtime_error("Error: invalid name");
        }
    }
};

#define MPOS    attrib<Vec3f>(GLattrib::Position)
#define MNORMAL attrib<Vec3f>(GLattrib::Normal)
#define MCOLOR  attrib<Vec3f>(GLattrib::Color)
#define MUV     attrib<Vec2f>(GLattrib::UV)
#define MPOSOFFSET attrib<Vec3f>(GLattrib::PosOffset)
#define MRADIUS attrib<float>(GLattrib::Radius)


#endif // _GLATTRIBS_H_
