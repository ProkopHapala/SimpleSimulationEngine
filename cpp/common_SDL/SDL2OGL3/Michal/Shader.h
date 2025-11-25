#ifndef _SHADER_H_
#define _SHADER_H_

#include "GLattribs.h"
#include "Mat4.h"
#include <cstring>
#include <iostream>
#include <string>
#include <unordered_set>
#include <utility>
#include "GLuniform.h"

GLuint linkProgram(GLuint vertexShader, GLuint fragmentShader);
GLuint compileShader(GLenum shaderType, const char* source);

class Shader{
private:
    GLint attrName2LocMap[GLattrib::ATTRIB_NAME_MAX];

    GLuint programId        = 0;
    GLuint vertexShaderId   = 0;
    GLuint fragmentShaderId = 0;

    char* __vertexShaderSource   = nullptr;
    char* __fragmentShaderSource = nullptr;

    unsigned int update_uniforms = 0;
    std::vector<GLuniform> uniforms;
    std::unordered_set<GLuint> tex_uniforms;

public:
    Shader(const char* vertexShaderSource, const char* fragmentShaderSource){
        __vertexShaderSource = (char*)malloc(strlen(vertexShaderSource) + 1);
        __fragmentShaderSource = (char*)malloc(strlen(fragmentShaderSource) + 1);
        strcpy(__vertexShaderSource, vertexShaderSource);
        strcpy(__fragmentShaderSource, fragmentShaderSource);
    }
    ~Shader(){
        if (programId)glDeleteProgram(programId);
        if (vertexShaderId)glDeleteShader(vertexShaderId);
        if (fragmentShaderId)glDeleteShader(fragmentShaderId);
    }
    

    inline GLuint getProgramId() const { return programId; }
    inline void ensure_handle(){
        if (programId) return;
        create_handle();
    }

    inline void use(){
        ensure_handle();
        GLES::useProgram(programId);
    
        // update uniforms
        while (update_uniforms){
            int i = __builtin_ctz(update_uniforms); // count trailing zeros
            GLuniform u = uniforms[i];
    
            switch (u.type) {
                case GLuniform::f1:   glUniform1f(i, u.data.f1); break;
                case GLuniform::f2:   glUniform2f(i, u.data.f2.x, u.data.f2.y); break;
                case GLuniform::f3:   glUniform3f(i, u.data.f3.x, u.data.f3.y, u.data.f3.z); break;
                case GLuniform::f4:   glUniform4f(i, u.data.f4.x, u.data.f4.y, u.data.f4.z, u.data.f4.w); break;
                case GLuniform::i1:   glUniform1i(i, u.data.i1); break;
                case GLuniform::i2:   glUniform2i(i, u.data.i2.x, u.data.i2.y); break;
                case GLuniform::i3:   glUniform3i(i, u.data.i3.x, u.data.i3.y, u.data.i3.z); break;
                case GLuniform::i4:   glUniform4i(i, u.data.i4.x, u.data.i4.y, u.data.i4.z, u.data.i4.w); break;
                case GLuniform::ui1:  glUniform1ui(i, u.data.ui1); break;
                case GLuniform::ui2:  glUniform2ui(i, u.data.ui2.x, u.data.ui2.y); break;
                case GLuniform::ui3:  glUniform3ui(i, u.data.ui3.x, u.data.ui3.y, u.data.ui3.z); break;
                case GLuniform::ui4:  glUniform4ui(i, u.data.ui4.x, u.data.ui4.y, u.data.ui4.z, u.data.ui4.w); break;
                case GLuniform::m3:   glUniformMatrix3fv(i, 1, GL_FALSE, u.data.m3.array); break;
                case GLuniform::m4:   glUniformMatrix4fv(i, 1, GL_FALSE, u.data.m4.array); break;
                case GLuniform::tex:  break;
                default: printf("ERROR: invalid uniform type\n");
            }
            update_uniforms &= ~(1 << i);
        }

        // update textures
        int texi = 0;
        for (auto loc : tex_uniforms) {
            glActiveTexture(GL_TEXTURE0 + texi);
            glBindTexture(GL_TEXTURE_2D, uniforms[loc].data.tex->getHandle());
            glUniform1i(loc, texi);
            texi++;
        }
    }
    


    inline GLint getUniformLocation(const char* name) {
        ensure_handle();
        GLint loc = glGetUniformLocation(programId, name);
        if (loc == -1) printf("Uniform '%s' not found\n", name);
        return loc;
    }

    inline void setUniformName(const char* name, GLuniform value){
        ensure_handle();
        GLuint loc = getUniformLocation(name);
        
        setUniformLoc(loc, value);
    }

    inline void setUniformLoc(GLuint loc, GLuniform value){
        // note: we don't need ensure_handle() here, because if you have loc, you must've used getUniformLocation()
        if (loc == -1) return;
        if (loc >= uniforms.size()){ uniforms.resize(loc + 1); }
        
        // TODO: optimise matrix writes to only update changed rows (if only one row changes)
        if (uniforms[loc] == value) return;
        update_uniforms |= (1 << loc); // TODO: check that loc is not > 32
        if (uniforms[loc].type == GLuniform::tex) tex_uniforms.erase(loc);
        uniforms[loc] = value;
        if (value.type == GLuniform::tex) tex_uniforms.insert(loc);
    }

    inline void setUniforms(const GLuniformSet& set){
        for (auto u : set.uniforms){
            setUniformName(u.first.c_str(), u.second);
        }
    }

    inline void setUniform1f(const char* name, GLfloat value)       { setUniformName(name, {.type=GLuniform::f1, .data={.f1=value}}); }
    inline void setUniform2f(const char* name, Vec2T<GLfloat> value){ setUniformName(name, {.type=GLuniform::f2, .data={.f2=value}}); }
    inline void setUniform3f(const char* name, Vec3T<GLfloat> value){ setUniformName(name, {.type=GLuniform::f3, .data={.f3=value}}); }
    inline void setUniform4f(const char* name, Vec4T<GLfloat> value){ setUniformName(name, {.type=GLuniform::f4, .data={.f4=value}}); }

    inline void setUniform1i(const char* name, GLint value)         { setUniformName(name, {.type=GLuniform::i1, .data={.i1=value}}); }
    inline void setUniform2i(const char* name, Vec2T<GLint> value)  { setUniformName(name, {.type=GLuniform::i2, .data={.i2=value}}); }
    inline void setUniform3i(const char* name, Vec3T<GLint> value)  { setUniformName(name, {.type=GLuniform::i3, .data={.i3=value}}); }
    inline void setUniform4i(const char* name, Vec4T<GLint> value)  { setUniformName(name, {.type=GLuniform::i4, .data={.i4=value}}); }

    inline void setUniform1ui(const char* name, GLuint value)       { setUniformName(name, {.type=GLuniform::ui1, .data={.ui1=value}}); }
    inline void setUniform2ui(const char* name, Vec2T<GLuint> value){ setUniformName(name, {.type=GLuniform::ui2, .data={.ui2=value}}); }
    inline void setUniform3ui(const char* name, Vec3T<GLuint> value){ setUniformName(name, {.type=GLuniform::ui3, .data={.ui3=value}}); }
    inline void setUniform4ui(const char* name, Vec4T<GLuint> value){ setUniformName(name, {.type=GLuniform::ui4, .data={.ui4=value}}); }

    inline void setUniform3m(const char* name, Mat3T<GLfloat> value){ setUniformName(name, {.type=GLuniform::m3, .data={.m3=value}}); }
    inline void setUniform4m(const char* name, Mat4T<GLfloat> value){ setUniformName(name, {.type=GLuniform::m4, .data={.m4=value}}); }


    inline GLint attrName2Loc(GLattrib::Name attrName){
        ensure_handle();
        if (attrName > GLattrib::ATTRIB_NAME_MAX){
            printf("ERROR: invalid attribute name\n");
            return -1;
        }
        GLint loc = attrName2LocMap[attrName];
        if (loc == -1){
            std::cout << "Warning: attribute name " << attrName << " was not found\n";
        }
        return loc;
    }

private:
    inline void create_handle(){
        if (programId){
            printf("ERROR: shader handle already exists!\n");
            return;
        }
    
        vertexShaderId   = compileShader(GL_VERTEX_SHADER,   __vertexShaderSource  );
        fragmentShaderId = compileShader(GL_FRAGMENT_SHADER, __fragmentShaderSource);
        programId = linkProgram(vertexShaderId, fragmentShaderId);
    
        for (int i = 0; i < GLattrib::ATTRIB_NAME_MAX; i++){
            attrName2LocMap[i] = glGetAttribLocation(programId, GLattrib::name2str(GLattrib::Name(i)));
            if (attrName2LocMap[i] != -1) glEnableVertexAttribArray(attrName2LocMap[i]); // TODO - do we ever need to disable VertexAttribArrays?
        }

        free(__vertexShaderSource  ); __vertexShaderSource   = nullptr;
        free(__fragmentShaderSource); __fragmentShaderSource = nullptr;
    }
};






// default shaders

template<typename T>
static void appendToSource(std::string& source, const attrib<T> attr, GLattrib::Name ifname, std::string append){
    if (attr.name == ifname){
        source += append;
    }
}

template<attrib...attribs>
const std::string buildDefaultVertexShaderSource(){
    std::string source = std::string("");
    source += "#version 100\n";

    // uniforms
    source += "uniform mat4 uMVPMatrix;\n";

    // attributes
    source += "attribute mediump vec4 vPosition;\n";
    (appendToSource(source, attribs, GLattrib::PosOffset, "attribute mediump vec3 vPosOffset;\n"), ...);
    (appendToSource(source, attribs, GLattrib::Normal   , "attribute mediump vec3 vNormal;\n"   ), ...);
    (appendToSource(source, attribs, GLattrib::Color    , "attribute mediump vec3 vColor;\n"    ), ...);
    (appendToSource(source, attribs, GLattrib::UV       , "attribute mediump vec2 vUV;\n"       ), ...);

    // varyings
    (appendToSource(source, attribs, GLattrib::Normal, "varying mediump vec3 fNormal;\n"), ...);
    (appendToSource(source, attribs, GLattrib::Color , "varying mediump vec3 fColor;\n" ), ...);
    (appendToSource(source, attribs, GLattrib::UV    , "varying mediump vec2 fUV;\n"    ), ...);

    // void main()
    source += "void main(){\n";
    source += "mediump vec4 world_pos = vPosition;\n";
    (appendToSource(source, attribs, GLattrib::PosOffset, "world_pos += vec4(vPosOffset, 0.0);\n"), ...);
    source += "gl_Position = uMVPMatrix * world_pos;\n";
    (appendToSource(source, attribs, GLattrib::Normal, "fNormal" " = vNormal;\n"), ...);
    (appendToSource(source, attribs, GLattrib::Color , "fColor"  " = vColor;\n" ), ...);
    (appendToSource(source, attribs, GLattrib::UV    , "fUV"     " = vUV;\n"    ), ...);
    source += "}";

    return source;
}

template<attrib...attribs>
const std::string buildDefaultFragmentShaderSource(bool tex, bool ucolor){
    std::string source = std::string("");
    source += "#version 100\n";

    // uniforms
    if (ucolor) source += "uniform mediump vec3 uColor;\n";
    if (tex) source += "uniform sampler2D uTexture;\n";

    // varyings
    (appendToSource(source, attribs, GLattrib::Normal, "varying mediump vec3 fNormal;\n"), ...);
    (appendToSource(source, attribs, GLattrib::Color , "varying mediump vec3 fColor;\n" ), ...);
    (appendToSource(source, attribs, GLattrib::UV    , "varying mediump vec2 fUV;\n"    ), ...);

    // void main()
    source += "void main(){\n";
    if (!ucolor) source += "gl_FragColor = vec4(1.0, 1.0, 1.0, 1.0);\n";
    if (ucolor) source += "gl_FragColor = vec4(uColor, 1.0);\n";
    if (tex) source += "gl_FragColor = gl_FragColor*texture2D(uTexture, fUV);\n";
    (appendToSource(source, attribs, GLattrib::Color, "gl_FragColor = gl_FragColor*vec4(fColor, 1.0);\n"), ...);

    (appendToSource(source, attribs, GLattrib::Normal, 
        "mediump float light = dot(fNormal, vec3(1.0, -1.0, 1.0));\n"
        "light = (light+1.0)/2.0;\n" // normalised to range (0; 1)
        "light = 0.3 + light*0.6;\n" // to range (0.3, 1.1)
        "gl_FragColor = gl_FragColor*vec4(light, light, light, 1.0);\n" // TODO: remove or implement lighting
    ), ...);


    source += "if (gl_FragColor.a == 0.0) discard;\n";
    source += "}";

    return source;
}

template<attrib...attribs>
Shader* defaultShader = new Shader(buildDefaultVertexShaderSource<attribs...>().c_str(), buildDefaultFragmentShaderSource<attribs...>(false, false).c_str());

template<attrib...attribs>
Shader* defaultcolorShader = new Shader(buildDefaultVertexShaderSource<attribs...>().c_str(), buildDefaultFragmentShaderSource<attribs...>(false, true).c_str());

template<attrib...attribs>
Shader* defaultTexShader = new Shader(buildDefaultVertexShaderSource<attribs...>().c_str(), buildDefaultFragmentShaderSource<attribs...>(true, false).c_str());

template<attrib...attribs>
Shader* defaultTexColorShader = new Shader(buildDefaultVertexShaderSource<attribs...>().c_str(), buildDefaultFragmentShaderSource<attribs...>(true, true).c_str());


#endif // _SHADER_H_
