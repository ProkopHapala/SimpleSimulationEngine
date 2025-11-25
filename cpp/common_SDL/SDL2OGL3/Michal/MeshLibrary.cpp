#include "MeshLibrary.h"
#include "MeshBuilder.h"

static GLMesh<MPOS> makePointCross() {
    GLMesh<MPOS> mesh(GL_LINES);
    mesh.addVertex({-1, 0, 0});
    mesh.addVertex({ 1, 0, 0});
    mesh.addVertex({ 0,-1, 0});
    mesh.addVertex({ 0, 1, 0});
    mesh.addVertex({ 0, 0,-1});
    mesh.addVertex({ 0, 0, 1});
    return mesh;
}

static GLMesh<MPOS> makeLine() {
    GLMesh<MPOS> mesh(GL_LINES);
    mesh.addVertex({0, 0, 0});
    mesh.addVertex({1, 1, 1});
    return mesh;
}

static GLMesh<MPOS> makeLine2D() {
    GLMesh<MPOS> mesh(GL_LINES);
    mesh.addVertex({0, 0, 0});
    mesh.addVertex({1, 1, 0});
    return mesh;
}

static GLMesh<MPOS> makePoint() {
    GLMesh<MPOS> mesh((GLenum)GL_POINTS);
    mesh.addVertex({0, 0, 0});
    return mesh;
}

static GLMesh<MPOS> makeWireCube() {
    GLMesh<MPOS> m = GLMesh<MPOS>(GL_LINES);
    // Bottom face
    m.addVertex({0, 0, 0}); m.addVertex({1, 0, 0});
    m.addVertex({1, 0, 0}); m.addVertex({1, 0, 1});
    m.addVertex({1, 0, 1}); m.addVertex({0, 0, 1});
    m.addVertex({0, 0, 1}); m.addVertex({0, 0, 0});
    
    // Top face
    m.addVertex({0, 1, 0}); m.addVertex({1, 1, 0});
    m.addVertex({1, 1, 0}); m.addVertex({1, 1, 1});
    m.addVertex({1, 1, 1}); m.addVertex({0, 1, 1});
    m.addVertex({0, 1, 1}); m.addVertex({0, 1, 0});
    
    // Connecting vertical lines
    m.addVertex({0, 0, 0}); m.addVertex({0, 1, 0});
    m.addVertex({1, 0, 0}); m.addVertex({1, 1, 0});
    m.addVertex({1, 0, 1}); m.addVertex({1, 1, 1});
    m.addVertex({0, 0, 1}); m.addVertex({0, 1, 1});
    
    return m;
}

static GLMesh<MPOS,MNORMAL> makeCubeWithNormals() {
    GLMesh<MPOS,MNORMAL> m = GLMesh<MPOS,MNORMAL>(GL_TRIANGLES);

    // positive x face
    m.addVertex({1, 0, 0}, {1, 0, 0});
    m.addVertex({1, 1, 0}, {1, 0, 0});
    m.addVertex({1, 1, 1}, {1, 0, 0});
    m.addVertex({1, 0, 0}, {1, 0, 0});
    m.addVertex({1, 1, 1}, {1, 0, 0});
    m.addVertex({1, 0, 1}, {1, 0, 0});

    // negative x face
    m.addVertex({0, 0, 0}, {-1, 0, 0});
    m.addVertex({0, 1, 0}, {-1, 0, 0});
    m.addVertex({0, 1, 1}, {-1, 0, 0});
    m.addVertex({0, 0, 0}, {-1, 0, 0});
    m.addVertex({0, 1, 1}, {-1, 0, 0});
    m.addVertex({0, 0, 1}, {-1, 0, 0});

    // positive y face
    m.addVertex({0, 1, 0}, {0, 1, 0});
    m.addVertex({1, 1, 0}, {0, 1, 0});
    m.addVertex({1, 1, 1}, {0, 1, 0});
    m.addVertex({0, 1, 0}, {0, 1, 0});
    m.addVertex({1, 1, 1}, {0, 1, 0});
    m.addVertex({0, 1, 1}, {0, 1, 0});

    // negative y face
    m.addVertex({0, 0, 0}, {0, -1, 0});
    m.addVertex({1, 0, 0}, {0, -1, 0});
    m.addVertex({1, 0, 1}, {0, -1, 0});
    m.addVertex({0, 0, 0}, {0, -1, 0});
    m.addVertex({1, 0, 1}, {0, -1, 0});
    m.addVertex({0, 0, 1}, {0, -1, 0});

    // positive z face
    m.addVertex({0, 0, 1}, {0, 0, 1});
    m.addVertex({1, 0, 1}, {0, 0, 1});
    m.addVertex({1, 1, 1}, {0, 0, 1});
    m.addVertex({0, 0, 1}, {0, 0, 1});
    m.addVertex({1, 1, 1}, {0, 0, 1});
    m.addVertex({0, 1, 1}, {0, 0, 1});

    // negative z face
    m.addVertex({0, 0, 0}, {0, 0, -1});
    m.addVertex({1, 0, 0}, {0, 0, -1});
    m.addVertex({1, 1, 0}, {0, 0, -1});
    m.addVertex({0, 0, 0}, {0, 0, -1});
    m.addVertex({1, 1, 0}, {0, 0, -1});
    m.addVertex({0, 1, 0}, {0, 0, -1});

    return m;
}

static GLMesh<MPOS> makeRectMesh(){
    GLMesh<MPOS> m = GLMesh<MPOS>(GL_TRIANGLE_FAN);
    m.addVertex( {0, 0, 0} );
    m.addVertex( {0, 1, 0} );
    m.addVertex( {1, 1, 0} );
    m.addVertex( {1, 0, 0} );
    return m;
}

static const GLMesh<MPOS> makeCircleMesh(){
	GLMesh<MPOS> m = GLMesh<MPOS>(GL_LINE_LOOP);
	const int n = 64;
	float dphi =  6.28318530718f / n;
	Vec2f drot; drot.fromAngle( dphi );
	Vec2f v = {1, 0};
	for ( int i=0; i<n; i++ ){
		m.addVertex( {v.x, v.y, 0} );
		v.mul_cmplx( drot );
	}
	return m;
}

template<bool instanced>
static std::string makeSphereVertexShader(){
    std::string sh = "#version 300 es\n";

    sh += "in mediump vec3 vPosition;\n";
    if constexpr (instanced){
        sh += "in mediump vec3  vPosOffset;\n";
        sh += "in mediump vec3  vColor;\n";
        sh += "in mediump float vRadius;\n";
    }

    sh += "out highp vec3 fPos;\n";
    sh += "out highp vec3 fDir;\n";
    if constexpr (instanced)
        sh += "out mediump vec3 fColor;\n";

    sh += "uniform mat4 uMVPinv;\n";
    sh += "uniform mat4 uMVPMatrix;\n";
    if constexpr (!instanced){
        sh += "uniform vec3  uPos;\n";
        sh += "uniform float uRadius;\n";
    }

    if constexpr (instanced){
        sh += "#define SpherePos vPosOffset\n";
        sh += "#define SphereRadius vRadius\n";
    }else{
        sh += "#define SpherePos uPos\n";
        sh += "#define SphereRadius uRadius\n";
    }

    sh += "vec4 model2clip(vec4 pos){";
    sh +=   "return uMVPMatrix * (pos*SphereRadius + vec4(SpherePos, 0.0)*pos.w);";
    sh += "}";

    sh += "vec4 clip2model(vec4 pos){";
    sh +=   "return ((uMVPinv * pos) - vec4(SpherePos, 0.0)*pos.w) / SphereRadius;";
    sh += "}";

    sh += "void main(){\n";
    sh +=   "mediump vec3 sphPos = model2clip(vec4(0.0, 0.0, 0.0, 1.0)).xyz;\n";
    sh +=   "mediump vec3 vecRight = normalize(clip2model(vec4(1.0, 0.0, 0.0, 0.0)).xyz);\n";
    sh +=   "mediump vec3 vecUp    = normalize(clip2model(vec4(0.0, 1.0, 0.0, 0.0)).xyz);\n";

    sh +=   "vecRight = model2clip(vec4(vecRight, 0.0)).xyz;\n";
    sh +=   "vecUp    = model2clip(vec4(vecUp   , 0.0)).xyz;\n";

    sh +=   "mediump vec3 newPos = vecRight*vPosition.x + vecUp*vPosition.y;\n";
    sh +=   "gl_Position = vec4(newPos + sphPos, 1.0);\n";

    sh +=   "fPos = clip2model(vec4(gl_Position.xy, 0.0, 1.0)).xyz;\n";
    sh +=   "fDir = clip2model(vec4(0.0, 0.0, 2.0, 0.0)).xyz;\n";
    if constexpr (instanced) sh += "fColor = vColor;\n";
    sh += "}";

    return sh;
}

template <bool instanced>
static std::string makeSphereFragmentShader(){
    std::string sh = "#version 300 es\n";

    sh += "in highp vec3 fPos;\n";
    sh += "in highp vec3 fDir;\n";
    if constexpr (instanced){
        sh += "in mediump vec3 fColor;\n";
        sh += "#define SphereColor fColor\n";
    }
    else{
        sh += "uniform mediump vec3 uColor;\n";
        sh += "#define SphereColor uColor\n";
    }

    sh += "layout(location=0) out mediump vec4 FragColor;\n";

    sh += "void main(){\n";
    sh +=   "highp float DirLen = length(fDir);\n";
    sh +=   "highp vec3 Dir = fDir / DirLen;\n";
    sh +=   "highp float Tc = -dot(fPos, Dir);\n";
    sh +=   "highp float dsqr = dot(fPos, fPos) - Tc*Tc;\n";
    sh +=   "if (dsqr >= 1.0) discard;\n";
    sh +=   "highp float T1c = sqrt(1.0 - dsqr);\n";
    sh +=   "highp float T1 = Tc - T1c;\n";

    sh +=   "highp vec3 normal = fPos + T1*Dir;\n";
    sh +=   "gl_FragDepth = clamp((T1 / DirLen) + .5, 0.0, 1.0);\n";

    sh +=   "mediump float light = dot(normal, vec3(1.0, -1.0, 1.0));\n";
    sh +=   "light = (light+1.0)/2.0;\n";
    sh +=   "light = .3 + light*.6;\n";
    sh +=   "FragColor = vec4(SphereColor, 1.0)*vec4(light, light, light, 1.0);\n";
    sh += "}";

    return sh;
}

static GLInstancedMeshBase<GLvbo<MPOS>, MPOSOFFSET, MRADIUS, MCOLOR> makeSphereInstanced(){
    GLInstancedMeshBase<GLvbo<MPOS>, MPOSOFFSET, MRADIUS, MCOLOR> m(GL_TRIANGLES, new Shader(makeSphereVertexShader<true>().c_str(), makeSphereFragmentShader<true>().c_str()));
    m.verts->push_back((Vec3f){-1, -1, -1});
    m.verts->push_back((Vec3f){ 1, -1, -1});
    m.verts->push_back((Vec3f){-1,  1, -1});

    m.verts->push_back((Vec3f){-1,  1, -1});
    m.verts->push_back((Vec3f){ 1, -1, -1});
    m.verts->push_back((Vec3f){ 1,  1, -1});
    return m;
}

static GLMeshBase<MPOS> makeSphere(){
    GLMeshBase<MPOS> m(GL_TRIANGLES, GL_STATIC_DRAW, new Shader(makeSphereVertexShader<false>().c_str(), makeSphereFragmentShader<false>().c_str()));
    m.verts->push_back((Vec3f){-1, -1, -1});
    m.verts->push_back((Vec3f){ 1, -1, -1});
    m.verts->push_back((Vec3f){-1,  1, -1});

    m.verts->push_back((Vec3f){-1,  1, -1});
    m.verts->push_back((Vec3f){ 1, -1, -1});
    m.verts->push_back((Vec3f){ 1,  1, -1});
    return m;
}

static GLMesh<MPOS> makeCross() {
    GLMesh<MPOS> m = GLMesh<MPOS>(GL_LINES);
    float sz = 0.5f;
    m.addVertex({-sz, 0, 0}); m.addVertex({sz, 0, 0});
    m.addVertex({0, -sz, 0}); m.addVertex({0, sz, 0});
    return m;
}

static GLMesh<MPOS> makeXMark() {
    GLMesh<MPOS> m = GLMesh<MPOS>(GL_LINES);
    float sz = 0.5f;
    m.addVertex({-sz, -sz, 0}); m.addVertex({sz, sz, 0});
    m.addVertex({-sz, sz, 0}); m.addVertex({sz, -sz, 0});
    return m;
}

static GLvbo<MPOS> makeOctSphere(){
    GLvbo<MPOS> m;
    MeshBuilder::addSphereOct(m, 8, 1, Vec3fZero);
    return m;
}

namespace MeshLibrary {
    GLMesh<MPOS> pointCross = makePointCross();
    GLMesh<MPOS> line = makeLine();
    GLMesh<MPOS> line2D = makeLine2D();
    GLMesh<MPOS> point = makePoint();
    GLMesh<MPOS> wireCube = makeWireCube();
    GLMesh<MPOS,MNORMAL> cubeWithNormals = makeCubeWithNormals();
    GLMesh<MPOS> rect = makeRectMesh();
    GLMesh<MPOS> circle = makeCircleMesh();
    GLMeshBase<MPOS> sphere = makeSphere();
    GLInstancedMeshBase<GLvbo<MPOS>, MPOSOFFSET, MRADIUS, MCOLOR> sphereInstanced = makeSphereInstanced();
    GLMesh<MPOS> cross = makeCross();
    GLMesh<MPOS> xmark = makeXMark();

    GLvbo<MPOS> octSphere = makeOctSphere();
    GLMesh<MPOS> octSphereMesh(&octSphere, GL_LINES);
    GLInstancedMeshBase<GLvbo<MPOS>, MPOSOFFSET> octSphereInstanced(&octSphere, GL_LINES);
}
