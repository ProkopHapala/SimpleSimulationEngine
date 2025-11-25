#include "GLES.h"
#include "Renderer.h"
#include "Vec2.h"
#include "Draw_o3.h"
#include "Vec3.h"
#include "quaternion.h"
#include "GLMesh.h"
#include "MeshLibrary.h"

#include "Draw3D_o3.h" // THE HEADER

// Dynamic meshes for operations that need to modify vertices
static GLMesh<MPOS>         tmpMesh1 = GLMesh<MPOS>        (GL_TRIANGLES, GL_STREAM_DRAW);
static GLMesh<MPOS,MNORMAL> tmpMesh2 = GLMesh<MPOS,MNORMAL>(GL_TRIANGLES, GL_STREAM_DRAW);
static GLMesh<MPOS,MCOLOR>  tmpMesh3 = GLMesh<MPOS,MCOLOR> (GL_TRIANGLES, GL_STREAM_DRAW);

void Draw3D::drawPoint( const Vec3f& vec ){
	MeshLibrary::point.setUniform3f("uColor", {1, 1, 1});
	MeshLibrary::point.draw(vec);
};

void Draw3D::drawPointCross( const Vec3f& vec, float sz, Vec3f color ){
    MeshLibrary::pointCross.setUniform3f("uColor", color);
	MeshLibrary::pointCross.draw(vec, (Vec3f){sz, sz, sz});
};

void Draw3D::drawLine( const Vec3f& p1, const Vec3f& p2, Vec3f color ){
    MeshLibrary::line.setUniform3f("uColor", color);
    MeshLibrary::line.draw(p1, p2-p1);
};

void Draw3D::drawVecInPos( const Vec3f& v, const Vec3f& pos, Vec3f color ){
    drawLine(pos, pos+v, color);
};

void Draw3D::drawVec( const Vec3f& vec, Vec3f color ){
	drawVecInPos( vec, Vec3fZero, color);
};

void Draw3D::drawArrow( const Vec3f& p1, const Vec3f& p2, float sz ){
    tmpMesh1.clear();
    Vec3f dir = p2-p1;
    float len = dir.norm();
    if(len < 1e-10f) return;
    dir.mul(1.0f/len);
    
    // Main line
    tmpMesh1.addVertex(p1);
    tmpMesh1.addVertex(p2);
    
    // Arrow head
    Vec3f up, right;
    dir.getSomeOrtho(up, right);
    up.mul(sz);
    right.mul(sz);
    
    Vec3f back = dir * (-sz);
    Vec3f tip = p2;
    Vec3f p;
    
    p = tip + back + up;    tmpMesh1.addVertex(tip); tmpMesh1.addVertex(p);
    p = tip + back - up;    tmpMesh1.addVertex(tip); tmpMesh1.addVertex(p);
    p = tip + back + right; tmpMesh1.addVertex(tip); tmpMesh1.addVertex(p);
    p = tip + back - right; tmpMesh1.addVertex(tip); tmpMesh1.addVertex(p);
    
    tmpMesh1.setUniform3f("uColor", opengl1renderer.color);
    tmpMesh1.draw(GL_LINES);
};

void Draw3D::vecsInPoss( int n, const Vec3d* vs, const Vec3d* ps, float sc, Vec3f color ){
    tmpMesh1.clear();
    for(int i=0; i<n; i++){
        Vec3f p1=(Vec3f)ps[i];
        Vec3f p2=(Vec3f)ps[i];
        p2.add_mul( (Vec3f)vs[i], sc);
        tmpMesh1.addVertex(p1);
        tmpMesh1.addVertex(p2);
    }
    tmpMesh1.setUniform3f("uColor", color);
    tmpMesh1.draw(GL_LINES);
};

void Draw3D::drawPolyLine( int n, Vec3d * ps, bool closed ){
    tmpMesh1.clear();
    for(int i=0; i<n; i++){
        tmpMesh1.addVertex((Vec3f)ps[i]);
    }
    tmpMesh1.setUniform3f("uColor", opengl1renderer.color);
    tmpMesh1.draw(closed ? GL_LINE_LOOP : GL_LINE_STRIP);
};

void Draw3D::drawTriangle( const Vec3f& p1, const Vec3f& p2, const Vec3f& p3 ){
    Vec3f d1 = p2 - p1;
    Vec3f d2 = p3 - p1;
    Vec3f nr;
    nr.set_cross(d1,d2);
    nr.normalize();
    tmpMesh2.clear();
    tmpMesh2.addVertex(p1, nr);
    tmpMesh2.addVertex(p2, nr);
    tmpMesh2.addVertex(p3, nr);
    tmpMesh2.setUniform3f("uColor", opengl1renderer.color);
    tmpMesh2.draw(GL_TRIANGLES);
}

void Draw3D::drawTriangle ( const Vec3f& p1,  const Vec3f& p2, const Vec3f& p3, bool filled ){
    if (filled) {
        drawTriangle(p1, p2, p3);
    } else {
        tmpMesh1.clear();
        tmpMesh1.addVertex(p1); tmpMesh1.addVertex(p2);
        tmpMesh1.addVertex(p2); tmpMesh1.addVertex(p3);
        tmpMesh1.addVertex(p3); tmpMesh1.addVertex(p1);
        tmpMesh1.setUniform3f("uColor", opengl1renderer.color);
        tmpMesh1.draw(GL_LINES);
    }
}

static void drawQuad( const Vec3f& p1, const Vec3f& p2, const Vec3f& p3, const Vec3f& p4 ){
    double r13=(p3-p1).norm2();
    double r24=(p4-p2).norm2();
    if(r13>r24){ 
        Draw3D::drawTriangle( p1, p2, p4 ); 
        Draw3D::drawTriangle( p2, p3, p4 ); 
    } else { 
        Draw3D::drawTriangle( p1, p2, p3 ); 
        Draw3D::drawTriangle( p3, p4, p1 ); 
    }
}

Vec3f lincomb(const Vec3f& p1, const Vec3f& p2, double v1, double v2 ){
    double f  = v1/(v1-v2);
    Vec3f p;
    p.set_lincomb( 1-f, p1, f, p2 );
    return p;
}

void Draw3D::drawTetraIso( Vec3f** ps, Quat4d vals ){
    bool b0 = vals.x>0;
    bool b1 = vals.y>0;
    bool b2 = vals.z>0;
    bool b3 = vals.w>0;
    int n01 = b0+b1;
    int n23 = b2+b3;
    int n   = n01+n23;
    if( n==0 || n==4 ) return;
    int i1,i2;
    int j1,j2,j3;
    if(n==1){        // triangle
        if(n01){
            if(b0){ i1=0; j1=1; j2=2; j3=3; } // b1
            else  { i1=1; j1=0; j2=3; j3=2; } // b2
        }else{
            if(b2){ i1=2; j1=3; j2=0; j3=1; } // b3
            else  { i1=3; j1=2; j2=1; j3=0; } // b4
        }
        drawTriangle(
            lincomb(*ps[i1],*ps[j1],vals.array[i1],vals.array[j1]),
            lincomb(*ps[i1],*ps[j2],vals.array[i1],vals.array[j2]),
            lincomb(*ps[i1],*ps[j3],vals.array[i1],vals.array[j3])
        );
    }else if (n==3){ // triangle
        if(n01==1){
            if(!b0){ i1=0; j1=1; j2=3; j3=2; } // b1
            else   { i1=1; j1=0; j2=2; j3=3; } // b2
        }else{
            if(!b2){ i1=2; j1=1; j2=0; j3=3; } // b3
            else   { i1=3; j1=1; j2=2; j3=0; } // b4
        }
        drawTriangle(
            lincomb(*ps[i1],*ps[j1],-vals.array[i1],-vals.array[j1]),
            lincomb(*ps[i1],*ps[j2],-vals.array[i1],-vals.array[j2]),
            lincomb(*ps[i1],*ps[j3],-vals.array[i1],-vals.array[j3])
        );
    }else if (n==2){ // quad
        if(n01==1){
            if(b0){
                if(b2){ i1=0; i2=2; j1=3; j2=1; } // b1-b3 | b2-b4
                else  { i1=0; i2=3; j1=1; j2=2; } // b1-b4 | b2-b3
            }else{
                if(b3){ j1=0; j2=2; i1=3; i2=1; } // b2-b3 | b1-b4
                else  { j1=0; j2=3; i1=1; i2=2; } // b2-b4 | b1-b3
            }
        }else{
            if(n01==2){ i1=0; i2=1; j1=2; j2=3; }
            else      { j1=0; j2=1; i1=2; i2=3; }
        }
        drawQuad(
            lincomb(*ps[i1],*ps[j1],vals.array[i1],vals.array[j1]),
            lincomb(*ps[i1],*ps[j2],vals.array[i1],vals.array[j2]),
            lincomb(*ps[i2],*ps[j2],vals.array[i2],vals.array[j2]),
            lincomb(*ps[i2],*ps[j1],vals.array[i2],vals.array[j1])
        );
    }
};

void Draw3D::drawMatInPos( const Mat3f& mat, const Vec3f& pos, const Vec3f& sc ){
    tmpMesh3.clear();
    // X axis (red)
    tmpMesh3.addVertex(pos, {1,0,0});
    tmpMesh3.addVertex({pos.x+mat.xx*sc.x, pos.y+mat.xy*sc.x, pos.z+mat.xz*sc.x}, {1,0,0});
    // Y axis (green)
    tmpMesh3.addVertex(pos, {0,1,0});
    tmpMesh3.addVertex({pos.x+mat.yx*sc.y, pos.y+mat.yy*sc.y, pos.z+mat.yz*sc.y}, {0,1,0});
    // Z axis (blue)
    tmpMesh3.addVertex(pos, {0,0,1});
    tmpMesh3.addVertex({pos.x+mat.zx*sc.z, pos.y+mat.zy*sc.z, pos.z+mat.zz*sc.z}, {0,0,1});

    tmpMesh3.draw(GL_LINES);
}

void Draw3D::drawShape( int shape, const Vec3f& pos, const Mat3f& rot, bool trasposed ){
    // TODO: This function needs to be reimplemented without using opengl1renderer.callList
    // For now, keeping the old implementation
    opengl1renderer.pushMatrix();
    float glMat[16];
    if( trasposed ){
        toGLMatT ( pos, rot, glMat );
    }else{
        toGLMat( pos, rot, glMat );
	}
	opengl1renderer.multMatrixf( glMat );
	opengl1renderer.callList( shape );
	opengl1renderer.popMatrix();
};

void Draw3D::drawShape    ( int shape, const Vec3f& pos, const Quat4f& qrot, const Vec3f& scale ){
	opengl1renderer.pushMatrix();
	float glMat[16];
	toGLMat ( pos, qrot, scale, glMat );
	opengl1renderer.multMatrixf( glMat );
	opengl1renderer.callList( shape );
	opengl1renderer.popMatrix();
};

int Draw3D::drawConeFan( int n, float r, const Vec3f& base, const Vec3f& tip ){
    int nvert=0;
    Vec3f a,b,c,c_hat;
    c.set_sub( tip, base );
    c_hat.set_mul( c, 1/c.norm() );
    c_hat.getSomeOrtho( a, b );
    a.normalize();
    b.normalize();
    float alfa = 2*M_PI/n;
    Vec2f rot,drot;
    rot .set(1.0f,0.0f);
    drot.set( cos( alfa ), sin( alfa ) );

    Vec3f q; q.set(c); q.add_mul( a, -r );
    float pnab =  c_hat.dot( q )/q.norm();
    float pnc  =  sqrt( 1 - pnab*pnab );

    tmpMesh2.clear();
    tmpMesh2.addVertex(tip, c_hat);

    nvert++;

    for(int i=0; i<=n; i++ ){
        Vec3f p,pn;
        p .set( rot.x*a.x +  rot.y*b.x, rot.x*a.y + rot.y*b.y, rot.x*a.z + rot.y*b.z    );
        pn.set( pnab*p.x + pnc*c_hat.x, pnab*p.y + pnc*c_hat.y, pnab*p.z + pnc*c_hat.z  );

        Vec3f vertex = {base.x + r*p.x, base.y + r*p.y, base.z + r*p.z};
        tmpMesh2.addVertex(vertex, pn);
        nvert++;
        rot.mul_cmplx( drot );
    }

    tmpMesh2.draw(GL_TRIANGLE_FAN);
    return nvert;
}

void Draw3D::drawSphere(Vec3f pos, float r, Vec3f color){
    MeshLibrary::sphere.uniforms.set3f("uColor", color);
    MeshLibrary::sphere.uniforms.set4m("uMVPMatrix", GLES::active_camera->viewProjectionMatrix());
    MeshLibrary::sphere.uniforms.set4m("uMVPinv", GLES::active_camera->inverseViewProjectionMatrix());
    MeshLibrary::sphere.uniforms.set3f("uPos", pos);
    MeshLibrary::sphere.uniforms.set1f("uRadius", r);

    MeshLibrary::sphere.draw();
}

int Draw3D::drawCircleAxis( int n, const Vec3f& pos, const Vec3f& v0, const Vec3f& uaxis, float R, float dca, float dsa ){
    int nvert=0;
    Vec3f v; v.set(v0);
    tmpMesh1.clear();
    for( int i=0; i<n; i++ ){
        tmpMesh1.addVertex(pos+v*R);
        nvert++;
        v.rotate_csa( dca, dsa, uaxis );
    }
    tmpMesh1.draw(GL_LINE_LOOP);
    return nvert;
}

int Draw3D::drawCircleAxis( int n, const Vec3f& pos, const Vec3f& v0, const Vec3f& uaxis, float R ){
    float dphi = 2*M_PI/n;
    float dca  = cos( dphi );
    float dsa  = sin( dphi );
    return drawCircleAxis( n, pos, v0, uaxis, R, dca, dsa );
}

int Draw3D::drawSphereOctLines( int n, float R, const Vec3f& pos, const Mat3f& rot, bool bRGB ){
    int nvert=0;
    float dphi = 2*M_PI/n;
    float dca  = cos( dphi );
    float dsa  = sin( dphi );

    if(bRGB) {
        // Red circle
        tmpMesh3.clear();
        Vec3f v = rot.b;
        for(int i=0; i<n; i++){
            tmpMesh3.addVertex({pos.x+v.x*R, pos.y+v.y*R, pos.z+v.z*R}, {1,0,0});
            v.rotate_csa(dca, dsa, rot.a);
        }
        tmpMesh3.addVertex({pos.x+rot.b.x*R, pos.y+rot.b.y*R, pos.z+rot.b.z*R}, {1,0,0}); // Close the loop
        tmpMesh3.draw(GL_LINE_STRIP);

        // Green circle
        tmpMesh3.clear();
        v = rot.c;
        for(int i=0; i<n; i++){
            tmpMesh3.addVertex({pos.x+v.x*R, pos.y+v.y*R, pos.z+v.z*R}, {0,1,0});
            v.rotate_csa(dca, dsa, rot.b);
        }
        tmpMesh3.addVertex({pos.x+rot.c.x*R, pos.y+rot.c.y*R, pos.z+rot.c.z*R}, {0,1,0}); // Close the loop
        tmpMesh3.draw(GL_LINE_STRIP);

        // Blue circle
        tmpMesh3.clear();
        v = rot.a;
        for(int i=0; i<n; i++){
            tmpMesh3.addVertex({pos.x+v.x*R, pos.y+v.y*R, pos.z+v.z*R}, {0,0,1});
            v.rotate_csa(dca, dsa, rot.c);
        }
        tmpMesh3.addVertex({pos.x+rot.a.x*R, pos.y+rot.a.y*R, pos.z+rot.a.z*R}, {0,0,1}); // Close the loop
        tmpMesh3.draw(GL_LINE_STRIP);
    } else {
        MeshLibrary::octSphereMesh.uniforms.set3f("uColor", COLOR_GREEN);
        MeshLibrary::octSphereMesh.draw(pos, R);
    }
    return nvert;
}

void Draw3D::drawSphereOctLinesInstanced( float r, const std::vector<Vec3f>& ps, Vec3f color ){
    MeshLibrary::octSphereInstanced.uniforms.set3f("uColor", color);
    MeshLibrary::octSphereInstanced.instances->clear();
    for( auto& p : ps ){
        MeshLibrary::octSphereInstanced.addInstance(p);
    }
    MeshLibrary::octSphereInstanced.draw(); // TODO: use r
}
void Draw3D::drawSphereOctLinesInstanced( float r, const Vec3d* ps, int n, Vec3f color ){
    MeshLibrary::octSphereInstanced.uniforms.set3f("uColor", color);
    MeshLibrary::octSphereInstanced.instances->clear();
    for(int i=0; i<n; i++){
        MeshLibrary::octSphereInstanced.addInstance((Vec3f)ps[i]);
    }
    MeshLibrary::octSphereInstanced.draw(); // TODO: use r
}

void Draw3D::drawPlanarPolygon( int n, const int * inds, const Vec3d * points ){
    if( n < 3 ) return;

    Vec3f a = (Vec3f)points[inds[0]];
    Vec3f b = (Vec3f)points[inds[1]];
    Vec3f c = (Vec3f)points[inds[2]];
    Vec3f normal;
    normal.set_cross( a-b, b-c );
    normal.normalize();

    tmpMesh2.clear();
    tmpMesh2.addVertex(a, normal);
    tmpMesh2.addVertex(b, normal);
    tmpMesh2.addVertex(c, normal);
    
    for( int i=3; i<n; i++ ){
        Vec3f vertex = (Vec3f)points[inds[i]];
        tmpMesh2.addVertex(vertex, normal);
    }
    tmpMesh2.draw(GL_TRIANGLE_FAN);
}

void Draw3D::drawPlanarPolygon( int ipl, Mesh& mesh ){
    Polygon * pl = mesh.polygons[ipl];
    Draw3D::drawPlanarPolygon( pl->ipoints.size(), &pl->ipoints.front(), &mesh.points.front() );
}

void Draw3D::drawPoints( int n, const Vec3d * points, float sz ){
    if(sz <= 0) {
        tmpMesh1.clear();
        for(int i=0; i<n; i++){
            tmpMesh1.addVertex((Vec3f)points[i]);
        }
        tmpMesh1.setUniform3f("uColor", opengl1renderer.color);
        tmpMesh1.draw(GL_POINTS);
    } else {
        tmpMesh1.clear();
        for(int i=0; i<n; i++){
            Vec3f vec = (Vec3f)points[i];
            tmpMesh1.addVertex({vec.x-sz, vec.y   , vec.z   }); tmpMesh1.addVertex({vec.x+sz, vec.y   , vec.z   });
            tmpMesh1.addVertex({vec.x   , vec.y-sz, vec.z   }); tmpMesh1.addVertex({vec.x   , vec.y+sz, vec.z   });
            tmpMesh1.addVertex({vec.x   , vec.y   , vec.z-sz}); tmpMesh1.addVertex({vec.x   , vec.y   , vec.z+sz});
        }
        tmpMesh1.setUniform3f("uColor", opengl1renderer.color);
        tmpMesh1.draw(GL_LINES);
    }
}

void Draw3D::drawLines( int nlinks, const int * links, const Vec3d * points ){
    tmpMesh1.clear();
    for(int i=0; i<nlinks*2; i+=2){
        Vec3f a = (Vec3f)points[links[i]];
        Vec3f b = (Vec3f)points[links[i+1]];
        tmpMesh1.addVertex(a);
        tmpMesh1.addVertex(b);
    }
    tmpMesh1.setUniform3f("uColor", opengl1renderer.color);
    tmpMesh1.draw(GL_LINES);
}

void Draw3D::drawTriangles( int nlinks, const int * links, const Vec3d * points, int mode ){
    if(mode == 2) { // Normal vectors
        tmpMesh1.clear();
        for(int i=0; i<nlinks*3; i+=3){
            Vec3f a = (Vec3f)points[links[i]];
            Vec3f b = (Vec3f)points[links[i+1]];
            Vec3f c = (Vec3f)points[links[i+2]];
            Vec3f nor = cross(a-b, b-c);
            nor.normalize();
            Vec3f cog = (a+b+c)*(1.0f/3.0f);
            tmpMesh1.addVertex(cog);
            tmpMesh1.addVertex(cog + nor);
        }
        tmpMesh1.setUniform3f("uColor", opengl1renderer.color);
        tmpMesh1.draw(GL_LINES);
    } else if(mode == 1) { // Wireframe
        tmpMesh1.clear();
        for(int i=0; i<nlinks*3; i+=3){
            Vec3f a = (Vec3f)points[links[i]];
            Vec3f b = (Vec3f)points[links[i+1]];
            Vec3f c = (Vec3f)points[links[i+2]];
            tmpMesh1.addVertex(a); tmpMesh1.addVertex(b);
            tmpMesh1.addVertex(b); tmpMesh1.addVertex(c);
            tmpMesh1.addVertex(c); tmpMesh1.addVertex(a);
        }
        tmpMesh1.setUniform3f("uColor", opengl1renderer.color);
        tmpMesh1.draw(GL_LINES);
    } else { // Solid triangles
        tmpMesh2.clear();
        for(int i=0; i<nlinks*3; i+=3){
            Vec3f a = (Vec3f)points[links[i]];
            Vec3f b = (Vec3f)points[links[i+1]];
            Vec3f c = (Vec3f)points[links[i+2]];
            Vec3f nor = cross(a-b, b-c);
            nor.normalize();
            tmpMesh2.addVertex(a, nor);
            tmpMesh2.addVertex(b, nor);
            tmpMesh2.addVertex(c, nor);
        }
        tmpMesh2.setUniform3f("uColor", opengl1renderer.color);
        tmpMesh2.draw(GL_TRIANGLES);
    }
}

void Draw3D::drawVectorArray(int n, const Vec3d* ps, const Vec3d* vs, double sc, double lmax){
    double l2max = sq(lmax/sc);
    tmpMesh1.clear();
    for(int i=0; i<n; i++){
        if(lmax>0){ if(vs[i].norm2()>l2max) continue; }
        Vec3d p1=ps[i];
        Vec3d p2=ps[i];
        p2.add_mul(vs[i], sc);
        tmpMesh1.addVertex((Vec3f)p1);
        tmpMesh1.addVertex((Vec3f)p2);
    }
    tmpMesh1.setUniform3f("uColor", opengl1renderer.color);
    tmpMesh1.draw(GL_LINES);
}

void Draw3D::drawVectorArray(int n, const Vec3d* ps, const Quat4f* qs, double sc, double lmax){
    double l2max = sq(lmax/sc);
    tmpMesh1.clear();
    for(int i=0; i<n; i++){
        if(lmax>0){ if(qs[i].f.norm2()>l2max) continue; }
        Vec3d p1=ps[i];
        Vec3d p2=ps[i];
        p2.add_mul((Vec3d)qs[i].f, sc);
        tmpMesh1.addVertex((Vec3f)p1);
        tmpMesh1.addVertex((Vec3f)p2);
    }
    tmpMesh1.setUniform3f("uColor", opengl1renderer.color);
    tmpMesh1.draw(GL_LINES);
}

void Draw3D::drawScalarArray(int n, const Vec3d* ps, const double* vs, double vmin, double vmax, const uint32_t* colors, int ncol ){
    tmpMesh3.clear();
    double sc = 1/(vmax-vmin);
    for(int i=0; i<n; i++){
        Vec3d p = ps[i];
        double c = (vs[i]-vmin)*sc;
        if(colors){
            uint32_t color = Draw::icolorScale(c, ncol, colors);
            tmpMesh3.addVertex((Vec3f)p, {((color>>16)&0xFF)/255.0f, ((color>>8)&0xFF)/255.0f, (color&0xFF)/255.0f});
        } else {
            tmpMesh3.addVertex((Vec3f)p, {(float)c,(float)c,(float)c});
        }
    }
    tmpMesh3.draw(GL_POINTS);
}

void Draw3D::drawScalarField(Vec2i ns, const Vec3d* ps, const double* data, double vmin, double vmax, const uint32_t* colors, int ncol){
    double clsc = 1/(vmax-vmin);
    tmpMesh3.clear();
    for(int iy=1; iy<ns.y; iy++){
        for(int ix=0; ix<ns.x; ix++){
            int i = (iy-1)*ns.x + ix;
            double c = clamp(clsc*(data[i]-vmin), 0, 1);
            Vec3f color;
            if(colors){ 
                uint32_t col = Draw::icolorScale(c, ncol, colors);
                color = {((col>>16)&0xFF)/255.0f, ((col>>8)&0xFF)/255.0f, (col&0xFF)/255.0f};
            } else {
                color = {(float)c,(float)c,(float)c};
            }
            tmpMesh3.addVertex((Vec3f)ps[i], color);

            i += ns.x;
            c = clamp(clsc*(data[i]-vmin), 0, 1);
            if(colors){
                uint32_t col = Draw::icolorScale(c, ncol, colors);
                color = {((col>>16)&0xFF)/255.0f, ((col>>8)&0xFF)/255.0f, (col&0xFF)/255.0f};
            } else {
                color = {(float)c,(float)c,(float)c};
            }
            tmpMesh3.addVertex((Vec3f)ps[i], color);
        }
    }
    tmpMesh3.draw(GL_TRIANGLE_STRIP);
}

void Draw3D::drawScalarField( Vec2i ns, const Quat4f* ps, const float* data, int pitch, int offset, double vmin, double vmax, const uint32_t* colors, int ncol ){
    printf( " Draw3D::drawScalarField() ns(%i,%i) vrange(%g,%g) @ps=%li @data=%li @colors=%li ncol=%i pitch=%i offset=%i \n", ns.x, ns.y, vmin, vmax, (long)ps, (long)data, (long)colors, ncol, pitch, offset );
    double clsc = 1/(vmax-vmin);
    if( (clsc<0)||(clsc>1e+300) ){ printf( "ERROR in drawScalarField() vrange(%g,%g) -> clsc=%g => exit() \n", vmin, vmax, clsc ); exit(0); }
    opengl1renderer.shadeModel(GL_SMOOTH);
    //opengl1renderer.enable( GL_POLYGON_SMOOTH);
    for(int iy=1;iy<ns.y;iy++){
        opengl1renderer.begin( GL_TRIANGLE_STRIP );
        for(int ix=0;ix<ns.x;ix++){
            Vec3f p;
            int i = (iy-1)*ns.x + ix;
            int ii = i*pitch+offset;
            //printf( "drawScalarField()[%i,%i] i=%i ii=%i \n", ix,iy, i, ii  );
            double c = clamp( clsc*(data[i]-vmin), 0., 1. );
            //printf( "c=%g, clsc=%g, data[i]=%g, vmin=%g \n", c, clsc, data[i], vmin );
            if(colors){ Draw::colorScale( c,ncol,colors); }else{ opengl1renderer.color3f(c,c,c); }
            //p = (gsh.dCell.a*(ix + (gsh.n.x*-0.5))) + (gsh.dCell.b*(iy-1 + (gsh.n.y*-0.5) ));
            p = ps[i].f;
            opengl1renderer.vertex3f(p.x,p.y,p.z);

            i += ns.x;
            ii = i*pitch+offset;
            c = clamp(  clsc*(data[ii]-vmin), 0., 1. );
            if(colors){ Draw::colorScale( c,ncol,colors); }else{ opengl1renderer.color3f(c,c,c); }
            p = ps[i].f;
            opengl1renderer.vertex3f(p.x,p.y,p.z);
        }
        opengl1renderer.end();
    }
}

void Draw3D::drawScalarGrid(Vec2i ns, const Vec3d& p0, const Vec3d& a, const Vec3d& b, const double* data, double vmin, double vmax, const uint32_t* colors, int ncol ){
    double clsc = 1/(vmax-vmin);
    tmpMesh3.clear();
    for(int iy=1; iy<ns.y; iy++){
        for(int ix=0; ix<ns.x; ix++){
            int i = (iy-1)*ns.x + ix;
            double c = clamp(clsc*(data[i]-vmin), 0, 1);
            Vec3f color;
            if(colors){
                uint32_t col = Draw::icolorScale(c, ncol, colors);
                color = {((col>>16)&0xFF)/255.0f, ((col>>8)&0xFF)/255.0f, (col&0xFF)/255.0f};
            } else {
                color = {(float)c,(float)c,(float)c};
            }
            Vec3d p = a*ix + b*(iy-1) + p0;
            tmpMesh3.addVertex((Vec3f)p, color);

            i += ns.x;
            c = clamp(clsc*(data[i]-vmin), 0, 1);
            if(colors){
                uint32_t col = Draw::icolorScale(c, ncol, colors);
                color = {((col>>16)&0xFF)/255.0f, ((col>>8)&0xFF)/255.0f, (col&0xFF)/255.0f};
            } else {
                color = {(float)c,(float)c,(float)c};
            }
            p.add(b);
            tmpMesh3.addVertex((Vec3f)p, color);
        }
    }
    tmpMesh3.draw(GL_TRIANGLE_STRIP);
}

void Draw3D::drawScalarGridLines(Vec2i ns, const Vec3d& p0, const Vec3d& a, const Vec3d& b, const Vec3d& up, const double* data, double sc, Vec2d vclamp ){
    for(int iy=1; iy<ns.y; iy++){
        tmpMesh1.clear();
        for(int ix=0; ix<ns.x; ix++){
            int i = iy*ns.x + ix;
            double val = _clamp(data[i], vclamp.x, vclamp.y);
            val *= sc;
            Vec3d p = p0 + a*ix + b*iy + up*val;
            tmpMesh1.addVertex((Vec3f)p);
        }
        tmpMesh1.draw(GL_LINE_STRIP);
    }

    for(int ix=0; ix<ns.x; ix++){
        tmpMesh1.clear();
        for(int iy=1; iy<ns.y; iy++){
            int i = iy*ns.x + ix;
            double val = _clamp(data[i], vclamp.x, vclamp.y) * sc;
            Vec3d p = p0 + a*ix + b*iy + up*val;
            tmpMesh1.addVertex((Vec3f)p);
        }
        tmpMesh1.draw(GL_LINE_STRIP);
    }
}

void Draw3D::drawColorScale( int n, const Vec3d& p0, const Vec3d& fw, const Vec3d& up, const uint32_t * colors, int ncol ){
    double step = 1./(n-1);
    tmpMesh3.clear();
    for(int iy=0; iy<n; iy++){
        double c = iy*step;
        uint32_t col = Draw::icolorScale(c, ncol, colors);
        Vec3f color = {((col>>16)&0xFF)/255.0f, ((col>>8)&0xFF)/255.0f, (col&0xFF)/255.0f};
        Vec3d p = fw*c + p0;
        tmpMesh3.addVertex((Vec3f)p, color);
        p.add(up);
        tmpMesh3.addVertex((Vec3f)p, color);
    }
    tmpMesh3.draw(GL_TRIANGLE_STRIP);
}

void Draw3D::drawSimplexGrid( int na, int nb, const Vec2d& da, const Vec2d& db,  const double * hs, const double * clrs, int ncolors, const uint32_t * cscale ){
    Vec2d pa; pa.set(0.0);
    if( !cscale ){ cscale=&Draw::colors_rainbow[0]; ncolors=Draw::ncolors; }
    int ii=0;
    
    for (int ia=0; ia<(na-1); ia++){
        Vec2d p; p.set(pa);
        if(clrs) {
            tmpMesh3.clear();
            for (int ib=0; ib<nb; ib++){
                double h=0.0;
                uint32_t col = Draw::icolorScale(clrs[ii], ncolors, cscale);
                Vec3f color = {((col>>16)&0xFF)/255.0f, ((col>>8)&0xFF)/255.0f, (col&0xFF)/255.0f};
                if(hs){ h=hs[ii]; }
                tmpMesh3.addVertex({(float)p.x, (float)p.y, (float)h}, color);

                col = Draw::icolorScale(clrs[ii+nb], ncolors, cscale);
                color = {((col>>16)&0xFF)/255.0f, ((col>>8)&0xFF)/255.0f, (col&0xFF)/255.0f};
                if(hs){ h=hs[ii+nb]; }
                tmpMesh3.addVertex({(float)(p.x+da.x), (float)(p.y+da.y), (float)h}, color);
                p.add(db);
                ii++;
            }
            tmpMesh3.draw(GL_TRIANGLE_STRIP);
        } else {
            tmpMesh2.clear();
            for (int ib=0; ib<nb; ib++){
                double h=0.0;
                if(hs){ h=hs[ii]; }
                Vec3f normal = {0.0f, 1.0f, 0.0f}; // Default normal pointing up
                tmpMesh2.addVertex({(float)p.x, (float)p.y, (float)h}, normal);
                if(hs){ h=hs[ii+nb]; }
                tmpMesh2.addVertex({(float)(p.x+da.x), (float)(p.y+da.y), (float)h}, normal);
                p.add(db);
                ii++;
            }
            tmpMesh2.draw(GL_TRIANGLE_STRIP);
        }
        pa.add(da);
    }
}

void Draw3D::drawSimplexGridLines( int na, int nb, const Vec2d& da, const Vec2d& db,  const double * hs ){
    Vec2d p,pa; pa.set(0.0);
    for (int ia=0; ia<(na-1); ia++){
        tmpMesh1.clear();
        p.set(pa);
        for (int ib=0; ib<nb; ib++){
            tmpMesh1.addVertex({(float)p.x, (float)p.y, (float)hs[ia*nb+ib]});
            p.add(db);
        }
        tmpMesh1.draw(GL_LINE_STRIP);

        tmpMesh1.clear();
        p.set(pa);
        for (int ib=0; ib<nb; ib++){
            int ii=ia*nb+ib;
            tmpMesh1.addVertex({(float)p.x, (float)p.y, (float)hs[ii]});
            tmpMesh1.addVertex({(float)(p.x+da.x), (float)(p.y+da.y), (float)hs[ii+nb]});
            p.add(db);
            ii++;
        }
        tmpMesh1.draw(GL_LINE_STRIP);
        pa.add(da);
    }
    p.set(pa);
    tmpMesh1.clear();
    for (int ib=0; ib<nb; ib++){
        tmpMesh1.addVertex({(float)p.x, (float)p.y, (float)hs[(na-1)*nb+ib]});
        p.add(db);
    }
    tmpMesh1.draw(GL_LINE_STRIP);
}

void Draw3D::drawSimplexGridLinesToned( int na, int nb, const Vec2d& da, const Vec2d& db,  const double * hs ){
    Vec2d p,pa; pa.set(0.0);
    for (int ia=0; ia<(na-1); ia++){
        tmpMesh3.clear();
        p.set(pa);
        for (int ib=0; ib<nb; ib++){
            float h = (float)hs[ia*nb+ib];
            tmpMesh3.addVertex({(float)p.x, (float)p.y, h}, {h,h*4,h*16});
            p.add(db);
        }
        tmpMesh3.draw(GL_LINE_STRIP);

        tmpMesh3.clear();
        p.set(pa);
        for (int ib=0; ib<nb; ib++){
            int ii=ia*nb+ib;
            float h = (float)hs[ii];
            tmpMesh3.addVertex({(float)p.x, (float)p.y, h}, {h,h*4,h*16});
            h = (float)hs[ii+nb];
            tmpMesh3.addVertex({(float)(p.x+da.x), (float)(p.y+da.y), h}, {h,h*4,h*16});
            p.add(db);
            ii++;
        }
        tmpMesh3.draw(GL_LINE_STRIP);
        pa.add(da);
    }
    p.set(pa);
    tmpMesh3.clear();
    for (int ib=0; ib<nb; ib++){
        float h = (float)hs[(na-1)*nb+ib];
        tmpMesh3.addVertex({(float)p.x, (float)p.y, h}, {h,h*4,h*16});
        p.add(db);
    }
    tmpMesh3.draw(GL_LINE_STRIP);
}

void Draw3D::drawRectGridLines( Vec2i n, const Vec3d& p0, const Vec3d& da, const Vec3d& db ){
    tmpMesh1.clear();
    Vec3d p  = p0;
    Vec3d dn = db*n.b;
    for (int ia=0; ia<n.a; ia++){
        tmpMesh1.addVertex((Vec3f)p);
        Vec3d p_ = p+dn;
        tmpMesh1.addVertex((Vec3f)p_);
        p.add(da);
    }

    p   = p0;
    dn  = da*n.a;
    for (int ib=0; ib<n.b; ib++){
        tmpMesh1.addVertex((Vec3f)p);
        Vec3d p_ = p+dn;
        tmpMesh1.addVertex((Vec3f)p_);
        p.add(db);
    }
    tmpMesh1.setUniform3f("uColor", opengl1renderer.color);
    tmpMesh1.draw(GL_LINES);
}

void Draw3D::drawTextBillboard( const char* str, Vec3f pos, float sz, bool ontop, int iend ){ // TODO: sz is unused
    Vec3f screen_pos = GLES::active_camera->world2Screen(pos);
    Draw::drawText(str, ontop ? (Vec3f){screen_pos.x, screen_pos.y, -1} : screen_pos, fontSizeDef, iend);
}

void Draw3D::drawInt( const Vec3d& pos, int i, int fontTex, float sz, const char* format ){
    char str[16]; sprintf(str,format,i); //printf("%s\n", str);
    Draw3D::drawText(str, pos, fontTex, sz, 0);
}
void Draw3D::drawDouble( const Vec3d& pos, double f, int fontTex, float sz, const char* format ){
    char str[24];  sprintf(str,format,f);
    Draw3D::drawText(str, pos, fontTex, sz, 0);
}
void Draw3D::pointLabels( int n, const Vec3d* ps, int fontTex, float sz ){
    for(int i=0; i<n; i++){ drawInt( ps[i], i, fontTex, sz ); }
}

void Draw3D::drawAxis3D( int n, Vec3d p0, Vec3d dp, double v0, double dval, int fontTex, float tickSz, float textSz, const char* format ){
    Vec3d a,b;
    dp.getSomeOrtho( a, b );
    Vec3d p = p0;
    // tick marks
    opengl1renderer.begin(GL_LINES);
    vertex(p); vertex(p+dp*n);
    for(int i=0; i<=n; i++){
        vertex(p-a*tickSz); vertex(p+a*tickSz);
        vertex(p+b*tickSz); vertex(p+b*tickSz);
        p.add(dp);
    }
    opengl1renderer.end();
    // labels
    p=p0;
    double val = v0;
    char str[64];
    for(int i=0; i<=n; i++){
        sprintf(str,format,val);
        //printf( "drawAxis3D()[%i] str(%s) \n", i, str );
        //drawText(p, a, b, str, fontTex, textSz );
        Draw3D::drawText(str, p, fontTex, textSz, 0);
        p.add(dp);
        val+=dval;
    }
}
void Draw3D::drawAxis3D( Vec3i ns, Vec3d p0, Vec3d ls, Vec3d v0s, Vec3d dvs, int fontTex, float tickSz, float textSz, const char* format ){
    //drawAxis3D( ns.x, p0, Vec3dX*ls.x, v0s.x, dvs.x, fontTex, tickSz, textSz, format );
    //drawAxis3D( ns.y, p0, Vec3dY*ls.y, v0s.y, dvs.y, fontTex, tickSz, textSz, format );
    //drawAxis3D( ns.z, p0, Vec3dZ*ls.z, v0s.z, dvs.z, fontTex, tickSz, textSz, format );
    drawAxis3D( ns.x, {p0.x,.0,.0}, Vec3dX*ls.x, v0s.x, dvs.x, fontTex, tickSz, textSz, format );
    drawAxis3D( ns.y, {0.,p0.x,.0}, Vec3dY*ls.y, v0s.y, dvs.y, fontTex, tickSz, textSz, format );
    drawAxis3D( ns.z, {0.,0.,p0.z}, Vec3dZ*ls.z, v0s.z, dvs.z, fontTex, tickSz, textSz, format );
}

void Draw3D::drawCurve( float tmin, float tmax, int n, Func1d3 func ){
    tmpMesh1.clear();
    float dt = (tmax-tmin)/n;
    for( float t=tmin; t<=tmax; t+=dt ){
        double x,y,z;
        func( t, x, y, z );
        tmpMesh1.addVertex({(float)x, (float)y, (float)z});
    }
    tmpMesh1.setUniform3f("uColor", opengl1renderer.color);
    tmpMesh1.draw(GL_LINE_STRIP);
}

void Draw3D::drawBox( float x0, float x1, float y0, float y1, float z0, float z1, float r, float g, float b ){
	MeshLibrary::cubeWithNormals.setUniform3f("uColor", {r, g, b});
    MeshLibrary::cubeWithNormals.draw((Vec3f){x0, y0, z0}, (Vec3f){x1-x0, y1-y0, z1-z0});
}

void Draw3D::drawBBox( const Vec3f& p0, const Vec3f& p1 ){
	MeshLibrary::wireCube.setUniform3f("uColor", opengl1renderer.color);
    MeshLibrary::wireCube.draw(p0, p1-p0);
}

void Draw3D::drawBBox( const Vec3f& p, float r ){ 
    drawBBox( Vec3f{p.x-r,p.y-r,p.z-r}, Vec3f{p.x+r,p.y+r,p.z+r} ); 
}

void Draw3D::drawTriclinicBox( const Mat3f& lvec, const Vec3f& c0, const Vec3f& c1 ){
    Vec3f p0,p1;
    tmpMesh1.clear();
           lvec.dot_to({c0.x,c0.y,c0.z},p0);
           lvec.dot_to({c0.x,c0.y,c1.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
           lvec.dot_to({c0.x,c1.y,c0.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
           lvec.dot_to({c1.x,c0.y,c0.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
           lvec.dot_to({c1.x,c1.y,c1.z},p0);
           lvec.dot_to({c0.x,c1.y,c1.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
           lvec.dot_to({c1.x,c0.y,c1.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
           lvec.dot_to({c1.x,c1.y,c0.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
    p0=p1; lvec.dot_to({c1.x,c0.y,c0.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
    p0=p1; lvec.dot_to({c1.x,c0.y,c1.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
    p0=p1; lvec.dot_to({c0.x,c0.y,c1.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
    p0=p1; lvec.dot_to({c0.x,c1.y,c1.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
    p0=p1; lvec.dot_to({c0.x,c1.y,c0.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
    p0=p1; lvec.dot_to({c1.x,c1.y,c0.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);

    tmpMesh1.draw(GL_LINES);
}

void Draw3D::drawTriclinicBoxT( const Mat3f& lvec, const Vec3f& c0, const Vec3f& c1 ){
    Vec3f p0,p1;
    tmpMesh1.clear();
           lvec.dot_to_T({c0.x,c0.y,c0.z},p0);
           lvec.dot_to_T({c0.x,c0.y,c1.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
           lvec.dot_to_T({c0.x,c1.y,c0.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
           lvec.dot_to_T({c1.x,c0.y,c0.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
           lvec.dot_to_T({c1.x,c1.y,c1.z},p0);
           lvec.dot_to_T({c0.x,c1.y,c1.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
           lvec.dot_to_T({c1.x,c0.y,c1.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
           lvec.dot_to_T({c1.x,c1.y,c0.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
    p0=p1; lvec.dot_to_T({c1.x,c0.y,c0.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
    p0=p1; lvec.dot_to_T({c1.x,c0.y,c1.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
    p0=p1; lvec.dot_to_T({c0.x,c0.y,c1.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
    p0=p1; lvec.dot_to_T({c0.x,c1.y,c1.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
    p0=p1; lvec.dot_to_T({c0.x,c1.y,c0.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);
    p0=p1; lvec.dot_to_T({c1.x,c1.y,c0.z},p1); tmpMesh1.addVertex(p0); tmpMesh1.addVertex(p1);

    tmpMesh1.draw(GL_LINES);
}

void Draw3D::drawAxis( float sc ){
    tmpMesh3.clear();
    // X axis (red)
    tmpMesh3.addVertex({ 0,0,0}, {1,0,0});
    tmpMesh3.addVertex({sc,0,0}, {1,0,0});
    // Y axis (green)
    tmpMesh3.addVertex({0, 0,0}, {0,1,0});
    tmpMesh3.addVertex({0,sc,0}, {0,1,0});
    // Z axis (blue)
    tmpMesh3.addVertex({0,0, 0}, {0,0,1});
    tmpMesh3.addVertex({0,0,sc}, {0,0,1});
    tmpMesh3.draw(GL_LINES);
}

