
#ifndef  DrawSphereMap_h
#define  DrawSphereMap_h

#include "SphereSampling.h"

void heightColor(double h){
    //glColor3f( h, h*0.3, 0.5+h*0.1 );
    //glColor3f( h*0.25, 0.5,1-h*0.25 );
    //double sh = sqrt(h);
    //if(h>0){ glColor3f( 0.0, 0.0, sqrt(h)*0.25 ); }else{   glColor3f( sqrt(-h)*0.25, 0.0, 0.0 ); };
    float c = h*0.05+0.5;
    if(h>0){ float f =sqrt(h)*0.25;  glColor3f( c, c, c+f); }else{ float f=sqrt(-h)*0.25;  glColor3f( c+f, c, c ); };
    //if(h>0){ glColor3f( 0.0, sin(h*15)*0.5+0.5, sqrt(h)*0.25 ); }else{   glColor3f( sqrt(-h)*0.25, sin(h*15)*0.5+0.5, 0.0 ); };
};

namespace SphereSampling{

bool bNormalize=false;
bool bRelief=false;

void drawDiTri( Vec2i n, const Vec3f& a, const Vec3f& b, const Vec3f& c, const Vec3f& d, float* hs, double hscale ){
    float da = 1.0/n.a;
    float db = 1.0/n.b;
    /*
    glBegin(GL_TRIANGLES);
    //glVertex3f(a.x,a.y,a.z);  glVertex3f(b.x,b.y,b.z);  glVertex3f(c.x,c.y,c.z);
    //glVertex3f(a.x,a.y,a.z);  glVertex3f(b.x,b.y,b.z);  glVertex3f(d.x,d.y,d.z);
    glColor3f(1.0,0.0,0.0); glVertex3f(a.x,a.y,a.z);
    glColor3f(0.0,0.0,1.0); glVertex3f(b.x,b.y,b.z);
    glColor3f(0.0,1.0,0.0); glVertex3f(c.x,c.y,c.z);
    glColor3f(1.0,0.0,0.0); glVertex3f(a.x,a.y,a.z);
    glColor3f(0.0,0.0,1.0); glVertex3f(b.x,b.y,b.z);
    glColor3f(0.0,0.0,0.0); glVertex3f(d.x,d.y,d.z);
    glEnd();
    return;
    */

    //glPointSize(2);
    //glPointSize(1);

    //int ioff = (hs-heights);
    //printf("ioff %i %i \n",ioff,ioff/(n.a*n.b));
    // Vec3f* sps = samplePs + ioff;  // a bit hack ... just debug

    for( int ia = 0; ia<n.a-1; ia++ ){
        float fa  =da*(ia);
        float fa_ =da*(ia+1);
        glBegin(GL_TRIANGLE_STRIP);
        //glBegin(GL_POINTS);
        for( int ib = 0; ib<n.b; ib++ ){
            Vec3f p;
            float h, fd;
            float fb  =(db*ib);
            int i;
            //fd = (fa_+fb)*0.5;
            //if( fa_>fb ){ float f=(fa_-fb); p = d*(fd-0.5*f) + c*(1-fd-0.5*f)   + a*f; }
            //else        { float f=(fb-fa_); p = d*(fd-0.5*f) + c*(1-fd-0.5*f)   + b*f; }
            diTri2cartes<float>( fa_, fb, a,b,c,d, p);

            i = (ia+1)*n.b+ib; //if(i>(n.a*n.b)) printf("i = %i", i);
            h = hs[i];
            //glColor3f(h,h,h);
            heightColor(h);
            //glColor3f(fa_,h,fb);
            if(bNormalize) p.normalize();
            if(bRelief)    p.mul(1+hscale*(h-0.5));
            glVertex3f( p.x,  p.y,  p.z  );
            //sps[i] = p;
            //fd = (fa+fb)*0.5;
            //if( fa_>fb ){ float f=(fa-fb); p = d*(fd-0.5*f) + c*(1-fd-0.5*f)   + a*f; }
            //else        { float f=(fb-fa); p = d*(fd-0.5*f) + c*(1-fd-0.5*f)   + b*f; }
            diTri2cartes<float>( fa, fb, a,b,c,d, p);

            i = ia*n.b+ib;
            h = hs[i];
            //glColor3f(h,h,h);
            heightColor(h);
            //glColor3f(fa,h,fb);
            if(bNormalize) p.normalize();
            if(bRelief)    p.mul(1+hscale*(h-0.5));
            glVertex3f( p.x,  p.y,  p.z  );
            //sps[i] = p;
            //sprintf( str,"%i,%i|%i", ia,ib,i );
            //Draw3D::drawText( str, (Vec3d)p, fontTex, 0.02, 0 );
            //glVertex3f( p_.x, p_.y, p_.z );
            //printf( "%i: %i %i %i | (%f,%f,%f) \n", i, ia, ib, iface, p.x, p.y, p.z );
        }
        glEnd();
    }
}

void drawDiTriWire( Vec2i n, const Vec3f& a, const Vec3f& b, const Vec3f& c, const Vec3f& d, float* hs, float hscale ){
    float da = 1.0/n.a;
    float db = 1.0/n.b;
    for( int ia = 0; ia<n.a; ia++ ){
        float fa  =da*(ia);
        float fa_ =da*(ia+1);

        //glBegin(GL_POINTS);

        //glColor3f(1.0,0.0,0.0);
        glBegin(GL_LINE_STRIP);
        for( int ib = 0; ib<n.b; ib++ ){
            Vec3f p;
            float h, fd;
            float fb  =(db*ib);
            int i;
            diTri2cartes<float>( fa, fb, a,b,c,d, p);
            i = ia*n.b+ib; //if(i>(n.a*n.b)) printf("i = %i", i);
            h = hs[i];
            if(bNormalize) p.normalize();
            if(bRelief)    p.mul(1+hscale*(h-0.5));
            glVertex3f( p.x,  p.y,  p.z  );
        }
        glEnd();

        if(ia>=(n.a-1)) continue;
        glBegin(GL_LINE_STRIP);
        //glColor3f(0.0,0.0,1.0);
        for( int ib = 0; ib<n.b; ib++ ){
            Vec3f p;
            float h, fd;
            float fb  =(db*ib);
            int i;
            diTri2cartes<float>( fa_, fb, a,b,c,d, p);
            i = (ia+1)*n.b+ib; //if(i>(n.a*n.b)) printf("i = %i", i);
            h = hs[i];
            if(bNormalize) p.normalize();
            if(bRelief)    p.mul(1+hscale*(h-0.5));
            glVertex3f( p.x,  p.y,  p.z  );
            diTri2cartes<float>( fa, fb, a,b,c,d, p);
            i = ia*n.b+ib;
            h = hs[i];
            if(bNormalize) p.normalize();
            if(bRelief)    p.mul(1+hscale*(h-0.5));
            glVertex3f( p.x,  p.y,  p.z  );
        }
        glEnd();

    }
}



/*
template<typename T>
inline T& arrayView( int i, Vec2i view, T * data ){
    return data[i*view.x+view.y];
}
*/

inline int index( int i, Vec2i view ){ return (i*view.x)+view.y; };

void drawDiTri_seam( int n, int n2, const Vec3f& a, const Vec3f& b, const Vec3f& c, float* hs, float* hs2, const Vec2i& view, const Vec2i& view2, float hn, float hscale ){

    Vec3f db = (b-a)*(1.0/(n ));
    Vec3f dc = (c-a)*(1.0/(n2));
    Vec3f p=a;
    glBegin(GL_TRIANGLE_STRIP);
    //glBegin(GL_LINES);
    //glBegin(GL_POINTS);
    for( int i = 0; i<(n+1); i++ ){
        float h;
        int ii;
        Vec3f p_;
        ii = index(n-i-1, view2);
        if(i==n){ h=hn; }else{ h=hs2 [ii]; };

        //glColor3f(h,h,h);
        heightColor(h);
        //glColor3f(1.0,0.0,0.0);
        p_=p;
        //p_=sps2[ii];
        if(bNormalize) p_.normalize();
        if(bRelief)    p_.mul(1+hscale*(h-0.5));
        glVertex3f( p_.x,  p_.y,  p_.z  );
        //glColor3f(1.0,1.0,0.0); glVertex3f( pp.x,  pp.y,  pp.z  );
        //printf( "(%f,%f,%f)  (%f,%f,%f) \n", p.x, p.y, p.z,  pp.x, pp.y, pp.z );

        if(i<n){
            //if(bT1){ h = hs[i*n]; }else{ h = hs[i]; };
            //h = hs2[ia*n.b+ib];
            ii = index(n-i-1, view);
            h=hs [ii]; //p_=sps[ii];
            //p=sps[ii];
            //printf( "sps  %i %i %f %f %f \n", i, ii, sps[ii].x,  sps[ii].y,  sps[ii].z );
            //glColor3f(fa,h,fb);
            //glColor3f(h,h,h);
            heightColor(h);
            //glColor3f(0.0,0.0,1.0);
            p_ = p + dc;
            //p_=sps[ii];
            if(bNormalize) p_.normalize();
            if(bRelief)    p_.mul(1+hscale*(h-0.5));
            glVertex3f( p_.x,  p_.y,  p_.z  );
            //glColor3f(1.0,1.0,0.0); glVertex3f( pp.x,  pp.y,  pp.z  );
            //glVertex3f( p.x,  p.y,  p.z  );
            //Draw3D::drawPointCross(p,0.05);
        }
        //if(bT1){ h = hs2[i*n]; }else{ h = hs2[i]; };

        //printf( "%i %f %f %f \n", i, p.x,  p.y,  p.z );
        //printf( "sps2 %i %i %f %f %f \n", i, ii, sps2[ii].x,  sps2[ii].y,  sps2[ii].z );
        p.add(db);
    }
    glEnd();
    //#undef index

}



void drawDiTriWire_seam( int n, int n2, const Vec3f& a, const Vec3f& b, const Vec3f& c, float* hs, float* hs2, const Vec2i& view, const Vec2i& view2, float hn, float hscale ){
    Vec3f db = (b-a)*(1.0/(n ));
    Vec3f dc = (c-a)*(1.0/(n2));
    Vec3f p=a;
    glBegin(GL_LINE_STRIP);
    for( int i = 0; i<(n+1); i++ ){
        float h;
        int ii;
        Vec3f p_;
        ii = index(n-i-1, view2);
        if(i==n){ h=hn; }else{ h=hs2[ii]; };
        p_=p;
        if(bNormalize) p_.normalize();
        if(bRelief)    p_.mul(1+hscale*(h-0.5));
        glVertex3f( p_.x,  p_.y,  p_.z  );
        if(i<n){
            ii = index(n-i-1, view);
            h=hs [ii]; //p_=sps[ii];
            p_ = p + dc;
            if(bNormalize) p_.normalize();
            if(bRelief)    p_.mul(1+hscale*(h-0.5));
            glVertex3f( p_.x,  p_.y,  p_.z  );
        }
        p.add(db);
    }
    glEnd();

}




/*
void drawDiTri_Inner( Vec2i n, const Vec3f& a, const Vec3f& b, const Vec3f& c, const Vec3f& d, float* height, float hscale ){
    float da = 1.0/(n.a);
    float db = 1.0/(n.b);
    for( int ia = 1; ia<n.a; ia++ ){
        float fa  =da*(ia-0.5);
        float fa_ =da*(ia+0.5);
        glBegin(GL_TRIANGLE_STRIP);
        for( int ib = 0; ib<n.b; ib++ ){
            float h;
            float fb  =(db*(ib+0.5));
            Vec3f p;

            if( (fa+fb)<1 ){ p = a*fa      + b*fb     + c*(1-fa-fb); }
            else           { p = a*(1-fb)  + b*(1-fa) + d*(fa+fb-1); }
            h = height[ia*n.b+ib];
            heightColor(h);
            //glColor3f(fa,0.0,fb);
            //p.normalize(); p.mul(1+0.1*(h-0.5));
            glVertex3f( p.x,  p.y,  p.z  );

            if( (fa_+fb)<1 ){ p = a*fa_    + b*fb     + c*(1-fa_-fb); }
            else            { p = a*(1-fb) + b*(1-fa_) + d*(fa_+fb-1); }
            h = height[(ia+1)*n.b+ib];
            //glColor3f(h,h,h);
            heightColor(h);
            //glColor3f(fa_,0.0,fb);
            //p.normalize(); p.mul(1+0.1*(h-0.5));
            glVertex3f( p.x,  p.y,  p.z  );
        }
        glEnd();
    }
}
*/

/*
void drawDiTri_edge( Vec2i n, const Vec3f& a, const Vec3f& b, const Vec3f& c, const Vec3f& d, bool bT1, bool bT2, float* hs2, float* hs1 ){
    float da = 1.0/(n.a);
    float db = 1.0/(n.b);
    glBegin(GL_TRIANGLE_STRIP);
    for( int ib = 0; ib<n.b; ib++ ){
        float fb  =(db*(ib+0.5));
        float h;
        Vec3f p;

        p = a*(1-fb) + b*fb + c*(da*0.5);
        if(bT1){ h = height[ib*n.a]; }else{ h = height[ib]; }; // beginnings
        glColor3f(h,h,h);
        //p.normalize(); p.mul(1+0.1*(h-0.5));
        glVertex3f( p.x,  p.y,  p.z  );

        p = a*(1-fb) + b*fb + d*(da*0.5);
        if(bT2){ h = height[ib*n.a+1]; }else{ h = height[(n.a-1)+ib]; }; // ends
        glColor3f(h,h,h);
        //p.normalize(); p.mul(1+0.1*(h-0.5));
        glVertex3f( p.x,  p.y,  p.z  );
    }
    glEnd();
}
*/

void drawIcosaMap( Vec2i n, float* heights, float hscale ){
    int nab = n.a*n.b;
    Quat4i* f2v = (Quat4i*)oct_fac2verts;
    Quat4i* f2n = (Quat4i*)oct_fac2neigh;
    Vec3d* vs = (Vec3d*)oct_polar_verts;
    bNormalize=true;
    bRelief   =true;
    for(int i=0; i<5; i++){
        int i2 = i+1; if(i2>=5) i2=0;
        Quat4i& iv  = f2v[i];
        Quat4i& iv2 = f2v[i+5];
        float* hs   = heights+(nab*i);
        float* hs2  = heights+(nab*i2);
        drawDiTri( n, (Vec3f)vs[iv.z], (Vec3f)vs[iv.w], (Vec3f)vs[iv.x], (Vec3f)vs[iv.y],      hs,       hscale  );
        drawDiTri( n, (Vec3f)vs[iv2.z], (Vec3f)vs[iv2.w], (Vec3f)vs[iv2.x], (Vec3f)vs[iv2.y],  hs+nab*5, hscale  );

        drawDiTri_seam( n.a, n.b, (Vec3f)vs[iv.y],  (Vec3f)vs[iv.w],  (Vec3f)vs[iv.x],  hs,       hs2,       {n.a,(n.b-1)},  {-1,n.b-1},          0.0           , hscale );
        drawDiTri_seam( n.a, n.b, (Vec3f)vs[iv.z],  (Vec3f)vs[iv.y],  (Vec3f)vs[iv.x],  hs,       hs2+5*nab, {-1,n.a*n.b-1}, {-n.a,n.a*(n.a-1)},  hs2[0]        , hscale );
        drawDiTri_seam( n.a, n.b, (Vec3f)vs[iv2.z], (Vec3f)vs[iv2.y], (Vec3f)vs[iv2.x], hs+5*nab, hs,        {-1,n.a*n.b-1}, {-n.a,n.a*(n.a-1)}, (hs2+5*nab)[0] , hscale );
        drawDiTri_seam( n.a, n.b, (Vec3f)vs[iv2.y], (Vec3f)vs[iv2.w], (Vec3f)vs[iv2.x], hs+5*nab, hs2+5*nab, {n.a,(n.b-1)},  {-1,n.b-1},          0.0           , hscale );
    }
}



void drawIcosaMapWire( Vec2i n, float* heights, float hscale ){
    int nab = n.a*n.b;
    Quat4i* f2v = (Quat4i*)oct_fac2verts;
    Quat4i* f2n = (Quat4i*)oct_fac2neigh;
    Vec3d* vs = (Vec3d*)oct_polar_verts;
    bNormalize=true;
    bRelief   =true;
    for(int i=0; i<5; i++){
        int i2 = i+1; if(i2>=5) i2=0;
        Quat4i& iv  = f2v[i];
        Quat4i& iv2 = f2v[i+5];
        float* hs   = heights+(nab*i);
        float* hs2  = heights+(nab*i2);
        drawDiTriWire( n, (Vec3f)vs[iv.z],  (Vec3f)vs[iv.w],  (Vec3f)vs[iv.x],  (Vec3f)vs[iv.y],      hs      , hscale );
        drawDiTriWire( n, (Vec3f)vs[iv2.z], (Vec3f)vs[iv2.w], (Vec3f)vs[iv2.x], (Vec3f)vs[iv2.y],  hs+nab*5, hscale );
        drawDiTriWire_seam( n.a, n.b, (Vec3f)vs[iv.y],  (Vec3f)vs[iv.w],  (Vec3f)vs[iv.x],  hs,       hs2,       {n.a,(n.b-1)},  {-1,n.b-1},          0.0           , hscale );
        drawDiTriWire_seam( n.a, n.b, (Vec3f)vs[iv.z],  (Vec3f)vs[iv.y],  (Vec3f)vs[iv.x],  hs,       hs2+5*nab, {-1,n.a*n.b-1}, {-n.a,n.a*(n.a-1)},  hs2[0]        , hscale );
        drawDiTriWire_seam( n.a, n.b, (Vec3f)vs[iv2.z], (Vec3f)vs[iv2.y], (Vec3f)vs[iv2.x], hs+5*nab, hs,        {-1,n.a*n.b-1}, {-n.a,n.a*(n.a-1)}, (hs2+5*nab)[0] , hscale );
        drawDiTriWire_seam( n.a, n.b, (Vec3f)vs[iv2.y], (Vec3f)vs[iv2.w], (Vec3f)vs[iv2.x], hs+5*nab, hs2+5*nab, {n.a,(n.b-1)},  {-1,n.b-1},          0.0           , hscale );
    }
}

}


#endif

