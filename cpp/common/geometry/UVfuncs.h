
#ifndef  UVfuncs_h
#define  UVfuncs_h


template<typename UVfunc>
Vec3f getUVFuncNormal( Vec2f uv, float eps,  UVfunc func ){
    Vec2f o;
    Vec3f nor,da,db;
    o=uv; o.a+=eps; da.set(func(o));
    o=uv; o.a-=eps; da.sub(func(o));
    o=uv; o.b+=eps; db.set(func(o));
    o=uv; o.b-=eps; db.sub(func(o));
    nor.set_cross(db,da); nor.normalize();
    return nor;
}

// =========  Quadric

inline Vec3f ConeUVfunc( Vec2f p, float R1, float R2, float L ){
    Vec2f csb; csb.fromAngle(p.b);
    float R = (1-p.a)*R1 + p.a*R2;
    return (Vec3f){csb.a*R,csb.b*R,L*p.a };
}

inline Vec3f SphereUVfunc( Vec2f p, float R ){
    Vec2f csa; csa.fromAngle(p.a);
    Vec2f csb; csb.fromAngle(p.b);
    return (Vec3f){csa.a*csb.a*R,csa.a*csb.b*R,csa.b*R };
}

inline Vec3f ParabolaUVfunc( Vec2f p, float K ){
    Vec2f csb; csb.fromAngle(p.b);
    float r = p.a;
    float l = p.a*p.a*K;
    return (Vec3f){csb.a*r,csb.b*r,l };
}

inline Vec3f HyperbolaRUVfunc( Vec2f p, float R, float K ){
    Vec2f csb; csb.fromAngle(p.b);
    float r = p.a;
    float l = K*sqrt( p.a*p.a + R*R); // - K*R;
    return (Vec3f){csb.a*r,csb.b*r,l };
}

inline Vec3f HyperbolaLUVfunc( Vec2f p, float R, float K ){
    Vec2f csb; csb.fromAngle(p.b);
    float l = p.a;
    float r = K*sqrt( p.a*p.a + R*R);
    return (Vec3f){csb.a*r,csb.b*r,l };
}


// =========  Quartic


inline Vec3f TorusUVfunc( Vec2f p, float r, float R ){
    Vec2f csa; csa.fromAngle(p.a);
    Vec2f csb; csb.fromAngle(p.b);
    return (Vec3f){csb.a*(R+r*csa.a),csb.b*(R+r*csa.a),r*csa.b};
}

inline Vec3f TeardropUVfunc( Vec2f p, float R1, float R2, float L ){
    Vec2f csa; csa.fromAngle(p.a*M_PI-M_PI*0.5);
    Vec2f csb; csb.fromAngle(p.b);
    float f =  0.5-csa.b*0.5;
    float R = (1-f)*R1 + f*R2;
    return (Vec3f){csa.a*csb.a*R,csa.a*csb.b*R,csa.b*R-L*f };
}

inline Vec3f HarmonicTubeUVfunc( Vec2f p, float R1, float R2, float L, float freq, float amp ){
    Vec2f csb; csb.fromAngle(p.b);
    float R = ((1-p.a)*R1 + p.a*R2)*(1.0+amp*cos(freq*p.a));
    //float R = (1-p.a)*R1 + p.a*R2;
    return (Vec3f){ csb.a*R, csb.b*R, L*p.a };
}





// ======== Special

inline Vec2f naca4digit( float u, float *coefs ){
    float c = coefs[0];
    float t = coefs[1];
    //float p  = coefs.b;
    //float m  = coefs.c;
    //if(u<0) t = -t;
    float x,y;
    if(u<0)t=-t;
    u  = fabs(u);
    x  =  u*u;
    //y = 5*t*u*(1-u);
    y =  5*t*( +0.2969*u+x*(-0.1260+x*(-0.3516+x*( +0.2843 - 0.1015*x ))));
    return {c*x,c*y};
}

Vec3f NACA4digitUVfunc( Vec2f p, float *coefs1, float *coefs2, float L ){
    Vec2f p1= naca4digit( p.a, coefs1 );
    Vec2f p2= naca4digit( p.a, coefs2 );
    float m = 1-p.b;
    //return (Vec3f){ m*p1.x+p.b*p2.x, m*p1.y+p.b*p2.y, L*p.b };
    return (Vec3f){ L*p.b, m*p1.y+p.b*p2.y, -(m*p1.x+p.b*p2.x)  };
}



#endif


