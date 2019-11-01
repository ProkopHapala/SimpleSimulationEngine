
#include "stdio.h"

//#include "VectorTypes.h"
#include "Vec3math.h"
#include "Mat3math.h"
#include "fastMath.h"

//#include "Camera.h"
//#include "GLView.h"

//#include "Vec4math.h"

//#include "Camera.h"
//#include "Camera.c"

//#include "Draw3D.c"
#include "GLView.c"

int main(){

    //float3x3 m1 = f3x3_Eye;

    float3  p  = {1.13,2.6,3.5};
    double3 pd = {2.13,3.6,4.5};
    //float3 p; p = (float3){1.13,2.6,3.5};
    //p.x = 1.13;
    //p.y = 2.6;
    //p.z = 3.2;

    //double3 pd0 = d3_Zero;

    float  fs[] = {1.1,3.2, 5.1,6.3,9,7,9,9.2};
    double ds[] = {1.5,3.3,-5.1,6.3,9,7,9,9.2};

    f3_read_array    ( &p, fs );
    p = f3_from_array( fs );

    d3_read_array     ( &pd, ds );
    pd = d3_from_array( ds );

    printf( " Hello World !!! \n" );
    printf( " p (%g,%g,%g) \n", p.x, p.y, p.z );
    float3 p_ = f3_add3f( p, 1.0, 2.3, 4.5 );
    printf( " p (%g,%g,%g) \n", p.x, p.y, p.z );
    printf( " p_(%g,%g,%g) \n", p_.x, p_.y, p_.z );

    printf( " pd(%g,%g,%g) \n", pd.x, pd.y, pd.z );
    double3 pd2 = d3_from_f3( p_ );
    double3 pd3 = d3_add( pd, pd2 );
    printf( " pd2(%g,%g,%g) \n", pd2.x, pd2.y, pd2.z );
    printf( " pd3(%g,%g,%g) \n", pd3.x, pd3.y, pd3.z );


    printf( " f3_norm2(p) %g f3_dot(p) %g \n", f3_norm2(p), f3_dot(p,p_) );

    GLView glview = GLView_( 1, 640, 480 );

    GLView_loop( &glview, 1000 );

}

