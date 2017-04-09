#ifndef  interpolate2D_h
#define  interpolate2D_h

//============== Globals

constexpr int nPoints   = 250;
Vec2f         points[nPoints];
float         vals  [nPoints];
float         vals_ [nPoints];

constexpr int nx   = 256;
constexpr int ny   = 256;
constexpr int ntot = nx*ny;
float hf[ntot];

void genXOR2D(int nx, int ny, float * data){
    float scale_val = 1.0/256.0;
    for(int iy=0;iy<ny;iy++){ for(int ix=0;ix<nx;ix++){
        data[iy*nx+ix] = ((ix^iy)&0xFF)*scale_val;
        //data[iy*nx+ix] = 0.5;
    } }
}

void genRCOS(int nx, int ny, float * data, float freq){
    Vec2f center = {0.5f*nx,0.5f*ny};
    for(int iy=0;iy<ny;iy++){ for(int ix=0;ix<nx;ix++){
        Vec2f p = {1.f*nx,1.0f*ny};
        float r = sqrt( sq(ix-center.x) + sq(iy-center.y) );
        data[iy*nx+ix] = sin(freq*r)*0.4f - 0.5f;
    } }
}

void genZero(int nx, int ny, float * data){
    float scale_val = 1.0/256.0;
    for(int iy=0;iy<ny;iy++){ for(int ix=0;ix<nx;ix++){  data[iy*nx+ix] = 0;  } }
}

void setArray(int n, float * arr, float f){ for(int i=0;i<n;i++){ arr[i] = f; } }

void pointsOnLine( const Vec2f& r0, const Vec2f& dr, int nPoints, Vec2f* points ){
    Vec2f p = r0;
    for(int i=0; i<nPoints; i++){
        points[i] = p;
        p.add(dr);
    }
}

void lerp(int nx, int ny, float * data, int nPoints, Vec2f* points, float * vals ){
    for(int i=0; i<nPoints; i++){
        Vec2f p = points[i];
        int ix  = (int)p.x; float dx = p.x - ix; float mx = 1-dx;
        int iy  = (int)p.y; float dy = p.y - iy;
        int ip  = nx*iy + ix;
        float val;
        val  = (1-dy)*( mx*data[ip] - dx*data[ip+1] );   ip+=nx;
        val +=    dy *( mx*data[ip] - dx*data[ip+1] );
        vals[i] = val;
    }
}

double checkDiff( int n, float * ps, float * p0s ){
    double errSum2 = 0;
    for(int i=0; i<n; i++){
        float err2 = sq(ps[i]-p0s[i]);
        errSum2   += err2;
        if ( err2 > 1e-8 ){
            printf( "%i is %g should be %g \n", i, ps[i], p0s[i] );
            exit(0);
            break;

        }
    }
    return errSum2;
}

#endif

