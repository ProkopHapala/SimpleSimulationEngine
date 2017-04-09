#ifndef  interpolate2D_h
#define  interpolate2D_h

//============== Globals

constexpr int nPoints   = 1<<12;
Vec2f         points[nPoints];
float         vals  [nPoints];
float         vals_ [nPoints];
Vec2f         Dvals [nPoints];

constexpr int nx   = 256;
constexpr int ny   = 256;
constexpr int ntot = nx*ny;
float hf [ntot  ];
float Dhf[ntot*2];

void arrayDerivs2D( int nx, int ny, float * F, float * DF, float dx, float dy ){
    float invdx = 1.0/dx;
    float invdy = 1.0/dy;
    for(int iy=1;iy<ny-1;iy++){ for(int ix=1;ix<nx-1;ix++){
        int i  = nx*iy + ix;
        int i2 = i<<1;
        DF[i2  ] = (F[i+1]  - F[i-1 ])*invdx;
        DF[i2+1] = (F[i+nx] - F[i-nx])*invdy;
    }}
}

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


void genPointsHash( int nPoints, Vec2f* points, Vec2f pmin, Vec2f pmax, int seed ){
    Vec2f span = (pmax-pmin)*(1.0f/256.0f);
    for(int i=0; i<nPoints; i++){
        int h = rand_hash2( i+seed );
        points[i].set( pmin.x + ((h   )&0xFF)*span.x,
                       pmin.y + ((h>>8)&0xFF)*span.y );
    }
}


void lerp(int nx, int ny, float * data, int nPoints, Vec2f* points, float * vals ){
    for(int i=0; i<nPoints; i++){
        Vec2f p = points[i];
        int ix  = (int)p.x; float dx = p.x - ix; float mx = 1-dx;
        int iy  = (int)p.y; float dy = p.y - iy;
        int ip  = nx*iy + ix;
        float val;
        val  = (1-dy)*( mx*data[ip] + dx*data[ip+1] );   ip+=nx;
        val +=    dy *( mx*data[ip] + dx*data[ip+1] );
        vals[i] = val;
    }
}

void lerpD(int nx, int ny, Vec2f* data, int nPoints, Vec2f* points, Vec2f * Dvals ){
    for(int i=0; i<nPoints; i++){
        Vec2f p = points[i];
        int ix  = (int)p.x; float dx = p.x - ix; float mx = 1-dx;
        int iy  = (int)p.y; float dy = p.y - iy; float my = 1-dy;
        int ip  = nx*iy + ix;
        Vec2f D;
        D.set_lincomb( my*mx, data[ip],  my*dx, data[ip+1] );   ip+=nx;
        D.add_lincomb( dy*mx, data[ip],  dy*dx, data[ip+1] );
        Dvals[i] = D;
        //printf( " %i (%g,%g) (%g,%g)  (%g,%g)\n", i, p.x, p.y, D.x, D.y,  data[ip].x, data[ip].y );
        //printf( " %i (%g,%g) (%g,%g)\n", i, p.x, p.y, Dvals[i].x, Dvals[i].y );
    }
    //printf(" ==== lerpD DONE \n");
    //exit(0);
}

void relaxPoints(int nx, int ny, Vec2f* data, int nPoints, Vec2f* points, Vec2f * Dvals, int niter, float dt, float damp ){
    for(int i=0; i<nPoints; i++){
        Vec2f p = points[i];
        Vec2f v = {0.0f,0.0f};
        for(int itr=0; itr<niter; itr++){
            int ix  = (int)p.x; float dx = p.x - ix; float mx = 1-dx;
            int iy  = (int)p.y; float dy = p.y - iy; float my = 1-dy;
            int ip  = nx*iy + ix;
            Vec2f f;
            f.set_lincomb( my*mx, data[ip],  my*dx, data[ip+1] );   ip+=nx;
            f.add_lincomb( dy*mx, data[ip],  dy*dx, data[ip+1] );

            v.mul    (damp);
            v.add_mul(f,-dt);
            p.add_mul(v,dt);

        }
        //printf("%i (%g,%g) (%g,%g) %i\n", i, points[i].x, points[i].y, p.x, p.y, niter );
        points[i] = p;
    }
    //printf(" ==== lerpD DONE \n");
    //exit(0);
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

