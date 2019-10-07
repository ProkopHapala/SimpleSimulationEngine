#ifndef SphereGaussSeidel_h
#define SphereGaussSeidel_h

#include "Vec3.h"
#include "fastmath.h"
#include "sweep.h"
//#include "kBoxes.h"
// ToDo : see /home/prokop/git/SimpleSimulationEngine/cpp/sketches_SDL/Shooter/test_BoxAndSweep.cpp









/*


To Do:

Icosahedral:
/home/prokop/git/SimpleSimulationEngine/cpp/common/maps/SphereSampling.h







Volume fraction : solid / circumscribed sphere

## octahedron: V = sqrt(2)/3 = 0.47140452079
circumscribed sphere:  4/3*pi*sqrt(0.5)^3 = 1.48096097939
ratio = 0.34906585039 ( 3.14159265359  )

## cube V = a^3 = 1
circumscribed sphere: 4/3*pi*sqrt(0.75)^3 = 2.72069904635
ratio = 0.36755259694 (2.72069904641)

## Rhombic Dodecahedron V = 2 (from unit cube)
a=1 V=3.07920144
circumscribed sphere: 4/3*pi * 1.154700538379251529^3 = 6.44906440617
ratio = 0.47746482927 ( 2.09439510239  ()

icosahedron volume = 2.18169  (a=1)
circumscribed sphere:  (r=0.95105651629515357211) V=3.60335944157
ratio 0.60545999792   ( 1.65163677771 )

dodecahedron  volume  7.66312 (a=1)
circumscribed sphere:  (r=1.401258538444073544676677) V = 11.5250661067
ratio 0.66490898438  (1.50396523958)

Truncated Octahedron
Volume (a=1) 11.3137085 a^3
circumscribed sphere: R=sqrt(5/2)   V = 16.5576471097
ratio 0.68329204174   1.46350306884

https://math.stackexchange.com/questions/1718062/which-archimedean-solid-takes-up-the-most-volume-in-its-circumscribed-sphere

https://en.wikipedia.org/wiki/Snub_cube

*/



//Icosahedral Face
//  /home/prokop/git/SimpleSimulationEngine/cpp/common/maps/SphereSampling.h

void getIcosaFace( Vec3d p, uint8_t& iface ){
//  http://www.kjmaclean.com/Geometry/Icosahedron.html
//  https://en.wikipedia.org/wiki/Regular_icosahedron
//  rc = sqrt(10+2*sqrt(5))/4 = sqrt(phi^2 + 1)/2 = 0.95105651629 * a
//  plate height = phi/(2*sqrt(phi^2+1)).a = 0.42532540417.a  = 0.4472135955 rc
//  top height   = 1/sqrt(phi^2+1).a       = 0.52573111211.a  = 0.5527864045 rc
    double phi10 = (atan2( p.z, p.x )+ M_PI) * 1.5915494309;  //  1.5915494309 = 10/(2*pi)
    int iphi     = (int)phi10;
    int ioff,i;
    Vec3d& d1=((Vec3d*)oct_edge_planes)[iphi];
    if( d1.dot(p)>0 ){ ioff=0; i=iphi/2; }else{ ioff=5; i=(iphi+1)/2;  if(i>=5) i=0; };
    int i2 = i+1; if(i2>=5) i2=0;
    iface=i+ioff;
    Vec3d& d2=((Vec3d*)oct_edge_planes)[iface+10];
    iface<<=1;
    if( d2.dot(p) > 0 ){ // uper half of rect
        iface++;
    }
}

static const double oct_edge_planes[] = {
// Equatorial edges [0..10]
 0.2628655560595669, 0.5257311121191337,-0.8090169943749476,
-0.6881909602355869, 0.5257311121191336, 0.5000000000000000,
 0.85065080835204,   0.5257311121191337, 0.0,
-0.6881909602355869, 0.5257311121191336,-0.5000000000000000,
 0.262865556059567,  0.5257311121191337, 0.8090169943749475,
 0.2628655560595667, 0.5257311121191337,-0.8090169943749476,
-0.6881909602355868, 0.5257311121191336, 0.5000000000000000,
 0.85065080835204,   0.5257311121191336, 0.0,
-0.688190960235587,  0.5257311121191337,-0.5000000000000000,
 0.2628655560595671, 0.5257311121191337, 0.8090169943749475,
// Upper Cap edges
 0.42532540417602,   0.85065080835204,  0.3090169943749475,
-0.1624598481164531, 0.85065080835204,  0.5,
-0.5257311121191336, 0.85065080835204,  0.0,
-0.1624598481164533, 0.85065080835204, -0.5,
 0.42532540417602,   0.85065080835204, -0.3090169943749476,
 // Bottom Cap edges
-0.5257311121191336, 0.85065080835204,  0.0,
-0.1624598481164533, 0.85065080835204, -0.5,
 0.4253254041760200, 0.85065080835204, -0.3090169943749476,
 0.4253254041760201, 0.85065080835204,  0.3090169943749473,
-0.1624598481164531, 0.85065080835204,  0.5,

};


// Octahedral Face

int assignOctahedronFace( Vec3d p ){
    return ( ((uint8_t)(p.x>0))<<2) | ( ((uint8_t)(p.y>0))<<1) | ( ((uint8_t)(p.z>0) ));
}
//  Truncating plane is in 2/3 of the face => coordinates ( (2/3)*v, (1/3)*sqrt(a) ) = normal





inline int assignCubicFace(Vec3d& d){
    int i0=0; if(d.x<0){ d.x=-d.x; i0+=3; }
    int i1=1; if(d.y<0){ d.y=-d.y; i1+=3; }
    int i2=2; if(d.z<0){ d.z=-d.z; i2+=3; }
    return (d.x>d.y)?((d.x>d.z)?i0:i2):((d.y>d.z)?i1:i2);
}

class SphereGaussSeidel{ public:

    constexpr static const int N_DIR   = 3;
    constexpr static const int N_NEIGH = N_DIR*2;

    int n;

    Vec3d  * pos  =0;
    double * Rs   =0;

    //Vec3d  * dpos =0; // estimate of error in position

    int      nneighs   =0;
    int    * neighs    =0;
    double * neighDist =0;

    inline int findBackNeigh(int i, int j){
        int ioff = i * N_NEIGH;
        int iend = ioff+N_NEIGH;
        for(int i=ioff;i<iend;i++){
            if( j==neighs[i] ) return i;
        }
        return -1;
    }

    void realloc_all(int n_){
        n=n_;
        _realloc( pos, n );
        _realloc( Rs,  n );
        realloc_internal();
    }

    void realloc_internal(){
        // -- these are external binds
        //Vec3d  * pos  =0;
        //double * Rs   =0;
        nneighs = n*N_NEIGH;
        _realloc( neighs,    nneighs );
        _realloc( neighDist, nneighs );
    }

    double eval_errors(){

        // brute force - improve this in future by some broad-space collision solver
        for(int i=0; i<nneighs; i++){
            neighs   [i]=-1;
            neighDist[i]=1e+300;
        }
        double errMax = 0;
        for(int i=0; i<n; i++){
            int ioff = i*N_NEIGH;
            double * ngl = neighDist + ioff;
            int    * ngi = neighs    + ioff;
            const Vec3d& pi = pos[i];
            double Ri       = Rs [i];
            for(int j=0; j<n; j++){
                if(j==i) continue;
                Vec3d d = pos[j]-pi;
                //double l2 = d.norm();
                int k    = assignCubicFace(d);
                //int kc   = i%3;
                //double dl = fabs(d.arr[kc]);
                double dk = d.array[k%3] - Ri - Rs[j];
                //printf( "(%i,%i):  k %i  %g  %g \n", i,j, k, d.array[k%3], dk );
                if(dk<errMax)errMax=dk;
                if( dk < ngl[k] ){
                    ngi[k]=j;
                    ngl[k]=dk;
                }
            }
        }
        return -errMax;
    }

    double update_pos( double treshhold, double bmix, double dmax ){
        //double tresh2 = sq(treshhold);
        const int NK = N_NEIGH/2;
        int nmove=0;
        double maxMove=0;
        for(int i=0; i<n; i++){
            //double dr = dpos.norm2();
            int ioff = i*N_NEIGH;
            double * ngl = neighDist + ioff;
            int    * ngi = neighs    + ioff;

            Vec3d dp;
            bool bmove = false;
            for(int k=0; k<N_DIR; k++){
                double l = ngl[k  ];
                double r = ngl[k+3];
                //printf( "%i %i %g %g  |  %g \n", i, k, r, l, treshhold  );
                if( (-l>treshhold)||(-r>treshhold) ) bmove=true;
                //double d = l-r;
                bool   br=r<-treshhold;
                bool   bl=l<-treshhold;
                double d=0;

                //if(i==49)printf("i==49 k==%i r,l %g,%g  pos[k] %g treshhold %g \n", k, r, l, pos[i].array[k], treshhold  );
                if(br){
                    //printf( "%i,%i  r %g \n", i,k, r );
                    if(bl){ // both r,l unsatisfied
                        d = (r-l)*bmix;
                        //printf( "%i,%i  r,l(%g,%g) :d %g \n", i,k, r,l, d );
                    }else{ // r unsitisfied - go negative
                        //d = -r*bmix; // d is positive
                        d = (-r + fmin(-r,l) )*0.5*bmix; // d is positive
                        //printf( " ->  d %g -r %g -r*k %g ", d, -r, -r*k );
                        //_setmin( d, l );
                        //printf( " ->  d %g l %g ", d, l );
                        //printf( "%i,%i  r %g  :d %g  | l  %g \n", i,k, r, d,     l );
                    }
                }else if(bl){ // l unsitisfied - go positive
                    //d = l*bmix;   // d is negative, r is positive
                    d = (l + fmax(l,-r) )*0.5*bmix; // d is positive
                    //_setmax( d, -r );
                    //printf( "%i,%i  l %g :d %g |  -r %g \n", i,k, l, d,    -r );
                }else{
                    d = -pos[i].array[k];
                    //double d_DEBUG=d;
                    if(d>0){ _setmin( d, l ); }else{ _setmax( d, -r ); };
                    //if( fabs(d)<treshhold )printf( "d %g d_DEBUG %g  l %g -r %g \n", d, d_DEBUG, l, -r );
                    d*=bmix;
                    if( fabs(d)>treshhold ) bmove=true;
                }
                //if(fabs(d)>treshhold) bmove=true;
                maxMove=fmax( maxMove, fabs(d) );
                dp.array[k] = clamp(d,-dmax,dmax);
            }
            //printf( "\n>> %i dp(%g,%g,%g)    l(%g,%g,%g)r(%g,%g,%g) \n", i, dp.x,dp.y,dp.z,   ngl[0],ngl[1],ngl[2],ngl[3],ngl[4],ngl[5]  );
            if(bmove){
                pos[i].add( dp );
                //printf( "move %i dp(%g,%g,%g) \n", i, dp.x,dp.y,dp.z );
                /*
                for(int k=0; k<N_NEIGH; k++){
                    ngi[k]=-1;
                    //ngl[k]=1+300;
                    int j = findBackNeigh(ngi[k],i);
                    if(j>=0){
                        neighs[j]=-1;
                        //neighDist[k]=1e+300;
                    }
                }
                */
                nmove++;
            }
        }
        printf( " nmove %i \n", nmove );
        //return nmove;
        return maxMove;
    }



};

#endif
