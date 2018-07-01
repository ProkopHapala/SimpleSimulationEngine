
#ifndef Radiosity_h
#define Radiosity_h

#include "Vec3.h"
#include "Mat3.h"
#include "geom3D.h"
#include "raytrace.h"

class SurfElement{ public:
    Vec3d pos;    //= Vec3dZero;
    Vec3d normal; //= Vec3dZ;
    double area;  //= 0.0;
    int isurf;//   = -1;   // helps with raytracing

    static inline double geomCoupling( const Vec3d& n, const Vec3d& p, const Vec3d& n0, const Vec3d& p0 ){
        Vec3d d = p - p0;
        double r = d.normalize();
        double c  = d.dot(n );
        double c0 = d.dot(n0);
        // Area is multiplied posteriori
        return c0*c/(r*r);
        //return c0*c/(r*r + area + area0); // correction for small size
    }

    inline double coupling( const SurfElement& b ){
        return geomCoupling( normal, pos, b.pos, b.normal );
    }

    SurfElement(){ isurf=-1; };
    SurfElement(Vec3d p, Vec3d nr, double area_,int isurf_):pos(p),normal(nr),area(area_),isurf(isurf_){};
};

class Radiosity{ public:
    //int nObstacles=0;
    //Triangle3D* obstacles=0;

    std::vector<Triangle3D> triangleObstacles;

    std::vector<SurfElement> elements;
    double* M=0;

    inline int trinagleToElements( int n, Triangle3D tri, int isurf ){
        const double off = 1.0/3;
        int ntot = n*n;
        //sprintf( " %i %i \n", n.a, n.b );
        Vec3d da = (tri.a - tri.c)*(1.0/n);
        Vec3d db = (tri.b - tri.c)*(1.0/n);
        Vec3d nr;
        double area = tri.normalArea(nr)/ntot;
        int i=0;
        for(int ia=0; ia<n; ia++){
            for(int ib=0; ib<(n-ia); ib++){
                Vec3d p = tri.c + da*(ia+off) + db*(ib+off);
                elements.push_back( SurfElement(p,nr,area,isurf) );
                i++;
                /*
                p.add( da*off + db*off );
                elements.push_back( SurfElement(p,nr,area,isurf) );
                i++;
                */
                //printf( "%i isurf %i \n", elements.size(), isurf );
            }
        }
    }

    inline void addTriangle( Triangle3D tri, double elemMaxSize ){
        // ToDo : automatic step bases on size (elemMaxSize)
        trinagleToElements( 5, tri, triangleObstacles.size() );
        triangleObstacles.push_back( tri );
    };

    inline double getOcclusion( const Vec3d& ray0, const Vec3d& hRay, double tmax, int ip1, int ip2 ){
        if( triangleObstacles.size()==0 ) return 0.0;
        // Check occlusion - TODO can be made better
        Vec3d hX,hY;
        hRay.getSomeOrtho(hX,hY);

        //glColor3f(0.0,0.0,1.0); Draw3D::drawVecInPos( hRay*100.0, ray0 );
        //glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( hX*100.0, ray0 );
        //glColor3f(0.0,1.0,0.0); Draw3D::drawVecInPos( hY*100.0, ray0 );

        //print( ray0 ); print( hRay ); print( hX ); print( hY ); printf("\n");


        //printf( "ip1 %i ip2 %i \n", ip1, ip2 );
        for( int i=0; i<triangleObstacles.size(); i++ ){
            if( (i==ip1) || (i==ip2) ) continue;
            //printf( "i %i \n", i, i );
            if( triangleObstacles[i].rayIn(ray0,hX,hY) ){ return 1.0; };
        }
        return 0.0;
    }

    inline int makeCouplingMatrix(){
        int n = elements.size();
        _realloc(M,n*n);
        for(int i=0; i<n; i++) M[i*n+i] = 0.0;
        for(int i=0; i<n; i++){
            SurfElement& eli = elements[i];
            for(int j=0; j<i; j++){
                SurfElement& elj = elements[j];
                if( eli.isurf == elj.isurf ){ M[i*n+j]=0.0; continue; }

                Vec3d  d = elj.pos - eli.pos;
                double r = d.normalize();
                double coupling = d.dot(eli.normal)*d.dot(elj.normal) / ( r*r + eli.area + elj.area );

                double occlusion = getOcclusion( eli.pos, d, r, eli.isurf, elj.isurf );

                coupling*=(1-occlusion);

                //glColor3f( occlusion, 0.0, 1-occlusion );
                //Draw::colorScale()
                //float c = sqrt(coupling)*0.3;
                float c = coupling*0.3;
                glColor3f( c, c, c );
                if(occlusion>0.5) glColor3f( 1.0, 0.0, 0.0 );
                Draw3D::drawLine( eli.pos, elj.pos );


                M[i*n+j] = coupling;
                M[j*n+i] = coupling;
            }
        }
    }

    inline void processTriangles( int ntri, Triangle3D* tris, double sz ){
        //nObstacles=ntri;
        //obstacles = tris;

        for(int i=0; i<ntri; i++){
            int nsamp = 5; // TODO :  automatially determine element size
            //trinagleToElements( {nsamp,nsamp}, tris[i], i );
            //trinagleToElements( nsamp, tris[i], i );
        }
        makeCouplingMatrix();
    }

};

#endif
