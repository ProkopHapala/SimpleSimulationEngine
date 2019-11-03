#ifndef TriangleRayTracer_h
#define TriangleRayTracer_h

#include "Vec3.h"
#include "Mat3.h"
#include "geom3D.h"
#include "raytrace.h"
#include "CMesh.h"
#include "Lingebra.h"

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

class TriangleRayTracer{ public:

    std::vector<Triangle3D>  triangleObstacles;
    std::vector<SurfElement> elements;

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
            }
        }
    }

    inline void addTriangle( Triangle3D tri, double elemMaxSize, bool active ){
        // ToDo : automatic step bases on size (elemMaxSize)
        if(active)trinagleToElements( 5, tri, triangleObstacles.size() );
        triangleObstacles.push_back( tri );
    };

    //inline fromMesh( CMesh& mesh, double elemMaxSize, bool justObstacle ){
    inline void fromMesh( int ntri, Vec3i* tris, Vec3d* verts, double elemMaxSize, bool active ){
        for(int i=0; i<ntri; i++){
            Triangle3D tri;
            tri.a = verts[tris[i].a];
            tri.b = verts[tris[i].b];
            tri.c = verts[tris[i].c];
            addTriangle( tri, elemMaxSize, active );
        }
    }

     inline void fromMesh( CMesh& mesh, double elemMaxSize, bool active ){
        fromMesh( mesh.ntri, mesh.tris, mesh.verts, elemMaxSize, active );
    }

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

};

#endif
