#ifndef TriangleRayTracer_h
#define TriangleRayTracer_h

#include "Vec3.h"
#include "Mat3.h"
#include "geom3D.h"
#include "raytrace.h"
#include "CMesh.h"
#include "Lingebra.h"

#include "Draw3D.h" // DEBUG

class SurfElement{ public:
    Vec3d pos;    //= Vec3dZero;
    Vec3d normal; //= Vec3dZ;
    double area;  //= 0.0;
    int isurf;//   = -1;   // helps with raytracing

    static inline double geomCoupling( const Vec3d& n, const Vec3d& p, const Vec3d& n0, const Vec3d& p0, double area0 ){
        //Vec3d  d = elj.pos - eli.pos;
        //double r2 = d.normalize();
        //double coupling = d.dot(eli.normal)*d.dot(elj.normal)/r2;
        //if( fabs(coupling) < couplingTrashold ){ M[i*n+j]=0.0; continue; };
        //coupling /= ( r2 + eli.area + elj.area );

        Vec3d d   = p - p0;
        double r2 = d.norm2();
        double c  = d.dot(n ); // not normalized
        double c0 = d.dot(n0); // not normalized
        // Area is multiplied posteriori
        return c0*c/(r2*(r2+area0)); // first r2 is normalization of c0*c; second is distance
        //return c0*c/(r*r + area + area0); // correction for small size
    }

    inline double geomCoupling( const SurfElement& b ){
        return geomCoupling( normal, pos, b.pos, b.normal, b.area + area );
    }

    SurfElement(){ isurf=-1; };
    SurfElement(Vec3d p, Vec3d nr, double area_,int isurf_):pos(p),normal(nr),area(area_),isurf(isurf_){};
};

class TriangleRayTracer{ public:

    std::vector<Triangle3D>  triangleObstacles;
    std::vector<SurfElement> elements;

    void clearTriangles(){
        triangleObstacles.clear();
        elements.clear();
    }

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

    inline double trapezElements(
        const Vec3d& nr, const Vec3d& p0, const Vec3d& a, const Vec3d& b, int isurf,
        double invStep, double xmax, double k, double y0
    ){
        double areaSum = 0;
        double na = (int)(xmax*invStep+0.5);  if(na==0)na=1;
        double dx = xmax/na;
        double kdx = k*dx;
        double x = dx*0.5;
        //double yo = y0;
        double ymax = y0  + k*x;
        //printf( "trapezElements xmax %g k %g y0 %g dx %g \n", xmax, k, y0, dx );
        for(int ia=0;ia<na; ia++){
            //double ymax = k*x + y0;
            int nb      = (int)(ymax*invStep+0.5); if(nb==0)nb=1;
            //double dy   = 0.5*(ymax+yo)/nb;
            double dy   = ymax/nb;
            double area = dy*dx;
            Vec3d pa = p0 + a*x;
            for(int ib=0;ib<nb; ib++){
                elements.push_back( SurfElement(pa+b*(dy*(ib+0.5)),nr,area,isurf) );
                areaSum+=area;
            }
            //yo = ymax;
            ymax += kdx;
            x+=dx;
        }
        return areaSum;
    }

    inline double trinagleToElements2( double maxSize, Triangle3D tri, int isurf ){
        Vec3d  ds []{tri.c-tri.b, tri.a-tri.c, tri.b-tri.a };
        double r2s[]{ds[0].norm2(),ds[1].norm2(),ds[2].norm2()};

        //printf( "r2s %g %g %g \n", r2s[0], r2s[1], r2s[2] );

        int i0,i1,i2;
        if(r2s[0]>r2s[1]){ if(r2s[0]>r2s[2]){i0=0;i1=1;i2=2;}else{i0=2;i1=0;i2=1;} }
        else             { if(r2s[1]>r2s[2]){i0=1;i1=2;i2=0;}else{i0=2;i1=0;i2=1;} }
        double r  =sqrt(r2s[i0]);
        double ir =1/r;
        Vec3d  h  =ds[i0]*ir;

        double la  = h.dot( ds[i1] );
        Vec3d  b   = ds[i1] - h*la;
        double lb  = b.normalize();
        if(la<0){ la=-la; h.mul(-1.); int i=i1;i1=i2;i2=i; }
        double la2 = r-la;

        //printf( ">r2s %g %g %g | %g %g %g \n", r2s[i0], r2s[i1], r2s[i2], la, la2, r );

        double invStep = 1/maxSize;

        Vec3d nr; nr.set_cross( h, b ); nr.normalize();

        /*
        glColor3f(1.,0.,0.); Draw3D::drawPointCross( tri.array[i1], maxSize );
        //glColor3f(0.,1.,0.); Draw3D::drawTriangle( tri.array[i2], false );
        glColor3f(0.,0.,1.); Draw3D::drawPointCross( tri.array[i2], maxSize );
        glColor3f(1.,1.,1.); Draw3D::drawTriangle( tri, false );
        glColor3f(1.,0.,0.); Draw3D::drawVecInPos( b*lb, tri.array[i1]      );
        glColor3f(0.,0.,1.); Draw3D::drawVecInPos( b*lb, tri.array[i2]      );
        glColor3f(.5,.5,.0); Draw3D::drawVecInPos( b*lb, tri.array[i1]+h*la );
        //glColor3f(.0,.5,.5); Draw3D::drawVecInPos( b*lb, tri.array[i2]+h*(-la2) );
        */


        //int na1 = (int)(     la *invStep + 0.5 ); if(na1==0)na1=1;
        //int na2 = (int)( ((r-la)*invStep + 0.5 ); if(na2==0)na2=1;
        //int nb  = (int)(     lb *invStep + 0.5 ); if(nb ==0)nb =1;
        return trapezElements( nr,tri.array[i1],h    ,b,isurf, invStep,  la,   lb/la , 0 )
             + trapezElements( nr,tri.array[i2],h*-1.,b,isurf, invStep,  la2,  lb/la2, 0 );
    }

    inline void addTriangle( Triangle3D tri, double elemMaxSize, bool active ){
        // ToDo : automatic step bases on size (elemMaxSize)
        //if(active)trinagleToElements( 5, tri, triangleObstacles.size() );
        if(active)trinagleToElements2( elemMaxSize, tri, triangleObstacles.size() );
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
