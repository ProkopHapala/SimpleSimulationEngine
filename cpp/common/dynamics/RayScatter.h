
#ifndef RayScatter_h
#define RayScatter_h

#include "Vec3.h"
#include "geom3D.h"
#include "raytrace.h"
#include "CMesh.h"

//#include "VecNm.h"

/*

ToDo:

For Optimal implementation of GPU / OpenCL we should use univorm datastructure => Triangle-surfaces
 * Each Triangle is a boundary between two material - Right-side / Left-side. Hiting the triangle means crossing this boundary
 * Triangles can be stored in acceleration Grid or Tree


*/


/*
enum class TRayObj{ sphere, cylinder, Tetrahedron, box, mesh };

struct RayedObject{
    int   kind;
    void* obj;
    int   imaterial;

    double ray( const Vec3d &ray0, const Vec3d &hRay ){
        // vtable vs switch :  https://stackoverflow.com/questions/4467485/vftable-performance-penalty-vs-switch-statement
        int   imin;
        Vec3d normal;
        switch(kind){
            case TRayObj::sphere :
                return ((Sphere3d*)obj)->ray(ray0,hRay);
                //inline double raySphere( const Vec3d& ray0, const Vec3d& hRay, double R, const Vec3d& center );
                break;
            case TRayObj::mesh :
                return ((Mesh*)obj)->ray( ray0, hRay, normal, imin );
                break;
            default:
        }
    }

}
*/

constexpr static const int nScatAng=16;

// ===================================
//       ScatterMatrial
// ===================================

struct ScatterMatrial{
    //double scatterProb;     // [1/m  ]  scattering probability
    double mfp;             // [m]  mean free path
    double Sang[nScatAng];  // [1/rad]  angular scattering distribution function ( cumulative distribution, i.e. integral of density )
    double Dang[nScatAng];  // [1/rad]  angular scattering density      function ( i.e. derivative of cummulative distribution )

    inline double getScatterCos()const{
        double rnd = randf();
        double oSi=0;
        double t;
        for(int i=0; i<nScatAng; i++){ // ToDo : This can be improved by bisection-search
            double Si = Sang[i];
            //printf( "S[%i] %g <? %g \n", i, Si, rnd );
            if(rnd<Si){
                t= i + (1-(Si-rnd)/(Si-oSi));
                break;
            }
            oSi=Si;
        }
        double alpha = M_PI*(t/nScatAng);
        //printf( "getScatterCos t %g alpha %g cos(a) %g  rnd %g \n", t, alpha, cos( alpha ), rnd );
        return cos( alpha ); // ToDo : in future we can avoid using cos() by clever choice of intervals
    }

    inline double interpolateDens(double x)const{
        x/=nScatAng;
        int    ix = (int)x;
        double dx = x-ix;
        double s0 = 0; if(ix>0)s0=Dang[ix-1];
        return s0*(1-dx)+Dang[ix]*dx;
    }

    void setLorenz(double width){
        double sum=0;
        double dx = 1./nScatAng;
        for(int i=0; i<nScatAng; i++){
            double val = 1/(1+sq(i*dx/width));
            //printf( "val[%i] %g \n", i , val );
            sum+=val;
            Dang[i]=val;
            Sang[i]=sum;
        }
        double renorm=1/sum;
        for(int i=0; i<nScatAng; i++){
            Sang[i]*=renorm;
            Dang[i]*=renorm*nScatAng;
        };
        double oSi=0;
        for(int i=0; i<nScatAng; i++){  printf("S[%i] %g %g \n", i, Sang[i], Sang[i]-oSi ); oSi=Sang[i]; };
    }

};

// ===================================
//       RayPath
// ===================================

struct RayPath{
    int n;
    int*   mats;  // ToDo : problem - there can be multiple materials along path between two points
    Vec3d* poss;
    // --- aux - for future acceleration
    //Vec3d* dir;
    //double* l;

    void eval( ScatterMatrial* materials, double pointSize ){
        double p = 1;
        Vec3d op=poss[0];
        Vec3d oh;
        for(int i=1; i<n; i++){
            Vec3d pi=poss[i];
            Vec3d d; d.set_sub( pi, op );
            double l = d.normalize();
            // Field of view
            double p_dist = 1 / ( 1 + sq(l/pointSize) );
            p*=p_dist;
            // Transparency Attenuation
            // https://en.wikipedia.org/wiki/Attenuation_length
            // https://en.wikipedia.org/wiki/Beer-Lambert_law
            // ToDo : problem - there can be multiple materials along path between two points
            int imat = mats[i];
            const ScatterMatrial& mat = materials[imat];
            double p_ate = exp( -l/mat.mfp ); // this can be done more efficiently
            p*=p_ate;
            // Scattering
            if(i>2){ // Kink - angular
                double c = d.dot( oh ); // cosine of scattering angle
                double p_ang = mat.interpolateDens( c );
                p*=p_ang;
            }
            op=pi;
            oh=d;
        }
    }

    void movePoint( int i, Vec2d pos ){
    }

};

// ===================================
//       RayPath
// ===================================

class RayScatterWorld{ public:
    std::vector<Vec3d>           points;
    //std::vector<RayTri>          triangles;
    std::vector<Vec3i>           triangles;
    std::vector<Vec2i>           tri2mat;
    std::vector<ScatterMatrial>  materials;
    //std::vector<RayedObject>     objs;
    //Rayt3d

    int nScatMax=64;
    bool bDraw = false;

    void addMesh( CMesh& mesh, Vec2i mat){
        int ip0=points.size();
        points.resize( ip0 + mesh.nvert );
        for(int i=0; i<mesh.nvert; i++){
            points[ip0+i] = mesh.verts[i];
        }
        int it0 = triangles.size();
        triangles.resize( it0 + mesh.ntri );
        tri2mat  .resize( it0 + mesh.ntri );
        for(int i=0; i<mesh.ntri; i++){
            triangles[it0+i] = ( mesh.tris[i] + ip0 );
            tri2mat  [it0+i] =  mat ;
        }
    }

    double ray( Vec3d& ray0, Vec3d& hRay, int& imat)const{
        Vec3d normal;
        int imin;
        double t = rayTriangles( ray0, hRay, triangles.size(), &triangles[0], &points[0], normal, imin );
        if( imin >= 0 ){
            double cdot = normal.dot( hRay );
            if( cdot>0 ){ imat = tri2mat[imin].b; }
            else        { imat = tri2mat[imin].a; }
        }else{
            imat = -1;
        }
        return t;
    }

    double scatterInMat(int imat, double tmax, Vec3d& hRay )const{
        const ScatterMatrial& mat = materials[imat];
        double prob  = tmax/mat.mfp;  // scattering probability
        int nscat    = (int)(prob*4)+1; // number of sub-steps - higher number here makes it more precise
        if(nscat>nScatMax) nscat=nScatMax;
        double dt    = tmax/nscat;
        double dprob = prob/nscat;
        double t     = 0;
        //printf( "scatter nscat %i prob %g \n", nscat, prob );
        for(int i=0; i<nscat; i++){
            double rnd=randf();
            //printf( "scatter[%i] dprob %g rnd %g \n", i, dprob, rnd );
            if( rnd<dprob ){
                t+=(1-rnd/dprob)*dt;
                break;
            }
            t+=dt;
        }
        // generate scattering angle and vector  - ToDo : this is perhaps inefficient way
        double ca = mat.getScatterCos();
        double sa = sqrt(1-ca*ca);
        Vec3d  a,b; hRay.getSomeOrtho(a,b);
        double phi = randf()*2*M_PI;
        double cph = cos(phi);
        double sph = sin(phi);
        hRay.mul(ca);
        hRay.add_mul(a,sa*cph);
        hRay.add_mul(b,sa*sph);
        //printf( "ca %g \n", ca );
        return t;
    }

    void scatter( int n, Vec3d& ray0, Vec3d& hRay )const{
        if(bDraw){
            glBegin(GL_LINE_STRIP);
            Draw3D::vertex( ray0 );
        }
        double tmax;
        for(int i=0; i<n; i++){
            int imat;
            tmax = ray( ray0, hRay, imat);
            if(tmax>1e+100) break;
            //printf( "scatter[%i] imat %i \n", i, imat );
            Vec3d ohRay = hRay;
            double t;
            if( imat>=0 ){ t = scatterInMat( imat, tmax, hRay ); }
            else         { t = tmax;                             }
            ray0.add_mul(ohRay, t*(1+1e-8) );
            if(bDraw)Draw3D::vertex( ray0 );
        }
        if(bDraw){
            if(tmax>1e+100){
                Draw3D::vertex( ray0+hRay*2.0 );
            }
            glEnd();
        }
    }

};


#endif
