
#ifndef  CityGeneration_h
#define  CityGeneration_h

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

#include "PushArray.h"


/*

TODO:
 - clear duplicated indexes

*/


// ConvexTreeMap ?
//class Quad{ public:
//    Vec2d* p00,p01,p10,p01;
//}
//inline void getButEdge( int i, int n, double* buff, def1, def2 ){}

template<typename T>
inline T getButEdge( int i, int n, T* ys, T y0, T yn ){
    if(i<0){ return y0; }else if(i>=n){ return yn; } else return ys[i];
}

inline int safe_mod(int i, int base){ if(base>0){ return i%base; }else{ return 0;} }



class QuadNode{ public:
    Vec2i n;
    //Vec2d *p00=NULL,*p01=NULL,*p10=NULL,*p11=NULL;
    int i00,i01,i10,i11;
    //Quat4i iqs;

    int isp0;
    double * as     = NULL;   // a split lines
    double * bs     = NULL;   // b split lines
    QuadNode * subs = NULL; // tiles

    inline void setCorners( int i00_, int i01_, int i10_, int i11_){ i00=i00_; i01=i01_; i10=i10_; i11=i11_; };

    QuadNode(){ setCorners(0,1,2,3); };
    //QuadNode(Vec2d *p00_,Vec2d *p01_,Vec2d *p10_,Vec2d *p11_){ p00=p00_; p01=p01_; p10=p10_; p11=p11_; };
    QuadNode( int i00_,int i01_,int i10_,int i11_){ setCorners(i00_,i01_,i10_,i11_); };

    inline int getQuadI ( int ia, int ib ){ return n.a*ib     + ia; }
    inline int getPointI( int ia, int ib ){ return (n.a+1)*ib + ia; }

    void split( const Vec2i& n_, const Vec2d& jitter ){
        n = n_;
        //printf( " (%i,%i) (%f,%f) \n", n.a, n.b, jitter.a, jitter.b );
        if(as) delete [] as; as = new double[n.a-1];
        if(bs) delete [] bs; bs = new double[n.b-1];
        double step;
        step=1.0d/n.a; for(int i=1; i<n.a; i++){ as[i-1] = ( i + randf(-jitter.a,jitter.a) )*step; };
        step=1.0d/n.b; for(int i=1; i<n.b; i++){ bs[i-1] = ( i + randf(-jitter.b,jitter.b) )*step; };
    }

    int makeSubPoints( std::vector<Vec2d>& ps ){
        Vec2d* ps0 = &ps[0];
        Vec2d c00=*(ps0+i00),c01=*(ps0+i01),c10=*(ps0+i10),c11=*(ps0+i11);
        isp0 = ps.size();
        for(int ib=0; ib<=n.b; ib++){
            double b,mb; b=getButEdge(ib-1,n.b-1,bs,0.0,1.0); mb=1-b;
            //printf( "ia %i %f %f \n", ia, a, ma  );
            Vec2d vi0 = c00*mb + c01*b;
            Vec2d vi1 = c10*mb + c11*b;
            for(int ia=0; ia<=n.a; ia++){
                double a,ma; a=getButEdge(ia-1,n.a-1,as,0.0,1.0); ma=1-a;
                Vec2d v = vi0*ma + vi1*a;
                //printf( "(%i,%i) (%f,%f) (%f,%f) \n", ia,ib,  a,b, v.x,v.y  );
                //if((ia>0)&&(ib>0))printf( "(%i,%i) %i \n", ia,ib, ps.size() );
                ps.push_back( v );
            }
        }
    }

    int makeSubPoints2( Vec2i n_, Vec2d jitter, std::vector<Vec2d>& ps ){
        n = n_;
        Vec2d* ps0 = &ps[0];
        Vec2d c00=*(ps0+i00),c01=*(ps0+i01),c10=*(ps0+i10),c11=*(ps0+i11);
        isp0 = ps.size();
        double step_a = 1.0d/n.a;
        double step_b = 1.0d/n.b;
        for(int ib=0; ib<=n.b; ib++){
            for(int ia=0; ia<=n.a; ia++){
                double a,b,ma,mb;
                a = randf(-jitter.a,jitter.a);
                b = randf(-jitter.b,jitter.b);
                if( (ia==0) || (ia==(n.a)) ) a=0;
                if( (ib==0) || (ib==(n.b)) ) b=0;
                a = ( ia + a )*step_a; ma=1-a;
                b = ( ib + b )*step_b; mb=1-b;
                Vec2d v = (c00*ma + c01*a)*mb  +  (c10*ma + c11*a)*b;
                ps.push_back( v );
            }
        }
    }

    int makeSubQuads( ){
        //printf("====\n");
        if(subs) delete [] subs; subs = new QuadNode[ n.totprod() ];
        //for(int ia=1; ia<=n.a; ia++){
        //    for(int ib=1; ib<=n.b; ib++){
        for(int ib=0; ib<n.b; ib++){
            for(int ia=0; ia<n.a; ia++){
                //int ij = getIndex(ia-1,ib-1);
                int iq   =      getQuadI (ia,ib);
                int ii00 = isp0+getPointI(ia,ib);
                int ii01 = ii00+1;
                int ii10 = ii00+n.a+1;
                int ii11 = ii10+1;
                //printf( "(%i,%i) %i \n", ia+1,ib+1, ii11 );
                //QuadNode* nd = new QuadNode( ii00, ii01, ii10, ii11 );
                subs[iq].setCorners(ii00, ii01, ii10, ii11);
            }
        }
    }

    void inset( Vec2d t00, Vec2d t11, std::vector<Vec2d>& ps ){
        //printf( "(%f,%f) (%f,%f) \n", t00.a,t00.b,  t11.a,t11.b );
        //Vec2d c00=*p00,c01=*p01,c10=*p10,c11=*p11;
        Vec2d* ps0 = &ps[0];
        Vec2d c00=*(ps0+i00),c01=*(ps0+i01),c10=*(ps0+i10),c11=*(ps0+i11);
        Vec2d c00_,c01_,c10_,c11_;
        Vec2d m00 = (Vec2d){1.0-t00.a,1.0-t00.b};
        Vec2d m11 = (Vec2d){1.0-t11.a,1.0-t11.b};
        c00_= (c00*m00.a + c01*t00.a)*m00.b  +  (c10*m00.a + c11*t00.a)*t00.b;
        c01_= (c00*m11.a + c01*t11.a)*m00.b  +  (c10*m11.a + c11*t11.a)*t00.b;
        c10_= (c00*m00.a + c01*t00.a)*m11.b  +  (c10*m00.a + c11*t00.a)*t11.b;
        c11_= (c00*m11.a + c01*t11.a)*m11.b  +  (c10*m11.a + c11*t11.a)*t11.b;
        //if( ps != NULL ){
        ps.push_back(c00_); i00=ps.size()-1;
        ps.push_back(c01_); i01=ps.size()-1;
        ps.push_back(c10_); i10=ps.size()-1;
        ps.push_back(c11_); i11=ps.size()-1;
        //}else{
        //    *p00=p00_; *p01=p01_; *p10=p10_; *p11=p11_;
        //}
    }

    void splitRecursive( int level, Vec2i nmax, Vec2i nmin, Vec2d jitter, double maxInset, std::vector<Vec2d>& ps ){
        //printf("splitRecursive level %i  %f\n", level, maxInset );
        int na = nmin.a+safe_mod( rand(), (nmax.a-nmin.a) );
        int nb = nmin.b+safe_mod( rand(), (nmax.b-nmin.b) );
        //printf( " (%i,%i) (%f,%f) \n", n.a, n.b, jitter.a, jitter.b );
        if(maxInset>1e-8){ inset( {randf(0,maxInset),randf(0,maxInset)}, {randf(1.0-maxInset,1.0),randf(1.0-maxInset,1.0)}, ps ); }

        //split( {na,nb}, jitter );
        //makeSubPoints( ps );
        makeSubPoints2( {na,nb}, jitter, ps );

        makeSubQuads();
        level--;
        if(level<=0) return;
        int ntot=n.totprod();
        for(int i=0; i<ntot; i++){
            subs[i].splitRecursive( level, nmax, nmin, jitter, maxInset, ps );
        };
    }

};




// Branch Roads by T-junctions starting from spline
//  - outsprings should not intersect previous roads

/*
class RoadBrancher{
    Vec2d ;

};
*/





//class SplineQuad(){};
//int splitSplineQuad(  ){}



#endif

