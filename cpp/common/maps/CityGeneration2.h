
#ifndef  CityGeneration2_h
#define  CityGeneration2_h

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

/*
class QuadNodeBin{ public:
    int kind=0;
    Vec2d s;
    QuadNode *a,b;

    void split(){

    }
};
*/

class Quad2d{ public:
    Vec2d p00,p01,p10,p11;
    inline void set( const Vec2d& p00_, const Vec2d& p01_, const Vec2d& p10_, const Vec2d& p11_ ){  p00=p00_; p01=p01_; p10=p10_; p11=p11_; };
};


template<typename Fcond, typename Fleaf>
class QuadSpliter{ public:

    double cmin=0.25,cmax=0.75;
    Fcond fcond;
    Fleaf fleaf;

    //QuadSpliter(Fcond fcond_, Fleaf fleaf_){ fcond=fcond_; fleaf=fleaf_; };
    QuadSpliter(Fcond fcond_, Fleaf fleaf_):fcond(fcond_),fleaf(fleaf_){};

    void splitBinRec( Quad2d& q, int level, bool side ){
        printf( "splitBinRec level: %i\n", level  );
        if( (level>0) && fcond(q) ){
            double c1 = randf(cmin,cmax); double m1=1-c1;
            double c2 = randf(cmin,cmax); double m2=1-c2;

            //if( rand()&1 ) side!=side;
            side=!side;

            Vec2d p1,p2;
            Quad2d q1,q2;
            if( side ){
                p1 = q.p00*m1 + q.p01*c1;
                p2 = q.p10*m2 + q.p11*c2;
                q1.set( q.p00, p1, q.p10, p2 );
                q2.set( p1, q.p01, p2, q.p11 );
            }else{
                p1 = q.p00*m1 + q.p10*c1;
                p2 = q.p01*m2 + q.p11*c2;
                q1.set( q.p00, q.p01, p1, p2 );
                q2.set( p1, p2,  q.p10, q.p11 );
            }
            level--;
            splitBinRec( q1, level, side );
            splitBinRec( q2, level, side );
        }else{
            printf( "to fleaf() \n" );
            fleaf(q);
        }
    }

};

class QuadNode2{ public:
    int n;
    int i00,i01,i10,i11;
    QuadNode2 * subs = NULL; // tiles

    inline void allocate( int n_ ){ subs; n=n_; if(subs) delete [] subs; subs = new QuadNode2[n_]; };
    inline void setCorners( int i00_, int i01_, int i10_, int i11_){ i00=i00_; i01=i01_; i10=i10_; i11=i11_; };

    QuadNode2(){ setCorners(0,1,2,3); };
    //QuadNode(Vec2d *p00_,Vec2d *p01_,Vec2d *p10_,Vec2d *p11_){ p00=p00_; p01=p01_; p10=p10_; p11=p11_; };
    QuadNode2( int i00_,int i01_,int i10_,int i11_){ setCorners(i00_,i01_,i10_,i11_); };

};

QuadNode2* makeQuad2( QuadNode2* qout, Vec2d p00, Vec2d p01, Vec2d p10, Vec2d p11, std::vector<Vec2d>& ps ){
    int i00=ps.size(); ps.push_back(p00);
    int i01=ps.size(); ps.push_back(p01);
    int i10=ps.size(); ps.push_back(p10);
    int i11=ps.size(); ps.push_back(p11);

    if( qout ){ qout    ->setCorners(i00,i01,i10,i11); }
    else      { qout = new QuadNode2(i00,i01,i10,i11); }
    return qout;
}

void splitOpen( QuadNode2& q, int n, uint8_t mask, std::vector<Vec2d>& ps ){
    int m = (mask&1) + ((mask>>1)&1) + ((mask>>2)&1);
    //q.n = n_*m;
    //q.subs = new QuadNode2[q.n];
    q.allocate( n*m );
    double step = 1.0/n;
    double ljit = 0.4;
    double s0   = 0.5;
    double sjit = 0.1;
    double wmin = 0.1;
    double wmax = 0.3;


    Vec2d c00 = ps[q.i00];
    Vec2d c01 = ps[q.i01];
    Vec2d c10 = ps[q.i10];
    Vec2d c11 = ps[q.i11];

    Vec2d p0l = c00;
    Vec2d p0r = c01;
    int i=0;
    for( int ii=1; ii<=n; ii++ ){

        double fx = step * ii;
        if(ii<n) fx += step*randf(-ljit,ljit);

        double mx = 1-fx;

        Vec2d p1l = c00*mx + c10*fx;
        Vec2d p1r = c01*mx + c11*fx;

        double s  = s0 + randf(-sjit,sjit);
        double w  = randf(-wmin,wmax);
        double s1 = s-w; double m1 = 1-s1;
        double s2 = s+w; double m2 = 1-s2;

        Vec2d p00 = p0l;
        Vec2d p01 = p0l*m1 + p0r*s1;
        Vec2d p02 = p0l*m2 + p0r*s2;
        Vec2d p03 = p0r;

        Vec2d p10 = p1l;
        Vec2d p11 = p1l*m1 + p1r*s1;
        Vec2d p12 = p1l*m2 + p1r*s2;
        Vec2d p13 = p1r;

        if(mask&1){ makeQuad2( &q.subs[i], p00, p10, p01, p11, ps ); i++; };
        if(mask&2){ makeQuad2( &q.subs[i], p01, p11, p02, p12, ps ); i++; };
        if(mask&4){ makeQuad2( &q.subs[i], p02, p12, p03, p13, ps ); i++; };

        p0l = p1l;
        p0r = p1r;
    }
}

void drawQuadNode2Rec( const QuadNode2& q, std::vector<Vec2d>& ps ){
    Vec2d& p00 = ps[q.i00];
    Vec2d& p01 = ps[q.i01];
    Vec2d& p10 = ps[q.i10];
    Vec2d& p11 = ps[q.i11];
    glBegin(GL_LINE_LOOP);
        glVertex3f(p00.x,p00.y, 10.0);
        glVertex3f(p01.x,p01.y, 10.0);
        glVertex3f(p11.x,p11.y, 10.0);
        glVertex3f(p10.x,p10.y, 10.0);
    glEnd();
    if(q.subs){
        for(int i=0; i<q.n; i++){
            drawQuadNode2Rec( q.subs[i], ps );
        }
    }
}




#endif

