
#ifndef  MeshBuilder2_h
#define  MeshBuilder2_h

// ==========================
//  Mesh::Builder2 aims to merge MeshBuilder.h, Truss.h, OMesh (in Mesh.h)  and GLMeshBuilder(in DrawOGL3.h)
// ==========================


#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include "datatypes.h"

#include "testUtils.h"
#include "MeshBuilder.h"

//#include "datatypes.h"

namespace Mesh{

template<typename T>
struct VertT{ // double8
    union{
	    struct{ 
            Vec3T<T> pos;
            Vec3T<T> nor;
            Vec2T<T> uv; 
        };
		struct{ 
            Quat4T<T> lo;
            Quat4T<T> hi;
        };
		T array[8];
	};
    VertT()=default;
    VertT( const Vec3T<T>& pos_, const Vec3T<T>& nor_=Vec3T<T>{0,0,0,0}, const Vec2T<T>& uv_=Vec2T<T>{0,0} ):pos(pos_),nor(nor_),uv(uv_){};
};
using Vert = VertT<double>; 

class Builder2{ public:
    //bool bnor = true;
    //bool bUVs = true;
    //bool bSplitSize = false;
    double max_size = 1.0;
    int   face_type=0;
    Vec2i edge_type=Vec2i{0,0}; // uv edges
    int ov,oov; // previous edge indexes;


    std::vector<Quat4i> blocks;  // {ivert,iedge,itri,ichunk}   // blocks of data, drawing of some complex object create one block
    std::vector<Vert>   verts;   // {pos,nor,uv}
    std::vector<Quat4i> edges;   // {a,b, ?, type}={f|e} or  {a,b, f1, f2} => we can store type of the edge in the 4rd component, or adjacent faces
    std::vector<Quat4i> tris;    // {a,b,c|type}={f|e}   => we can store type of the face in the 4th component
    std::vector<Quat4i> chunks;  // {istrip0, by_type, } can represent polygons, lines_strips, triangles_strips 
    std::vector<int>    strips;  // indexes of primitives in the chunks
    
    int draw_mode = TRIANGLES;
    Vec3f  penColor;
    
    void clear(){ blocks.clear(); verts.clear(); edges.clear(); tris.clear(); chunks.clear(); strips.clear(); }
    void printSizes(){printf( "MeshBuilder::printSizes() blocks=%i verts=%i edges=%i tris=%i chunks=%i strips=%i \n", blocks.size(), verts.size(), edges.size(), tris.size(), chunks.size(), strips.size() );}

    int findClosestVert(const Vec3d& p0,int i0=0,int n=-1){
        if(n==-1) n=verts.size();
        double r2min=1e+300;
        int imin=-1;
        for(int j=0;j<n;j++){ double r2=(verts[i0+j].pos-p0).norm2(); if(r2<r2min){r2min=r2; imin=j;} }
        return imin;
    }


    inline Quat4i latsBlock()const{ return Quat4i{(int)verts.size(),(int)edges.size(),(int)tris.size(),(int)chunks.size()}; }
    inline int block(){ int i=blocks.size(); blocks.push_back( latsBlock() ); return i; };
    inline int vert( const Vec3d& pos, const Vec3d& nor=Vec3dZero, const Vec2d& uv=Vec2dZero ){ verts.push_back(Vert(pos,nor,uv)); return verts.size()-1; }
    inline int edge( int a, int b, int t=-1, int t2=-1 ){ edges.push_back(Quat4i{a,b,t2,t}); return edges.size()-1; }
    inline int tri ( int a, int b, int c,    int t=-1  ){ tris .push_back(Quat4i{a,b,c,t});  return tris .size()-1; }

    inline int stick( Vec3d a, Vec3d b, int t=-1 ){ 
        int ia = vert(a); 
        int ib = vert(b);
        return edge(ia,ib,t); 
    }

    inline int stickTo( int ia, Vec3d b, int t=-1 ){ 
        int ib = vert(b);
        return edge(ia,ib,t); 
    }



    inline int bondsBetweenVertRanges( Vec2i v1s, Vec2i v2s, double Rmax, int et=-1 ){
        double R2max = Rmax*Rmax;
        int n1 = v1s.b-v1s.a;
        int n2 = v2s.b-v2s.a;
        int nb =0;
        for(int i=0; i<n1; i++){
            Vec3d& p1 = verts[v1s.a+i].pos;
            for(int j=0; j<n2; j++){
                Vec3d& p2 = verts[v2s.a+j].pos;
                double r2 = (p1-p2).norm2();
                //printf( "bondsBetweenVertRanges[%i,%i] r=%g Rmax=%g \n", i, j, sqrt(r2), Rmax );
                if( r2<R2max ){ 
                    edge( v1s.a+i, v2s.a+j, et ); 
                    nb++;
                }
            }
        }
        return nb;
    }

    inline int vstrip(Vec3d p0, Vec3d p1, int n, int et=-1 ){
        Vec3d d=p1-p0; 
        d.mul(1./n);
        Vec3d p = p0;
        int i0 = verts.size();
        //printf("vstri(%i){", i0 );
        for(int ii=0; ii<n+1; ii++){
            int i = vert(p);  // printf("%i ", i );
            if((et>-1)&&(ii>0))edge(i-1,i,et);  // long
            p.add(d);
        }
        //printf("}END\n" );
        return i0;
    }

    inline int fstrip( int ip0, int ip1, int n, int ft=-1, Vec2i et={-1,-1} ){ 
        int i0 = tris.size();
        //printf("fstrip(%i,%i)\n", ip0, ip1 );
        for(int ii=0; ii<n+1; ii++){
            int it = ip0*1000+ii;
            //printf("it %i \n", it );
            int i = ip0+ii;
            int j = ip1+ii;
            if(et.x>-1)edge(j,i,et.x);         // perp
            //if(et.x>-1)edge(j,i,it);
            if(ii>0){
                if(et.y>-1)edge(i,j-1,et.y);   // diag
                //if(et.y>-1)edge(j-1,i,it);
                if(ft>-1){ 
                    //tri(j-1,i,j,ft); 
                    //printf("tri(%i,%i,%i)\n", i,j,j-1 );
                    //printf("tri(%i,%i,%i)\n", i,j,j-1 );
                    tri(j-1,i  ,j  ,ft);  //printf("tri(%i,%i,%i)\n", i,j,j-1 );
                    //tri(j-1,i  ,i-1,ft);
                    tri(i  ,j-1,i-1,ft);
                    //tri(ip0,ip1,ip1-1,ft); printf("tri(%i,%i,%i)\n", ip0,ip1,ip1-1 );
                    //tri(ip0,ip1,ip1-1, it ); printf("tri(%i,%i,%i)\n", ip0,ip1,ip1-1 );
                }
            }
        }
        return i0;
    }

    /* Too Complicated logic inside
    inline int vstrip(Vec3d p0, Vec3d p1, int n, Quat4i t=Quat4i{-1,-1,-1,-1}){
        Vec3d d=p1-p0; 
        d.mul(1./n);
        Vec3d p = p0;
        int oi=-1; 
        int i0 = verts.size();
        for(int ii=0; ii<n+1; ii++){
            int i = vert(p);
            int j = i-n;
            if(t.y>-1)edge(j,i,t.y);         // perp
            if(oi>0){
                if(t.x>-1)edge(i-1,i,t.x);   // long
                if(t.y>-1)edge(j-n,i,t.z);   // diag
                if(t.w>-1){ tri(j-1,j,i,t.w); tri(i,i-1,j,t.w); }
            p.add(d);
            oi=i;
        }
        return i0;
    }
    */


    inline int edgst(int v,int t=-1){  int i=edge(ov,v,t); ov=v;         return i; };
    inline int trist(int v,int t=-1){  int i=tri (ov,v,t); oov=ov; ov=v; return i; };

    inline int chunk( int t, int n, int* inds=0 ){ 
        int i0=strips.size(); 
        chunks.push_back(Quat4i{i0,n,-1,t}); 
        if(inds){
            //for(int i=0;i<n;i++){ strips.push_back(inds[i]); } return chunks.size()-1;  // ToDo: can we copy array to vectro faster ?
            // https://stackoverflow.com/questions/259297/how-do-you-copy-the-contents-of-an-array-to-a-stdvector-in-c-without-looping
            strips.insert(strips.end(), &inds[0], &inds[n]);
            //strips.resize(strips.size()+n); memcpy(&strips[strips.size()-n], &inds[0], n*sizeof(int) );
        }
        return chunks.size()-1;
    }  

    inline void quad( Quat4i q, int face_type=-2, int edge_type=-3 ){
        if( face_type>-2 ){ tri( q.x, q.y, q.z, face_type ); tri( q.x, q.w, q.y, face_type ); }
        if( edge_type>-3 ){ edge( q.y, q.z, edge_type );  }
        if( edge_type>-2 ){ edge( q.x, q.y, edge_type ); edge( q.y, q.w, edge_type ); edge( q.w, q.z, edge_type ); edge( q.z, q.x, edge_type ); }
    }

    int rope( int ip0, int ip1, int typ=-1, int n=1 ){   
        //printf( "MeshBuilder2::rope(%i,%i,n=%i,t=%i)\n", ip0,ip1, n,typ );
        int i0 = verts.size();
        Vec3d p0,d;
        if( n!=1 ){
            p0 = verts[ip0].pos;
            d  = verts[ip1].pos-p0;
            if( n<0 ){ 
                double l = d.norm(); n=(int)(l/max_size)+1;   
                //printf( "rope() l=%g n=%i ]n", l, n ); 
            }         
            //printf( "rope() n=%i d(%g,%g,%g)  p0(%g,%g,%g) \n", n,  d.x,d.y,d.z,  p0.x,p0.y,p0.z );   
            d.mul(1./(double)n);  // we need this only it n>1
            //printf( "rope() n=%i d(%g,%g,%g)  p0(%g,%g,%g) \n", n,  d.x,d.y,d.z,  p0.x,p0.y,p0.z );   
        }
        int oi=ip0;
        for(int ii=1; ii<n; ii++){
            Vec3d p = p0+d*ii;
            int i=vert(p);
            edge(oi,i, typ);
            //Vec3d& v = verts.back().pos;
            //printf( "rope(%i,%i)[%i] p(%g,%g,%g)\n",  ip0,ip1,ii, v.x,v.y,v.z );
            oi=i;
        }
        edge( oi,ip1, typ );
        // ToDo: implementation by LINE_STRIP ?
        return i0;
    };

    inline int ring( Vec3d p, Vec3d a, Vec3d b, Vec2d cs, int n=4 ){
        Vec2d rot = Vec2d{1.0,0.0};
        int i0  = verts.size();
        int ip0 = ip0 + n;
        for(int i=0; i<n; i++){
            int ip = vert(p + a*rot.x + b*rot.y );
            edge( ip0, ip, edge_type.x );
            rot.mul_cmplx( cs );
            ip0 = ip;
        }
        return i0;
    }
    inline int ring( Vec3d p, Vec3d ax, Vec3d up, double R, int n=4 ){
        double dphi = 2*M_PI/n;
        if( n<0 ){ double dL   = R*dphi; int n2=(int)(dL/max_size)+1; }
        Vec2d cs; cs.fromAngle( dphi );
        Vec3d a,b; 
        b.set_cross(ax,up); b.normalize();
        a.set_cross(b ,ax); a.normalize();
        return ring(p,a*R,b*R,cs,n );
    }

    void tube( Vec3d p0, Vec3d p1, Vec3d up, Vec2d R, Vec2i n={4,1} ){
        Vec3d ax = p1-p0; double L = ax.normalize();
        if( n.x<0 ){ n.x=(int)(L/max_size)+1; }
        if( n.y<0 ){ double r = fmax(R.x,R.y); n.y=(int)(r/max_size)+1; }
        double dL = L/n.y;
        Vec2d cs; cs.fromAngle( 2*M_PI/n.x );
        Vec3d a,b; 
        b.set_cross(ax,up); b.normalize();
        a.set_cross(b ,ax); a.normalize();
        int oe;
        double dR = (R.y-R.x)/n.y;
        for(int iy=0; iy<n.y; iy++){
            int e = edges.size();
            double r = R.x + dR*iy;
            ring( p0+ax*(dL*iy), a*r, b*r, cs, n.x );
            if(iy>0){
                for(int ix=0; ix<n.x; ix++){
                    Vec2i e1 = edges[e+ix].f.xy();
                    Vec2i e2 = edges[oe+ix].f.xy(); 
                    tri( e1.x, e1.y, e2.x, face_type ); 
                    tri( e1.y, e2.y, e2.x, face_type );
                    edge( e1.x, e2.x, edge_type.y );        
                    //edge( e1.y, e2.y, edge_type.y );  // it would do it twice
                    edge( e1.y, e2.x, edge_type.y );    // diagonal edge
                }
            }
            oe=e;
        }
    }

    //int strip( Vec3d p0, Vec3d p1, Vec3d d, int n,  ){
    //    for(int i=0;i<n;i++){
    //    }
    //}


    //inline int vstrip(Vec3d p0, Vec3d p1, int n, int et=-1 ){
    //inline int fstrip( int ip1, int ip0, int n, int ft=-1, Vec2i et={-1,-1} ){ 
    

    int plate( Vec3d p00, Vec3d p01, Vec3d p10, Vec3d p11, Quat4i t={-1,-1,-1,-1}, Vec2i n={1,1}, int fillType=1 ){
        Vec3d dx0=p01-p00; double lx0=dx0.norm();
        Vec3d dx1=p11-p10; double lx1=dx1.norm();
        Vec3d dy0=p10-p00; double ly0=dy0.norm();
        Vec3d dy1=p11-p01; double ly1=dy1.norm();
        if(n.x<0){ int n1=(int)(lx0/max_size)+1; int n2=(int)(lx1/max_size)+1;  n.x=_max(n1,n2); }
        if(n.y<0){ int n1=(int)(ly0/max_size)+1; int n2=(int)(ly1/max_size)+1;  n.y=_max(n1,n2); }
        //n.x=5;
        //n.y=5;
        //Quat4i t_;
        int i0=verts.size();
        int oiv;
        for(int iy=0;iy<n.y+1;iy++){
            double cy=iy/(double)n.y;  double my=1-cy; 
            Vec3d p0y = p00*my+p10*cy;
            Vec3d p1y = p01*my+p11*cy;
            //if(iy==0){ t_=Quat4i{t.x,-2,-2,-2}; }else{t=t;}
            //strip(p0,p1, n.x, t);
            int iv = vstrip( p0y, p1y, n.x, t.x );
            if(iy>0){
                fstrip( oiv, iv, n.x, t.w, {t.y,t.z} );
                //fstrip( oiv, iv, n.x, -2, {t.y,t.z} );

            }
            oiv=iv;
        }
        return i0;
    }

    // ToDo: This is too complicated, put we should remove it or move it elsewhere
    int plate_quad( int ip00, int ip01, int ip10, int ip11, Quat4i typs={-1,-1,-1,-1}, Vec2i n={1,1}, int fillType=1 ){ 
        printf( "plate_quad(%i,%i;%i,%i) typs{%i,%i,%i,%i} n{%i,%i} fillType=%i \n", ip00,ip01,ip10,ip11, typs.x,typs.y,typs.z,typs.w, n.x, n.y, fillType ); 
        int i0 = verts.size();
        Vec3d p00,p01,p10,p11;
        //DEBUG
        p00 = verts[ip00].pos; p01 = verts[ip01].pos; p10 = verts[ip10].pos; p11 = verts[ip11].pos;
        if(n.x<0){ int n1=(int)((p00-p01).norm()/max_size)+1; int n2=(int)((p10-p11).norm()/max_size)+1;  n.x=_max(n1,n2); }
        if(n.y<0){ int n1=(int)((p10-p00).norm()/max_size)+1; int n2=(int)((p11-p01).norm()/max_size)+1;  n.y=_max(n1,n2); }
        printf( "p00(%g,%g,%g) p01(%g,%g,%g) p10(%g,%g,%g) p11(%g,%g,%g)\n", p00.x,p00.y,p00.z,  p01.x,p01.y,p01.z,  p10.x,p10.y,p10.z,  p11.x,p11.y,p11.z );
        n.x=3; n.y=3;
        printf( "plate_quad n(%i,%i)\n", n.x, n.y);
        int ex[n.x];
        int ey[n.y]; 
        //DEBUG
        // --- edges
        //typs.x = 234234;
        ex[0    ]=edges.size(); int i0x = rope( ip00, ip01, typs.x, n.y );
        ex[n.x-1]=edges.size(); int i1x = rope( ip10, ip11, typs.x, n.y );
        ey[0    ]=edges.size(); int i0y = rope( ip00, ip10, typs.x, n.x );
        ey[n.y-1]=edges.size(); int i1y = rope( ip01, ip11, typs.x, n.x );
        //DEBUG

        //rope( i0y+5, i1y+5, typs.y, n.x );
        // WARRNING: This will not work properly because the points generated here for x-strips and y-strips are generated twice, they are not colapsed => we should generate points first, only then add edges
        for(int iy=1; iy<n.y-1; iy++){ ex[iy]=rope( i0y+iy, i1y+iy, typs.y, n.x ); }
        for(int ix=1; ix<n.x-1; ix++){ ex[ix]=rope( i0x+ix, i1x+ix, typs.y, n.y ); }
        
        /*
        //DEBUG
        for(int iy=1; iy<n.y; iy++){
            for(int ix=1; ix<n.x; ix++){
                Vec2i e1= edges[ex[ix+1]].f.xy();
                Vec2i e2= edges[ex[ix  ]].f.xy();
                if(fillType>0){
                    tri( e1.x, e1.y, e2.x, typs.z );
                    tri( e1.y, e2.y, e2.x, typs.z );
                    edge( e1.y, e2.x, typs.y );        // diagonal edge
                }else{
                    tri( e2.x, e2.y, e1.x, typs.z );
                    tri( e2.y, e1.y, e1.x, typs.z );
                    edge( e2.y, e1.x, typs.y );
                }
            }
        }
        */
        //DEBUG
        // ToDo: implementation by TRIANGLE_STRIP ?
        return i0;
    };

    int export_pos( Vec3d* ps, int i0=0, int i1=-1 ){   if(i1<0){ i1=verts.size()-i1; };
        for(int i=i0; i<=i1; i++){ ps[i]=verts[i].pos; }
        return i1-i0+1;
    }

    int export_pos( float4* ps, int i0=0, int i1=-1 ){     if(i1<0){ i1=verts.size()-i1; };
        for(int i=i0; i<=i1; i++){ ps[i]=*(float4*)&verts[i].lo;}
        return i1-i0+1;
    }

    int export_edges( Vec2i* eds, int i0=0, int i1=-1 ){   if(i1<0){ i1=edges.size()-i1; };
        for(int i=i0; i<=i1; i++){ eds[i]=edges[i].lo; }
        return i1-i0+1;
    }

    int export_tris( Quat4i* tri, int i0=0, int i1=-1 ){  if(i1<0){ i1=edges.size()-i1; };
        for(int i=i0; i<=i1; i++){ tri[i]=tris[i]; }
        return i1-i0+1;
    }



/*

    void   newSub( int mode = TRIANGLES ){ subs.push_back( Vec3i{(int)vpos.size(),(int)inds.size(),(int)mode} ); }
    inline Vec2i subVertRange(int i){ int i0=0; if(i<0)i=subs.size()+i; i0=subs[i-1].a; return {i0,subs[i].a}; }
    inline Vec2i subIndRange (int i){ int i0=0; if(i<0)i=subs.size()+i; i0=subs[i-1].b; return {i0,subs[i].b}; }

    void move ( Vec2i iv, Vec3f shift ){ for(int i=iv.a; i<iv.b; i++){ vpos[i].add(shift); } }
    void scale( Vec2i iv, Vec3f sc    ){
        Vec3f invSc = {1/sc.x,1/sc.y,1/sc.z};
        for(int i=iv.a; i<iv.b; i++){
            vpos[i].mul(sc);
            if(bnor) vnor[i].mul(invSc); vnor[i].normalize();
        }
    }
    void rotate( Vec2i iv, Vec3f p0, Vec3f p1, float angle ){
        Vec3f uax=p1-p0; uax.normalize();
        Vec2f cs; cs.fromAngle(angle);
        for(int i=iv.a; i<iv.b; i++){
            Vec3f v = vpos[i]-p0;  v.rotate_csa(cs.a,cs.b,uax); vpos[i]=v+p0;
            if(bnor) vnor[i].rotate_csa(cs.a,cs.b,uax);
        }
    }

    void applyMatrix( Vec2i iv, Mat3f mat ){
        for(int i=iv.a; i<iv.b; i++){
            mat.dot_to(vpos[i], vpos[i]);
            if(bnor) mat.dot_to(vnor[i],vnor[i]);
        }
    }

    void applyMatrixT( Vec2i iv, Mat3f mat ){
        for(int i=iv.a; i<iv.b; i++){
            mat.dot_to_T(vpos[i], vpos[i]);
            if(bnor) mat.dot_to_T(vnor[i],vnor[i]);
        }
    }

    void duplicateSub( int i ){
        Vec2i ivs = subVertRange(i);
        Vec2i iis = subIndRange (i);
        int di    = vpos.size() - ivs.a;
        for(int i=iis.a; i<iis.b; i++){
            inds.push_back( inds[i]+di );
        }
        for(int i=ivs.a; i<ivs.b; i++){
            vpos         .push_back( vpos[i] );
            if(bnor) vnor.push_back( vnor[i] );
            if(bUVs) vUVs.push_back( vUVs[i] );
        }
        newSub();
    }

    void addLine( Vec3f p1, Vec3f p2 ){
        vpos.push_back(p1);  vnor.push_back(penColor);
        vpos.push_back(p2);  vnor.push_back(penColor);
    };

    void addLine( Vec3d p1, Vec3d p2 ){
        vpos.push_back((Vec3f)p1);  vnor.push_back(penColor);
        vpos.push_back((Vec3f)p2);  vnor.push_back(penColor);
    };

    void addPointCross( Vec3f p, float d ){
        addLine( p+(Vec3f){d,0,0}, p+(Vec3f){-d, 0, 0} );
        addLine( p+(Vec3f){0,d,0}, p+(Vec3f){ 0,-d, 0} );
        addLine( p+(Vec3f){0,0,d}, p+(Vec3f){ 0, 0,-d} );
    };

    void addArrow( Vec3f p1, Vec3f p2, float d ){

    };

    void addLines( int n, int * inds, Vec3f* verts ){
        for(int i=0; i<n; i++){
            int i2=i*2;
            addLine( verts[inds[i2]], verts[inds[i2+1]] );
        }
    }

    void addTriangle( const Vec3f& a, const Vec3f& b, const Vec3f& c ){
        Vec3f nor; nor.set_cross( a-c, b-c ); nor.normalize();
        vpos.push_back( a ); vnor.push_back( nor );
        vpos.push_back( b ); vnor.push_back( nor );
        vpos.push_back( c ); vnor.push_back( nor );
    }

    void addQuad( const Vec3f& p00, const Vec3f& p01, const Vec3f& p10, const Vec3f& p11 ){
        addTriangle( p00, p01, p10 );
        addTriangle( p11, p10, p01 );
    }

    void addTube4( Vec3f p1, Vec3f p2, Vec3f up, float r1, float r2 ){
        up.normalize();
        Vec3f lf; lf.set_cross( p2-p1, up ); lf.normalize();
        addQuad( p1+up*r1, p1+lf*r1, p2+up*r2, p2+lf*r2 );
        addQuad( p1+lf*r1, p1-up*r1, p2+lf*r2, p2-up*r2 );
        addQuad( p1-up*r1, p1-lf*r1, p2-up*r2, p2-lf*r2 );
        addQuad( p1-lf*r1, p1+up*r1, p2-lf*r2, p2+up*r2 );
    };

    int addCircleAxis( int n, const Vec3f& pos, const Vec3f& v0, const Vec3f& uaxis, float R ){
        printf( "MeshBuilder::addCircleAxis() n=%i pos(%6.3f,%6.3f,%6.3f) v0(%6.3f,%6.3f,%6.3f) uaxis(%6.3f,%6.3f,%6.3f) R=%3.3f )\n", n, pos.x,pos.y,pos.z, v0.x,v0.y,v0.z, uaxis.x,uaxis.y,uaxis.z, R );
        float dphi = 2*M_PI/n;
        float dca  = cos( dphi );
        float dsa  = sin( dphi );
        int nvert=0;
        Vec3f v; v.set(v0);
        //glBegin( GL_LINE_LOOP );
        for( int i=0; i<n; i++ ){
            //glVertex3f( pos.x+v.x*R, pos.y+v.y*R, pos.z+v.z*R ); nvert++;
            //printf( " drawCircleAxis %i (%3.3f,%3.3f,%3.3f) \n", i, v.x, v.y, v.z );
            Vec3f ov = v;
            v.rotate_csa( dca, dsa, uaxis );
            addLine( ov, v );
        }
        //glEnd();
        return nvert;
    }

    void moveSub     ( int i, Vec3f shift ){ move ( subVertRange(i), shift); }
    void scaleSub    ( int i, Vec3f sc    ){ scale( subVertRange(i), sc   ); }
    void rotateSub   ( int i, Vec3f p0, Vec3f p1, float angle ){ rotate( subVertRange(i), p0, p1, angle ); }


    void printSizes(){
        printf( "MeshBuilder::printSizes() nsubs=%i nvpos=%i nvnor=%i nvUVs=%i ninds=%i \n", subs.size(), vpos.size(), vnor.size(), vUVs.size(), inds.size() );
    }

  void write_obj( char* fname ){
        printf( "MeshBuilder::write_obj(%s)\n", fname );
        int nsubs  = subs.size();
        Vec3i osub = subs[0];

        FILE * pFile;
        pFile = fopen (fname,"w");

        int nobj =0;
        int nvert=0;
        int nnor =0;
        for( int i=1; i<nsubs; i++ ){
            Vec3i sub = subs[i];
            int mode  = osub.z;
            if      (mode == TRIANGLES ){      // un-indexed triangles
                //printf(         "o OBJ_TRIANGLES.%i  [ %i ... %i ] \n", nobj, osub.x, sub.x );
                fprintf( pFile, "o OBJ_TRIANGLES.%i \n", nobj ); nobj++;
                int iii = 0;
                for(int j=osub.x; j<sub.x; j++){
                    Vec3f vp = vpos[j];
                    Vec3f vn = vnor[j];
                    fprintf(pFile, "v   %f %f %f\n", vp.x, vp.y, vp.z ); nvert++;
                    fprintf(pFile, "vn  %f %f %f\n", vn.x, vn.y, vn.z ); nnor ++;
                    if(iii%3==2) fprintf( pFile, "f %i//%i %i//%i %i//%i \n", nvert-2,nnor-2,  nvert-1,nnor-1,   nvert,nnor );
                    iii++;
                }
            }else if(mode == TRIANGLE_STRIP ) {  // Indexed Triangles
                //printf(         "o OBJ_TRIANGLE_STRIP.%i  [ %i ... %i ] \n", nobj, osub.x, sub.x );
                fprintf( pFile, "o OBJ_TRIANGLE_STRIP.%i \n", nobj ); nobj++;
                
                //  //  ToDo - This does not work for some reason => brute force polygonization
                // int iv0 = nvert-osub.x;
                // int in0 = nnor -osub.x;
                // for(int j=osub.x; j<sub.x; j++){
                //     Vec3f vp = vpos[j];
                //     Vec3f vn = vnor[j];
                //     fprintf(pFile, "v   %f %f %f\n", vp.x, vp.y, vp.z ); nvert++;
                //     fprintf(pFile, "vn  %f %f %f\n", vn.x, vn.y, vn.z ); nnor ++;
                // }
                // int iii = 0;
                // for(int j=osub.y; j<sub.y; j+=3 ){
                //     int i0=inds[j  ];
                //     int i1=inds[j+1];
                //     int i2=inds[j+2];
                //     fprintf( pFile, "f %i//%i %i//%i %i//%i \n",    i0+iv0, i0+in0,       i1+iv0, i1+in0,     i2+iv0, i2+in0 );
                //     iii++;
                // }
                
                int iii = 0;
                for(int j=osub.y; j<sub.y; j++ ){
                    int ii = inds[j  ];
                    Vec3f vp = vpos[ii];
                    Vec3f vn = vnor[ii];
                    fprintf(pFile, "vn  %f %f %f\n", vn.x, vn.y, vn.z ); nnor ++;
                    fprintf(pFile, "v   %f %f %f\n", vp.x, vp.y, vp.z ); nvert++;
                    if(iii%3==2) fprintf( pFile, "f %i//%i %i//%i %i//%i \n", nvert-2,nnor-2,  nvert-1,nnor-1,   nvert,nnor );
                    iii++;
                }

            }else if(mode == LINES ) {                
                // printf(         "o OBJ_LINES.%i  [ %i ... %i ] \n", nobj, osub.x, sub.x );
                // fprintf( pFile, "o OBJ_LINES.%i \n", nobj ); nobj++;
                // int iii = 0;
                // for(int j=osub.x; j<sub.x; j++){
                //     //Vec3f vnor = mesh.vnor[j];
                //     Vec3f vp = vpos[j];
                //     fprintf(pFile, "v   %f %f %f\n", vp.x, vp.y, vp.z ); nvert++;
                //     //if(iii%2==1) printf( "l %i %i \n", nvert-1, nvert );
                //     if(iii%2==1) printf( "f %i %i \n", nvert-1, nvert );
                //     iii++;
                // }
            }
            osub=sub;
        }
        fclose(pFile);
    }
*/

}; // class Mesh::Builder2



}; // namespace Mesh

#endif


