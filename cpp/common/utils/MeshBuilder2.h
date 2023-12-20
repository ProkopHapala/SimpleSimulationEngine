
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

#include "MeshBuilder.h"

namespace Mesh{

struct Vert{ // double8
    Vec3d pos;
    Vec3d nor;  // normal ( can also store color )
    Vec2d uv;  
    Vert()=default;
    Vert( const Vec3d& pos_, const Vec3d& nor_=Vec3dZero, const Vec2d& uv_=Vec2dZero ):pos(pos_),nor(nor_),uv(uv_){};
};

class Builder2{ public:
    bool bnor = true;
    bool bUVs = true;
    std::vector<Quat4i> subs;    // {ivert,iedge,itri,ichunk}   // blocks of data, drawing of some complex object create one block
    std::vector<Vert>   verts;   // {pos,nor,uv}
    std::vector<Quat4i> edges;   // {a,b, ?, type}={f|e} or  {a,b, f1, f2} => we can store type of the edge in the 4rd component, or adjacent faces
    std::vector<Quat4i> tris;    // {a,b,c|type}={f|e}   => we can store type of the face in the 4th component
    std::vector<Quat4i> chunks;  // {istrip0, by_type, } can represent polygons, lines_strips, triangles_strips 
    std::vector<int>    strips;  // indexes of primitives in the chunks
    
    int draw_mode = TRIANGLES;
    Vec3f  penColor;
    
    void clear(){ subs.clear(); verts.clear(); edges.clear(); tris.clear(); chunks.clear(); strips.clear(); }

    int vert( const Vec3d& pos, const Vec3d& nor, const Vec2d& uv=Vec2dZero ){ verts.push_back(Vert(pos,nor,uv)); return verts.size()-1; }
    int edge( int a, int b, int t=-1, int t2=-1 ){ edges.push_back(Quat4i{a,b,t2,t}); return edges.size()-1; }
    int tri ( int a, int b, int c,    int t=-1  ){ tris .push_back(Quat4i{a,b,c,t});  return tris .size()-1; }

    int chunk( int t, int n, int* inds ){ 
        int i0=strips.size(); 
        chunks.push_back(Quat4i{i0,n,-1,t}); 
        //for(int i=0;i<n;i++){ strips.push_back(inds[i]); } return chunks.size()-1;  // ToDo: can we copy array to vectro faster ?
        // https://stackoverflow.com/questions/259297/how-do-you-copy-the-contents-of-an-array-to-a-stdvector-in-c-without-looping
        strips.insert(strips.end(), &inds[0], &inds[n]);
        //strips.resize(strips.size()+n); memcpy(&strips[strips.size()-n], &inds[0], n*sizeof(int) );
        return chunks.size()-1;
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


