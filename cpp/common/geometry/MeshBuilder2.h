
#ifndef  MeshBuilder2_h
#define  MeshBuilder2_h

// ==========================
//  Mesh::Builder2 aims to merge MeshBuilder.h, Truss.h, OMesh (in Mesh.h)  and GLMeshBuilder(in DrawOGL3.h)
// ==========================


#include <vector>
#include <unordered_map>
#include <unordered_set>


#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "datatypes.h"
#include "Slots.h"

//#include "MeshBuilder.h"

using LoopDict = std::unordered_map<int,Slots<int,2>>;

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
    std::vector<Quat4i> chunks;  // {iv0,ie0,n,type} can represent polygons, lines_strips, triangles_strips 
    std::vector<int>    strips;  // indexes of primitives in the chunks

    bool use_vert2edge  = false;
    bool bPolygonToTris = true;
    std::unordered_map<uint64_t,int> vert2edge; // map from pair of vert index to vert index in verts

    // Edit settings
    bool bExitError    = true;
    int ngon_max       = 16;
    double R_snapVert  = 0.1;   // tolerance for snapping vertices used in findVert()
    bool bAdditiveSelect = false;
    int selection_mode = 3;  // 0-none, 1-vert, 2-edge, 3-face
    enum class ChunkType{ face=0, edgestrip=1, trianglestrip=2 };
    enum class SelectionMode{ vert=1, edge=2, face=3 };
    std::vector<int>        selection; //  vector is orderd - usefull for e.g. edge-loops    indices of selected vertices (eventually edges, faces ? )
    std::unordered_set<int> selset;    // edge index for vert

    //int draw_mode = TRIANGLES;
    Vec3f  penColor;
    

    // ======= Inline Functions

    inline Quat4i latsBlock()const{ return Quat4i{(int)verts.size(),(int)edges.size(),(int)tris.size(),(int)chunks.size()}; }
    inline int block(){ int i=blocks.size(); blocks.push_back( latsBlock() ); return i; };
    inline int vert( const Vec3d& pos, const Vec3d& nor=Vec3dZero, const Vec2d& uv=Vec2dZero ){ 
        //printf( "Mesh::Builder2::vert() %3i pos: %16.10f %16.10f %16.10f \n", verts.size(), pos.x,pos.y,pos.z );
        verts.push_back(Vert(pos,nor,uv)); return verts.size()-1; 
    }
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

    inline int getOtherEdgeVert(int ie, int iv){
        Quat4i& e = edges[ie];
        return e.x==iv ? e.y : e.x;
    }
    
    inline int edgst(int v,int t=-1){  int i=edge(ov,v,t); ov=v;         return i; };
    inline int trist(int v,int t=-1){  int i=tri (ov,v,t); oov=ov; ov=v; return i; };
    inline int chunk( const Quat4i ch ){
        chunks.push_back( ch );
        return chunks.size()-1;
    }
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

    inline int selectVertRange(int i0, int iend){
        for(int i=i0;i<iend;i++){ selection.push_back(i); }
        return selection.size();
    }

    inline int* getChunkStrip( int ich ){
        Quat4i ch = chunks[ich];
        return strips.data() + ch.x;
    };


    // ======= Functions


    int selectVertsAlongLine( Vec3d p0, Vec3d p1, double r=0.1, bool bSort=true );

    int selectVertsAlongPolyline( double r=0.1, bool bSort=true );

    int conected_vertex( const Vec3d& p, int stickType, int n, int* iverts );
    int select_in_cylinder( const Vec3d& p0, const Vec3d& fw, double r, double l );
    int select_in_box( const Vec3d& p0, const Vec3d& fw, const Vec3d& up, const Vec3d& Lmin, const Vec3d& Lmax );
    int make_anchor_point( const Vec3d& p, int stickType, const Vec3d& fw, double r, double l );

    int extrudeFace( int ich, double L, Quat4i stickTypes=Quat4i{-1,-1,-1,-1}, Quat4i maks={1,1,1,1} );

    // ======= Functions

    Vec3d getCOG(int n, const int* ivs) const;
    void alling_polygons( int n, const int* ivs1, int* ivs2, int ipiv=0 );
    int bridge_quads( Quat4i q1, Quat4i q2, int n, Quat4i stickTypes, Quat4i mask, bool bAlling=false );
    int extrudeVertLoop( int n, int* iverts, Vec3d d, bool bEdges, bool bFace, bool bTris, bool bSort );
    int loadChunk( int ich, int* iedges=0, int* iverts=0 );
    Vec3d getChunkNormal( int ich );

    void clear();


    bool sortPotentialEdgeLoop( int n, Vec2i* edges, int* iverts );
    bool sortEdgeLoop( int n, int* iedges, int* iverts=0 );
    int findEdgeByVerts_brute( Vec2i verts );
    int findEdgeByVerts_map( const Vec2i verts );
    int findEdgeByVerts( const Vec2i verts );
    int findOrAddEdges( const Vec2i verts, int t=-1, int t2=-1 );
    void buildVerts2Edge();


    int plateBetweenVertStrips( int n, int* ivs1, int* ivs2, int nsub );
    int plateBetweenEdges( int nsub=1, double r=0.1, bool bSort=true );

    int polygonChunk( int n, int* iedges, int* ivs, bool bPolygonToTris );
    int polygon( int n, int* iedges );
    int polygonToTris( int i );
    int selectionToFace();
    int clearSelection();
    int pickVertex( const Vec3d& ray0, const Vec3d& hRay, double R );
    int pickEdge( const Vec3d& ro, const Vec3d& rh, double Rmax );
    int toggleSelSet(  int i );
    int pickTriangle( const Vec3d& ro, const Vec3d& rh, bool bReturnFace=false );
    int pickEdgeSelect( const Vec3d& ro, const Vec3d& rh, double Rmax );
    int pickSelect( const Vec3d& ro, const Vec3d& rh, double Rmax );
    int selectRectEdge( const Vec3d& p0, const Vec3d& p1, const Mat3d& rot );
    int selectRectVert( const Vec3d& p0, const Vec3d& p1, const Mat3d& rot );
    int selectRect( const Vec3d& p0, const Vec3d& p1, const Mat3d& rot  );
    void makeSelectrionUnique();
    int findClosestVert(const Vec3d& p0,int i0=0,int n=-1);
    int findVert(const Vec3d& p0, double Rmax, int n=-1, int* sel=0 );
    int bondsBetweenVertRanges( Vec2i v1s, Vec2i v2s, double Rmax, int et=-1 );
    int vstrip(Vec3d p0, Vec3d p1, int n, int et=-1 );
    int fstrip( int ip0, int ip1, int n, int ft=-1, Vec2i et={-1,-1} );
    void box( Vec3d p, Vec3d ls, Mat3d rot );
    void snapBoxFace( const Vec3d& p0, const Mat3d& rot, double La, double Lb );
    void frustrumFace( const Vec3d& p0, const Mat3d& rot, double La, double Lb, double h, double Lbh, double Lah  );
    void snapFrustrumFace( const Vec3d& p0, const Mat3d& rot, double La, double Lb, double h, double Lbh, double Lah, bool bFace=true );
    void prismFace( const Vec3d& p0, const Mat3d& rot, double La, double Lb, double h, double Lbh );
    void snapPrismFace( const Vec3d& p0, const Mat3d& rot, double La, double Lb, double h, double Lbh, bool bFace=true );

    void quad( Quat4i q, int face_type=-2, int edge_type=-3 );
    int  rope( int ip0,  int ip1,  int typ=-1, int n=1 );
    int  ring( Vec3d p,  Vec3d a,  Vec3d b, Vec2d cs, int n=4 );
    int  ring( Vec3d p,  Vec3d ax, Vec3d up, double R, int n=4 );
    void tube( Vec3d p0, Vec3d p1, Vec3d up, Vec2d R, Vec2i n={4,1} );
    //int strip( Vec3d p0, Vec3d p1, Vec3d d, int n,  ){
    //    for(int i=0;i<n;i++){
    //    }
    //}
    //inline int vstrip(Vec3d p0, Vec3d p1, int n, int et=-1 ){
    //inline int fstrip( int ip1, int ip0, int n, int ft=-1, Vec2i et={-1,-1} ){ 
    int plate( Vec3d p00, Vec3d p01, Vec3d p10, Vec3d p11, Quat4i t={-1,-1,-1,-1}, Vec2i n={1,1}, int fillType=1 );

    // ToDo: This is too complicated, put we should remove it or move it elsewhere
    int plate_quad( int ip00, int ip01, int ip10, int ip11, Quat4i typs={-1,-1,-1,-1}, Vec2i n={1,1}, int fillType=1 );
    int export_pos ( Vec3d* ps, int i0=0, int i1=-1 );
    int  export_pos( float4* ps, int i0=0, int i1=-1 );
    int  export_edges( Vec2i* eds, int i0=0, int i1=-1 );
    int  export_tris( Quat4i* tri, int i0=0, int i1=-1 );
    
    void printSelection();
    void printSelectedVerts();

    void printSizes();
    void printVerts();
    void printEdges();



    int girder1( Vec3d p0, Vec3d p1, Vec3d up, int n, double width, Quat4i stickTypes, bool bCaps=false );
    int triangle_strip( Vec3d p0, Vec3d p1, Vec3d up, int n, double width, int stickType, bool bCaps=false );
    int plateOnGriders( Vec2i ns, Vec2i prange1, Vec2i prange2, Vec2i byN, Vec2i offs, Vec2d span1, Vec2d span2, Quat4i stickTypes );
    int girder1_caps( int ip0, int ip1, int kind );
    int girder1( int ip0, int ip1, Vec3d up, int n, double width, Quat4i stickTypes );
    int wheel( Vec3d p0, Vec3d p1, Vec3d ax, int n, Vec2d wh, Quat4i stickTypes );
    int ngon( Vec3d p0, Vec3d p1, Vec3d ax, int n,  int stickType );
    int rope( Vec3d p0, Vec3d p1, int n,  int stickType );
    int panel( Vec3d p00, Vec3d p01, Vec3d p10, Vec3d p11, Vec2i n, double width, Quat4i stickTypes );

}; // class Mesh::Builder2



}; // namespace Mesh

#endif


