
#ifndef  MeshBuilder2_h
#define  MeshBuilder2_h

// ==========================
//  Mesh::Builder2 aims to merge MeshBuilder.h, Truss.h, OMesh (in Mesh.h)  and GLMeshBuilder(in DrawOGL3.h)
// ==========================


#include <vector>
#include <unordered_map>
#include <unordered_set>

// #include <ranges>
// #include <functional>
// #include <algorithm>



#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "datatypes.h"
#include "Slots.h"
#include "CMesh.h"

#include "globals.h"
#include "Buckets.h"

#include "Selection.h"





//#include "MeshBuilder.h"

using LoopDict = std::unordered_map<int,Slots<int,2>>;

//#include "datatypes.h"

namespace Mesh{

// Flags controlling Builder2::bevel options
namespace BevelFlags{ enum{ EdgeWedge=1<<0, DiagWedge=1<<1, DiagFlat=1<<2, FacesFlat=1<<3, FacesWedge=1<<4 }; }
//inline BevelFlags operator|(BevelFlags a, BevelFlags b){ return static_cast<BevelFlags>(static_cast<uint16_t>(a)|static_cast<uint16_t>(b)); }
//inline BevelFlags operator&(BevelFlags a, BevelFlags b){ return static_cast<BevelFlags>(static_cast<uint16_t>(a)&static_cast<uint16_t>(b)); }

// #define _unpack_BevelFlags(flags) \
//     bool bEdgeWedge  = (bool)(flags & BevelFlags::EdgeWedge); \
//     bool bDiagWedge  = (bool)(flags & BevelFlags::DiagWedge); \
//     bool bDiagFlat   = (bool)(flags & BevelFlags::DiagFlat); \
//     bool bFacesFlat  = (bool)(flags & BevelFlags::FacesFlat); \
//     bool bFacesWedge = (bool)(flags & BevelFlags::FacesWedge);

#define _unpack_BevelFlag(i,name) bool b##name=(bool)(i& BevelFlags::name); 
#define _unpack_BevelFlags(i) \
_unpack_BevelFlag(i,EdgeWedge )\
_unpack_BevelFlag(i,DiagWedge )\
_unpack_BevelFlag(i,DiagFlat  )\
_unpack_BevelFlag(i,FacesFlat )\
_unpack_BevelFlag(i,FacesWedge)


static double RvertCollapse = 1e-3;

    struct ObjMask {
        enum : uint8_t {
            Verts    = 1 << 0,
            Normals  = 1 << 1,
            UVs      = 1 << 2,
            Tris     = 1 << 3,
            Polygons = 1 << 4,
            Edges    = 1 << 5,
        };
    };


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

class Builder2 : public SelectionBanks { public:
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
    // polygon chunks structure: {iv0,ie0,n,type} where type is ChunkType::face
    

    

    bool use_vert2edge  = false;
    bool bPolygonToTris = true;
    std::unordered_map<uint64_t,int> vert2edge;        // map from pair of vert index to vert index in edges
    
    //std::vector<std::unordered_set<int>> edgesOfVerts; // edges of each vertex
    Buckets edgesOfVerts;

    // Edit settings
    bool bExitError    = true;
    int ngon_max       = 16;
    double R_snapVert  = 0.1;   // tolerance for snapping vertices used in findVert()
    bool bAdditiveSelect = false;
    int selection_mode = 3;  // 0-none, 1-vert, 2-edge, 3-face
    enum class ChunkType{ face=0, edgestrip=1, trianglestrip=2 };
    enum class SelectionMode{ vert=1, edge=2, face=3 };

    double default_node_size = 1.0;

    // bool bSelectionSet = true;
    // std::vector<int>        selection; // vector is orderd - usefull for e.g. edge-loops    indices of selected vertices (eventually edges, faces ? )
    // std::unordered_set<int> selset;    // edge index for vert

    // --- from SelectionBanks
    // int icurSelection    = 0;
    // Selection* curSelection = 0;
    // std::vector<Selection> selections;

    //int draw_mode = TRIANGLES;
    Vec3f  penColor;
    

    // ======= Inline Functions

    Builder2(int nSel=10) : SelectionBanks(nSel){}


    inline int   _number_of_points()     { return verts.size(); };
    inline Vec3d _get_point       (int i){ return verts[i].pos; };
    inline int   _number_of_edges ()     { return edges.size(); };
    inline Vec2i _get_edge        (int i){ return edges[i].lo;  };
    inline bool  _add_to_selection(int i){ return curSelection->add(i); };
    // template<typename MESH, typename DistFunc, typename Selection> int _selectRectVerts( MESHr& mesh, DistFunc& distFunc){ 
    int selectRectVerts( Vec3d& p0, double Rmax ){    
        return _selectRectVerts( *this, [&](Vec3d p){ return (p-p0).norm2();},Rmax);
    }

    // int selectRectVerts( Vec3d p0, Vec3d p1, const Mat3d& rot=Mat3dIdentity ){
    //     _order(p0.x,p1.x);
    //     _order(p0.y,p1.y);
    //     return _selectRectVerts( *this, [p0,p1,rot](Vec3d p){  
    //         Vec3d Tp0,Tp1;
    //         rot.dot_to(p0,Tp0);
    //         rot.dot_to(p1,Tp1);
    //         Tp0.z=-1e+300;
    //         Tp1.z=+1e+300;
    //         return (p-p0).norm2();}
    //     );
    // }
        
    //auto edges_range() { return std::views::all(edges); }
    //auto verts_range() { return std::views::all(verts); }

    
    void selectEdgesBySDF(const std::function<double(const Vec3d&)>& sdf, double threshold = 0.0) {
        curSelection->selectByPredicate(
            std::views::all(edges), 
            [&](const Quat4i& edge) -> bool {
                const Vec3d& pA = verts[edge.lo.i].pos;
                const Vec3d& pB = verts[edge.lo.j].pos;
                return (sdf(pA) < threshold) && (sdf(pB) < threshold);
            }
        );
    }

    void selectVertsBySDF(const std::function<double(const Vec3d&)>& sdf, double threshold = 0.0) {
        curSelection->selectByPredicate(std::views::all(verts), [&](const Vert& vert) -> bool {return (sdf(vert.pos) < threshold); });
    }

    

    inline Quat4i latsBlock()const{ return Quat4i{(int)verts.size(),(int)edges.size(),(int)tris.size(),(int)chunks.size()}; }
    inline int block(){ int i=blocks.size(); blocks.push_back( latsBlock() ); return i; };
    inline int vert( const Vec3d& pos, const Vec3d& nor=Vec3dZero, const Vec2d& uv=Vec2dZero ){ 
        //printf( "Mesh::Builder2::vert() %3i pos: %16.10f %16.10f %16.10f \n", verts.size(), pos.x,pos.y,pos.z );
        // _assert( // check vertex min distance
        //     double Rmin=1e-3;
        //     for(int i=0;i<verts.size();i++){
        //         Vec3d p = verts[i].pos;
        //         double r2 = (p-pos).norm2();
        //         if( r2<Rmin*Rmin ){
        //             printf( "Mesh::Builder2::vert() ERROR [%3i] iverts(%3i,%3i) rij(%g)<%g p0(%g,%g,%g) p1(%g,%g,%g) \n", verts.size(), i,verts.size(), sqrt(r2),Rmin, p.x,p.y,p.z, pos.x,pos.y,pos.z );
        //             exit(0);
        //         }
        //     }
        // }
        // Assert that no vertex is found nearby. If one is found (iv>=0), the assertion fails and the action is executed.
        // The action is now safe because it only runs when 'iv' is a valid index.
        _assert( int iv=findVert(pos,RvertCollapse), iv<0, { Vec3d p=verts[iv].pos; printf( "Mesh::Builder2::vert() ERROR [%i] vertex already exists: new_vert_ind=%i is too close to old_vert_ind=%i, rij(%g)<RvertCollapse(%g) p_new(%g,%g,%g) p_old(%g,%g,%g) \n", (int)verts.size(), (int)verts.size(), iv, sqrt((p-pos).norm2()), RvertCollapse, pos.x,pos.y,pos.z, p.x,p.y,p.z ); } );
        verts.push_back(Vert(pos,nor,uv)); return verts.size()-1; 
    }
    inline int edge( int a, int b, int t=-1, int t2=-1 ){ 
        //printf( "Mesh::Builder2::edge() [%3i] (%3i,%3i) t: %i t2: %i \n", edges.size(), a,b,t,t2 );
        // { // check vertex min distance
        //     if(a==b){ printf( "Mesh::Builder2::edge() ERROR [%3i] iverts(%3i,%3i) are the same! \n", edges.size(), a,b ); exit(0); }
        //     double Rmin=1e-3;
        //     Vec3d p0 = verts[a].pos;
        //     Vec3d p1 = verts[b].pos;
        //     double r2 = (p0-p1).norm2();
        //     if( r2<(Rmin*Rmin) ){
        //         printf( "Mesh::Builder2::edge() ERROR [%3i] iverts(%3i,%3i) rij(%g)<%g p0(%g,%g,%g) p1(%g,%g,%g) \n", edges.size(), a,b, sqrt(r2),Rmin, p0.x,p0.y,p0.z, p1.x,p1.y,p1.z );
        //         exit(0);
        //     }
        // }
        //_assert( {},(a<verts.size()) && (b<verts.size()) ,    printf( "Mesh::Builder2::edge() ERROR [%3i] iverts(%3i,%3i) out of bounds 0..vert.size(%i)\n", edges.size(), a,b, verts.size() ) );
        _assert(                         {}          , a!=b , printf( "Mesh::Builder2::edge() ERROR [%3i] iverts(%3i,%3i) are the same\n",  edges.size(), a,b ) );
        _assert( int ie=findEdgeByVerts_brute({a,b}) , ie<0 , printf( "Mesh::Builder2::edge() ERROR [%3i] iverts(%3i,%3i) already exists ie=%i \n", (int)edges.size(), a,b, ie )  );
        edges.push_back(Quat4i{a,b,t2,t}); return edges.size()-1; 
    }

    inline int tri ( int a, int b, int c,    int t=-1  ){ tris .push_back(Quat4i{a,b,c,t});  return tris .size()-1; }


    inline void add_verts(int n, Vec3d* ps, Vec3d* nors=0, Vec2d* uvs=0){
        for(int i=0;i<n;i++){ 
            Vec3d nor = nors ? nors[i] : Vec3dZero;
            Vec2d uv  = uvs  ? uvs[i]  : Vec2dZero;
            vert( ps[i], nor, uv ); 
        } 
    }
    inline void add_edges(int n, Vec2i* es, int* types=0, int* types2=0){ 
        for(int i=0;i<n;i++){ 
            int t = types ? types[i]  : -1;
            int t2= types2? types2[i] : -1;
            edge( es[i].x, es[i].y, t, t2 ); 
        } 
    }

    inline void add_tris(int n, Vec3i* ts, int* types=0){ 
        for(int i=0;i<n;i++){ 
            int t = types ? types[i] : -1;
            tri( ts[i].x, ts[i].y, ts[i].z, t ); 
        } 
    }
    
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
        std::vector<int>& selection = curSelection->vec;
        for(int i=i0;i<iend;i++){ selection.push_back(i); }
        return selection.size();
    }

    inline int* getChunkStrip( int ich ){
        if(ich>=chunks.size()){ printf("ERROR in Builder2::getChunkStrip() ich=%i our of range [0,chunks.size()=%i] \n", ich, chunks.size() ); exit(0); }
        Quat4i ch = chunks[ich];
        return strips.data() + ch.x;
    };


    // ======= Functions



    int selectVertsAlongLine( Vec3d p0, Vec3d p1, double r=0.1, bool bSort=true );
    int selectVertsAlongPolyline( double r, bool bSort, int n, int* edges );
    int selectVertEdgeAngle( int iv, Vec3d hdir, double cosMin=0.9, int ie_ignore=-1 );
    int selectEdgeStrip( int ie, int iv, double cosMin=0.9, int nmax=10, bool bUpdateDir=true );
    int selectEdgeStrip2( int ie, double cosMin=0.9, Vec3i nmaxs=Vec3i{10,-1,-1}, bool bUpdateDir=true );
    
    int select_in_box     ( const Vec3d& p0, const Vec3d& fw, const Vec3d& up, const Vec3d& Lmin, const Vec3d& Lmax );
    
    Vec2i conect_vertex  ( int iv, int stickType, int n, int* iverts );
    int conected_vertex  ( const Vec3d& p, int stickType, int n, int* iverts );
    int make_anchor_point( const Vec3d& p, int stickType, double Rcolapse=0.1, double r=1.0, const Vec3d* fw=0, double l=1.0, int i0=0, int n=0 );
    int make_anchor_points( int nv, Vec3d* vs, int* ivrts, int anchorType, double Rcolapse, double r=1.0, const Vec3d* fw=0, double l=1.0 );

    int extrudeFace( int ich, double L, Quat4i stickTypes=Quat4i{-1,-1,-1,-1}, Quat4i maks={1,1,1,1} );

    // ======= Functions

    Vec3d getCOG(int n, const int* ivs) const;
    void alling_polygons( int n, const int* ivs1, int* ivs2, int ipiv=0 );
    int bridge_quads( Quat4i q1, Quat4i q2, int n, Quat4i stickTypes, Quat4i mask, bool bAlling=true );

    // ======= Mesh Editing & Transformations (from legacy MeshBuilder)
    void move_verts( const std::vector<int>& indices, const Vec3d& shift );
    void scale_verts( const std::vector<int>& indices, const Vec3d& p, const Vec3d& s );
    void rotate_verts( const std::vector<int>& indices, const Vec3d& p, const Mat3d& rot );
    int duplicateBlock( int iblock );

    int extrudeVertLoop( int n, int* iverts, Vec3d d, bool bEdges, bool bFace, bool bTris, bool bSort );
    
    int   loadChunk( int ich, int* iedges=0, int* iverts=0 );
    int   polygonChunk( int n, int* iedges, const int* ivs, bool bPolygonToTris );
    int   polygon( int n, int* iedges );
    int   polygonToTris( int i );
    Vec3d polygonNormal ( int ich ) const;
    //Vec3d getChunkNormal( int ich ) const;
    Vec3d getChunkCOG   ( int ich ) const;

    int   findMostFacingNormal(Vec3d hray, int nch, int* chs, double cosMin=0.0, bool bTwoSide=false, double distWeight=0.0, Vec3d ray0=Vec3dZero )const;
    int   findMostFacingNormal(Vec3d hray, Vec2i chrange,     double cosMin=0.0, bool bTwoSide=false, double distWeight=0.0, Vec3d ray0=Vec3dZero )const;

    void clear();


    bool sortPotentialEdgeLoop( int n, Vec2i* edges, int* iverts );
    bool sortEdgeLoop( int n, int* iedges, int* iverts=0 );
    
    int findEdgeByVerts_brute( Vec2i verts );
    int findEdgeByVerts_map( const Vec2i verts );
    int findEdgeByVerts( const Vec2i verts );
    int findOrAddEdges( const Vec2i verts, int t=-1, int t2=-1 );
    void buildVerts2Edge();
    void build_edgesOfVerts(bool bClear=true);
    int loadNeighbours( int iv, int* ivs, int* ies, int n=-1 );
    Vec3d vertNormalByEdges( int iv, bool bNormalizeEach=false);
    void sortVertEdgesByNormal( Vec3d p, Vec3d nor, int n, int* ies );
    int bevel_vert(int iv, double L, double h, bool bPoly=true, bool bEdgeWedge=false, int* ies=0, Vec3d nor=Vec3dZero );
    int bevel( int ne, int* ies, double L, double h, int nseg=1, uint16_t bevel_flags = BevelFlags::EdgeWedge|BevelFlags::FacesWedge );

    int select_edge_by_verts( int iv, int n, int* ies );
    int select_verts_of_edge( int ne, int* ies, std::vector<int>* ivs=0 );
    void normalsTowardPoint( int nv, int* ivs, Vec3d p, double sc=1.0 );

    int plateBetweenVertStrips( int n, int* ivs1, int* ivs2, int nsub );
    int plateBetweenEdges( int nsub=1, double r=0.1, bool bSort=true );


    Vec2i addVerts( int n, const Vec3d* ps );
    Vec2i addEdges( int n, const Vec2i* iedges, const int* types,  const int* types2, int iv0=0 );
    Vec2i addFaces( int n, const int* nVerts,   const int* iverts, bool bPolygonToTris, int iv0=0 );
    Quat4i addCMesh(const CMesh& cmesh, bool bFaces, Vec3d p0=Vec3dZero, Vec3d sc=Vec3dOne, Mat3d* rot=0, int edge_type=-1 );
    int selectionToFace();

    // ---- Selection
    int clearSelection();
    int toggleSelSet( int i );
    void makeSelectrionUnique();
    // find vertices by distance from point 
    int findClosestVert(const Vec3d& p0,int i0=0,int n=-1);
    int findVert(const Vec3d& p0, double Rmax, int n=-1, int* sel=0 );
    int closestInSelection(const Vec3d& p0, double Rmax, int n, int* sel );
    int select_in_sphere  ( const Vec3d& p0, double r, int i0=0, int imax=-1 );
    int select_in_cylinder( const Vec3d& p0, const Vec3d& fw, double r, double l, int i0=0, int imax=-1 );
    // pick by ray intersection
    int pickVertex( const Vec3d& ray0, const Vec3d& hRay, double R );
    int pickEdge( const Vec3d& ro, const Vec3d& rh, double Rmax );
    int pickTriangle( const Vec3d& ro, const Vec3d& rh, bool bReturnFace=false );
    int pickEdgeSelect( const Vec3d& ro, const Vec3d& rh, double Rmax );
    int pickSelect( const Vec3d& ro, const Vec3d& rh, double Rmax );
    // select by rectangle (box)
    int selectRectEdge( const Vec3d& p0, const Vec3d& p1, const Mat3d& rot=Mat3dIdentity );
    int selectRectVert( const Vec3d& p0, const Vec3d& p1, const Mat3d& rot=Mat3dIdentity );
    int selectRect( const Vec3d& p0, const Vec3d& p1,     const Mat3d& rot=Mat3dIdentity );

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
    int  rope( int ip0,  int ip1,  int typ=-1, int nseg=1 );
    void ropes( int n, int nseg, const Vec2i* ends, int typ );
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
    int export_pos ( Vec3d* ps, int i0=0, int i1=-1 )const;
    int export_pos( float4* ps, int i0=0, int i1=-1 )const;
    int export_edges( Vec2i* eds, int i0=0, int i1=-1 )const;
    int export_tris( Quat4i* tri, int i0=0, int i1=-1 )const;

    // ======= Debug Drawing (from legacy MeshBuilder)
    void addPointCross( const Vec3d& p, double d );
    void addArrow( const Vec3d& p1, const Vec3d& p2, double d );

    void write_obj( const char* fname, uint8_t mask = 0xFF )const;
    void read_obj( const char* fname, uint8_t mask = 0xFF );
    
    void printSelection( bool bDetail=false )const;
    void printSelectedVerts()const;

    void printVert(int iv)const;
    void printEdge(int ie)const;

    void printSizes()const;
    void printVerts()const;
    void printEdges()const;
    void printChunkIndices( int ich )const;
    void printChunkRange( int ich, int ich2=-1 )const;
    void printFaces( bool bNormal=true, bool bCOG=true )const;
    int checkAllPointsConnected(bool bExit=true, bool bPrint=true) const;

    int girder1( Vec3d p0, Vec3d p1, Vec3d up, int n, double width, Quat4i stickTypes, bool bCaps=false );
    int triangle_strip( Vec3d p0, Vec3d p1, Vec3d up, int n, double width, int stickType, bool bCaps=false );
    int plateOnGriders( Vec2i ns, Vec2i prange1, Vec2i prange2, Vec2i byN, Vec2i offs, Vec2d span1, Vec2d span2, Quat4i stickTypes );
    int girder1_caps( int ip0, int ip1, int kind );
    int girder1( int ip0, int ip1, Vec3d up, int n, double width, Quat4i stickTypes );
    int wheel( Vec3d p0, Vec3d p1, Vec3d ax, int n, Vec2d wh, Quat4i stickTypes );
    int ngon( Vec3d p0, Vec3d p1, Vec3d ax, int n,  int stickType );
    //int rope( Vec3d p0, Vec3d p1, int n,  int stickType = -1 );

    int rope ( Vec3d p0,  Vec3d p1, int nseg, int ropeType, int anchorType, double Rcolapse=0.1, double r=-1.0 );
    int ropes( int nv, Vec3d* vs, int ne, int nseg, const Vec2i* ends, int ropeType, int anchorType, double Rcolapse=0.1, double r=-1.0 );
    int panel( Vec3d p00, Vec3d p01, Vec3d p10, Vec3d p11, Vec2i n, double width, Quat4i stickTypes );



    void facingNodes( const CMesh& cmesh, int nnod, const Vec3d* points, Vec2i* out_chs, double* node_sizes=0, const int nplane=0, const int* planes=0, const int* planeVs=0 );

    void bridgeFacingPolygons( Vec3d p0, Vec3d p1, const Vec2i ch1, const Vec2i ch2, int nseg=4, Quat4i stickTypes=Quat4i{-1,-1,-1,-1}, Quat4i maks={1,1,1,1} );
    void bridgeFacingPolygons( int nrod, const Vec2i* edges, const Vec3d* points, int nseg, const Vec2i* chs,  Quat4i stickTypes=Quat4i{-1,-1,-1,-1}, Quat4i maks={1,1,1,1} );

}; // class Mesh::Builder2



}; // namespace Mesh

#endif
