
#ifndef  ConstructionBlock_h
#define  ConstructionBlock_h

/// @brief This file defines the ConstructionBlock class used for building and eddition 3D structures like trusses from platonic solids in a way similar to ZomeTool construction set 


#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"
#include "Solids.h"



inline int cubeIndex3D( uint8_t code ){
    /// encode index of cube face, vertex and edge of 4D cube into 8 bits, 2bits for each of 4 axes xyzw, low bits indicate if axis is non-zero, high bits indicate if axis is positive or negative
    /// bitmap: [x|y|z|w][sx|sy|sz|sw]
    bool x = code & 0b001;
    bool y = code & 0b010;
    bool z = code & 0b100;
    int axBitSum = (int)x + (int)y + (int)z; // sum of non-zero axes
    if      ( axBitSum == 3 ){ // vertex
        return (code >> 4) & 0x111; // top sign bits encode the vertex index
    }else{
        bool sx = code & 0b0010000;
        bool sy = code & 0b0100000;
        bool sz = code & 0b1000000;
        if( axBitSum == 2 ){ // edge
            constexpr int off = 8;  // after 8 vertices
            if     (!x){ return ( (int)sz+(int)sy*2) + (off  ); } // along x axis
            else if(!y){ return ( (int)sx+(int)sy*2) + (off+4); } // along y axis
            else       { return ( (int)sy+(int)sx*2) + (off+8); } // along z axis
        }else if( axBitSum == 1 ){ // face
            constexpr int off = 8+12; // after 8 vertices and 12 edges
            if     (x){ return (int)sx + (off  ); }
            else if(y){ return (int)sy + (off+2); }
            else      { return (int)sz + (off+4); }
        }
    }
    return -1; // invalid code
}

inline int iface( Vec3d p ){
    double xx = p.x*p.x;
    double yy = p.y*p.y;
    double zz = p.z*p.z;
    // which axis is the longest? - triple sort
    //int iax;
    if( xx > yy ){
        if( xx > zz ){ return 0 + (p.x>0); }  // x is longest
        else         { return 4 + (p.z>0); }  // z is longest
    }else{
        if( yy > zz ){ return 2 + (p.y>0); }  // y is longest
        else         { return 4 + (p.z>0); }  // z is longest
    }
}


/// @TODO: Geometry-topological connectivity on Platonic Solids (Convex Polyhedron):   We can take any Convex Polyhedron as a node and take its faces as bond interfaces. 
//         - We only need to make sure that the faces on both ends of n-fold bond are the same (n-gons )
//         - To properly connect the faces if they are not perfectly aligned we need to compute cosine of angle between the princcipal directions for the faces

inline Vec3d getFaceCOG( int n, const int* inds, const Vec3d* points ){
    Vec3d cog=Vec3d{0,0,0};
    for( int i=0; i<n; i++ ){
        int iv  = inds[i];
        cog.add( points[iv] );
    }
    cog.mul(1./n);
    return cog;
}

inline Mat3d getFaceRot( int n, const int* inds, const Vec3d* points, Vec3d* cog_=0 ){
    Vec3d cog;
    if (cog_==0){ cog = getFaceCOG(n, inds, points); } 
    else        { cog = *cog_; }
    Mat3d db=Mat3dZero;
    for( int i=0; i<n; i++ ){
        int iv  = inds[i];
        Vec3d d = points[iv] - cog;
        db.add_outer( d, d );
    }
    Vec3d evals; db.eigenvals(evals);
    Mat3d rot=Mat3dZero;
    db.eigenvec(evals.x, rot.a);
    db.eigenvec(evals.y, rot.b);
    db.eigenvec(evals.z, rot.c);
    return rot;
}

// WTF does this function do?
inline int mapFaces( Vec2i fc1, Vec2i fc2, const int* inds, const Vec3d* points, Vec2i* out_edges ){
    // fc1, fc2 are faces ranges in the index buffer {i0,n} 
    if( fc1.y!=fc2.y ){
        printf("mapFaces: fc1.y(%i)!=fc2.y(%i) => cannot match\n", fc1.y,fc2.y);
        return 0;
    }
    int n = fc1.y;
    Vec3d g1 = getFaceCOG( n, inds+fc1.x, points );
    Vec3d g2 = getFaceCOG( n, inds+fc2.x, points );
    Vec3d ax = (g2-g1).normalized();
    //Mat3d rot; rot.c = g2-g1;
    //rot.c.normalize();
    //rot.c.getSomeOrtho( rot.a, rot.b );
    Vec3d h1=Vec3dZero; // handness vector 1
    Vec3d h2=Vec3dZero; // handness vector 2
    //Vec2d urot=Vec3dZero;
    Vec3d d1s[n];
    Vec3d d2s[n];
    bool invert = h1.dot(h2) < 0;
    for(int i=0; i<n; i++){
        int iv1  = inds[fc1.x+i];
        int iv2  = inds[fc2.x+i];
        Vec3d d1 = (points[iv1] - g1).normalized();
        Vec3d d2 = (points[iv2] - g1).normalized();
        d1s[i] = d1;
        d2s[i] = d2;
        h1.add_cross( d1, ax );
        h2.add_cross( d2, ax );
        //Vec3d urot.add_udiv_cmplx(d1,d2);
    }
    bool bInv = h1.dot(h2) < 0;
    double cbest  = 0;
    int    ishift = 0;
    for(int i=0; i<n; i++){
        double c = 0.0;
        if( bInv   ){ for(int j=0; j<n; j++){ c += d1s[i].dot(d2s[i  ]); } } 
        else        { for(int j=0; j<n; j++){ c += d1s[i].dot(d2s[n-i]); } } // reverse order
        if(c>cbest ){ cbest = c; ishift = i; }
    }
    for(int i=0; i<n; i++){
        int iv1 = inds[fc1.x+i];
        int i2;
        if( bInv   ){ i2 = (i+ishift  )%n; }
        else        { i2 = (i-ishift+n)%n; }
        out_edges[i] = Vec2i{ iv1,  inds[fc2.x+i2] };
    }

}


struct BlockFace{
    static const int nid_max = 9;
    int typ;
    int rot; // is the face rotated by 90 degrees?
    int nid;
    int ids[nid_max];
    Vec3d Lhs;
    BlockFace( int typ_=0, int rot_=0 ): typ(typ_), rot(rot_), nid(0), ids{-1,-1,-1, -1,-1,-1, -1,-1,-1}, Lhs{0.5,0.5,0.5} {};

    void clean( int typ_=0, int rot_=0 ){
        typ=typ_;
        rot=rot_;
        nid=0;
        for(int i=0; i<nid_max; i++){ ids[i] = -1; }
        Lhs = Vec3d{0.5,0.5,0.5};
    }

    bool addId( int id, int i ){
        if( i>=nid_max ) return false;
        if( ids[i]>=0  ) return false;
        ids[i] = id;
        nid++;
        return true;
    }

    int removeIdAt( int i ){
        if( i>=nid_max ) return -1;
        int id = ids[i];
        ids[i] = -1;
        nid--;
        return id;
    }

    int findId( int id ){
        for(int i=0; i<nid_max; i++){
            if( ids[i]==id ){ return i; }
        }
        return -1;
    }

    int removeId( int id ){
        int i = findId(id);
        if( i<0 ) return -1;
        return removeIdAt(i);
    }

    /// TODO: find edge by dir


};

/// @brief class ConstructionBlock use case cube do describe topology and geometry of attached edges into truss  
class ConstructionBlock{ public:
    static const int nfaces = 6;
    Vec3d pos=Vec3dZero;
    Vec3d Ls=Vec3dOne;
    BlockFace faces[nfaces]; 
    // block orientation is 

    ConstructionBlock() = default;
    ConstructionBlock( Vec3d p, Vec3d L , int ftyp=0){ pos=p; Ls=L; cleanFaces(ftyp); };
    ConstructionBlock( Vec3d p, double L, int ftyp=0 ){ pos=p; Ls={L,L,L}; cleanFaces(ftyp); };

    void cleanFaces( int ftyp=0 ){ for(int i=0; i<nfaces; i++){ faces[i].clean(ftyp, 0); } }

    Vec2i findId( int id ){
        for(int i=0; i<nfaces; i++){
            int j = faces[i].findId(id);
            if( j>=0 ){ return Vec2i{i,j}; }
        }
        return Vec2i{-1,-1};
    }

    Vec2i findFace( Vec3d d ){
        // we fine face with minimum cosinus
        double cbest = -1.0;
        Vec2i  ibest = Vec2i{-1,-1}; 
        for(int i=0; i<nfaces; i++){
            Vec3d nr = Solids::Cube_normals[i];
            double c = nr.dot(d);
            //double c = faces[i].bestFace(d); // ToDo: later we need to consider sub-faces inside each face
            if( c>cbest ){
                cbest = c;
                ibest = Vec2i{i,0};
            }
        }
        return ibest;
    }

    Vec2i addId( int id, Vec3d d){
        Vec2i where = findFace(d);
        if( where.i<0 ) return Vec2i{-1,-1};
        if( !faces[where.i].addId(id,where.j) ) return Vec2i{-1,-1};
        return where;
    }

    bool addId( int id, Vec2i where ){
        if( where.i>=nfaces ) return false;
        if( !faces[where.i].addId(id,where.j) ) return false;
        return true;
    }

    // Replaces the logical edge ID in a slot with a new ID (typically a geometric chunk ID)
    // and returns the old logical ID that was stored there.
    int assignChunkToSlot( int new_id, Vec2i where ){
        int oid = faces[where.i].ids[where.j];
        faces[where.i].ids[where.j] = new_id;
        return oid;
    }

    int removeIdAt( Vec2i where ){
        return faces[where.i].removeIdAt(where.j);
    }

    Vec2i removeId( int id ){
        Vec2i where = findId(id);
        if( where.x<0 ) return Vec2i{-1,-1};
        removeIdAt( where );
        return where;
    }

}; 


class BlockBuilder{ public:
    static const int edge_end_offset = 100000000;
    std::vector<ConstructionBlock> blocks;
    std::vector<Quat4i> edges;

    void clear(){
        blocks.clear();
        edges.clear();
    }

    void addBlock( Vec3d p, Vec3d L  ){ 
        //blocks.push_back(ConstructionBlock(p,L)); 
        blocks.emplace_back(p,L);
    }
    int addBlock( Vec3d p, double L, int ftyp=1 ){ 
        //blocks.push_back(ConstructionBlock(p,L)); 
        blocks.emplace_back(p,L,ftyp);
        return blocks.size()-1;
    }
    // can we use emplace ?


    // ToDo: There is serious problem how to identify (index) faces and sub-faces by unique id systematically
    int connectBlocks( int i, int j ){
        Vec3d d = blocks[i].pos - blocks[j].pos; 
        d.normalize();
        int id = edges.size();
        Vec2i f1 = blocks[i].addId( id                 , d*-1.0 );
        Vec2i f2 = blocks[j].addId( id+edge_end_offset , d      );
        printf("connectBlocks() edges.push_back(Quat4i{%i,%i,%i,%i});  id: %i \n", i, j, id, f1.x, id );
        edges.push_back( Quat4i{i,j,f1.x,f2.x}  );
        return id;
    }
   
};

//#ifdef MeshBuilder2_h
#include "MeshBuilder2.h"
//namespace Mesh{
namespace Mesh{

    class ConstructionBlockToMeshBuilder{ public:
        
        Builder2* mesh=0;
        std::unordered_map<int,int> edge2chunk;  // todo: later we perhaps need to use slots
        Quat4i stickTypes{-1,-1,-1,-1};
        Quat4i stickMaks{1,1,1,1};

        void printEdge2chunk(){
            printf("printEdge2chunk() n=%i \n", edge2chunk.size() );
            for(auto& e: edge2chunk){
                printf(" %3i -> %3i \n", e.first, e.second);
            }
        }


    int replace_chunk( ConstructionBlock& block, Vec2i where, int ich=-1 ){
        if(ich<0){ ich = mesh->chunks.size()+ich; };
        int id = block.assignChunkToSlot( ich, where );
        //printf("replace_chunk() id %i ich %i where %i %i \n", id, ich, where.x, where.y );
        if( id>=0 ) edge2chunk[id] = ich;
        return id;
    }

    void drawFace( ConstructionBlock& block, int iface, const Vec3d& p0, const Mat3d& rot, Vec2d Ls, bool bStoreFaceIds=false ){
        Builder2& mesh = *this->mesh;
        const BlockFace& f = block.faces[iface];
        //printf("drawFace: iface=%i p0(%g,%g,%g) f.typ %i f.rot %i \n", iface, p0.x, p0.y, p0.z, f.typ, f.rot );
        switch(f.typ){
            case 1:{ // single edge
                if( f.rot ){ mesh.snapBoxFace     ( p0, Mat3d{rot.b,rot.a,rot.c}, Ls.y, Ls.x ); } // is the face rotated by 90 degrees?
                else       { mesh.snapBoxFace     ( p0, rot,                      Ls.x, Ls.y ); }
                if(bStoreFaceIds){ replace_chunk( block, {iface,0}, -1 ); }
            }break;  
            case 2:{ // fork 2 edges
                if( f.rot ){ mesh.snapPrismFace   ( p0, Mat3d{rot.b,rot.a,rot.c}, Ls.y, Ls.x, f.Lhs.z, f.Lhs.x ); }
                else       { mesh.snapPrismFace   ( p0, rot,                      Ls.x, Ls.y, f.Lhs.z, f.Lhs.x ); }
                if(bStoreFaceIds){ 
                    replace_chunk( block, {iface,1}, -2 ); 
                    replace_chunk( block, {iface,2}, -1 ); 
                    //block.addId( mesh.chunks.size()-1, {iface,0} ); 
                }
            }break;
            case 3: {
                if( f.rot ){ mesh.snapFrustrumFace( p0, Mat3d{rot.b,rot.a,rot.c}, Ls.y, Ls.x, f.Lhs.z, f.Lhs.x, f.Lhs.y ); }
                else       { mesh.snapFrustrumFace( p0, rot,                      Ls.x, Ls.y, f.Lhs.z, f.Lhs.x, f.Lhs.y ); }
                if(bStoreFaceIds){ 
                    replace_chunk( block, {iface,0}, -5 ); // Front face (quad)
                    replace_chunk( block, {iface,1}, -4 ); // top
                    replace_chunk( block, {iface,3}, -3 ); // left
                    replace_chunk( block, {iface,4}, -2 ); // right 
                    replace_chunk( block, {iface,2}, -1 ); // botton
                }
            }break;// fork 3 edges

            case 5: {}break;// fork 5 edges
        }
    }

    void drawBlock( ConstructionBlock& block, const Mat3d& rot=Mat3dIdentity, bool bStoreFaceIds=false ){
        //printf("drawBlockBuilder()\n");
        Builder2& mesh = *this->mesh;
        //for(int i=0; i<6; i++){
        const Vec3d& L = block.Ls;
        int i0 = mesh.verts.size();
        mesh.box( block.pos, L, rot );
        mesh.selectVertRange( i0, mesh.verts.size() );
        //mesh.printSelectedVerts();
        
        // drawFace( block, 0, block.pos+rot.a* L.a, Mat3d{ rot.b    ,rot.c    ,rot.a    }, {L.y,L.z} );
        // drawFace( block, 1, block.pos+rot.a*-L.a, Mat3d{ rot.b*-1.,rot.c*-1.,rot.a*-1.}, {L.y,L.z} );
        // drawFace( block, 2, block.pos+rot.b* L.b, Mat3d{ rot.c    ,rot.a    ,rot.b    }, {L.z,L.x} );
        // drawFace( block, 3, block.pos+rot.b*-L.b, Mat3d{ rot.c*-1.,rot.a*-1.,rot.b*-1.}, {L.z,L.x} );
        // drawFace( block, 4, block.pos+rot.c* L.c, Mat3d{ rot.a    ,rot.b    ,rot.c    }, {L.x,L.y} );
        // drawFace( block, 5, block.pos+rot.c*-L.c, Mat3d{ rot.a*-1.,rot.b*-1.,rot.c*-1.}, {L.x,L.y} );

        Mat3d rot_; 
        for(int i=0; i<6; i++){
            rot_.fromDirUp( Solids::Cube_normals[i], Solids::Cube_ups[i]  );
            drawFace( block, i, block.pos+rot_.c*L.a, rot_, {L.y,L.z}, bStoreFaceIds );
        }

        mesh.selection.clear();
    }

    void drawBlockBuilder( BlockBuilder& skelet, int nseg=4, bool bAllign=true ){
        //printf("drawBlockBuilder() nblock=%i nedge=%i \n", skelet.blocks.size(), skelet.edges.size() );
        Builder2& mesh = *this->mesh;
        edge2chunk.clear();
        for(int i=0; i<skelet.blocks.size(); i++){
            //printf("\n drawBlockBuilder() block[ %i ]\n", i);
            drawBlock( skelet.blocks[i], Mat3dIdentity, true );
        }
        //printEdge2chunk();
        for(int i=0; i<skelet.edges.size(); i++){
            //printf("drawEdge[ %i ]\n", i);
            Quat4i e = skelet.edges[i];
            ConstructionBlock& b1 = skelet.blocks[e.x];
            ConstructionBlock& b2 = skelet.blocks[e.y];
            int ich1   = edge2chunk[ i ];
            int ich2   = edge2chunk[ i + BlockBuilder::edge_end_offset ];
            //Vec2i fe1 = b1.findId( if1 );
            //Vec2i fe2 = b2.findId( if2 );
            Quat4i q1 = *(Quat4i*)mesh.getChunkStrip( ich1 );
            Quat4i q2 = *(Quat4i*)mesh.getChunkStrip( ich2 );
            //printf("drawBlockBuilder() edge[ %i ] e %i %i ich %i %i  q1 %i %i %i %i q2 %i %i %i %i\n", i, e.x, e.y, ich1, ich2, q1.x,q1.y,q1.z,q1.w,   q2.z,q2.w,q2.x,q2.y );
            mesh.bridge_quads( q1, q2, nseg, stickTypes, stickMaks, bAllign );
        }
        //printf("drawBlockBuilder() END\n");
    }

    }; // ConstructionBlockToMeshBuilder

};

//#endif

// Draw Blocks
//#ifdef Draw3D_h

namespace Draw3D{

// void drawPrismFace( const Vec3d& p0, const Mat3d& rot, double w, double h ){
//     Vec3d p1  = p0+rot.a* w;
//     Vec3d p2  = p0+rot.a*-w;
//     Vec3d p1c = p1 + rot.c*h;
//     Vec3d p2c = p2 + rot.c*h;
//     //Draw3D::drawPointCross( p1, 0.1 );
//     //Draw3D::drawPointCross( p2, 0.1 );
//     glBegin(GL_LINES);
//     Draw3D::vertex( p1c        ); Draw3D::vertex( p2c );
//     Draw3D::vertex( p1+rot.b*-w); Draw3D::vertex( p1c );
//     Draw3D::vertex( p1+rot.b* w); Draw3D::vertex( p1c );
//     Draw3D::vertex( p2+rot.b*-w); Draw3D::vertex( p2c );
//     Draw3D::vertex( p2+rot.b* w); Draw3D::vertex( p2c );
//     glEnd();
    
// }

void drawFrustrumFace( const Vec3d& p0, const Mat3d& rot, double La, double Lb, double h, double Lbh, double Lah ){
    //Lbh=0;
    Quat4d ps[4];
    ps[0].w=Lb;     ps[0].f = p0+rot.a* La;
    ps[1].w=Lb-Lbh; ps[1].f = p0+rot.a* (La-Lah) + rot.c*h;
    ps[2].w=Lb-Lbh; ps[2].f = p0+rot.a*-(La-Lah) + rot.c*h;
    ps[3].w=Lb;     ps[3].f = p0+rot.a*-La;
    glBegin(GL_LINES);
    for(int i=0; i<4; i++){
        const Quat4d& q = ps[i];
        Vec3d p1 = q.f + rot.b* q.w;
        Vec3d p2 = q.f + rot.b*-q.w;
        Draw3D::vertex( p1 ); Draw3D::vertex( p2 );
        if(i<3){
            const Quat4d& q2 = ps[i+1];
            Draw3D::vertex( p1 ); Draw3D::vertex( q2.f+rot.b* q2.w );
            Draw3D::vertex( p2 ); Draw3D::vertex( q2.f+rot.b*-q2.w );
        }
    }
    glEnd();
}

void drawPrismFace( const Vec3d& p0, const Mat3d& rot, double La, double Lb, double h, double Lbh ){
    //printf("drawPrismFace: La: %f Lb: %f h: %f Lbh: %f  \n", La,Lb,h,Lbh);
    Quat4d ps[3];
    ps[0].w=Lb;       ps[0].f = p0+rot.a* La;
    ps[1].w=(Lb-Lbh); ps[1].f = p0+rot.c* h ;
    ps[2].w=Lb;       ps[2].f = p0+rot.a*-La;
    //Draw3D::drawPointCross( p1, 0.1 );
    //Draw3D::drawPointCross( p2, 0.1 );
    glBegin(GL_LINES);
    for(int i=0; i<3; i++){
        const Quat4d& q = ps[i];
        Vec3d p1 = q.f + rot.b* q.w;
        Vec3d p2 = q.f + rot.b*-q.w;
        Draw3D::vertex( p1 ); Draw3D::vertex( p2 );
        if(i<2){
            const Quat4d& q2 = ps[i+1];
            Draw3D::vertex( p1 ); Draw3D::vertex( q2.f+rot.b* q2.w );
            Draw3D::vertex( p2 ); Draw3D::vertex( q2.f+rot.b*-q2.w );
        }
    }
    glEnd();
    
}

void drawFace(const ConstructionBlock& block, int iface, const Vec3d& p0, const Mat3d& rot, Vec2d Ls ){
    const BlockFace& f = block.faces[iface];
    switch(f.typ){
        case 1:{ // single edge
            Draw3D::drawLine( p0, p0+rot.c*Ls.x      );
            //Draw3D::drawLine( p0, p0+rot.b*Ls.y*1.0  );
            //Draw3D::drawLine( p0+rot.b*Ls.y*1.0, p0+rot.b*Ls.y*1.0+rot.c*Ls.z*-1.0  );
            //Draw3D::drawLine( p0, p0+rot.a*Ls.x*0.25 );
        }break;  
        case 2:{ // fork 2 edges
            if( f.rot ){ drawPrismFace   ( p0, Mat3d{rot.b,rot.a,rot.c}, Ls.y, Ls.x, f.Lhs.z, f.Lhs.x); }
            else       { drawPrismFace   ( p0, rot,                      Ls.x, Ls.y, f.Lhs.z, f.Lhs.x ); }
        }break;
        case 3: {
            if( f.rot ){ drawFrustrumFace( p0, Mat3d{rot.b,rot.a,rot.c}, Ls.y, Ls.x, f.Lhs.z, f.Lhs.x, f.Lhs.y ); }
            else       { drawFrustrumFace( p0, rot,                      Ls.x, Ls.y, f.Lhs.z, f.Lhs.x, f.Lhs.y ); }
        }break;// fork 3 edges

        case 5: {}break;// fork 5 edges
    }
}

void drawBlock(const ConstructionBlock& block, const Mat3d& rot=Mat3dIdentity ){
    //for(int i=0; i<6; i++){
    const Vec3d& L = block.Ls;
    Draw3D::drawBox( block.pos, L, rot );
    
    drawFace(block, 0, block.pos+rot.a* L.a, Mat3d{ rot.b    ,rot.c    ,rot.a    }, {L.y,L.z} );
    drawFace(block, 1, block.pos+rot.a*-L.a, Mat3d{ rot.b*-1.,rot.c*-1.,rot.a*-1.}, {L.y,L.z} );

    drawFace(block, 2, block.pos+rot.b* L.b, Mat3d{ rot.c    ,rot.a    ,rot.b    }, {L.z,L.x} );
    drawFace(block, 3, block.pos+rot.b*-L.b, Mat3d{ rot.c*-1.,rot.a*-1.,rot.b*-1.}, {L.z,L.x} );
    
    drawFace(block, 4, block.pos+rot.c* L.c, Mat3d{ rot.a    ,rot.b    ,rot.c    }, {L.x,L.y} );
    drawFace(block, 5, block.pos+rot.c*-L.c, Mat3d{ rot.a*-1.,rot.b*-1.,rot.c*-1.}, {L.x,L.y} );
    


    // drawFace(block, 0, block.pos+rot.a* L, Mat3d{ rot.b    ,rot.c*-1,rot.a    }, L );
    // drawFace(block, 1, block.pos+rot.a*-L, Mat3d{ rot.b*-1.,rot.c   ,rot.a*-1.}, L );

    // drawFace(block, 2, block.pos+rot.b* L, Mat3d{ rot.c    ,rot.a*-1,rot.b    }, L );
    // drawFace(block, 3, block.pos+rot.b*-L, Mat3d{ rot.c*-1.,rot.a   ,rot.b*-1.}, L );
    
    // drawFace(block, 4, block.pos+rot.c* L, Mat3d{ rot.a    ,rot.b*-1,rot.c    }, L );
    // drawFace(block, 5, block.pos+rot.c*-L, Mat3d{ rot.a*-1.,rot.b   ,rot.c*-1.}, L );
    //}

}

    void drawBlockBuilder( BlockBuilder& skelet ){
        // for(int i=0; i<skelet.blocks.size(); i++){
        //     drawBlock( skelet.blocks[i], Mat3dIdentity, true );
        // }
        glBegin(GL_LINES);
        for(int i=0; i<skelet.edges.size(); i++){
            Quat4i e = skelet.edges[i];
            ConstructionBlock& b1 = skelet.blocks[e.x];
            ConstructionBlock& b2 = skelet.blocks[e.y];
            //mesh.bridge_quads( q1, q2, nseg, stickTypes, stickMaks );
            Draw3D::vertex( b1.pos ); Draw3D::vertex( b2.pos );
        }
        glEnd();
    }

}; // Draw3D
//#endif // Draw3D_h

#endif
