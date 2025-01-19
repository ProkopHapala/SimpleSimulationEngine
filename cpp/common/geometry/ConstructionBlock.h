
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
    int typ;
    int rot;
    int nedge;
    int edges[9];
    Vec3d Lhs;
    BlockFace(): typ(0), rot(0), nedge(0), edges{-1,-1,-1, -1,-1,-1, -1,-1,-1}, Lhs{0.5,0.5,0.5} {};
};

/// @brief class ConstructionBlock use case cube do describe topology and geometry of attached edges into truss  
class ConstructionBlock{ public:
   Vec3d pos=Vec3dZero;
   Vec3d Ls=Vec3dOne;
   BlockFace faces[6]; // 
}; 


class BlockBuilder{ public:
    std::vector<ConstructionBlock> blocks;
    std::vector<Vec2i> edges;

};

//#ifdef MeshBuilder2_h
#include "MeshBuilder2.h"
//namespace Mesh{
namespace Mesh{

    void drawFace( Builder2& mesh, const ConstructionBlock& block, int iface, const Vec3d& p0, const Mat3d& rot, Vec2d Ls ){
        const BlockFace& f = block.faces[iface];
        switch(f.typ){
            case 1:{ // single edge
                //Draw3D::drawLine( p0, p0+rot.c*Ls.x      );
            }break;  
            case 2:{ // fork 2 edges
                //if( f.rot ){ mesh.prismFace   ( p0, Mat3d{rot.b,rot.a,rot.c}, Ls.y, Ls.x, f.Lhs.z, f.Lhs.x); }
                //else       { mesh.prismFace   ( p0, rot,                      Ls.x, Ls.y, f.Lhs.z, f.Lhs.x ); }
                if( f.rot ){ mesh.snapPrismFace   ( p0, Mat3d{rot.b,rot.a,rot.c}, Ls.y, Ls.x, f.Lhs.z, f.Lhs.x); }
                else       { mesh.snapPrismFace   ( p0, rot,                      Ls.x, Ls.y, f.Lhs.z, f.Lhs.x ); }
            }break;
            case 3: {
                if( f.rot ){ mesh.snapFrustrumFace( p0, Mat3d{rot.b,rot.a,rot.c}, Ls.y, Ls.x, f.Lhs.z, f.Lhs.x, f.Lhs.y ); }
                else       { mesh.snapFrustrumFace( p0, rot,                      Ls.x, Ls.y, f.Lhs.z, f.Lhs.x, f.Lhs.y ); }
            }break;// fork 3 edges

            case 5: {}break;// fork 5 edges
        }
    }

    void drawBlock( Builder2& mesh, const ConstructionBlock& block, const Mat3d& rot=Mat3dIdentity ){
        //for(int i=0; i<6; i++){
        const Vec3d& L = block.Ls;
        int i0 = mesh.verts.size();
        mesh.box( block.pos, L, rot );
        mesh.selectVertRange( i0, mesh.verts.size() );
        
        drawFace( mesh, block, 0, block.pos+rot.a* L.a, Mat3d{ rot.b    ,rot.c    ,rot.a    }, {L.y,L.z} );
        drawFace( mesh, block, 1, block.pos+rot.a*-L.a, Mat3d{ rot.b*-1.,rot.c*-1.,rot.a*-1.}, {L.y,L.z} );

        drawFace( mesh, block, 2, block.pos+rot.b* L.b, Mat3d{ rot.c    ,rot.a    ,rot.b    }, {L.z,L.x} );
        drawFace( mesh, block, 3, block.pos+rot.b*-L.b, Mat3d{ rot.c*-1.,rot.a*-1.,rot.b*-1.}, {L.z,L.x} );
        
        drawFace( mesh, block, 4, block.pos+rot.c* L.c, Mat3d{ rot.a    ,rot.b    ,rot.c    }, {L.x,L.y} );
        drawFace( mesh, block, 5, block.pos+rot.c*-L.c, Mat3d{ rot.a*-1.,rot.b*-1.,rot.c*-1.}, {L.x,L.y} );
        mesh.selection.clear();
    }

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

}; // Draw3D
//#endif // Draw3D_h

#endif


