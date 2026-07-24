/// @file test_Zome.cpp
/// @brief Interactive ZomeTool-like 3D construction editor using icosahedral 62-slot nodes.
///
/// Implements a Zometool-style construction set where nodes are 62-slot connectors
/// based on the 31 symmetry-unique axes of an icosahedron (6 vertex + 10 face + 15 edge).
/// Three connector families map to the three axis classes:
/// - **Red** (5-fold, pentagonal): vertex axes, slots 0–11
/// - **Yellow** (3-fold, triangular): face-normal axes, slots 12–31
/// - **Blue** (2-fold, rectangular): edge-midpoint axes, slots 32–61
///
/// Each family has struts in 3 golden-ratio lengths (1, φ, φ²), giving 9 strut types total.
/// The 60 rotational symmetries of the icosahedral group are generated and used for
/// symmetry-aware grow operations.
///
/// Key interaction features:
/// - **Discrete angular shell grow**: mouse wheel cycles through discrete angular layers
///   around the selected slot's direction. Each layer corresponds to a distinct angle
///   between the selected slot axis and other free slot axes — determined by the
///   icosahedral symmetry, not arbitrary increments. Only struts **on** the shell (at
///   exactly that angle) are previewed and grown, not everything inside the cone.
/// - **Two-node connection**: Shift+click selects a second node; green/red line indicates
///   whether a strut can connect them. T key creates the strut if compatible.
/// - **Picking**: ray-sphere for nodes, ray-point distance for edges and slot arrows.
///   Ray origin pushed behind scene to avoid missing close objects.
///
/// Controls:
///   Left-click:  select node → select free slot → press G to grow
///   Right-click: delete node or edge
///   Mouse wheel: cycle discrete angular shell layers around selected strut direction
///   G: grow strut(s) at current angular shell layer
///   S: jump to layer 1 (first shell)   A: jump to max layer (all slots)
///   Shift+click: select second node for connection check
///   T: create strut between two selected nodes (if compatible)
///   1/2/3: strut length (short/medium/long)
///   C: clear selection   R: reset

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <algorithm>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec2.h"
#include "CMesh.h"
#include "Solids.h"
#include "raytrace.h"

#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"
#include "GLUtils.h"
#include "Zome.h"

inline Vec3d toVec3d( const rigidbuild::Vec3& v ){ return (Vec3d){ v.x, v.y, v.z }; }
inline Vec3f toVec3f( const rigidbuild::Vec3& v ){ return (Vec3f){ (float)v.x, (float)v.y, (float)v.z }; }
inline rigidbuild::Vec3 toRB( const Vec3d& v ){ return { v.x, v.y, v.z }; }

// Draw a prism with n-sided polygonal cross-section between two points.
void drawPrism( const Vec3f& base, const Vec3f& tip, float rx, float ry, int nsides ){
    Vec3f axis = tip - base; axis.normalize();
    Vec3f a, b; axis.getSomeOrtho( a, b ); a.normalize(); b.normalize();
    std::vector<Vec3f> bv(nsides), tv(nsides);
    for( int i=0; i<nsides; i++ ){
        float ang = 2.0f*(float)M_PI*i/nsides;
        Vec3f off = a*(rx*cosf(ang)) + b*(ry*sinf(ang));
        bv[i] = base + off;
        tv[i] = tip  + off;
    }
    glBegin( GL_QUADS );
    for( int i=0; i<nsides; i++ ){
        int j = (i+1)%nsides;
        Vec3f e1 = bv[j]-bv[i], e2 = tv[i]-bv[i];
        Vec3f n; n.set_cross(e1,e2); n.normalize();
        glNormal3f(n.x,n.y,n.z);
        glVertex3f(bv[i].x,bv[i].y,bv[i].z);
        glVertex3f(bv[j].x,bv[j].y,bv[j].z);
        glVertex3f(tv[j].x,tv[j].y,tv[j].z);
        glVertex3f(tv[i].x,tv[i].y,tv[i].z);
    }
    glEnd();
}

using namespace rigidbuild;
constexpr uint32_t NO_ID    = INVALID_ID;
constexpr uint16_t NO_SLOT  = 0xFFFF;
constexpr uint16_t DEL_TYPE = 0xFFFF;

// ======================  TestApp

class TestAppZome : public AppSDL2OGL_3D {
public:
    Construction construction;

    NodeTypeId nodeTypeId = 0;
    ConnectorId connRed=0, connYellow=0, connBlue=0;
    std::vector<StrutTypeId> strutTypeIds;   // 9: 3 families × 3 lengths
    std::vector<double> strutLengths;         // {1.0, φ, φ²}

    /// @brief Selection state — which node/slot the user clicked, plus hover state for highlighting.
    NodeId selectedNode = NO_ID;
    SlotId selectedSlot = NO_SLOT;
    NodeId selectedNode2= NO_ID;  ///< second node for two-node connection check
    NodeId hoveredNode  = NO_ID;
    EdgeId hoveredEdge  = NO_ID;
    SlotId hoveredSlot  = NO_SLOT;
    int selectedLength  = 1;  ///< strut length index: 0=1.0, 1=φ, 2=φ²

    /// @brief Discrete angular shell layers for directional symmetry grow.
    /// Mouse wheel cycles through layers. Layer 0 = only selected slot.
    /// Each layer N corresponds to the N-th distinct angle between the selected
    /// slot direction and other free slot directions — these are discrete values
    /// determined by icosahedral symmetry (e.g. 20.9°, 31.7°, 36.0°, 54.7°, ...).
    /// Only struts **at** that exact angle are selected (shell, not filled cone).
    int symLayer = 0;
    std::vector<double> layerAngles;  ///< sorted unique angular distances from selected slot to all free slots
    double currentAngle = 0.0;        ///< angle (radians) of the current shell layer

    /// @brief Grow preview state — previewed strut endpoints for live visual feedback.
    bool hasPreview = false;
    GrowCandidate preview;                 ///< preview for the selected slot itself
    std::vector<GrowCandidate> previewOrbit;  ///< previews for other slots on the current shell

    /// @brief Two-node connection check state.
    bool canConnect = false;
    GrowCandidate connectCandidate;  ///< valid grow request if canConnect is true

    virtual void draw();
    virtual void camera();
    virtual void eventHandling( const SDL_Event& event );

    void setupTypes();          ///< build 62-slot node type, 9 strut types, 60 symmetry ops
    void buildInitial();        ///< reset construction to single root node
    void renderConstruction();  ///< draw all nodes (spheres) and edges (prisms) with family colors
    void renderSlotArrows( NodeId nid );  ///< draw colored arrows for each free slot on selected node
    void renderPreview();       ///< draw preview struts (blue lines) + cone outline + connection indicator

    void getMouseRay( Vec3d& ro, Vec3d& rd );  ///< convert mouse position to world-space ray (origin pushed behind scene)
    NodeId pickNode ( const Vec3d& ro, const Vec3d& rd );  ///< ray-sphere intersection with node spheres
    EdgeId pickEdge ( const Vec3d& ro, const Vec3d& rd );  ///< ray-point distance to edge segments
    SlotId pickSlot ( NodeId nid, const Vec3d& ro, const Vec3d& rd );  ///< ray-point distance to slot arrows

    void updateHover();         ///< pick node/edge/slot under mouse for highlight
    void updatePreview();       ///< recompute preview struts for current shell layer
    void computeLayers();       ///< compute sorted list of discrete angles from selected slot to all free slots
    void growFromSlot( NodeId nid, SlotId sid );  ///< grow single strut from one slot
    void growCone();            ///< grow all struts on the current angular shell layer
    void checkConnection();    ///< test if any strut type/length can bridge selectedNode → selectedNode2
    void connectTwoNodes();    ///< commit connection between two selected nodes via connectExisting
    void deleteEdge  ( EdgeId eid );
    void removeNode  ( NodeId nid );

    int  getFamilyIndex( ConnectorId c );
    Vec3f familyColor  ( int fam );
    const char* familyName( int fam );

    bool isNodeDeleted( NodeId nid );
    bool isEdgeDeleted( EdgeId eid );
    int  countActiveNodes();
    int  countActiveEdges();

    TestAppZome( int& id, int WIDTH_, int HEIGHT_ );
};

// ======================  Setup

/// @brief Build the 62-slot icosahedral node type, 9 strut types, and 60 symmetry ops.
/// Extracts 6 vertex axes, 10 face-normal axes, 15 edge-midpoint axes from an icosahedron,
/// then generates all 60 rotational symmetries by composing rotations about those axes.
void TestAppZome::setupTypes(){
    // Three connectors: red (5-fold), yellow (3-fold), blue (2-fold)
    connRed    = construction.addConnector( {5, 0, false} );
    connYellow = construction.addConnector( {3, 1, false} );
    connBlue   = construction.addConnector( {2, 2, false} );

    // Icosahedron vertices (normalized) — from Solids.h with a=1, b=φ
    const double phi = GOLDEN_RATIO;
    std::vector<Vec3> verts = {
        {0, -1, -phi}, {0, -1, +phi}, {0, +1, -phi}, {0, +1, +phi},
        {-1, -phi, 0}, {-1, +phi, 0}, {+1, -phi, 0}, {+1, +phi, 0},
        {-phi, 0, -1}, {+phi, 0, -1}, {-phi, 0, +1}, {+phi, 0, +1}
    };
    for( auto& v : verts ) v = normalized(v);

    // --- 1. Vertex axes (6 pairs → 12 slots, 5-fold red) ---
    std::vector<Vec3> vertexAxes;
    std::vector<bool> usedV(12, false);
    for( int i=0; i<12; i++ ){
        if( usedV[i] ) continue;
        for( int j=i+1; j<12; j++ ){
            if( usedV[j] ) continue;
            if( dot(verts[i], verts[j]) < -0.99 ){
                vertexAxes.push_back(verts[i]);
                usedV[i] = usedV[j] = true;
                break;
            }
        }
    }

    // --- 2. Face normals (10 pairs → 20 slots, 3-fold yellow) ---
    int faces[20][3] = {
        {0,8,2}, {0,4,8}, {0,6,4}, {0,9,6}, {0,2,9},
        {3,10,1}, {3,5,10}, {3,7,5}, {3,11,7}, {3,1,11},
        {4,1,10}, {8,10,5}, {2,5,7}, {9,7,11}, {6,11,1},
        {1,4,6}, {10,8,4}, {5,2,8}, {7,9,2}, {11,6,9}
    };
    std::vector<Vec3> faceNorms;
    for( int f=0; f<20; f++ ){
        Vec3 e1 = verts[faces[f][1]] - verts[faces[f][0]];
        Vec3 e2 = verts[faces[f][2]] - verts[faces[f][0]];
        faceNorms.push_back( normalized(cross(e1, e2)) );
    }
    std::vector<Vec3> faceAxes;
    std::vector<bool> usedF(20, false);
    for( int i=0; i<20; i++ ){
        if( usedF[i] ) continue;
        for( int j=i+1; j<20; j++ ){
            if( usedF[j] ) continue;
            if( dot(faceNorms[i], faceNorms[j]) < -0.99 ){
                faceAxes.push_back(faceNorms[i]);
                usedF[i] = usedF[j] = true;
                break;
            }
        }
    }

    // --- 3. Edge midpoint axes (15 pairs → 30 slots, 2-fold blue) ---
    int edges[30][2] = {
        {0,2}, {0,8}, {0,4}, {0,6}, {0,9},
        {3,1}, {3,10}, {3,5}, {3,7}, {3,11},
        {2,8}, {8,4}, {4,6}, {6,9}, {9,2},
        {1,10}, {10,5}, {5,7}, {7,11}, {11,1},
        {1,4}, {10,8}, {5,2}, {7,9}, {11,6},
        {1,6}, {10,4}, {5,8}, {7,2}, {11,9}
    };
    std::vector<Vec3> edgeMids;
    for( int e=0; e<30; e++ ){
        Vec3 mid = verts[edges[e][0]] + verts[edges[e][1]];
        edgeMids.push_back( normalized(mid) );
    }
    std::vector<Vec3> edgeAxes;
    std::vector<bool> usedE(30, false);
    for( int i=0; i<30; i++ ){
        if( usedE[i] ) continue;
        for( int j=i+1; j<30; j++ ){
            if( usedE[j] ) continue;
            if( dot(edgeMids[i], edgeMids[j]) < -0.99 ){
                edgeAxes.push_back(edgeMids[i]);
                usedE[i] = usedE[j] = true;
                break;
            }
        }
    }

    printf( "Icosahedral axes: %zu vertex, %zu face, %zu edge\n",
            vertexAxes.size(), faceAxes.size(), edgeAxes.size() );

    // --- Build 62-slot node type ---
    NodeTypeBuilder builder(0.3);
    Vec3 upHint{0, 0, 1};

    // Red: 6 vertex axes (slots 0-11)
    for( int i=0; i<6; i++ )
        builder.addAxisPair( connRed, 0, i, 0, vertexAxes[i], upHint );
    // Yellow: 10 face axes (slots 12-31)
    for( int i=0; i<10; i++ )
        builder.addAxisPair( connYellow, 1, i, 1, faceAxes[i], upHint );
    // Blue: 15 edge axes (slots 32-61)
    for( int i=0; i<15; i++ )
        builder.addAxisPair( connBlue, 2, i, 2, edgeAxes[i], upHint );

    nodeTypeId = construction.addNodeType( builder.finish() );

    // --- Strut types: 3 lengths × 3 families = 9 ---
    strutLengths = {1.0, phi, phi*phi};
    for( int fam=0; fam<3; fam++ ){
        ConnectorId c = (fam==0) ? connRed : (fam==1) ? connYellow : connBlue;
        for( int len=0; len<3; len++ )
            strutTypeIds.push_back(
                construction.addStrutType( makeStraightStrut( strutLengths[len], c, c ) ) );
    }

    printf( "Built icosahedral node: %zu slots, %zu strut types\n",
            construction.nodeTypes[nodeTypeId].slots.size(), strutTypeIds.size() );

    // --- Generate 60 icosahedral symmetry operations ---
    auto& nt = construction.nodeTypes[nodeTypeId];
    std::vector<Mat3> rots;
    Mat3 id{{1,0,0},{0,1,0},{0,0,1}};
    rots.push_back(id);
    for( auto& ax : vertexAxes ) for( int k=1; k<=4; k++ ) rots.push_back( rotAxis(ax, 2.0*PI*k/5) );
    for( auto& ax : faceAxes   ) for( int k=1; k<=2; k++ ) rots.push_back( rotAxis(ax, 2.0*PI*k/3) );
    for( auto& ax : edgeAxes   ) rots.push_back( rotAxis(ax, PI) );
    // Deduplicate
    std::vector<Mat3> uniq;
    for( auto& r : rots ){
        bool dup=false;
        for( auto& u : uniq ){
            if( fabs(r.c0.x-u.c0.x)<1e-10 && fabs(r.c0.y-u.c0.y)<1e-10 && fabs(r.c0.z-u.c0.z)<1e-10 &&
                fabs(r.c1.x-u.c1.x)<1e-10 && fabs(r.c1.y-u.c1.y)<1e-10 && fabs(r.c1.z-u.c1.z)<1e-10 &&
                fabs(r.c2.x-u.c2.x)<1e-10 && fabs(r.c2.y-u.c2.y)<1e-10 && fabs(r.c2.z-u.c2.z)<1e-10 ){ dup=true; break; }
        }
        if( !dup ) uniq.push_back(r);
    }
    for( auto& r : uniq ){
        try { nt.symmetry.push_back( makeSymmetryOp( nt, construction.connectors, r, 1e-6, 2.0 ) ); }
        catch( ... ) {}
    }
    printf( "Icosahedral symmetry ops: %zu (expected 60)\n", nt.symmetry.size() );
}

void TestAppZome::buildInitial(){
    construction.nodes.clear();
    construction.edges.clear();
    selectedNode = NO_ID; selectedSlot = NO_SLOT;
    hoveredNode  = NO_ID; hoveredEdge  = NO_ID; hoveredSlot = NO_SLOT;
    construction.addNode( nodeTypeId, Pose{} );
    printf( "Zome: reset. Nodes=%d Edges=%d\n", countActiveNodes(), countActiveEdges() );
}

// ======================  Family helpers

int TestAppZome::getFamilyIndex( ConnectorId c ){
    if( c == connRed )    return 0;
    if( c == connYellow ) return 1;
    if( c == connBlue )   return 2;
    return -1;
}

Vec3f TestAppZome::familyColor( int fam ){
    if( fam == 0 ) return {0.9f, 0.2f, 0.2f};  // red
    if( fam == 1 ) return {0.9f, 0.8f, 0.2f};  // yellow
    return {0.2f, 0.4f, 0.9f};                  // blue
}

const char* TestAppZome::familyName( int fam ){
    if( fam == 0 ) return "Red";
    if( fam == 1 ) return "Yellow";
    return "Blue";
}

// ======================  Helpers

bool TestAppZome::isNodeDeleted( NodeId nid ){ return construction.nodes[nid].type == DEL_TYPE; }
bool TestAppZome::isEdgeDeleted( EdgeId eid ){ return construction.edges[eid].node[0] == NO_ID; }

int TestAppZome::countActiveNodes(){
    int n=0; for( NodeId i=0; i<construction.nodes.size(); i++ ) if(!isNodeDeleted(i)) n++; return n;
}
int TestAppZome::countActiveEdges(){
    int n=0; for( EdgeId i=0; i<construction.edges.size(); i++ ) if(!isEdgeDeleted(i)) n++; return n;
}

// ======================  Picking

void TestAppZome::getMouseRay( Vec3d& ro, Vec3d& rd ){
    ro = (Vec3d)ray0 - (Vec3d)cam.rot.c * 100.0;
    rd = (Vec3d)cam.rot.c;
    rd.normalize();
}

NodeId TestAppZome::pickNode( const Vec3d& ro, const Vec3d& rd ){
    NodeId best = NO_ID; double bestT = 1e30;
    for( NodeId i=0; i<construction.nodes.size(); i++ ){
        if( isNodeDeleted(i) ) continue;
        Vec3d pos = toVec3d( construction.nodes[i].world.p );
        double t = raySphere( ro, rd, 0.25, pos );
        if( t>0 && t<bestT ){ bestT=t; best=i; }
    }
    return best;
}

EdgeId TestAppZome::pickEdge( const Vec3d& ro, const Vec3d& rd ){
    EdgeId best = NO_ID; double bestD2 = 0.04;
    for( EdgeId i=0; i<construction.edges.size(); i++ ){
        if( isEdgeDeleted(i) ) continue;
        auto& e = construction.edges[i];
        Vec3d p0 = toVec3d( construction.nodes[e.node[0]].world.p );
        Vec3d p1 = toVec3d( construction.nodes[e.node[1]].world.p );
        for( int j=0; j<=8; j++ ){
            double s = j/8.0;
            Vec3d p = p0*(1-s) + p1*s;
            double t; double d2 = rayPointDistance2( ro, rd, p, t );
            if( d2 < bestD2 ){ bestD2 = d2; best = i; break; }
        }
    }
    return best;
}

SlotId TestAppZome::pickSlot( NodeId nid, const Vec3d& ro, const Vec3d& rd ){
    if( nid==NO_ID || isNodeDeleted(nid) ) return NO_SLOT;
    auto& node = construction.nodes[nid];
    auto& nt   = construction.nodeTypes[node.type];
    SlotId best = NO_SLOT; double bestD2 = 0.025; // tighter for 62 slots
    for( SlotId s=0; s<nt.slots.size(); s++ ){
        if( !construction.slotFree(nid, s) ) continue;
        Pose worldSlot = node.world * nt.slots[s].local;
        Vec3d arrowEnd = toVec3d( worldSlot.p ) + toVec3d( worldSlot.R.c2 ) * 0.4;
        Vec3d p0 = toVec3d( worldSlot.p );
        for( int j=0; j<=6; j++ ){
            double t = j/6.0;
            Vec3d p = p0*(1-t) + arrowEnd*t;
            double tRay; double d2 = rayPointDistance2( ro, rd, p, tRay );
            if( d2 < bestD2 ){ bestD2 = d2; best = s; break; }
        }
    }
    return best;
}

void TestAppZome::updateHover(){
    Vec3d ro, rd; getMouseRay( ro, rd );
    hoveredNode = pickNode( ro, rd );
    hoveredEdge = (hoveredNode==NO_ID) ? pickEdge( ro, rd ) : NO_ID;
    hoveredSlot = NO_SLOT;
    if( selectedNode != NO_ID && hoveredNode == NO_ID && hoveredEdge == NO_ID )
        hoveredSlot = pickSlot( selectedNode, ro, rd );
}

/// @brief Compute sorted list of discrete angles from selected slot direction to all other free slots.
/// These are the only angles at which struts can appear relative to the selected strut —
/// they are determined by icosahedral symmetry, not arbitrary increments.
/// For a fresh node this yields ~18 distinct angles: 20.9°, 31.7°, 36.0°, 54.7°, ..., 180°.
void TestAppZome::computeLayers(){
    layerAngles.clear();
    if( selectedNode == NO_ID || selectedSlot == NO_SLOT ) return;
    if( isNodeDeleted(selectedNode) || !construction.slotFree(selectedNode, selectedSlot) ) return;
    auto& nt = construction.nodeTypes[nodeTypeId];
    auto& node = construction.nodes[selectedNode];
    Pose wsSel = node.world * nt.slots[selectedSlot].local;
    rigidbuild::Vec3 dirSel = normalized(wsSel.R.c2);
    for( SlotId s=0; s<nt.slots.size(); s++ ){
        if( s == selectedSlot ) continue;
        if( !construction.slotFree(selectedNode, s) ) continue;
        Pose ws = node.world * nt.slots[s].local;
        rigidbuild::Vec3 dir = normalized(ws.R.c2);
        double ang = acos( std::clamp(dot(dirSel, dir), -1.0, 1.0) );
        // Check if this angle is already in the list (within tolerance)
        bool found = false;
        for( double a : layerAngles ) if( fabs(a - ang) < 1e-4 ) { found = true; break; }
        if( !found ) layerAngles.push_back( ang );
    }
    std::sort( layerAngles.begin(), layerAngles.end() );
    printf( "Discrete angles: %zu layers:", layerAngles.size() );
    for( double a : layerAngles ) printf( " %.1f°", a*180/M_PI );
    printf( "\n" );
}

/// @brief Recompute preview struts for the current angular shell layer.
/// Layer 0: preview only the selected slot. Layer N>0: preview all free slots
/// whose angle from the selected slot direction equals layerAngles[N-1] (within tolerance).
/// This selects struts **on** the shell, not inside the cone.
void TestAppZome::updatePreview(){
    hasPreview = false;
    previewOrbit.clear();
    if( selectedNode == NO_ID || selectedSlot == NO_SLOT ) return;
    if( isNodeDeleted(selectedNode) || !construction.slotFree(selectedNode, selectedSlot) ) return;
    auto& nt = construction.nodeTypes[nodeTypeId];
    auto& node = construction.nodes[selectedNode];
    Pose wsSel = node.world * nt.slots[selectedSlot].local;
    rigidbuild::Vec3 dirSel = normalized(wsSel.R.c2);
    ConnectorId c = nt.slots[selectedSlot].connector;
    int fam = getFamilyIndex(c);
    if( fam < 0 ) return;
    // Determine current angle threshold
    currentAngle = 0.0;
    if( symLayer > 0 && symLayer <= (int)layerAngles.size() )
        currentAngle = layerAngles[symLayer - 1];
    // Always include the selected slot itself
    GrowRequest q0;
    q0.sourceNode     = selectedNode;
    q0.sourceSlot     = selectedSlot;
    q0.strutType      = strutTypeIds[fam * 3 + selectedLength];
    q0.targetNodeType = nodeTypeId;
    q0.targetSlot     = nt.slots[selectedSlot].opposite;
    try { preview = construction.previewGrow( q0 ); hasPreview = true; }
    catch( const std::exception& e ) { printf( "Preview failed: %s\n", e.what() ); }
    // If symLayer > 0, also preview slots within currentAngle
    if( symLayer > 0 ){
        for( SlotId s=0; s<nt.slots.size(); s++ ){
            if( s == selectedSlot ) continue;
            if( !construction.slotFree(selectedNode, s) ) continue;
            Pose ws = node.world * nt.slots[s].local;
            rigidbuild::Vec3 dir = normalized(ws.R.c2);
            double ang = acos( std::clamp(dot(dirSel, dir), -1.0, 1.0) );
            if( fabs(ang - currentAngle) > 1e-4 ) continue;
            ConnectorId sc = nt.slots[s].connector;
            int sfam = getFamilyIndex(sc);
            if( sfam < 0 ) continue;
            GrowRequest q;
            q.sourceNode     = selectedNode;
            q.sourceSlot     = s;
            q.strutType      = strutTypeIds[sfam * 3 + selectedLength];
            q.targetNodeType = nodeTypeId;
            q.targetSlot     = nt.slots[s].opposite;
            try { previewOrbit.push_back( construction.previewGrow( q ) ); }
            catch( ... ) {}
        }
        printf( "Layer %d (%.1f°): %zu struts in preview\n", symLayer, currentAngle*180/M_PI, previewOrbit.size()+1 );
    }
}

/// @brief Check if any strut type/length can bridge selectedNode to selectedNode2.
/// Tries all 9 strut types from all free slots on selectedNode, checking if the
/// resulting target node position matches selectedNode2's position (within 0.01 tolerance).
void TestAppZome::checkConnection(){
    canConnect = false;
    if( selectedNode==NO_ID || selectedNode2==NO_ID ) return;
    if( isNodeDeleted(selectedNode) || isNodeDeleted(selectedNode2) ) return;
    auto& nt = construction.nodeTypes[nodeTypeId];
    for( SlotId s1=0; s1<nt.slots.size() && !canConnect; s1++ ){
        if( !construction.slotFree(selectedNode, s1) ) continue;
        ConnectorId c1 = nt.slots[s1].connector;
        int fam1 = getFamilyIndex(c1);
        if( fam1 < 0 ) continue;
        for( int len=0; len<3 && !canConnect; len++ ){
            GrowRequest q;
            q.sourceNode     = selectedNode;
            q.sourceSlot     = s1;
            q.strutType      = strutTypeIds[fam1*3 + len];
            q.targetNodeType = nodeTypeId;
            q.targetSlot     = nt.slots[s1].opposite;
            try {
                GrowCandidate gc = construction.previewGrow( q );
                Vec3d p2 = toVec3d( construction.nodes[selectedNode2].world.p );
                Vec3d dp = toVec3d( gc.targetNodeWorld.p ) - p2;
                if( dp.norm() < 0.01 ){
                    canConnect = true;
                    connectCandidate = gc;
                }
            } catch( ... ) {}
        }
    }
}

/// @brief Commit a strut between two selected nodes using connectExisting (no duplicate node).
/// Finds the matching slot on selectedNode2 and calls construction.connectExisting().
void TestAppZome::connectTwoNodes(){
    if( !canConnect ) return;
    auto& nt = construction.nodeTypes[nodeTypeId];
    // Find the target slot on node2 that matches
    Vec3d targetPos = toVec3d( connectCandidate.targetNodeWorld.p );
    for( SlotId s2=0; s2<nt.slots.size(); s2++ ){
        if( !construction.slotFree(selectedNode2, s2) ) continue;
        Pose ws2 = construction.nodes[selectedNode2].world * nt.slots[s2].local;
        if( (toVec3d(ws2.p) - targetPos).norm() < 0.01 ){
            try {
                construction.connectExisting( connectCandidate.request, selectedNode2 );
                printf( "Connected nodes %u-%u. Edges=%d\n", selectedNode, selectedNode2, countActiveEdges() );
                selectedNode2 = NO_ID;
                canConnect = false;
            } catch( const std::exception& e ){
                printf( "Connect failed: %s\n", e.what() );
            }
            return;
        }
    }
    printf( "Connect: no matching slot found on node2\n" );
}

// ======================  Editing

void TestAppZome::growFromSlot( NodeId nid, SlotId sid ){
    if( nid==NO_ID || sid==NO_SLOT || isNodeDeleted(nid) ) return;
    if( !construction.slotFree(nid, sid) ){ printf("Slot occupied\n"); return; }
    auto& nt = construction.nodeTypes[nodeTypeId];
    ConnectorId c = nt.slots[sid].connector;
    int fam = getFamilyIndex(c);
    if( fam < 0 ){ printf("Unknown connector\n"); return; }
    GrowRequest q;
    q.sourceNode     = nid;
    q.sourceSlot     = sid;
    q.strutType      = strutTypeIds[fam * 3 + selectedLength];
    q.targetNodeType = nodeTypeId;
    q.targetSlot     = nt.slots[sid].opposite;
    try {
        auto [newNode, newEdge] = construction.grow( q );
        printf( "Grew: %s L=%.3f  node=%u edge=%u  Nodes=%d Edges=%d\n",
                familyName(fam), strutLengths[selectedLength],
                newNode, newEdge, countActiveNodes(), countActiveEdges() );
        selectedNode = newNode; selectedSlot = NO_SLOT;
        updatePreview();
    } catch( const std::exception& e ) {
        printf( "Grow failed: %s\n", e.what() );
    }
}

/// @brief Grow all struts on the current angular shell layer.
/// Iterates all free slots, grows those whose angle from the selected slot direction
/// matches currentAngle (within tolerance). Resets symLayer and selectedSlot after growing.
void TestAppZome::growCone(){
    if( selectedNode==NO_ID || selectedSlot==NO_SLOT ) return;
    if( isNodeDeleted(selectedNode) ) return;
    auto& nt = construction.nodeTypes[nodeTypeId];
    auto& node = construction.nodes[selectedNode];
    Pose wsSel = node.world * nt.slots[selectedSlot].local;
    rigidbuild::Vec3 dirSel = normalized(wsSel.R.c2);
    int count = 0;
    for( SlotId s=0; s<nt.slots.size(); s++ ){
        if( !construction.slotFree(selectedNode, s) ) continue;
        Pose ws = node.world * nt.slots[s].local;
        rigidbuild::Vec3 dir = normalized(ws.R.c2);
        double ang = acos( std::clamp(dot(dirSel, dir), -1.0, 1.0) );
        if( fabs(ang - currentAngle) > 1e-4 ) continue;
        ConnectorId c = nt.slots[s].connector;
        int sfam = getFamilyIndex(c);
        if( sfam < 0 ) continue;
        GrowRequest q;
        q.sourceNode     = selectedNode;
        q.sourceSlot     = s;
        q.strutType      = strutTypeIds[sfam * 3 + selectedLength];
        q.targetNodeType = nodeTypeId;
        q.targetSlot     = nt.slots[s].opposite;
        try { construction.grow( q ); count++; }
        catch( ... ) {}
    }
    printf( "Grew %d struts (layer %d, %.1f°). Nodes=%d Edges=%d\n",
            count, symLayer, currentAngle*180/M_PI, countActiveNodes(), countActiveEdges() );
    selectedSlot = NO_SLOT;
    symLayer = 0; currentAngle = 0;
    updatePreview();
}

void TestAppZome::deleteEdge( EdgeId eid ){
    if( eid==NO_ID || isEdgeDeleted(eid) ) return;
    auto& e = construction.edges[eid];
    construction.nodes[e.node[0]].edgeAtSlot[e.slot[0]] = NO_ID;
    construction.nodes[e.node[1]].edgeAtSlot[e.slot[1]] = NO_ID;
    e.node[0] = NO_ID; e.node[1] = NO_ID;
    printf( "Deleted edge %u. Nodes=%d Edges=%d\n", eid, countActiveNodes(), countActiveEdges() );
}

void TestAppZome::removeNode( NodeId nid ){
    if( nid==NO_ID || isNodeDeleted(nid) ) return;
    auto& node = construction.nodes[nid];
    for( SlotId s=0; s<node.edgeAtSlot.size(); s++ ){
        EdgeId eid = node.edgeAtSlot[s];
        if( eid==NO_ID ) continue;
        auto& e = construction.edges[eid];
        NodeId other = (e.node[0]==nid) ? e.node[1] : e.node[0];
        SlotId oslot = (e.node[0]==nid) ? e.slot[1] : e.slot[0];
        construction.nodes[other].edgeAtSlot[oslot] = NO_ID;
        e.node[0] = NO_ID; e.node[1] = NO_ID;
    }
    node.type = DEL_TYPE;
    if( selectedNode==nid ){ selectedNode=NO_ID; selectedSlot=NO_SLOT; }
    printf( "Deleted node %u. Nodes=%d Edges=%d\n", nid, countActiveNodes(), countActiveEdges() );
}

// ======================  Rendering

/// @brief Draw all nodes as spheres and edges as family-colored polygonal prisms.
/// Cross-section shape encodes family: pentagon (red/5-fold), triangle (yellow/3-fold),
/// rectangle (blue/2-fold). Selected/hovered nodes and edges are highlighted.
void TestAppZome::renderConstruction(){
    // Draw edges (struts) as polygonal prisms, cross-section by family
    for( EdgeId i=0; i<construction.edges.size(); i++ ){
        if( isEdgeDeleted(i) ) continue;
        auto& e = construction.edges[i];
        Vec3f p0 = toVec3f( construction.nodes[e.node[0]].world.p );
        Vec3f p1 = toVec3f( construction.nodes[e.node[1]].world.p );
        Vec3f d = p1-p0; d.normalize();
        Vec3f a = p0 + d*0.15f, b = p1 - d*0.15f;
        int fam = getFamilyIndex( construction.nodeTypes[nodeTypeId].slots[e.slot[0]].connector );
        Vec3f col = familyColor(fam);
        if( i == hoveredEdge ) col = {1.0f, 0.3f, 0.3f};
        glColor3f( col.x, col.y, col.z );
        // Cross-section: pentagon(5) red, triangle(3) yellow, rectangle(4) blue
        if( fam==0 )      drawPrism( a, b, 0.06f, 0.06f, 5 );
        else if( fam==1 ) drawPrism( a, b, 0.06f, 0.06f, 3 );
        else              drawPrism( a, b, 0.08f, 0.04f, 4 );
    }

    // Draw nodes as spheres
    for( NodeId i=0; i<construction.nodes.size(); i++ ){
        if( isNodeDeleted(i) ) continue;
        Vec3f pos = toVec3f( construction.nodes[i].world.p );
        if      ( i==selectedNode  ) glColor3f( 0.9f, 0.7f, 0.2f );
        else if ( i==selectedNode2 ) glColor3f( 0.2f, 1.0f, 0.2f );
        else if ( i==hoveredNode   ) glColor3f( 0.3f, 0.9f, 0.3f );
        else                        glColor3f( 0.8f, 0.8f, 0.8f );
        Draw3D::drawSphere_oct( 3, 0.15f, pos );
    }

    // Draw slot arrows on selected node
    if( selectedNode != NO_ID && !isNodeDeleted(selectedNode) )
        renderSlotArrows( selectedNode );
}

/// @brief Draw colored arrows for each slot on the selected node.
/// Free slots show as family-colored arrows; occupied slots as dark crosses.
/// Selected slot is cyan, hovered slot is yellow, others dimmed for clarity.
void TestAppZome::renderSlotArrows( NodeId nid ){
    auto& node = construction.nodes[nid];
    auto& nt   = construction.nodeTypes[node.type];
    glDisable( GL_LIGHTING );
    for( SlotId s=0; s<nt.slots.size(); s++ ){
        Pose ws = node.world * nt.slots[s].local;
        Vec3f pos = toVec3f( ws.p );
        Vec3f dir = toVec3f( ws.R.c2 ) * 0.35f;
        if( !construction.slotFree(nid, s) ){
            // Occupied: small dark dot
            glColor3f( 0.3f, 0.3f, 0.3f );
            Draw3D::drawPointCross( pos, 0.03f );
        } else {
            int fam = getFamilyIndex( nt.slots[s].connector );
            Vec3f col = familyColor(fam);
            if( s == selectedSlot )      col = {0.2f, 0.9f, 1.0f};  // cyan
            else if( s == hoveredSlot )  col = {1.0f, 1.0f, 0.0f }; // yellow hover
            // Dim non-hovered, non-selected arrows for clarity
            if( s != selectedSlot && s != hoveredSlot ) col = col * 0.6f;
            glColor3f( col.x, col.y, col.z );
            Draw3D::drawArrow( pos, pos + dir, 0.02f );
        }
    }
    glEnable( GL_LIGHTING );
}

/// @brief Draw preview struts, cone outline, and two-node connection indicator.
/// Layer 0: single blue line + sphere for the selected slot preview.
/// Layer N>0: green line for selected strut, blue lines for shell struts,
/// cyan circle outline at currentAngle around the selected slot direction.
/// Two-node mode: green line if connectable, red if not.
void TestAppZome::renderPreview(){
    glDisable( GL_LIGHTING );
    if( hasPreview && symLayer == 0 ){
        Vec3f p0 = toVec3f( construction.nodes[selectedNode].world.p );
        Vec3f p1 = toVec3f( preview.targetNodeWorld.p );
        glColor3f( 0.5f, 0.8f, 1.0f );
        Draw3D::drawLine( p0, p1 );
        glColor3f( 0.4f, 0.6f, 0.8f );
        Draw3D::drawSphere_oct( 2, 0.12f, p1 );
    }
    if( hasPreview && symLayer > 0 ){
        Vec3f p0 = toVec3f( construction.nodes[selectedNode].world.p );
        // Draw selected strut brighter
        if( hasPreview ){
            Vec3f p1 = toVec3f( preview.targetNodeWorld.p );
            glColor3f( 0.8f, 1.0f, 0.8f );
            Draw3D::drawLine( p0, p1 );
            Draw3D::drawSphere_oct( 2, 0.10f, p1 );
        }
        glColor3f( 0.5f, 0.8f, 1.0f );
        for( auto& gc : previewOrbit ){
            Vec3f p1 = toVec3f( gc.targetNodeWorld.p );
            Draw3D::drawLine( p0, p1 );
            Draw3D::drawSphere_oct( 2, 0.08f, p1 );
        }
        // Draw cone outline at currentAngle
        auto& node = construction.nodes[selectedNode];
        auto& nt = construction.nodeTypes[nodeTypeId];
        Pose wsSel = node.world * nt.slots[selectedSlot].local;
        Vec3f axis = toVec3f( wsSel.R.c2 ); axis.normalize();
        Vec3f a, b; axis.getSomeOrtho( a, b ); a.normalize(); b.normalize();
        float coneLen = 1.5f;
        float coneR = coneLen * tanf( (float)currentAngle );
        glColor3f( 0.3f, 0.5f, 0.7f );
        glBegin( GL_LINE_LOOP );
        for( int i=0; i<32; i++ ){
            float ang = 2.0f*(float)M_PI*i/32;
            Vec3f p = toVec3f( wsSel.p ) + axis*coneLen + a*(coneR*cosf(ang)) + b*(coneR*sinf(ang));
            glVertex3f( p.x, p.y, p.z );
        }
        glEnd();
    }
    // Two-node connection indicator
    if( selectedNode2 != NO_ID && !isNodeDeleted(selectedNode2) ){
        Vec3f p0 = toVec3f( construction.nodes[selectedNode].world.p );
        Vec3f p1 = toVec3f( construction.nodes[selectedNode2].world.p );
        if( canConnect ) glColor3f( 0.2f, 1.0f, 0.2f );  // green
        else             glColor3f( 1.0f, 0.2f, 0.2f );  // red
        Draw3D::drawLine( p0, p1 );
    }
    glEnable( GL_LIGHTING );
}

// ======================  Event handling

void TestAppZome::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_KEYDOWN:
            switch( event.key.keysym.sym ){
                case SDLK_g:
                    if( symLayer == 0 ) growFromSlot( selectedNode, selectedSlot );
                    else                growCone();
                    break;
                case SDLK_s: symLayer = 1; updatePreview(); printf("Shell layer 1\n"); break;
                case SDLK_a: symLayer = (int)layerAngles.size(); updatePreview(); printf("Shell layer %d (all)\n", symLayer); break;
                case SDLK_t: connectTwoNodes(); break;
                case SDLK_c: selectedNode=NO_ID; selectedSlot=NO_SLOT; selectedNode2=NO_ID; canConnect=false; updatePreview(); break;
                case SDLK_r: buildInitial(); break;
                case SDLK_1: selectedLength=0; updatePreview(); printf("Length=%.3f\n", strutLengths[0]); break;
                case SDLK_2: selectedLength=1; updatePreview(); printf("Length=%.3f\n", strutLengths[1]); break;
                case SDLK_3: selectedLength=2; updatePreview(); printf("Length=%.3f\n", strutLengths[2]); break;
            }
            break;
        case SDL_MOUSEWHEEL:
            if( selectedSlot != NO_SLOT && !layerAngles.empty() ){
                int maxLayer = (int)layerAngles.size();
                if( event.wheel.y > 0 ) symLayer = std::min( symLayer + 1, maxLayer );
                else                    symLayer = std::max( symLayer - 1, 0 );
                updatePreview();
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT: {
                    Vec3d ro, rd; getMouseRay( ro, rd );
                    // Shift+click: select second node for connection
                    if( SDL_GetModState() & KMOD_SHIFT ){
                        NodeId n = pickNode( ro, rd );
                        if( n != NO_ID && n != selectedNode ){
                            selectedNode2 = n;
                            checkConnection();
                            printf( "Selected node2 %u (connectable: %s)\n", n, canConnect?"YES":"NO" );
                        } else {
                            selectedNode2 = NO_ID; canConnect = false;
                        }
                        break;
                    }
                    if( selectedNode != NO_ID ){
                        SlotId s = pickSlot( selectedNode, ro, rd );
                        if( s != NO_SLOT ){
                            selectedSlot = s;
                            symLayer = 0;
                            computeLayers();
                            updatePreview();
                            auto& nt = construction.nodeTypes[nodeTypeId];
                            int fam = getFamilyIndex( nt.slots[s].connector );
                            printf( "Selected slot %u (%s) on node %u  layer=%d\n", s, familyName(fam), selectedNode, symLayer );
                            break;
                        }
                    }
                    NodeId n = pickNode( ro, rd );
                    if( n != NO_ID ){
                        selectedNode = n; selectedSlot = NO_SLOT; selectedNode2 = NO_ID; canConnect = false;
                        updatePreview();
                        printf( "Selected node %u\n", n );
                    } else {
                        selectedNode = NO_ID; selectedSlot = NO_SLOT; selectedNode2 = NO_ID; canConnect = false;
                        updatePreview();
                    }
                } break;
                case SDL_BUTTON_RIGHT: {
                    Vec3d ro, rd; getMouseRay( ro, rd );
                    NodeId n = pickNode( ro, rd );
                    if( n != NO_ID && n != selectedNode ){
                        removeNode( n );
                    } else {
                        EdgeId e = pickEdge( ro, rd );
                        if( e != NO_ID ) deleteEdge( e );
                    }
                    updatePreview();
                } break;
            }
            break;
    }
    AppSDL2OGL_3D::eventHandling( event );
}

// ======================  Camera

void TestAppZome::camera(){
    ((Quat4f)qCamera).toMatrix(cam.rot);
    Vec3f target = {0.0f, 0.0f, 0.0f};
    if( selectedNode != NO_ID && !isNodeDeleted(selectedNode) )
        target = toVec3f( construction.nodes[selectedNode].world.p );
    cam.pos = target;
    cam.zoom   = zoom;
    cam.aspect = ASPECT_RATIO;
    Cam::ortho( cam, true );
}

// ======================  Draw

void TestAppZome::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    GLfloat lightPos[] = { -cam.rot.c.x, -cam.rot.c.y, -cam.rot.c.z, 0.0 };
    glLightfv( GL_LIGHT0, GL_POSITION, lightPos );

    glEnable( GL_LIGHTING );
    glShadeModel( GL_SMOOTH );

    updateHover();
    renderConstruction();
    renderPreview();

    glDisable( GL_LIGHTING );
    Draw3D::drawAxis( 3.0f );
}

// ======================  Constructor

TestAppZome::TestAppZome( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ){
    setupTypes();
    buildInitial();
    cam.pos  = (Vec3f){0.0f, 0.0f, 0.0f};
    cam.zmin = -1000.0f;
    cam.zmax =  1000.0f;
    zoom = 3.0f;
    printf( "\n=== ZomeTool Interactive Demo (62-slot icosahedral) ===\n" );
    printf( "Left-click: select node, then select free slot\n" );
    printf( "Shift+click: select second node for connection check\n" );
    printf( "Mouse wheel: cycle discrete angular shell layers around selected strut\n" );
    printf( "G: grow strut(s) at current shell layer\n" );
    printf( "S: jump to layer 1 (first shell)  A: jump to max layer (all)\n" );
    printf( "T: connect two selected nodes (if compatible)\n" );
    printf( "Right-click: delete node or edge\n" );
    printf( "1/2/3: strut length (%.3f / %.3f / %.3f)\n", strutLengths[0], strutLengths[1], strutLengths[2] );
    printf( "C: clear  R: reset\n" );
    printf( "Colors: Red=pentagon(5-fold) Yellow=triangle(3-fold) Blue=rectangle(2-fold)\n" );
    printf( "========================================================\n\n" );
}

// ===================== MAIN

TestAppZome * testApp;

int main( int argc, char *argv[] ){
    SDL_Init( SDL_INIT_VIDEO );
    SDL_GL_SetAttribute( SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1 );
    int junk;
    testApp = new TestAppZome( junk, 800, 600 );
    testApp->loop( 1000000 );
    return 0;
}
