#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

namespace rigidbuild {

constexpr double PI = 3.1415926535897932384626433832795;
constexpr uint32_t INVALID_ID = std::numeric_limits<uint32_t>::max();

using ConnectorId = uint16_t;
using NodeTypeId  = uint16_t;
using StrutTypeId = uint16_t;
using SlotId      = uint16_t;
using NodeId      = uint32_t;
using EdgeId      = uint32_t;

// -----------------------------------------------------------------------------
// Minimal rigid-transform math. Mat3 is column-major: c0,c1,c2 are local axes
// expressed in the parent coordinate system.
// -----------------------------------------------------------------------------

struct Vec3 {
    double x=0.0, y=0.0, z=0.0;
};

inline Vec3 operator+(Vec3 a, Vec3 b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
inline Vec3 operator-(Vec3 a, Vec3 b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
inline Vec3 operator-(Vec3 a)         { return {-a.x,-a.y,-a.z}; }
inline Vec3 operator*(Vec3 a, double s){ return {a.x*s,a.y*s,a.z*s}; }
inline Vec3 operator*(double s, Vec3 a){ return a*s; }
inline Vec3 operator/(Vec3 a, double s){ return {a.x/s,a.y/s,a.z/s}; }

inline double dot(Vec3 a, Vec3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
inline Vec3 cross(Vec3 a, Vec3 b) {
    return {
        a.y*b.z-a.z*b.y,
        a.z*b.x-a.x*b.z,
        a.x*b.y-a.y*b.x
    };
}
inline double norm2(Vec3 a) { return dot(a,a); }
inline double norm(Vec3 a)  { return std::sqrt(norm2(a)); }

inline Vec3 normalized(Vec3 a) {
    const double n = norm(a);
    if (n < 1e-14) throw std::runtime_error("Cannot normalize a zero vector");
    return a/n;
}

struct Mat3 {
    Vec3 c0{1,0,0};
    Vec3 c1{0,1,0};
    Vec3 c2{0,0,1};
};

inline Vec3 operator*(const Mat3& A, Vec3 v) {
    return A.c0*v.x + A.c1*v.y + A.c2*v.z;
}

inline Mat3 operator*(const Mat3& A, const Mat3& B) {
    return { A*B.c0, A*B.c1, A*B.c2 };
}

inline Mat3 transpose(const Mat3& A) {
    return {
        {A.c0.x, A.c1.x, A.c2.x},
        {A.c0.y, A.c1.y, A.c2.y},
        {A.c0.z, A.c1.z, A.c2.z}
    };
}

inline Mat3 rotZ(double a) {
    const double c = std::cos(a);
    const double s = std::sin(a);
    return {
        { c, s, 0},
        {-s, c, 0},
        { 0, 0, 1}
    };
}

inline Mat3 rotAxis(Vec3 axis, double a) {
    axis = normalized(axis);
    const double x=axis.x, y=axis.y, z=axis.z;
    const double c=std::cos(a), s=std::sin(a), t=1.0-c;
    // Construct rows, then transpose into our column representation.
    const Vec3 r0{t*x*x+c,   t*x*y-s*z, t*x*z+s*y};
    const Vec3 r1{t*x*y+s*z, t*y*y+c,   t*y*z-s*x};
    const Vec3 r2{t*x*z-s*y, t*y*z+s*x, t*z*z+c};
    return {
        {r0.x,r1.x,r2.x},
        {r0.y,r1.y,r2.y},
        {r0.z,r1.z,r2.z}
    };
}

// Builds an orthonormal interface frame.
// c2 = outward/forward direction, c1 is as close as possible to upHint.
inline Mat3 makeFrame(Vec3 forward, Vec3 upHint) {
    const Vec3 z = normalized(forward);
    Vec3 y = upHint - z*dot(upHint,z);
    if (norm2(y) < 1e-12) {
        const Vec3 fallback = (std::abs(z.x) < 0.7) ? Vec3{1,0,0}
                                                   : Vec3{0,1,0};
        y = fallback - z*dot(fallback,z);
    }
    y = normalized(y);
    const Vec3 x = normalized(cross(y,z));
    y = cross(z,x);
    return {x,y,z};
}

struct Pose {
    Mat3 R{};
    Vec3 p{};
};

inline Pose operator*(const Pose& A, const Pose& B) {
    return { A.R*B.R, A.p + A.R*B.p };
}

inline Pose inverse(const Pose& A) {
    const Mat3 Rt = transpose(A.R);
    return {Rt, Rt*(-A.p)};
}

// Interface convention:
//   local +Z points OUT of the owning rigid part.
// Two interfaces mate with opposite Z axes. A discrete axial phase may be added.
inline Pose mateOf(const Pose& interfaceA, double phase) {
    const Mat3 flipX{{1,0,0},{0,-1,0},{0,0,-1}}; // 180 deg around local X
    return interfaceA * Pose{rotZ(phase)*flipX, {0,0,0}};
}

inline bool near(Vec3 a, Vec3 b, double eps=1e-8) {
    return norm2(a-b) <= eps*eps;
}

inline double rotationError(const Mat3& A, const Mat3& B) {
    // Frobenius error; adequate for editor validation.
    return std::sqrt(
        norm2(A.c0-B.c0) + norm2(A.c1-B.c1) + norm2(A.c2-B.c2)
    );
}

// -----------------------------------------------------------------------------
// Immutable catalogue data
// -----------------------------------------------------------------------------

struct ConnectorType {
    uint8_t fold = 1;       // frame phase repeats after 2*pi/fold
    uint16_t family = 0;    // semantic family: blue/yellow/red/custom...
    bool chiral = false;    // reserved for future handed matching rules
};

// Stable semantic identity of a slot. Dense SlotId is still used for arrays.
struct SlotKey {
    uint16_t family = 0;
    uint16_t axis   = 0;    // unoriented axis number within the family
    int8_t sign     = +1;   // +1 or -1 along that axis
};

struct SlotDef {
    ConnectorId connector = 0;
    SlotKey key{};
    Pose local{};                 // interface frame in node-local coordinates
    SlotId opposite = 0;
    uint16_t orbit = 0;           // symmetry orbit, often one per family
    uint16_t matchTag = 0;        // optional Ammann/arrow/color rule
    int8_t polarity = 0;          // -1,0,+1 for female/neutral/male-like rules
};

struct SlotMap {
    SlotId slot = 0;
    int16_t phaseStep = 0;  // transformed frame = mapped frame * Rz(step*2pi/fold)
};

struct SymmetryOp {
    Mat3 localRotation{};
    std::vector<SlotMap> map; // one entry per source slot
};

struct NodeType {
    double radius = 0.0;
    std::vector<SlotDef> slots;
    std::vector<SymmetryOp> symmetry;
};

struct StrutEnd {
    ConnectorId connector = 0;
    Pose local{};          // interface frame in strut-local coordinates
    uint16_t matchTag = 0;
    int8_t polarity = 0;
};

struct StrutType {
    double length = 1.0;
    StrutEnd end[2];
};

// Convenient builder for opposite slot pairs.
class NodeTypeBuilder {
public:
    explicit NodeTypeBuilder(double radius) { out_.radius = radius; }

    // Adds two opposite slots. axisIndex should be stable within a family.
    std::pair<SlotId,SlotId> addAxisPair(
        ConnectorId connector,
        uint16_t family,
        uint16_t axisIndex,
        uint16_t orbit,
        Vec3 axis,
        Vec3 upHint,
        uint16_t matchTag=0,
        int8_t polarity=0
    ) {
        axis = normalized(axis);
        const SlotId plus  = static_cast<SlotId>(out_.slots.size());
        const SlotId minus = static_cast<SlotId>(plus+1);

        SlotDef a;
        a.connector = connector;
        a.key = {family,axisIndex,+1};
        a.local = {makeFrame(axis,upHint), axis*out_.radius};
        a.opposite = minus;
        a.orbit = orbit;
        a.matchTag = matchTag;
        a.polarity = polarity;

        SlotDef b;
        b.connector = connector;
        b.key = {family,axisIndex,-1};
        b.local = {makeFrame(-axis,upHint), -axis*out_.radius};
        b.opposite = plus;
        b.orbit = orbit;
        b.matchTag = matchTag;
        b.polarity = polarity;

        out_.slots.push_back(a);
        out_.slots.push_back(b);
        return {plus,minus};
    }

    NodeType finish() { return std::move(out_); }

private:
    NodeType out_;
};

// Build the slot permutation and axial phase induced by a local symmetry rotation.
// This is useful after constructing the slot table: the group operation becomes
// an exact table lookup during editing.
inline SymmetryOp makeSymmetryOp(
    const NodeType& type,
    const std::vector<ConnectorType>& connectors,
    const Mat3& rotation,
    double directionTolerance=1e-7,
    double phaseTolerance=1e-5
) {
    SymmetryOp op;
    op.localRotation = rotation;
    op.map.resize(type.slots.size());

    for (SlotId i=0; i<type.slots.size(); ++i) {
        const SlotDef& src = type.slots[i];
        const Vec3 zt = rotation * src.local.R.c2;
        const Vec3 xt = rotation * src.local.R.c0;

        SlotId best = static_cast<SlotId>(INVALID_ID);
        double bestDot = -2.0;

        for (SlotId j=0; j<type.slots.size(); ++j) {
            const SlotDef& dst = type.slots[j];
            if (dst.connector != src.connector) continue;
            const double d = dot(zt,dst.local.R.c2);
            if (d > bestDot) { bestDot=d; best=j; }
        }

        if (best==static_cast<SlotId>(INVALID_ID) ||
            1.0-bestDot > directionTolerance) {
            throw std::runtime_error("Rotation does not permute this node's slots");
        }

        const SlotDef& dst = type.slots[best];
        const Vec3 xd = dst.local.R.c0;
        const Vec3 yd = dst.local.R.c1;
        const double angle = std::atan2(dot(xt,yd), dot(xt,xd));

        const uint8_t fold = std::max<uint8_t>(1, connectors[src.connector].fold);
        const double stepAngle = 2.0*PI/double(fold);
        int step = int(std::llround(angle/stepAngle));
        const double quantized = step*stepAngle;

        // Compare modulo 2*pi.
        double err = std::remainder(angle-quantized, 2.0*PI);
        if (std::abs(err) > phaseTolerance) {
            throw std::runtime_error(
                "Rotation maps slot direction but not its keyed interface frame"
            );
        }
        step %= fold;
        if (step<0) step += fold;
        op.map[i] = {best,static_cast<int16_t>(step)};
    }
    return op;
}

inline StrutType makeStraightStrut(
    double length,
    ConnectorId connector0,
    ConnectorId connector1
) {
    StrutType s;
    s.length = length;

    // Strut local axis is +Z. End 0 points outward toward -Z; end 1 toward +Z.
    const Mat3 flipX{{1,0,0},{0,-1,0},{0,0,-1}};
    s.end[0].connector = connector0;
    s.end[0].local = {flipX,{0,0,0}};
    s.end[1].connector = connector1;
    s.end[1].local = {Mat3{}, {0,0,length}};
    return s;
}

// -----------------------------------------------------------------------------
// Mutable construction graph
// -----------------------------------------------------------------------------

struct NodeInstance {
    NodeTypeId type = 0;
    Pose world{};
    std::vector<EdgeId> edgeAtSlot; // clear authoring layout; optimize later
};

struct EdgeInstance {
    NodeId node[2]{INVALID_ID,INVALID_ID};
    SlotId slot[2]{0,0};
    StrutTypeId type = 0;
    uint16_t phaseStep[2]{0,0};
    Pose world{}; // strut pose
};

struct OrientedEdge {
    EdgeId edge = INVALID_ID;
    bool reversed = false;
};

struct FaceInstance {
    uint16_t type = 0;
    std::vector<OrientedEdge> boundary;
};

struct OrientedFace {
    uint32_t face = INVALID_ID;
    bool reversed = false;
};

struct CellInstance {
    uint16_t type = 0;
    std::vector<OrientedFace> boundary;
};

struct GrowRequest {
    NodeId sourceNode = INVALID_ID;
    SlotId sourceSlot = 0;
    StrutTypeId strutType = 0;

    NodeTypeId targetNodeType = 0;
    SlotId targetSlot = 0;

    uint16_t sourcePhaseStep = 0;
    uint16_t targetPhaseStep = 0;
};

struct GrowCandidate {
    GrowRequest request{};
    Pose strutWorld{};
    Pose targetNodeWorld{};
};

class Construction {
public:
    std::vector<ConnectorType> connectors;
    std::vector<NodeType> nodeTypes;
    std::vector<StrutType> strutTypes;

    std::vector<NodeInstance> nodes;
    std::vector<EdgeInstance> edges;
    std::vector<FaceInstance> faces;
    std::vector<CellInstance> cells;

    ConnectorId addConnector(const ConnectorType& c) {
        connectors.push_back(c);
        return static_cast<ConnectorId>(connectors.size()-1);
    }

    NodeTypeId addNodeType(NodeType t) {
        nodeTypes.push_back(std::move(t));
        return static_cast<NodeTypeId>(nodeTypes.size()-1);
    }

    StrutTypeId addStrutType(StrutType s) {
        strutTypes.push_back(std::move(s));
        return static_cast<StrutTypeId>(strutTypes.size()-1);
    }

    NodeId addNode(NodeTypeId type, const Pose& world) {
        if (type>=nodeTypes.size()) throw std::out_of_range("Bad NodeTypeId");
        NodeInstance n;
        n.type = type;
        n.world = world;
        n.edgeAtSlot.assign(nodeTypes[type].slots.size(), INVALID_ID);
        nodes.push_back(std::move(n));
        return static_cast<NodeId>(nodes.size()-1);
    }

    bool slotFree(NodeId n, SlotId s) const {
        return n<nodes.size() &&
               s<nodes[n].edgeAtSlot.size() &&
               nodes[n].edgeAtSlot[s]==INVALID_ID;
    }

    bool connectorCompatible(const SlotDef& slot, const StrutEnd& end) const {
        if (slot.connector != end.connector) return false;
        if (slot.matchTag && end.matchTag && slot.matchTag != end.matchTag) return false;
        if (slot.polarity && end.polarity && slot.polarity == end.polarity) return false;
        return true;
    }

    double phaseAngle(ConnectorId c, uint16_t phaseStep) const {
        const uint8_t fold = std::max<uint8_t>(1,connectors.at(c).fold);
        return 2.0*PI*double(phaseStep%fold)/double(fold);
    }

    GrowCandidate previewGrow(const GrowRequest& q) const {
        if (q.sourceNode>=nodes.size()) throw std::out_of_range("Bad source node");
        if (q.strutType>=strutTypes.size()) throw std::out_of_range("Bad strut type");
        if (q.targetNodeType>=nodeTypes.size()) throw std::out_of_range("Bad target type");

        const NodeInstance& A = nodes[q.sourceNode];
        const NodeType& typeA = nodeTypes[A.type];
        const NodeType& typeB = nodeTypes[q.targetNodeType];
        if (q.sourceSlot>=typeA.slots.size()) throw std::out_of_range("Bad source slot");
        if (q.targetSlot>=typeB.slots.size()) throw std::out_of_range("Bad target slot");
        if (!slotFree(q.sourceNode,q.sourceSlot))
            throw std::runtime_error("Source slot is occupied");

        const SlotDef& slotA = typeA.slots[q.sourceSlot];
        const SlotDef& slotB = typeB.slots[q.targetSlot];
        const StrutType& strut = strutTypes[q.strutType];

        if (!connectorCompatible(slotA,strut.end[0]))
            throw std::runtime_error("Source slot and strut end are incompatible");
        if (!connectorCompatible(slotB,strut.end[1]))
            throw std::runtime_error("Target slot and strut end are incompatible");

        const Pose worldSlotA = A.world * slotA.local;
        const Pose desiredEnd0 =
            mateOf(worldSlotA, phaseAngle(slotA.connector,q.sourcePhaseStep));

        const Pose strutWorld = desiredEnd0 * inverse(strut.end[0].local);
        const Pose worldEnd1 = strutWorld * strut.end[1].local;

        const Pose desiredSlotB =
            mateOf(worldEnd1, phaseAngle(slotB.connector,q.targetPhaseStep));
        const Pose nodeBWorld = desiredSlotB * inverse(slotB.local);

        return {q,strutWorld,nodeBWorld};
    }

    // Adds a new node and one connecting strut. Collision tests and exact-coordinate
    // canonicalization are deliberately separate policy layers.
    std::pair<NodeId,EdgeId> commitGrow(const GrowCandidate& c) {
        // Re-check occupancy because another edit may have happened after preview.
        if (!slotFree(c.request.sourceNode,c.request.sourceSlot))
            throw std::runtime_error("Source slot became occupied");

        const NodeId b = addNode(c.request.targetNodeType,c.targetNodeWorld);
        const EdgeId e = static_cast<EdgeId>(edges.size());

        EdgeInstance edge;
        edge.node[0] = c.request.sourceNode;
        edge.node[1] = b;
        edge.slot[0] = c.request.sourceSlot;
        edge.slot[1] = c.request.targetSlot;
        edge.type = c.request.strutType;
        edge.phaseStep[0] = c.request.sourcePhaseStep;
        edge.phaseStep[1] = c.request.targetPhaseStep;
        edge.world = c.strutWorld;
        edges.push_back(edge);

        nodes[edge.node[0]].edgeAtSlot[edge.slot[0]] = e;
        nodes[edge.node[1]].edgeAtSlot[edge.slot[1]] = e;
        return {b,e};
    }

    std::pair<NodeId,EdgeId> grow(const GrowRequest& q) {
        return commitGrow(previewGrow(q));
    }

    // Connect to an already-existing node if its target slot has the pose predicted
    // by previewGrow. Useful when a polygonal loop closes.
    EdgeId connectExisting(
        const GrowRequest& q,
        NodeId targetNode,
        double positionTolerance=1e-7,
        double rotationTolerance=1e-6
    ) {
        if (targetNode>=nodes.size()) throw std::out_of_range("Bad target node");
        const GrowCandidate c = previewGrow(q);
        const NodeInstance& B = nodes[targetNode];

        if (B.type != q.targetNodeType)
            throw std::runtime_error("Existing target has wrong node type");
        if (!slotFree(targetNode,q.targetSlot))
            throw std::runtime_error("Existing target slot is occupied");
        if (!near(B.world.p,c.targetNodeWorld.p,positionTolerance) ||
            rotationError(B.world.R,c.targetNodeWorld.R)>rotationTolerance) {
            throw std::runtime_error("Existing target does not geometrically match");
        }

        const EdgeId e = static_cast<EdgeId>(edges.size());
        EdgeInstance edge;
        edge.node[0]=q.sourceNode; edge.node[1]=targetNode;
        edge.slot[0]=q.sourceSlot; edge.slot[1]=q.targetSlot;
        edge.type=q.strutType;
        edge.phaseStep[0]=q.sourcePhaseStep;
        edge.phaseStep[1]=q.targetPhaseStep;
        edge.world=c.strutWorld;
        edges.push_back(edge);

        nodes[edge.node[0]].edgeAtSlot[edge.slot[0]]=e;
        nodes[edge.node[1]].edgeAtSlot[edge.slot[1]]=e;
        return e;
    }

    // Grow the symmetry orbit of one slot on a central node.
    // opIndices refer to NodeType::symmetry and should normally include identity.
    // All candidates are validated before anything is committed.
    std::vector<std::pair<NodeId,EdgeId>> growOrbit(
        NodeId center,
        SlotId seedSlot,
        const std::vector<uint16_t>& opIndices,
        StrutTypeId strut,
        NodeTypeId targetType,
        SlotId targetSlot,
        uint16_t sourcePhaseStep=0,
        uint16_t targetPhaseStep=0
    ) {
        if (center>=nodes.size()) throw std::out_of_range("Bad center node");
        const NodeType& nt = nodeTypes[nodes[center].type];
        if (seedSlot>=nt.slots.size()) throw std::out_of_range("Bad seed slot");

        std::vector<GrowCandidate> planned;
        std::unordered_set<uint32_t> seenSlots;

        for (uint16_t oi : opIndices) {
            if (oi>=nt.symmetry.size()) throw std::out_of_range("Bad symmetry op");
            const SlotMap sm = nt.symmetry[oi].map.at(seedSlot);
            if (!seenSlots.insert(sm.slot).second) continue; // stabilizer duplicate

            const ConnectorId c = nt.slots[sm.slot].connector;
            const uint8_t fold = std::max<uint8_t>(1,connectors[c].fold);
            const uint16_t phase =
                uint16_t((sourcePhaseStep + sm.phaseStep) % fold);

            GrowRequest q;
            q.sourceNode=center;
            q.sourceSlot=sm.slot;
            q.strutType=strut;
            q.targetNodeType=targetType;
            q.targetSlot=targetSlot;
            q.sourcePhaseStep=phase;
            q.targetPhaseStep=targetPhaseStep;
            planned.push_back(previewGrow(q));
        }

        std::vector<std::pair<NodeId,EdgeId>> result;
        result.reserve(planned.size());
        for (const GrowCandidate& c : planned) result.push_back(commitGrow(c));
        return result;
    }
};

} // namespace rigidbuild