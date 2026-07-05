
#ifndef SpaceCraftFromOBJ_h
#define SpaceCraftFromOBJ_h

/// @file SpaceCraftFromOBJ.h
/// @brief Low-res OBJ sketch ↔ **SpaceCraft** bridge — partial topology import and sketch-LOD mesh export.
///
/// OBJ is a *layout* format here, not a full ship spec: `v`/`l`/`f` plus `usemtl kind:WorkshopType`
/// (e.g. `girder:GS1_long`, `radiator:LithiumHeatPipe`). Girder cross-section, `nseg`, plate spans, etc.
/// stay **deferred** (`DEFERRED_I` / `DEFERRED_F`) until Lua calls `fulfillTag` / `setGirderParams`.
///
/// Design choices:
/// - One OBJ vertex → one **Node** (no welding at sketch LOD) — preserves editability in Lua
/// - `l` polylines split into 2-node segment girders; use `rope:`/`path:` prefix for intentional chains
/// - Quad `f` faces → **Shield**/**Radiator** by matching *opposite* girder edges; default spans `{0,1}`
/// - **BuildCraft_sketch** emits bare edges + `plate2quad` tris — visualization only, not truss dynamics
/// - **writeSpaceCraftSketchOBJ** round-trips topology with workshop tags (not full FEM mesh)
///
/// Used by `EditSpaceCraft.h` (`fromObj`, `fromObjString`) and `spaceCraftMeshExport -lod sketch`.
/// Path resolution: relative paths use `g_luaScriptDir` (set from `-s` script directory).
///
/// Open issues / caveats:
/// - **No OBJ import for rings, sliders, tanks, guns** — edges + quad plates only; rest stays Lua/procedural.
/// - **Plate inference is fragile**: quad `f` must have two *opposite* edges that are existing girders; ngons warn-skip.
///   Default `g1span`/`g2span` = `{0,1}` — correct for full-edge plates, wrong if face corners ≠ girder endpoints.
/// - **`fulfillTag` / `setGirderParams` are girder-only** — no `fulfillTagRopes`; deferred ropes skipped in blocks LOD.
/// - **Workshop tags ≠ materials**: `usemtl girder:GS1_long` does not register **StickMaterial** — Lua must define it.
/// - **Polylines** (`l 1 2 3 4`) become chained segment girders unless tagged `rope:`/`path:`/`powerline:`.
/// - **No vertex welding** — re-import or `writeSpaceCraftSketchOBJ` can duplicate nodes (plates add extra `v`).
/// - **Sketch LOD is not sim-ready**: no `.truss`, no truss dynamics; TODO raycast mesh source.
/// - **`-svg` export** uses edge-only projection (tris dropped) — shield/radiator fans invisible in SVG.
/// - **Mixing `fromObj` + `BoundNode`/`Slider`**: sliders need host `pointRange` + valid `nseg`/`mseg` for
///   `sideToPath` (see **SpaceCraftComponents.h**) — fulfill girders and run **blocks** mesh build first.

#include <stdio.h>
#include <string.h>
#include <cctype>
#include <vector>
#include <string>

#include "SpaceCraft.h"
#include "MeshBuilder2.h"

namespace SpaceCrafting{

// Directory of the active Lua script (`-s` path); used to resolve relative `fromObj("hull.obj")` paths.
inline char g_luaScriptDir[512] = "";

/// Optional rigid transform applied to every imported vertex (before Lua edits).
struct ObjImportOpts {
    Vec3d  pos   = Vec3dZero;
    Mat3d  rot   = Mat3dIdentity;
    double scale = 1.0;
    char   basePath[512] = "";
};

/// Component IDs created by one import pass — returned to Lua as `{nodes={}, girders={}, ...}`.
struct ObjImportResult {
    std::vector<int> nodes;
    std::vector<int> girders;
    std::vector<int> ropes;
    std::vector<int> shields;
    std::vector<int> radiators;
};

inline void resolveObjPath( const char* path, const char* baseOverride, char* out, int outLen ){
    if(!path || !out || outLen<1){ return; }
    if(path[0]=='/' || (strlen(path)>1 && path[1]==':')){
        strncpy(out, path, outLen-1); out[outLen-1]=0; return;
    }
    const char* base = (baseOverride && baseOverride[0]) ? baseOverride : g_luaScriptDir;
    if(!base || !base[0]){ strncpy(out, path, outLen-1); out[outLen-1]=0; return; }
    snprintf(out, outLen, "%s/%s", base, path);
}

inline Vec3d transformImportPos( const Vec3d& p, const ObjImportOpts& o ){
    Vec3d ps = p; ps.mul(o.scale);
    return o.pos + o.rot.dot(ps);
}

inline void parseUseMtlKind( const char* mtl, char* kind, int kindLen, char* tag, int tagLen ){
    kind[0]=0; tag[0]=0;
    if(!mtl) return;
    const char* colon = strchr(mtl, ':');
    if(colon){
        int nk = (int)(colon - mtl); if(nk>=kindLen) nk=kindLen-1;
        strncpy(kind, mtl, nk); kind[nk]=0;
        strncpy(tag, mtl, tagLen-1); tag[tagLen-1]=0;
    }else{
        strncpy(kind, mtl, kindLen-1); kind[kindLen-1]=0;
        snprintf(tag, tagLen, "%s", mtl);
    }
}

inline bool readLine( FILE* f, char* buf, int buflen, std::string& lineStore ){
    if(!fgets(buf, buflen, f)) return false;
    lineStore = buf;
    return true;
}

inline bool readLineFromString( const char*& p, const char* end, char* buf, int buflen, std::string& lineStore ){
    const char* eol = p;
    while(eol<end && *eol!='\n' && *eol!='\r') eol++;
    int n = (int)(eol - p);
    if(n>=buflen) n=buflen-1;
    if(n<=0) return false;
    memcpy(buf, p, n); buf[n]=0; p = eol;
    while(p<end && (*p=='\n'||*p=='\r')) p++;
    lineStore = buf;
    return true;
}

inline void parseFaceIndices( const char* line, std::vector<int>& face_v ){
    face_v.clear();
    const char* p = line + 1;
    while(*p){
        while(*p && isspace((unsigned char)*p)) p++;
        if(*p=='\0') break;
        int v=0, vt=0, vn=0;
        if(sscanf(p, "%d/%d/%d", &v, &vt, &vn)==3) face_v.push_back(v);
        else if(sscanf(p, "%d//%d", &v, &vn)==2) face_v.push_back(v);
        else if(sscanf(p, "%d/%d", &v, &vt)==2) face_v.push_back(v);
        else if(sscanf(p, "%d", &v)==1) face_v.push_back(v);
        while(*p && !isspace((unsigned char)*p)) p++;
    }
}

inline void parseLineIndices( const char* line, std::vector<int>& line_v ){
    line_v.clear();
    const char* p = line + 1;
    while(*p){
        while(*p && isspace((unsigned char)*p)) p++;
        if(*p=='\0') break;
        int v=0;
        if(sscanf(p, "%d", &v)==1) line_v.push_back(v);
        while(*p && !isspace((unsigned char)*p)) p++;
    }
}

/// @brief Attach shield/radiator from a quad face by pairing opposite girder edges (g01+g23 or g12+g30).
inline int importPlateFromQuad( SpaceCraft& craft, const int* nv, int n, const char* tag, ObjImportResult& res ){
    if(n!=4){ printf("WARNING importPlateFromQuad: face n=%i (need 4)\n", n); return -1; }
    int nid[4];
    for(int i=0;i<4;i++){
        if(nv[i]<0 || nv[i]>=(int)craft.nodes.size()){ printf("WARNING importPlateFromQuad: bad node %i\n", nv[i]); return -1; }
        nid[i]=nv[i];
    }
    int g01 = craft.findGirderByNodes(nid[0], nid[1]);
    int g12 = craft.findGirderByNodes(nid[1], nid[2]);
    int g23 = craft.findGirderByNodes(nid[2], nid[3]);
    int g30 = craft.findGirderByNodes(nid[3], nid[0]);
    int g1=-1, g2=-1;
    if(g01>=0 && g23>=0){ g1=g01; g2=g23; }
    else if(g12>=0 && g30>=0){ g1=g12; g2=g30; }
    else{ printf("WARNING importPlateFromQuad: no opposite girder pair for face tag '%s'\n", tag); return -1; }
    char kind[32]; char fullTag[SKETCH_TAG_LEN];
    parseUseMtlKind(tag, kind, 32, fullTag, SKETCH_TAG_LEN);
    Vec2d span{0,1};
    if(strncmp(kind, "radiator", 8)==0){
        int id = craft.add_RadiatorFromSketch(g1, g2, fullTag, span, span);
        res.radiators.push_back(id); return id;
    }
    int id = craft.add_ShieldFromSketch(g1, g2, fullTag, span, span);
    res.shields.push_back(id); return id;
}

/// @brief Parse one OBJ line into SpaceCraft components; `usemtl` prefix selects girder vs rope vs plate kind.
inline void processObjLine( SpaceCraft& craft, const std::string& line, std::vector<Vec3d>& objVerts,
    std::vector<int>& objToNode, char* curMtl, ObjImportOpts& opts, ObjImportResult& res ){
    if(line.empty() || line[0]=='#') return;
    char type[8]={0};
    sscanf(line.c_str(), "%7s", type);
    if(type[0]=='v' && type[1]==0){
        Vec3d p; if(sscanf(line.c_str(), "v %lf %lf %lf", &p.x,&p.y,&p.z)!=3) return;
        p = transformImportPos(p, opts);
        objVerts.push_back(p);
        int nid = craft.add_Node(p);
        objToNode.push_back(nid);
        res.nodes.push_back(nid);
    }else if(strncmp(type, "usemtl", 6)==0){
        const char* p = line.c_str() + 6;
        while(*p && isspace((unsigned char)*p)) p++;
        strncpy(curMtl, p, SKETCH_TAG_LEN-1); curMtl[SKETCH_TAG_LEN-1]=0;
        size_t L=strlen(curMtl);
        while(L>0 && (curMtl[L-1]=='\n'||curMtl[L-1]=='\r')) curMtl[--L]=0;
    }else if(type[0]=='l'){
        std::vector<int> lv; parseLineIndices(line.c_str(), lv);
        if(lv.size()<2) return;
        char kind[32]; char tag[SKETCH_TAG_LEN];
        parseUseMtlKind(curMtl, kind, 32, tag, SKETCH_TAG_LEN);
        bool bGirder = (strncmp(kind, "girder", 6)==0);
        bool bRope   = (strncmp(kind, "rope", 4)==0) || (strncmp(kind, "path", 4)==0) || (strncmp(kind, "powerline", 9)==0);
        if(!bGirder && !bRope){ bGirder=true; } // default edge → girder
        for(size_t i=0; i+1<lv.size(); i++){
            int ia = lv[i]-1, ib = lv[i+1]-1;
            if(ia<0||(size_t)ia>=objToNode.size()||ib<0||(size_t)ib>=objToNode.size()) continue;
            int na = objToNode[ia], nb = objToNode[ib];
            if(bRope){
                int id = craft.add_RopeDeferred(na, nb, tag);
                res.ropes.push_back(id);
            }else{
                if(lv.size()>2 && i==0) printf("NOTE: girder polyline split into segments (use rope:/path: for chains)\n");
                int id = craft.add_GirderDeferred(na, nb, tag);
                res.girders.push_back(id);
            }
        }
    }else if(type[0]=='f'){
        std::vector<int> fv; parseFaceIndices(line.c_str(), fv);
        if(fv.size()<3) return;
        if(fv.size()==4){
            int nv[4];
            for(int i=0;i<4;i++){
                int oi = fv[i]-1;
                if(oi<0||(size_t)oi>=objToNode.size()) return;
                nv[i]=objToNode[oi];
            }
            importPlateFromQuad(craft, nv, 4, curMtl, res);
        }else{
            printf("WARNING: ngon face n=%zu skipped (only quads for shield/radiator)\n", fv.size());
        }
    }
}

/// @brief Load edge-only OBJ file into **SpaceCraft**; leaves girder/rope params deferred for Lua.
inline bool importOBJToSpaceCraft( SpaceCraft& craft, const char* path, const ObjImportOpts& opts, ObjImportResult& res ){
    FILE* f = fopen(path, "r");
    if(!f){ printf("ERROR importOBJToSpaceCraft: cannot open '%s'\n", path); return false; }
    printf("importOBJToSpaceCraft('%s')\n", path);
    res = ObjImportResult{};
    std::vector<Vec3d> objVerts; (void)objVerts;
    std::vector<int> objToNode;
    char curMtl[SKETCH_TAG_LEN] = "girder:default";
    char buf[4096]; std::string line;
    ObjImportOpts o = opts;
    while(readLine(f, buf, sizeof(buf), line)) processObjLine(craft, line, objVerts, objToNode, curMtl, o, res);
    fclose(f);
    printf("importOBJToSpaceCraft DONE nodes=%zu girders=%zu ropes=%zu shields=%zu radiators=%zu\n",
        res.nodes.size(), res.girders.size(), res.ropes.size(), res.shields.size(), res.radiators.size());
    return true;
}

/// @brief Same as importOBJToSpaceCraft but from an in-memory string (Lua `fromObjString([[...]])`).
inline bool importOBJStringToSpaceCraft( SpaceCraft& craft, const char* text, const ObjImportOpts& opts, ObjImportResult& res ){
    if(!text){ return false; }
    res = ObjImportResult{};
    std::vector<Vec3d> objVerts; (void)objVerts;
    std::vector<int> objToNode;
    char curMtl[SKETCH_TAG_LEN] = "girder:default";
    char buf[4096]; std::string line;
    const char* p = text; const char* end = text + strlen(text);
    ObjImportOpts o = opts;
    while(p<end && readLineFromString(p, end, buf, sizeof(buf), line)) processObjLine(craft, line, objVerts, objToNode, curMtl, o, res);
    printf("importOBJStringToSpaceCraft DONE nodes=%zu girders=%zu\n", res.nodes.size(), res.girders.size());
    return true;
}

// --- sketch LOD mesh build (visualization / print; TODO: raycast source, NOT truss dynamics) ---

inline void sketch_plate_to_mesh( const SpaceCraft& craft, const Plate& pl, Mesh::Builder2& mesh ){
    if(pl.g1<0 || pl.g2<0 || pl.g1>=(int)craft.girders.size() || pl.g2>=(int)craft.girders.size()) return;
    Quad3d qd;
    craft.plate2quad(pl, qd);
    int i0 = mesh.vert(qd.p00); int i1 = mesh.vert(qd.p01); int i2 = mesh.vert(qd.p11); int i3 = mesh.vert(qd.p10);
    mesh.tri(i0,i1,i2); mesh.tri(i0,i2,i3);
}

inline void sketch_wire_box( Mesh::Builder2& mesh, const Vec3d& c, const Mat3d& R, const Vec3d& half ){
    Vec3d p[8];
    for(int i=0;i<8;i++){
        Vec3d q{ ((i&1)?1:-1)*half.x, ((i&2)?1:-1)*half.y, ((i&4)?1:-1)*half.z };
        p[i] = c + R.dot(q);
    }
    int iv[8]; for(int i=0;i<8;i++) iv[i]=mesh.vert(p[i]);
    int e[12][2]={{0,1},{1,3},{3,2},{2,0},{4,5},{5,7},{7,6},{6,4},{0,4},{1,5},{2,6},{3,7}};
    for(int i=0;i<12;i++) mesh.edge(iv[e[i][0]], iv[e[i][1]]);
}

/// @brief Sketch-LOD mesh builder — topology wireframe + plate fans; ignores deferred girder params.
inline void BuildCraft_sketch( Mesh::Builder2& mesh, SpaceCraft& craft ){
    printf("### BuildCraft_sketch()\n");
    for(Node* o: craft.nodes){
        if(o->boundTo==0 && o->component_kind()!=(int)ComponetKind::Slider){
            o->ivert = (int)mesh.verts.size();
            mesh.vert(o->pos);
        }
    }
    for(Girder* o: craft.girders){
        o->update_nodes();
        int i0=o->nodes.x->ivert, i1=o->nodes.y->ivert;
        if(i0>=0 && i1>=0) mesh.edge(i0, i1, 0);
    }
    for(Rope* o: craft.ropes){
        o->update_nodes();
        int i0=o->nodes.x->ivert, i1=o->nodes.y->ivert;
        if(i0>=0 && i1>=0) mesh.edge(i0, i1, 1);
    }
    const int ringSeg=12;
    for(Ring* o: craft.rings){
        Vec3d c=o->pose.pos, ax=o->pose.rot.b, rad=o->pose.rot.c;
        int iv0=-1, iv1=-1;
        for(int i=0;i<=ringSeg;i++){
            double t = (2*M_PI*i)/ringSeg;
            Vec3d p = c + ax*(cos(t)*o->R) + rad*(sin(t)*o->R);
            iv1 = mesh.vert(p);
            if(iv0>=0) mesh.edge(iv0, iv1, 0);
            iv0=iv1;
        }
    }
    for(Radiator* o: craft.radiators) sketch_plate_to_mesh(craft, *o, mesh);
    for(Shield* o: craft.shields) sketch_plate_to_mesh(craft, *o, mesh);
    for(Tank* o: craft.tanks){
        Vec3d h{ o->span.a*0.5, o->span.b*0.5, o->span.c*0.5 };
        sketch_wire_box(mesh, o->pose.pos, o->pose.rot, h);
    }
    for(Thruster* o: craft.thrusters){
        Vec3d h{ o->span.a*0.5, o->span.b, o->span.c*0.5 };
        sketch_wire_box(mesh, o->pose.pos, o->pose.rot, h);
    }
    printf("BuildCraft_sketch() DONE nvert=%zu nedge=%zu ntri=%zu\n", mesh.verts.size(), mesh.edges.size(), mesh.tris.size());
}

/// @brief Export sketch topology OBJ with `usemtl kind:WorkshopType` tokens (not tessellated truss).
inline void writeSpaceCraftSketchOBJ( const char* fname, const SpaceCraft& craft ){
    printf("writeSpaceCraftSketchOBJ(%s)\n", fname);
    FILE* f = fopen(fname, "w");
    if(!f){ printf("ERROR writeSpaceCraftSketchOBJ cannot open '%s'\n", fname); return; }
    fprintf(f, "# SpaceCraft sketch LOD export (usemtl = workshop type tokens)\n");
    for(Node* o: craft.nodes){
        if(o->boundTo==0 && o->component_kind()!=(int)ComponetKind::Slider)
            fprintf(f, "v %g %g %g\n", o->pos.x, o->pos.y, o->pos.z);
    }
    // node id → obj index (1-based, free nodes only)
    std::vector<int> nodeObj(craft.nodes.size(), -1);
    int oi=1;
    for(size_t i=0;i<craft.nodes.size();i++){
        Node* o=craft.nodes[i];
        if(o->boundTo==0 && o->component_kind()!=(int)ComponetKind::Slider) nodeObj[i]=oi++;
    }
    auto emitEdge=[&](int na, int nb, const char* tag){
        if(na<0||(size_t)na>=nodeObj.size()||(size_t)nb>=nodeObj.size()) return;
        if(nodeObj[na]<0||nodeObj[nb]<0) return;
        if(tag && tag[0]) fprintf(f, "usemtl %s\n", tag);
        fprintf(f, "l %i %i\n", nodeObj[na], nodeObj[nb]);
    };
    for(Girder* o: craft.girders){
        const char* t = o->sketchTag[0] ? o->sketchTag : "girder:unknown";
        emitEdge(o->nodes.x->id, o->nodes.y->id, t);
    }
    for(Rope* o: craft.ropes){
        const char* t = o->sketchTag[0] ? o->sketchTag : "rope:unknown";
        emitEdge(o->nodes.x->id, o->nodes.y->id, t);
    }
    for(Shield* o: craft.shields){
        if(o->g1<0||o->g2<0) continue;
        Girder* g1=craft.girders[o->g1]; Girder* g2=craft.girders[o->g2];
        Quad3d qd; craft.plate2quad(*o, qd);
        int iv[4]; Vec3d ps[4]={qd.p00,qd.p01,qd.p11,qd.p10};
        for(int i=0;i<4;i++){ fprintf(f, "v %g %g %g\n", ps[i].x,ps[i].y,ps[i].z); iv[i]=oi++; }
        const char* t = o->sketchTag[0] ? o->sketchTag : "shield:unknown";
        fprintf(f, "usemtl %s\n", t);
        fprintf(f, "f %i %i %i %i\n", iv[0],iv[1],iv[2],iv[3]);
        (void)g1;(void)g2;
    }
    for(Radiator* o: craft.radiators){
        if(o->g1<0||o->g2<0) continue;
        Quad3d qd; craft.plate2quad(*o, qd);
        int iv[4]; Vec3d ps[4]={qd.p00,qd.p01,qd.p11,qd.p10};
        for(int i=0;i<4;i++){ fprintf(f, "v %g %g %g\n", ps[i].x,ps[i].y,ps[i].z); iv[i]=oi++; }
        const char* t = o->sketchTag[0] ? o->sketchTag : "radiator:unknown";
        fprintf(f, "usemtl %s\n", t);
        fprintf(f, "f %i %i %i %i\n", iv[0],iv[1],iv[2],iv[3]);
    }
    fclose(f);
    printf("writeSpaceCraftSketchOBJ DONE\n");
}

} // namespace SpaceCrafting

#endif
