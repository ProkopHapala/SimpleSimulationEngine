#ifndef SpaceCraft2Mesh_blocks_h
#define SpaceCraft2Mesh_blocks_h

#include <vector>

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "MeshBuilder2.h"
#include "SpaceCraft.h"

namespace Mesh {

using namespace SpaceCrafting;

// Simplified helpers mirroring BuildCraft_truss (no CMesh/nodeMeshes)

inline void node_to_mesh(Node* o, Builder2& mesh) {
    // Bound nodes should not create new verts; free nodes become verts
    if (o->boundTo == 0) {
        o->ivert = mesh.verts.size();
        mesh.vert(o->pos);
    }
}

inline void girder_to_mesh(Girder* o, Builder2& mesh) {
    o->update_nodes();
    mesh.block();
    mesh.girder1(o->nodes.x->pos, o->nodes.y->pos, o->up, o->nseg, o->wh.a, o->st);
    mesh.girder1_caps(o->nodes.x->ivert, o->nodes.y->ivert, o->st.x);
    Quat4i& b = mesh.blocks.back();
    o->pointRange = {b.x, (int)mesh.verts.size()};
    o->stickRange = {b.y, (int)mesh.edges.size()};
}

inline void ring_to_mesh(Ring* o, Builder2& mesh) {
    o->update_nodes();
    mesh.block();
    mesh.wheel(o->pose.pos, o->pose.pos + o->pose.rot.b * o->R, o->pose.rot.c, o->nseg, o->wh, o->st);
    Quat4i& b = mesh.blocks.back();
    o->pointRange = {b.x, (int)mesh.verts.size()};
    o->stickRange = {b.y, (int)mesh.edges.size()};
}

inline void slider_to_mesh(Slider* o, Builder2& mesh) {
    static const double g_Rcollapse = 0.1;
    static const double g_Ranchor   = 10.0;
    static const int    g_AnchorStickType = 0;
    int pre_verts_size = mesh.verts.size();
    int pre_edges_size = mesh.edges.size();
    Vec2i host_point_range = o->boundTo->pointRange;
    o->ivert = mesh.make_anchor_point(o->pos, g_AnchorStickType, g_Rcollapse, g_Ranchor, 0, 0, host_point_range.x, host_point_range.y);
    o->pointRange = {pre_verts_size, (int)mesh.verts.size()};
    o->stickRange = {pre_edges_size, (int)mesh.edges.size()};
    printf("slider_to_mesh() o->id %i o->ivert=%i o->pointRange(%i,%i) o->stickRange(%i,%i)\n", o->id, o->ivert, o->pointRange.x, o->pointRange.y, o->stickRange.x, o->stickRange.y);
}

inline void rope_to_mesh(Rope* o, Builder2& mesh) {
    o->update_nodes();
    mesh.block();
    Quat4i& b = mesh.blocks.back();
    mesh.rope(o->nodes.x->ivert, o->nodes.y->ivert, o->face_mat, o->nseg);
    o->pointRange = {b.x, (int)mesh.verts.size()};
    o->stickRange = {b.y, (int)mesh.edges.size()};
}

inline int checkSpaceCraftMesh(const Builder2& mesh, const SpaceCraft& craft, bool bPrint=true, bool bExit=true ){
    std::vector<bool> connected(mesh.verts.size(), false);
    for (const auto& edge : mesh.edges) {
        if (edge.x < mesh.verts.size()) connected[edge.x] = true;
        if (edge.y < mesh.verts.size()) connected[edge.y] = true;
    }
    int unconnected_count = 0;
    for (size_t i = 0; i < mesh.verts.size(); ++i) {
        if (!connected[i]) {
            if (bPrint){
                printf("WARNING: checkSpaceCraftMesh() vertex %3zu at (%10.3f, %10.3f, %10.3f) is not connected to any edge.\n",  i, mesh.verts[i].pos.x, mesh.verts[i].pos.y, mesh.verts[i].pos.z);
                printf("        - vertex %zu belong to SpaceCraft component: \n", i);
            }
            int j = craft.find_mesh_element((int)i, 'v', bPrint );
            if(j>=0){ unconnected_count++; }
        }
    }
    if (unconnected_count > 0 && bExit) {
        if(bPrint) printf("ERROR: checkSpaceCraftMesh() found %i unconnected vertices.\n", unconnected_count);
    }
    return unconnected_count;
}

inline void BuildCraft_blocks( Builder2& mesh, SpaceCraft& craft, double max_size=-1, double node_scale=1.0 ){
    printf( "### BuildCraft_blocks() [standalone-truss]\n" );
    if(max_size>0){ mesh.max_size=max_size; };
    (void)node_scale;

    mesh.block();
    int ip0 = mesh.verts.size();
    (void)ip0;

    printf("BuildCraft_blocks().nodes n=%zu\n", craft.nodes.size());
    for(Node* o: craft.nodes){ node_to_mesh(o, mesh); }

    printf("BuildCraft_blocks().girders n=%zu\n", craft.girders.size() );
    for(Girder* o: craft.girders){ girder_to_mesh(o, mesh); }

    printf("BuildCraft_blocks().ropes n=%zu\n", craft.ropes.size() );
    for(Rope* o: craft.ropes){ rope_to_mesh(o, mesh); }

    printf("BuildCraft_blocks().rings n=%zu\n", craft.rings.size() );
    for(Ring* o: craft.rings){ ring_to_mesh(o, mesh); }

    // sliders, radiators, welds etc. are ignored for now in pure truss view

    checkSpaceCraftMesh(mesh, craft);
    printf("BuildCraft_blocks() DONE! : npoint %zu nstick %zu nblocks %zu\n", mesh.verts.size(), mesh.edges.size(), mesh.blocks.size());
}


// --- File export helper (mesh + workshop -> text truss file, no TrussDynamics dependency) ----

// Text format (whitespace-separated, '#' starts comment):
//   # TrussSim v1
//   meta <npoint> <nbond>
//   p <id>  <x> <y> <z>  <mass> <nneigh>
//   ... (npoint lines)
//   s <id>  <i> <j>  <kpull> <kpush> <l0> <damping>
//   ... (nbond lines)

inline void exportSimToFile( const char* fname, const Builder2& mesh, const SpaceCraftWorkshop& shop ){
    if(!fname){ printf("exportSimToFile(): fname is null\n"); return; }
    FILE* f = fopen(fname, "w");
    if(!f){ printf("exportSimToFile(): cannot open '%s' for writing\n", fname); return; }

    const int np = (int)mesh.verts.size();
    const int nb = (int)mesh.edges.size();

    // per-point accumulators
    std::vector<double> mass(np, 0.0);
    std::vector<int>    nneigh(np, 0);

    struct StickRec{ int i,j; double kpull,kpush,l0,damping; };
    std::vector<StickRec> sticks(nb);

    // derive masses and spring parameters directly from mesh edges + stick materials
    for(int ib=0; ib<nb; ib++){
        const Quat4i& e = mesh.edges[ib];
        const int ia = e.x;
        const int ibp = e.y;
        if(ia<0 || ia>=np || ibp<0 || ibp>=np){
            printf("exportSimToFile(): edge[%d] has invalid verts (%d,%d) for np=%d\n", ib, ia, ibp, np);
            fclose(f);
            return;
        }
        nneigh[ia]++;
        nneigh[ibp]++;

        if(e.w < 0 || e.w >= (int)shop.stickMaterials.vec.size()){
            printf("exportSimToFile(): edge[%d].type=%d out of range stickMaterials.size()=%zu\n", ib, e.w, shop.stickMaterials.vec.size());
            fclose(f);
            return;
        }
        const StickMaterial& mat = *shop.stickMaterials.vec[e.w];

        const Vec3d& pa = mesh.verts[ia].pos;
        const Vec3d& pb = mesh.verts[ibp].pos;
        double l_geom = (pb - pa).norm();
        double l0 = l_geom * (1.0 - mat.preStrain);
        double m  = l0 * mat.linearDensity;

        mass[ia] += 0.5 * m;
        mass[ibp] += 0.5 * m;

        StickRec rec;
        rec.i = ia;
        rec.j = ibp;
        rec.l0 = l0;
        rec.kpull   = (l0 > 0.0) ? (mat.Kpull / l0) : 0.0;
        rec.kpush   = (l0 > 0.0) ? (mat.Kpush / l0) : 0.0;
        rec.damping = mat.damping;
        sticks[ib] = rec;
    }

    // header
    fprintf(f, "# TrussSim v1\n");
    fprintf(f, "meta %d %d\n", np, nb);

    // points
    for(int i=0; i<np; i++){
        const Vec3d& p = mesh.verts[i].pos;
        fprintf(f, "p %d %.16g %.16g %.16g %.16g %d\n", i, p.x, p.y, p.z, mass[i], nneigh[i]);
    }

    // sticks
    for(int ib=0; ib<nb; ib++){
        const StickRec& s = sticks[ib];
        fprintf(f, "s %d %d %d %.16g %.16g %.16g %.16g\n", ib, s.i, s.j, s.kpull, s.kpush, s.l0, s.damping);
    }

    fclose(f);
    printf("exportSimToFile(): written np=%d nb=%d to '%s'\n", np, nb, fname);
}


} // namespace Mesh

#endif
