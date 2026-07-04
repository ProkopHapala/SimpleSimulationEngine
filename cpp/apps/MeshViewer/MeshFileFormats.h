
#ifndef MeshFileFormats_h
#define MeshFileFormats_h

/// @file MeshFileFormats.h
/// @brief Mesh snapshot I/O: .obj+ (extended OBJ with header counts) and .npz (packed binary arrays).
///
/// A MeshSnapshot is a single time-step of a mesh animation: vertices, edges, and faces.
/// Each snapshot can have different counts (topology may change between frames, like .xyz
/// for molecules where atom count/bonds change). This is the key difference from standard
/// OBJ which assumes a single static mesh.
///
/// Two formats:
/// - **.obj+**: Human-readable OBJ with optional header comment `# nverts nedges nfaces`.
///   If the header is present, the reader can pre-allocate arrays in one pass (no double scan).
///   Without the header, it falls back to the standard two-pass OBJ reader.
///   Multiple snapshots in one file are separated by `# --- snapshot N ---` lines.
/// - **.npz**: Binary packed format. Header: magic(4) + nverts(4) + nedges(4) + nfaces(4) + ngons(4).
///   Then raw arrays: verts[nverts*3] (float64), edges[nedges*2] (int32), tris[nfaces*3] (int32),
///   ngons[nfaces] (int32). One file per snapshot, or concatenated with snapshot separators.
///
/// The .npz format is designed for fast loading: one fread for the header, one fread per array.
/// No parsing, no text conversion. For animation sequences, files are named frame_0000.npz, etc.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <dirent.h>
#include <sys/stat.h>
#include <algorithm>
#include <cmath>

#include "Vec3.h"
#include "datatypes.h"
#include "CMesh.h"

// ========== Free functions on CMesh ==========

inline void cmesh_computeNormals(const CMesh* m, Vec3d* nors) {
    for (int i = 0; i < m->nvert; i++) nors[i] = Vec3dZero;
    for (int i = 0; i < m->ntri; i++) {
        Vec3i iv = m->tris[i];
        Vec3d a = m->verts[iv.a], b = m->verts[iv.b], c = m->verts[iv.c];
        Vec3d nor; nor.set_cross(b - a, c - a);
        nors[iv.a].add(nor); nors[iv.b].add(nor); nors[iv.c].add(nor);
    }
    for (int i = 0; i < m->nvert; i++) nors[i].normalize();
}

inline Vec3d cmesh_center(const CMesh* m) {
    Vec3d c = Vec3dZero;
    for (int i = 0; i < m->nvert; i++) c.add(m->verts[i]);
    if (m->nvert > 0) c.mul(1.0 / m->nvert);
    return c;
}

inline double cmesh_boundingRadius(const CMesh* m, const Vec3d& center) {
    double rmax = 0;
    for (int i = 0; i < m->nvert; i++) {
        double r = (m->verts[i] - center).norm();
        if (r > rmax) rmax = r;
    }
    return rmax;
}

/// MeshSnapshot: owning wrapper around CMesh — adds alloc/dealloc/dtor to the non-owning C-style CMesh.
/// Inheritance so it's implicitly convertible to CMesh* for passing to C APIs, GPU upload, etc.
struct MeshSnapshot : public CMesh {

    void alloc(int nv, int ne, int nf) {
        dealloc();
        nvert = nv; nedge = ne; ntri = nf; nfaces = nf;
        verts = new Vec3d[nv];
        if (ne > 0) edges = new Vec2i[ne]; else edges = nullptr;
        if (nf > 0) { tris = new Vec3i[nf]; ngons = new int[nf]; } else { tris = nullptr; ngons = nullptr; }
    }

    void dealloc() {
        delete[] verts; verts = nullptr;
        delete[] edges; edges = nullptr;
        delete[] tris;  tris  = nullptr;
        delete[] ngons; ngons = nullptr;
        nvert = nedge = ntri = nfaces = 0;
    }

    MeshSnapshot() { nvert = 0; nedge = 0; ntri = 0; nfaces = 0; verts = nullptr; edges = nullptr; tris = nullptr; ngons = nullptr; }
    ~MeshSnapshot() { dealloc(); }
};

struct MeshAnimation {
    std::vector<CMesh*> snapshots;
    int current = 0;

    ~MeshAnimation() { for (auto* s : snapshots) delete (MeshSnapshot*)s; }

    int size() const { return snapshots.size(); }
    CMesh* get(int i) const { return snapshots[((i % size()) + size()) % size()]; }
    CMesh* currentSnap() const { return size() > 0 ? get(current) : nullptr; }
    void next() { if (size() > 0) current = (current + 1) % size(); }
    void prev() { if (size() > 0) current = (current - 1 + size()) % size(); }
};

// ========== File scanning ==========

#include "GUI.h" // for TreeViewTree, TreeViewItem

inline bool hasSuffix(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

inline bool isMeshFile(const std::string& name) {
    return hasSuffix(name, ".obj") || hasSuffix(name, ".obj+") || hasSuffix(name, ".npz") || hasSuffix(name, ".npy");
}

/// Populate a TreeViewTree from a directory. Recursively builds subdirectories.
/// Files are included only if they are mesh files (.obj, .obj+, .npz, .npy) or if bShowAll=true.
/// Sets path and bDir on each TreeViewItem.
inline void dir2tree(TreeViewTree& node, const char* name, const std::string& prefix, bool bShowAll=false) {
    node.content.caption = name;
    std::string path = (prefix.length() == 0) ? name : (prefix + "/" + name);
    node.content.path = path;
    DIR* dir = opendir(path.c_str());
    if (dir) {
        node.content.bDir = true;
        std::vector<std::pair<std::string, bool>> entries; // name, isDir
        struct dirent* ent;
        while ((ent = readdir(dir)) != nullptr) {
            if (ent->d_name[0] == '.') continue;
            std::string childName = ent->d_name;
            std::string childPath = path + "/" + childName;
            struct stat st;
            if (stat(childPath.c_str(), &st) != 0) continue;
            bool isDir = S_ISDIR(st.st_mode);
            if (!isDir && !bShowAll && !isMeshFile(childName)) continue;
            entries.push_back({childName, isDir});
        }
        closedir(dir);
        std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
            if (a.second != b.second) return a.second > b.second; // dirs first
            return a.first < b.first;
        });
        for (auto& [childName, isDir] : entries) {
            TreeViewTree* child = new TreeViewTree();
            child->parrent = &node;
            node.branches.push_back(child);
            dir2tree(*child, childName.c_str(), path, bShowAll);
        }
    } else {
        node.content.bDir = false;
    }
}

// ========== OBJ+ reader ==========

inline bool readOBJPlus(const std::string& path, MeshAnimation& anim) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) { printf("Cannot open %s\n", path.c_str()); return false; }

    // Check for header: "# nverts nedges nfaces"
    char line[1024];
    int hdr_nv = -1, hdr_ne = -1, hdr_nf = -1;

    // Peek at first line
    long pos0 = ftell(f);
    if (fgets(line, sizeof(line), f)) {
        if (sscanf(line, "# %d %d %d", &hdr_nv, &hdr_ne, &hdr_nf) == 3) {
            // header present — single-pass read
        } else {
            // no header, rewind
            fseek(f, pos0, SEEK_SET);
            hdr_nv = hdr_ne = hdr_nf = -1;
        }
    }

    // Read all lines, collecting snapshots
    // Snapshots are separated by "# --- snapshot N ---"
    std::vector<Vec3d> verts;
    std::vector<Vec2i> edges;
    std::vector<Vec3i> tris;
    std::vector<int>   ngons;

    auto flushSnapshot = [&]() {
        if (verts.empty()) return;
        MeshSnapshot* snap = new MeshSnapshot();
        snap->alloc(verts.size(), edges.size(), tris.size());
        memcpy(snap->verts, verts.data(), verts.size() * sizeof(Vec3d));
        if (edges.size() > 0) memcpy(snap->edges, edges.data(), edges.size() * sizeof(Vec2i));
        if (tris.size() > 0) {
            memcpy(snap->tris, tris.data(), tris.size() * sizeof(Vec3i));
            for (int i = 0; i < (int)tris.size(); i++) snap->ngons[i] = 3;
        }
        anim.snapshots.push_back((CMesh*)snap);
        verts.clear(); edges.clear(); tris.clear(); ngons.clear();
    };

    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#') {
            // check for snapshot separator
            if (strstr(line, "--- snapshot")) { flushSnapshot(); continue; }
            // check for header counts (already handled, but skip)
            continue;
        }
        if (line[0] == 'v' && line[1] == ' ') {
            double x, y, z;
            if (sscanf(line, "v %lf %lf %lf", &x, &y, &z) == 3) {
                if (std::isnan(x) || std::isnan(y) || std::isnan(z)) {
                    static int nanWarn = 0;
                    if (nanWarn < 5) printf("WARNING: NaN vertex at line %ld, replacing with (0,0,0)\n", ftell(f));
                    nanWarn++;
                    x = y = z = 0.0;
                }
                verts.push_back({x, y, z});
            }
        } else if (line[0] == 'l' && line[1] == ' ') {
            int a, b;
            if (sscanf(line, "l %d %d", &a, &b) == 2) edges.push_back({a - 1, b - 1}); // OBJ is 1-indexed
        } else if (line[0] == 'f' && line[1] == ' ') {
            // face: f a b c  (or f a//n b//n c//n, f a/b/c ...)
            int a, b, c;
            if (sscanf(line, "f %d %*s %d %*s %d", &a, &b, &c) == 3 ||
                sscanf(line, "f %d/%*s %d/%*s %d/%*s", &a, &b, &c) == 3 ||
                sscanf(line, "f %d %d %d", &a, &b, &c) == 3) {
                tris.push_back({a - 1, b - 1, c - 1});
                ngons.push_back(3);
            }
        }
    }
    fclose(f);
    flushSnapshot(); // flush last snapshot
    return anim.size() > 0;
}

// ========== NPZ binary reader ==========

// NPZ format (simple packed binary, not numpy .npz):
// Header: magic[4] = "MNPZ", nverts(int32), nedges(int32), nfaces(int32)
// Data: verts[nverts*3] (float64), edges[nedges*2] (int32), tris[nfaces*3] (int32)
// For animation: concatenated snapshots, each starting with magic

inline bool readNPZ(const std::string& path, MeshAnimation& anim) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) { printf("Cannot open %s\n", path.c_str()); return false; }

    const char MAGIC[4] = {'M', 'N', 'P', 'Z'};
    while (true) {
        char magic[4];
        if (fread(magic, 1, 4, f) != 4) break;
        if (memcmp(magic, MAGIC, 4) != 0) { printf("Bad NPZ magic\n"); break; }
        int32_t nv, ne, nf;
        if (fread(&nv, 4, 1, f) != 1) break;
        if (fread(&ne, 4, 1, f) != 1) break;
        if (fread(&nf, 4, 1, f) != 1) break;
        MeshSnapshot* snap = new MeshSnapshot();
        snap->alloc(nv, ne, nf);
        fread(snap->verts, sizeof(Vec3d), nv, f);
        if (ne > 0) fread(snap->edges, sizeof(Vec2i), ne, f);
        if (nf > 0) {
            fread(snap->tris, sizeof(Vec3i), nf, f);
            for (int i = 0; i < nf; i++) snap->ngons[i] = 3;
        }
        anim.snapshots.push_back((CMesh*)snap);
    }
    fclose(f);
    return anim.size() > 0;
}

// ========== NPZ binary writer ==========

inline bool writeNPZ(const std::string& path, const MeshAnimation& anim) {
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) { printf("Cannot write %s\n", path.c_str()); return false; }
    const char MAGIC[4] = {'M', 'N', 'P', 'Z'};
    for (int i = 0; i < anim.size(); i++) {
        const CMesh* s = anim.snapshots[i];
        fwrite(MAGIC, 1, 4, f);
        int32_t nv = s->nvert, ne = s->nedge, nf = s->ntri;
        fwrite(&nv, 4, 1, f); fwrite(&ne, 4, 1, f); fwrite(&nf, 4, 1, f);
        fwrite(s->verts, sizeof(Vec3d), nv, f);
        if (ne > 0) fwrite(s->edges, sizeof(Vec2i), ne, f);
        if (nf > 0) fwrite(s->tris, sizeof(Vec3i), nf, f);
    }
    fclose(f);
    return true;
}

// ========== SVG writer ==========

struct SVGExportOpts {
    bool   showFaces      = true;
    double faceOpacity    = 0.5;    // 0=transparent, 1=opaque
    bool   showBackfaces  = false;  // cull backfaces if false
    bool   showEdges      = true;
    bool   showPoints     = false;  // off by default
    double pointRadius    = 1.5;    // small points
    bool   showVertNums   = false;  // vertex index labels
    bool   showEdgeNums   = false;  // edge index labels
    int    vertNumOffset  = 0;      // starting index for vertex labels (for multi-mesh)
    int    edgeNumOffset  = 0;      // starting index for edge labels
    int    fontSize       = 10;
    int    projAxis       = 2;      // 0=YZ plane, 1=XZ plane, 2=XY plane (drop axis)
    int    width          = 800;
    int    height         = 600;
    int    margin         = 10;     // tight fit
    const char* bgColor   = "#fefefe";
    const char* faceColor = "#4488cc";
    const char* edgeColor = "#226633";
    const char* pointColor= "#cc6600";
    bool   showNormals    = false;  // draw face normals as arrows
    double normalLength   = 0.1;    // fraction of bbox size
    const char* normalColor= "#ff0000";
    const char* vertNumColor = "#0000aa";
    const char* edgeNumColor = "#aa0000";
};

// Compute 2D bbox of a single mesh in projected coordinates
inline void svg_meshBBox(const CMesh& mesh, int ix, int iy, double& minx, double& miny, double& maxx, double& maxy) {
    for (int i = 0; i < mesh.nvert; i++) {
        double px = mesh.verts[i].array[ix];
        double py = mesh.verts[i].array[iy];
        if (px < minx) minx = px; if (px > maxx) maxx = px;
        if (py < miny) miny = py; if (py > maxy) maxy = py;
    }
}

// Write one mesh into an already-open SVG file using pre-computed scale/offset.
// vertOffset/edgeOffset allow continuous numbering across multiple meshes.
inline void svgWriteMesh(FILE* f, const CMesh& mesh, const SVGExportOpts& opts,
                         int ix, int iy, double scale, double ox, double oy,
                         int vertOffset, int edgeOffset) {
    auto sx = [&](double px) { return ox + px * scale; };
    auto sy = [&](double py) { return opts.height - (oy + py * scale); };

    // Faces first (so edges/labels draw on top)
    if (opts.showFaces && mesh.ntri > 0) {
        char opacityStr[16];
        snprintf(opacityStr, sizeof(opacityStr), "%.2f", opts.faceOpacity);
        fprintf(f, "<g fill=\"%s\" fill-opacity=\"%s\" stroke=\"none\">\n", opts.faceColor, opacityStr);
        for (int i = 0; i < mesh.ntri; i++) {
            Vec3i iv = mesh.tris[i];
            if (iv.a >= mesh.nvert || iv.b >= mesh.nvert || iv.c >= mesh.nvert) continue;
            double ax_ = mesh.verts[iv.a].array[ix], ay_ = mesh.verts[iv.a].array[iy];
            double bx_ = mesh.verts[iv.b].array[ix], by_ = mesh.verts[iv.b].array[iy];
            double cx_ = mesh.verts[iv.c].array[ix], cy_ = mesh.verts[iv.c].array[iy];
            double signedArea = (bx_ - ax_) * (cy_ - ay_) - (cx_ - ax_) * (by_ - ay_);
            if (!opts.showBackfaces && signedArea < 0) continue;
            fprintf(f, "  <polygon points=\"%.2f,%.2f %.2f,%.2f %.2f,%.2f\"/>\n",
                sx(ax_), sy(ay_), sx(bx_), sy(by_), sx(cx_), sy(cy_));
        }
        fprintf(f, "</g>\n");
    }

    // Edges (from explicit edge list, or derived from triangles)
    if (opts.showEdges) {
        fprintf(f, "<g stroke=\"%s\" stroke-width=\"1\" fill=\"none\">\n", opts.edgeColor);
        if (mesh.nedge > 0) {
            for (int i = 0; i < mesh.nedge; i++) {
                int a = mesh.edges[i].x, b = mesh.edges[i].y;
                if (a >= mesh.nvert || b >= mesh.nvert) continue;
                fprintf(f, "  <line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\"/>\n",
                    sx(mesh.verts[a].array[ix]), sy(mesh.verts[a].array[iy]),
                    sx(mesh.verts[b].array[ix]), sy(mesh.verts[b].array[iy]));
            }
        } else if (mesh.ntri > 0) {
            for (int i = 0; i < mesh.ntri; i++) {
                Vec3i iv = mesh.tris[i];
                if (iv.a >= mesh.nvert || iv.b >= mesh.nvert || iv.c >= mesh.nvert) continue;
                fprintf(f, "  <line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\"/>\n",
                    sx(mesh.verts[iv.a].array[ix]), sy(mesh.verts[iv.a].array[iy]),
                    sx(mesh.verts[iv.b].array[ix]), sy(mesh.verts[iv.b].array[iy]));
                fprintf(f, "  <line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\"/>\n",
                    sx(mesh.verts[iv.b].array[ix]), sy(mesh.verts[iv.b].array[iy]),
                    sx(mesh.verts[iv.c].array[ix]), sy(mesh.verts[iv.c].array[iy]));
                fprintf(f, "  <line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\"/>\n",
                    sx(mesh.verts[iv.c].array[ix]), sy(mesh.verts[iv.c].array[iy]),
                    sx(mesh.verts[iv.a].array[ix]), sy(mesh.verts[iv.a].array[iy]));
            }
        }
        fprintf(f, "</g>\n");
    }

    // Points
    if (opts.showPoints) {
        fprintf(f, "<g fill=\"%s\">\n", opts.pointColor);
        for (int i = 0; i < mesh.nvert; i++)
            fprintf(f, "  <circle cx=\"%.2f\" cy=\"%.2f\" r=\"%.1f\"/>\n",
                sx(mesh.verts[i].array[ix]), sy(mesh.verts[i].array[iy]), opts.pointRadius);
        fprintf(f, "</g>\n");
    }

    // Vertex numbers
    if (opts.showVertNums) {
        fprintf(f, "<g fill=\"%s\" font-size=\"%d\" font-family=\"monospace\">\n", opts.vertNumColor, opts.fontSize);
        for (int i = 0; i < mesh.nvert; i++)
            fprintf(f, "  <text x=\"%.2f\" y=\"%.2f\">%d</text>\n",
                sx(mesh.verts[i].array[ix]) + 2, sy(mesh.verts[i].array[iy]) - 2, i + vertOffset);
        fprintf(f, "</g>\n");
    }

    // Edge numbers (at edge midpoints)
    if (opts.showEdgeNums) {
        fprintf(f, "<g fill=\"%s\" font-size=\"%d\" font-family=\"monospace\">\n", opts.edgeNumColor, opts.fontSize);
        if (mesh.nedge > 0) {
            for (int i = 0; i < mesh.nedge; i++) {
                int a = mesh.edges[i].x, b = mesh.edges[i].y;
                if (a >= mesh.nvert || b >= mesh.nvert) continue;
                double mx = (mesh.verts[a].array[ix] + mesh.verts[b].array[ix]) * 0.5;
                double my = (mesh.verts[a].array[iy] + mesh.verts[b].array[iy]) * 0.5;
                fprintf(f, "  <text x=\"%.2f\" y=\"%.2f\">%d</text>\n", sx(mx) + 1, sy(my) - 1, i + edgeOffset);
            }
        } else if (mesh.ntri > 0) {
            int ei = 0;
            for (int i = 0; i < mesh.ntri; i++) {
                Vec3i iv = mesh.tris[i];
                if (iv.a >= mesh.nvert || iv.b >= mesh.nvert || iv.c >= mesh.nvert) continue;
                double mx, my;
                mx = (mesh.verts[iv.a].array[ix] + mesh.verts[iv.b].array[ix]) * 0.5;
                my = (mesh.verts[iv.a].array[iy] + mesh.verts[iv.b].array[iy]) * 0.5;
                fprintf(f, "  <text x=\"%.2f\" y=\"%.2f\">%d</text>\n", sx(mx) + 1, sy(my) - 1, ei++ + edgeOffset);
                mx = (mesh.verts[iv.b].array[ix] + mesh.verts[iv.c].array[ix]) * 0.5;
                my = (mesh.verts[iv.b].array[iy] + mesh.verts[iv.c].array[iy]) * 0.5;
                fprintf(f, "  <text x=\"%.2f\" y=\"%.2f\">%d</text>\n", sx(mx) + 1, sy(my) - 1, ei++ + edgeOffset);
                mx = (mesh.verts[iv.c].array[ix] + mesh.verts[iv.a].array[ix]) * 0.5;
                my = (mesh.verts[iv.c].array[iy] + mesh.verts[iv.a].array[iy]) * 0.5;
                fprintf(f, "  <text x=\"%.2f\" y=\"%.2f\">%d</text>\n", sx(mx) + 1, sy(my) - 1, ei++ + edgeOffset);
            }
        }
        fprintf(f, "</g>\n");
    }

    // Face normals (lines from centroid)
    if (opts.showNormals && mesh.ntri > 0) {
        double bminx = 1e30, bminy = 1e30, bmaxx = -1e30, bmaxy = -1e30;
        svg_meshBBox(mesh, ix, iy, bminx, bminy, bmaxx, bmaxy);
        double diag = sqrt((bmaxx-bminx)*(bmaxx-bminx) + (bmaxy-bminy)*(bmaxy-bminy));
        if (diag < 1e-12) diag = 1.0;
        double nLen = diag * opts.normalLength;

        fprintf(f, "<g stroke=\"%s\" stroke-width=\"0.8\" fill=\"none\">\n", opts.normalColor);
        for (int i = 0; i < mesh.ntri; i++) {
            Vec3i iv = mesh.tris[i];
            if (iv.a >= mesh.nvert || iv.b >= mesh.nvert || iv.c >= mesh.nvert) continue;
            Vec3d a = mesh.verts[iv.a], b = mesh.verts[iv.b], c = mesh.verts[iv.c];
            Vec3d nor; nor.set_cross(b - a, c - a); nor.normalize();
            Vec3d cen = (a + b + c) * (1.0/3.0);
            fprintf(f, "  <line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\"/>\n",
                sx(cen.array[ix]), sy(cen.array[iy]),
                sx(cen.array[ix] + nor.array[ix]*nLen), sy(cen.array[iy] + nor.array[iy]*nLen));
        }
        fprintf(f, "</g>\n");
    }
}

// Write multiple meshes into one SVG (tight bbox across all meshes)
inline bool writeSVGMulti(const std::string& path, const CMesh** meshes, int nMeshes, const SVGExportOpts& opts = SVGExportOpts()) {
    if (nMeshes == 0) return false;
    int ax = opts.projAxis;
    int ix = (ax == 0) ? 1 : 0;
    int iy = (ax == 2) ? 1 : 2;

    // Compute combined 2D bbox across all meshes
    double minx = 1e30, miny = 1e30, maxx = -1e30, maxy = -1e30;
    for (int m = 0; m < nMeshes; m++) {
        if (meshes[m] && meshes[m]->nvert > 0)
            svg_meshBBox(*meshes[m], ix, iy, minx, miny, maxx, maxy);
    }
    double rangex = maxx - minx; if (rangex < 1e-12) rangex = 1.0;
    double rangey = maxy - miny; if (rangey < 1e-12) rangey = 1.0;

    // Scale to fit within max width×height, then shrink canvas to actual geometry
    double availW = opts.width  - 2 * opts.margin;
    double availH = opts.height - 2 * opts.margin;
    double scale = (availW / rangex < availH / rangey) ? availW / rangex : availH / rangey;
    double ox = opts.margin - minx * scale;
    double oy = opts.margin - miny * scale;

    // Actual canvas size = geometry extent + margin on each side
    int svgW = (int)(rangex * scale + 2 * opts.margin + 0.5);
    int svgH = (int)(rangey * scale + 2 * opts.margin + 0.5);

    FILE* f = fopen(path.c_str(), "w");
    if (!f) { printf("Cannot write %s\n", path.c_str()); return false; }
    fprintf(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\">\n", svgW, svgH);
    fprintf(f, "<rect width=\"100%%%%\" height=\"100%%%%\" fill=\"%s\"/>\n", opts.bgColor);

    int vertOff = 0, edgeOff = 0;
    for (int m = 0; m < nMeshes; m++) {
        if (!meshes[m] || meshes[m]->nvert == 0) continue;
        SVGExportOpts mopts = opts;
        mopts.vertNumOffset = vertOff;
        mopts.edgeNumOffset = edgeOff;
        mopts.height = svgH;  // svgWriteMesh uses opts.height for Y-flip
        svgWriteMesh(f, *meshes[m], mopts, ix, iy, scale, ox, oy, vertOff, edgeOff);
        vertOff += meshes[m]->nvert;
        edgeOff += (meshes[m]->nedge > 0) ? meshes[m]->nedge : meshes[m]->ntri * 3;
    }

    fprintf(f, "</svg>\n");
    fclose(f);
    int totalV = 0, totalE = 0, totalF = 0;
    for (int m = 0; m < nMeshes; m++) { if (meshes[m]) { totalV += meshes[m]->nvert; totalE += meshes[m]->nedge; totalF += meshes[m]->ntri; } }
    printf("Wrote SVG: %s (%d meshes, %d verts, %d edges, %d faces)\n", path.c_str(), nMeshes, totalV, totalE, totalF);
    return true;
}

// Convenience: single mesh
inline bool writeSVG(const std::string& path, const CMesh& mesh, const SVGExportOpts& opts = SVGExportOpts()) {
    const CMesh* meshes[1] = { &mesh };
    return writeSVGMulti(path, meshes, 1, opts);
}

// ========== OBJ+ writer ==========

inline bool writeOBJPlus(const std::string& path, const MeshAnimation& anim) {
    FILE* f = fopen(path.c_str(), "w");
    if (!f) { printf("Cannot write %s\n", path.c_str()); return false; }
    for (int i = 0; i < anim.size(); i++) {
        const CMesh* s = anim.snapshots[i];
        fprintf(f, "# %d %d %d\n", s->nvert, s->nedge, s->ntri);
        if (anim.size() > 1) fprintf(f, "# --- snapshot %d ---\n", i);
        for (int j = 0; j < s->nvert; j++)
            fprintf(f, "v %.6g %.6g %.6g\n", s->verts[j].x, s->verts[j].y, s->verts[j].z);
        for (int j = 0; j < s->nedge; j++)
            fprintf(f, "l %d %d\n", s->edges[j].x + 1, s->edges[j].y + 1);
        for (int j = 0; j < s->ntri; j++)
            fprintf(f, "f %d %d %d\n", s->tris[j].x + 1, s->tris[j].y + 1, s->tris[j].z + 1);
    }
    fclose(f);
    return true;
}

// ========== Universal loader ==========

inline bool loadMeshFile(const std::string& path, MeshAnimation& anim) {
    if (hasSuffix(path, ".npz") || hasSuffix(path, ".npy")) return readNPZ(path, anim);
    return readOBJPlus(path, anim); // .obj and .obj+ both handled by same reader
}

#endif
