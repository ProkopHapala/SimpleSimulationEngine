// Headless test: load OBJ files and export SVG with various options
// Usage: ./test_svg_export <file.obj> [output_dir] [file2.obj for multi-mesh test]

#include <string>
#include <dirent.h>
#include <sys/stat.h>
#include "MeshFileFormats.h"

int main(int argc, char** argv) {
    if (argc < 2) { printf("Usage: %s <file.obj> [output_dir]\n", argv[0]); return 1; }
    std::string inputFile = argv[1];
    std::string outDir = (argc >= 3) ? argv[2] : ".";
    if (outDir[outDir.size()-1] != '/') outDir += "/";

    // Ensure output dir exists
    mkdir(outDir.c_str(), 0755);

    MeshAnimation anim;
    if (!loadMeshFile(inputFile, anim)) { printf("Failed to load %s\n", inputFile.c_str()); return 1; }
    CMesh* snap = anim.currentSnap();
    if (!snap) { printf("No snapshot\n"); return 1; }
    printf("Loaded: %d verts, %d edges, %d faces\n", snap->nvert, snap->nedge, snap->ntri);

    // Derive base name
    std::string base = inputFile;
    size_t slash = base.find_last_of('/');
    if (slash != std::string::npos) base = base.substr(slash + 1);
    size_t dot = base.find_last_of('.');
    if (dot != std::string::npos) base = base.substr(0, dot);

    // Export 1: default (faces on, opacity 0.5, backface cull, edges on, XY plane)
    SVGExportOpts opts1;
    writeSVG(outDir + base + "_default.svg", *snap, opts1);

    // Export 2: faces off, edges only
    SVGExportOpts opts2;
    opts2.showFaces = false;
    opts2.showEdges = true;
    writeSVG(outDir + base + "_edges_only.svg", *snap, opts2);

    // Export 3: faces on, transparent (opacity 0.2), backfaces on
    SVGExportOpts opts3;
    opts3.showFaces = true;
    opts3.faceOpacity = 0.2;
    opts3.showBackfaces = true;
    opts3.showEdges = true;
    writeSVG(outDir + base + "_transparent_bothfaces.svg", *snap, opts3);

    // Export 4: faces on, opaque, backface cull, no edges, XZ plane (side view)
    SVGExportOpts opts4;
    opts4.showFaces = true;
    opts4.faceOpacity = 1.0;
    opts4.showBackfaces = false;
    opts4.showEdges = false;
    opts4.projAxis = 1; // drop Y, show XZ
    writeSVG(outDir + base + "_faces_opaque_side.svg", *snap, opts4);

    // Export 5: everything on, backfaces on, YZ plane (front view)
    SVGExportOpts opts5;
    opts5.showFaces = true;
    opts5.faceOpacity = 0.4;
    opts5.showBackfaces = true;
    opts5.showEdges = true;
    opts5.showPoints = true;
    opts5.projAxis = 0; // drop X, show YZ
    writeSVG(outDir + base + "_all_front.svg", *snap, opts5);

    // Export 6: edges + vertex numbers + edge numbers (topology debug view)
    SVGExportOpts opts6;
    opts6.showFaces = false;
    opts6.showEdges = true;
    opts6.showVertNums = true;
    opts6.showEdgeNums = true;
    opts6.fontSize = 8;
    writeSVG(outDir + base + "_numbers.svg", *snap, opts6);

    // Export 6b: faces + normals
    SVGExportOpts opts6b;
    opts6b.showFaces = true;
    opts6b.faceOpacity = 0.3;
    opts6b.showEdges = false;
    opts6b.showNormals = true;
    opts6b.normalLength = 0.15;
    writeSVG(outDir + base + "_normals.svg", *snap, opts6b);

    // Export 7: multi-mesh — load second file if provided, compose both
    if (argc >= 4) {
        std::string inputFile2 = argv[3];
        MeshAnimation anim2;
        if (loadMeshFile(inputFile2, anim2)) {
            CMesh* snap2 = anim2.currentSnap();
            if (snap2) {
                const CMesh* meshes[2] = { snap, snap2 };
                SVGExportOpts optsM;
                optsM.showFaces = false;
                optsM.showEdges = true;
                optsM.showVertNums = true;
                writeSVGMulti(outDir + base + "_multi.svg", meshes, 2, optsM);
            }
        }
    }

    printf("Done. SVGs in %s\n", outDir.c_str());
    return 0;
}
