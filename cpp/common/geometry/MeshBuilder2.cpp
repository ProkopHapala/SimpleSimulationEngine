#include "MeshBuilder2.h"

#include <algorithm>

#include "arrayAlgs.h"
#include "raytrace.h"
#include "testUtils.h"

namespace Mesh {

void Builder2::clear(){ blocks.clear(); verts.clear(); edges.clear(); tris.clear(); chunks.clear(); strips.clear(); }
void Builder2::printSizes(){printf( "MeshBuilder::printSizes() blocks=%i verts=%i edges=%i tris=%i chunks=%i strips=%i \n", blocks.size(), verts.size(), edges.size(), tris.size(), chunks.size(), strips.size() );}

// int Builder2::getOtherEdgeVert(int ie, int iv){
//     Quat4i& e = edges[ie];
//     return e.x==iv ? e.y : e.x;
// }

bool Builder2::sortPotentialEdgeLoop( int n, Vec2i* edges, int* iverts ){
    // sort edges by shared points so that they for edge loop
    LoopDict point2edge;
    for(int i=0; i<n; i++){
        //edges[i].order();    // for comparison
        printf( "sortPotentialEdgeLoop() edges[%i]=%i,%i \n", i, edges[i].x, edges[i].y );
        Vec2i& e = edges[i];
        point2edge[e.x].add(i);
        point2edge[e.y].add(i);
    }
    for( auto & p : point2edge ){ printf( "sortPotentialEdgeLoop() point2edge[%i]=(%i,%i) \n", p.first, p.second.data[0], p.second.data[1] ); }
    //for( auto& p : point2edge ){ printf( "sortEdgeLoop() point2edge[%i]: (%i,%i) \n", p.first, p.second.data[0], p.second.data[1] ); }
    int iedges[n];
    int oe    = 0; 
    Vec2i& e0 = edges[oe];
    int iv    = e0.y;
    if( iverts ){ 
        printf( "sortPotentialEdgeLoop() i: %i iv=%i ie=%i \n", 0, iv, oe );
        iverts[0] = iv; 
    }
    for(int i=1; i<n; i++){
        auto& v = point2edge[iv];
        int ie;
        if( v.data[0]==oe ){ ie = v.data[1]; }
        else               { ie = v.data[0]; }
        if( ie==-1 ){
            printf( "sortPotentialEdgeLoop() no edge found for iv=%i \n", iv ); 
            return false; 
        } // edge loop ends
        iedges[i] = ie;
        iv = edges[ie].other(iv);
        oe = ie;
        printf( "sortPotentialEdgeLoop() i: %i iv=%i ie=%i \n", i, iv, ie );
        if( iverts ){ iverts[i] = iv; }
    }

    return true;
}

bool Builder2::sortEdgeLoop( int n, int* iedges, int* iverts ){
    // sort edges by shared points so that they for edge loop
    LoopDict point2edge;
    for(int i=0; i<n; i++){
        int ie = iedges[i];
        Quat4i& e = edges[ie];
        point2edge[e.x].add(ie);
        point2edge[e.y].add(ie);
    }
    //for( auto& p : point2edge ){ printf( "sortEdgeLoop() point2edge[%i]: (%i,%i) \n", p.first, p.second.data[0], p.second.data[1] ); }
    int oe    = iedges[0];
    Quat4i& e0 = edges[oe];
    int iv    = e0.y;
    if( iverts ){ iverts[0] = iv; }
    for(int i=1; i<n; i++){
        //point2edge[iv].remove(ie0);
        auto& v = point2edge[iv];
        int ie;
        if( v.data[0]==oe ){ ie = v.data[1]; }
        else               { ie = v.data[0]; }
        if( ie==-1 ){ return false; } // edge loop ends
        iedges[i] = ie;
        oe = ie;
        iv = getOtherEdgeVert(ie,iv);
        if( iverts ){ iverts[i] = iv; }
    }
    return true;
}

int Builder2::findEdgeByVerts_brute( Vec2i verts ){
    // todo: we can speed this up by using unordered_map<uint64_t,int> point2edge;
    verts.order();
    for( int i=0; i<edges.size(); i++ ){
        Vec2i e = edges[i].lo;
        e.order();
        if( e==verts ){ return i; }
    }
    return -1;
}

int Builder2::findEdgeByVerts_map( const Vec2i verts ){
    uint64_t key = symetric_id( verts );
    auto it = vert2edge.find(key);
    if( it==vert2edge.end() ){ return -1; }
    return it->second;
}
int Builder2::findEdgeByVerts( const Vec2i verts ){
    int ie= (use_vert2edge) ? findEdgeByVerts_map(verts) : findEdgeByVerts_brute(verts);
    if( bExitError) if( ie==-1 ){ printf( "ERROR: findEdgeByVerts() edge(%i,%i) not found => exit() \n", verts.x, verts.y ); exit(0); }
    return ie;
}

int Builder2::findOrAddEdges( const Vec2i verts, int t, int t2 ){
    int ie = findEdgeByVerts( verts );
    if( ie==-1 ){ ie = edge(verts.x, verts.y, t, t2); }
    return ie;
}

void Builder2::buildVerts2Edge(){
    vert2edge.clear();
    for( int i=0; i<edges.size(); i++ ){
        Vec2i e = edges[i].lo;
        uint64_t key = symetric_id( e );
        vert2edge[key] = i;
    }
}

int Builder2::polygonChunk( int n, int* iedges, const int* ivs, bool bPolygonToTris ){
    int i0 = strips.size();
    int ich = chunk( Quat4i{i0, i0+n, n, (int)ChunkType::face} );
    for(int i=0; i<n; i++){ strips.push_back( ivs[i]    ); } // store verts
    for(int i=0; i<n; i++){ strips.push_back( iedges[i] ); } // store edges
    if( bPolygonToTris){ polygonToTris( ich ); }
    return ich;
}

int Builder2::polygon( int n, int* iedges ){
    int ivs[n];
    if( !sortEdgeLoop( n, iedges, ivs ) ) return false;
    return polygonChunk( n, iedges, ivs, bPolygonToTris );
}

int Builder2::polygonToTris( int i ){
    Quat4i ch = chunks[i];
    int n    = ch.z;
    int i0   = ch.x;
    int* ivs = &strips[i0];
    //int* ies = &strips[i0+n];
    //int i1 = strips.size();
    //int i2 = i1 + n-2;
    //chunks[i] = Quat4i{i1, i2, n-2, (int)ChunkType::face};
    for(int j=1; j<n-1; j++){
        tris.push_back( Quat4i{ivs[0], ivs[j], ivs[j+1], i} );
    }
    return n;
}

Vec3d Builder2::polygonNormal( int ich ){
    Quat4i ch = chunks[ich];
    int n   = ch.z;
    int i0  = ch.x;
    int* ivs = &strips[i0];
    //printf("Builder2::polygonNormal(ich=%d): type: %d n=%d i0=%d\n", ich, ch.w, n, i0);
    Vec3d nrm = Vec3dZero;
    Vec3d a = verts[ivs[0]].pos;
    Vec3d b = verts[ivs[1]].pos;
    for(int j=1; j<n-1; j++){
        Vec3d c = verts[ivs[j+1]].pos;
        Vec3d nr; nr.set_cross( c-b, c-a);
        //nr.normalize();
        //printf("polygonNormal[%i]: %d ivs(%d, %d, %d) Normal: (%g,%g,%g) A(%g,%g,%g) B(%g,%g,%g) C(%g,%g,%g)\n", j, ivs[0], ivs[j], ivs[j+1], nr.x,nr.y,nr.z, a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z);
        nrm.add(nr);
        a = b;
        b = c;
    }
    nrm.normalize();
    //printf("polygonNormal(ich=%d): Final normalized normal: (%g,%g,%g)\n", ich, nrm.x,nrm.y,nrm.z);
    return nrm;
}

int Builder2::findMostFacingNormal(Vec3d hray, int nch, int* chs, double cosMin, bool bTwoSide ){
    //printf("Builder2::findMostFacingNormal() nch=%i cosMin=%g bTwoSide=%i\n", nch, cosMin, bTwoSide );
    int ibest=-1;
    double cmax = -1.0;
    for(int i=0; i<nch; i++){
        int ich = chs[i];
        //printf("findMostFacingNormal[%3i]: ich %3i \n", i, ich );
        Vec3d nr = polygonNormal(ich); // This will now print its own debug info
        double c = hray.dot(nr);
        //printf("  [DEBUG] findMostFacingNormal: hray=(%g,%g,%g), normal=(%g,%g,%g), dot_product (c)=%g\n", hray.x,hray.y,hray.z, nr.x,nr.y,nr.z, c);
        //printf("  [DEBUG] findMostFacingNormal: cosMin=%g, bTwoSide=%d, current_cmax=%g\n", cosMin, bTwoSide, cmax);
        if(bTwoSide) c=fabs(c);
        //printf("  [DEBUG] findMostFacingNormal: After bTwoSide adjustment, c=%g\n", c);
        if( c>cosMin && c>cmax){ ibest=chs[i]; cmax=c; }
    }
    if(bExitError){
        if(ibest<0){ printf("ERROR in Builder2::findMostFacingNormal(): ibest %3i \n", ibest ); 
            for(int i=0; i<nch; i++){
                int ich = chs[i];
                Quat4i ch = chunks[ich];
                Vec3d nr = polygonNormal(ich);
                //printf("chunk[%i] id %i ch(%3i,%3i,%3i,%3i) nr(%g,%g,%g) \n", i, ich,    ch.x, ch.y, ch.z, ch.w, nr.x, nr.y, nr.z );
            }
            exit(0); 
        }
    }

    return ibest;
}

int Builder2::findMostFacingNormal(Vec3d hray, Vec2i chrange, double cosMin, bool bTwoSide ){
    int nch = chrange.y - chrange.x;
    //printf("Builder2::findMostFacingNormal() chrange(%i,%i) nch=%i \n", chrange.x, chrange.y, nch );
    int chs[nch];
    for(int i=0; i<nch; i++){ chs[i] = chrange.x+i; }
    return findMostFacingNormal(hray, nch, chs, cosMin, bTwoSide );
}

Vec2i Builder2::addVerts( int n, const Vec3d* ps ){
    int i0 = verts.size();
    for(int i=0; i<n; i++){ vert(ps[i]); }
    return Vec2i{i0, (int)verts.size()-1};
}

Vec2i Builder2::addEdges( int n, const Vec2i* iedges, const int* types, const int* types2, int iv0 ){
    //int ich0 = chunks.size();
    int ie0 = edges.size();
    for (int i = 0; i < n; ++i) {
        int t1=(types  )? types [i]:-1;
        int t2=(types2 )? types2[i]:-1;
        edge(iedges[i].x+iv0,iedges[i].y+iv0,t1,t2);
        //edges.push_back(Quat4i{iedges[i].x,iedges[i].y,t1,t2});
    }
    return Vec2i{ie0, (int)edges.size()-1};
};

Vec2i Builder2::addFaces( int nf, const int* nVerts, const int* verts, bool bPolygonToTris, int iv0 ){
    //std::vector<int> plane_chunks;
    //const int* iverts = Solids::Octahedron_planeVs;
    //printf("Adding octahedron planes as polygons...\n");
    int ich0 = chunks.size();
    for (int i = 0; i < nf; ++i) {
        int nvi = nVerts[i];
        int verts_[nvi]; for(int i=0; i<nvi; i++){ verts_[i] = verts[i]+iv0; }
        int iedges[nvi];
        for (int j = 0; j < nvi; ++j) {
            Vec2i vpair = {verts_[j], verts_[(j + 1) % nvi]};
            //printf( "Builder2::addFaces() vpair(%3i,%3i) iv0 %3i\n", vpair.a, vpair.b, iv0 );
            iedges[j]   = findEdgeByVerts(vpair);
            if (iedges[j] == -1) {printf("ERROR: -extrude_octahedron could not find edge between %d and %d. This should not happen for a well-defined CMesh.\n", vpair.a, vpair.b);exit(0); }
        }
        int chunk_id = polygonChunk(nvi, iedges, verts_, true); // true to triangulate
        //plane_chunks.push_back(chunk_id);
        //printf("  Added plane %i with %i vertices as chunk %i\n", i, nvi, chunk_id);
        verts += nvi;
    }
    return Vec2i{ich0, (int)chunks.size()-1};
};

Quat4i Builder2::addCMesh(const CMesh& cmesh, bool bFaces, Vec3d p0, Vec3d sc, Mat3d* rot, int edge_type ) {
    //b.clear();
    // Add vertices
    printf("Builder2::addCMesh() cmesh.nvert=%i cmesh.nedge=%i cmesh.nfaces=%i p0(%g,%g,%g) sc(%g,%g,%g)\n", cmesh.nvert, cmesh.nedge, cmesh.nfaces, p0.x,p0.y,p0.z, sc.x,sc.y,sc.z);
    int iv0  = verts.size();
    int ie0  = edges.size();
    int it0  = tris.size();
    int ich0 = chunks.size();
    for (int i = 0; i < cmesh.nvert; ++i) {
        Vec3d p_ = cmesh.verts[i];
        Vec3d p = p_*sc + p0;
        if(rot){ p = rot->dot(p); }
        int iv=vert(p);
        printf("Builder2::addCMesh() i=%3i iv=%3i p(%g,%g,%g) p_(%g,%g,%g)\n", i,iv,p.x,p.y,p.z, p_.x,p_.y,p_.z );
    }
    // Add edges
    for (int i = 0; i < cmesh.nedge; ++i) {
        edge(cmesh.edges[i].a + iv0, cmesh.edges[i].b + iv0, edge_type );
    }
    buildVerts2Edge(); // Important for findEdgeByVerts to be fast
    // Add faces (as polygons)
    if (bFaces) {
        const int* ifaces = cmesh.faces;
        for (int i = 0; i < cmesh.nfaces; ++i) {
            int n = cmesh.ngons[i];
            int iverts[n];
            int iedges[n];
            printf("Builder2::addCMesh() iface=%i n=%i\n", i,n );
            for(int j=0; j<n; ++j) {
                iverts[j] = ifaces[j] + iv0;
                {
                    Vec3d p = verts[iverts[j]].pos;
                    printf("Builder2::addCMesh() i iverts[%i]=%i p(%g,%g,%g)\n", j, iverts[j], p.x,p.y,p.z);
                }
            }
            for (int j = 0; j < n; ++j) {
                Vec2i vpair = {iverts[j], iverts[(j + 1) % n]};
                iedges[j] = findEdgeByVerts(vpair);
                if (iedges[j] == -1) {
                    // This should not happen if CMesh is well-defined.
                    printf("WARNING: CMesh2Builder could not find edge between %d and %d. Adding it.\n", vpair.a, vpair.b);
                    iedges[j] = edge(vpair.a, vpair.b);
                }
            }
            polygonChunk(n, iedges, iverts, false); // true for triangulating
            ifaces += n;
        }
    }
    return Quat4i{iv0,ie0,it0,ich0};
}

int Builder2::selectionToFace(){
    //printf( "selToFace() n=%i \n", selection.size() );
    int n = selset.size();
    if(n>ngon_max){
        printf( "WARRNING: selectionToFace() n(%i) > ngon_max(%i) => return \n", n, ngon_max );
        return 0;
    }
    selection.clear();
    for( int i : selset ){ selection.push_back(i);}
    polygon( n, selection.data() );
    return 0;
}

int Builder2::clearSelection(){
    int n = selection.size();
    selection.clear();
    selset.clear();
    return n;
}

int Builder2::pickVertex( const Vec3d& ray0, const Vec3d& hRay, double R ){
    printf( "pickVertex() ray0(%g,%g,%g) hRay(%g,%g,%g) R=%g \n", ray0.x, ray0.y, ray0.z, hRay.x, hRay.y, hRay.z, R );
    double tmin =  1e+300;
    int imin    = -1;
    for(int i=0; i<verts.size(); i++){
        //if(ignore)if(ignore[i])continue;
        double ti = raySphere( ray0, hRay, R, verts[i].pos );
        if(ti<tmin){ imin=i; tmin=ti; }
    }
    printf( "pickVertex() DONE ipick %i \n", imin );
    return imin;
}

int Builder2::pickEdge( const Vec3d& ro, const Vec3d& rh, double Rmax ){
    //printf( "pickEdge() ro(%g, %g, %g) rh(%g, %g, %g) \n", ro.x, ro.y, ro.z, rh.x, rh.y, rh.z );
    return rayPickBond( ro, rh, edges.size(), [&](int ib,Vec3d&pa,Vec3d&pb){ 
        Vec2i b = edges[ib].lo; pa=verts[b.i].pos; pb=verts[b.j].pos;
    }, Rmax, false );
}

int Builder2::toggleSelSet(  int i ){
    if( selset.find(i)!=selset.end() ){ selset.erase(i); return -1; }else{ selset.insert(i); return 1; }; return 0;
}

int Builder2::pickTriangle( const Vec3d& ro, const Vec3d& rh, bool bReturnFace ){
    //printf( "pickTriangle() ro(%g,%g,%g) rh(%g,%g,%g) n", ro.x, ro.y, ro.z, rh.x, rh.y, rh.z );
    Vec3d hX,hY;
    rh.getSomeOrtho(hX,hY);
    //printf( "pickTriangle() hX(%g,%g,%g) hY(%g,%g,%g) \n", hX.x, hX.y, hX.z, hY.x, hY.y, hY.z );
    double Lmin = 1e+300;
    int    imin = -1;
    for( int i=0; i<tris.size(); i++ ){
        Vec3d normal;
        Quat4i t = tris[i];
        double L = rayTriangle2( ro, rh, hX, hY, verts[t.x].pos, verts[t.y].pos, verts[t.z].pos, normal );
        //printf( "pickTriangle() i=%i L=%g imin=%i Lmin=%g \n", i, L, imin, Lmin );
        if( L<Lmin ){ 
            //printf( "pickTriangle() i=%i L=%g imin=%i Lmin=%g \n", i, L, imin, Lmin );
            Lmin=L; imin=i; 
        }
    }
    if(bReturnFace && (imin>=0) ){ return tris[imin].w; }
    return imin;
}

int Builder2::pickEdgeSelect( const Vec3d& ro, const Vec3d& rh, double Rmax ){ int i=pickEdge( ro, rh, Rmax ); if(i>=0){  toggleSelSet( i ); } return i; }

int Builder2::pickSelect( const Vec3d& ro, const Vec3d& rh, double Rmax ){
    //printf( "pickSelect() selection_mode %i  \n", selection_mode);
    int ipick=-1;
    switch( (SelectionMode)selection_mode ){
        case SelectionMode::vert: ipick= pickVertex    ( ro, rh, Rmax ); break;
        case SelectionMode::edge: ipick= pickEdgeSelect( ro, rh, Rmax ); break;
        case SelectionMode::face: ipick= pickTriangle  ( ro, rh, true ); break;
    }
    if( bAdditiveSelect && (ipick>=0) ){ 
        // is ipick in selset?
        if( selset.find(ipick)!=selset.end() ){ 
            selset.erase(ipick); 
            selection.erase( std::remove( selection.begin(), selection.end(), ipick ), selection.end() );
        }else{ 
            selset.insert(ipick);
            selection.push_back(ipick); 
        }
    }
    return -1;
}

int Builder2::selectRectEdge( const Vec3d& p0, const Vec3d& p1, const Mat3d& rot ){ 
    //printf( "Mesh::Builder2::selectRectEdge() p0(%g,%g,%g) p1(%g,%g,%g) \n", p0.x,p0.y,p0.z, p1.x,p1.y,p1.z );
    Vec3d Tp0,Tp1;
    //Mat3d rot = (Mat3d)cam.rot;
    rot.dot_to(p0,Tp0);
    rot.dot_to(p1,Tp1);
    _order(Tp0.x,Tp1.x);
    _order(Tp0.y,Tp1.y);
    Tp0.z=-1e+300;
    Tp1.z=+1e+300;
    int nfound=0;
    for(int i=0; i<edges.size(); i++ ){
        Vec2i b = edges[i].lo;
        Vec3d pa,pb;
        rot.dot_to( verts[b.i].pos,pa);
        rot.dot_to( verts[b.j].pos,pb);
        if( pa.isBetween(Tp0,Tp1) && pb.isBetween(Tp0,Tp1) ){
            selset.insert( i );
            nfound++;
        }
    }
    return nfound;
}

int Builder2::selectRectVert( const Vec3d& p0, const Vec3d& p1, const Mat3d& rot ){ 
    //printf( "Mesh::Builder2::selectRectEdge() p0(%g,%g,%g) p1(%g,%g,%g) \n", p0.x,p0.y,p0.z, p1.x,p1.y,p1.z );
    Vec3d Tp0,Tp1;
    //Mat3d rot = (Mat3d)cam.rot;
    rot.dot_to(p0,Tp0);
    rot.dot_to(p1,Tp1);
    _order(Tp0.x,Tp1.x);
    _order(Tp0.y,Tp1.y);
    Tp0.z=-1e+300;
    Tp1.z=+1e+300;
    int nfound=0;
    for(int i=0; i<verts.size(); i++ ){
        Vec3d p;
        rot.dot_to( verts[i].pos,p);
        if( p.isBetween(Tp0,Tp1) ){
            selset.insert( i );
            nfound++;
        }
    }
    return nfound;
}

int Builder2::selectRect( const Vec3d& p0, const Vec3d& p1, const Mat3d& rot  ){
    selset.clear();
    //printf( "pickSelect() selection_mode %i  \n", selection_mode);
    switch( (SelectionMode)selection_mode ){
        case SelectionMode::edge: return selectRectEdge( p0, p1, rot ); break;
        case SelectionMode::vert: return selectRectVert( p0, p1, rot ); break;
    }
    return 0;
}

int Builder2::select_in_sphere( const Vec3d& p0, double r ){
    int n=0;
    for(int i=0;i<verts.size();i++){
        Vec3d  d =  verts[i].pos-p0;
        if( d.norm2() > r*r ) continue; // check bounds along the axis
        selection.push_back(i);
        n++;
    }
    return n;
}

int Builder2::select_in_cylinder( const Vec3d& p0, const Vec3d& fw, double r, double l ){
    int n=0;
    for(int i=0;i<verts.size();i++){
        Vec3d  d =  verts[i].pos-p0;
        double cH = fw.dot(d);
        if( (cH<0) || (cH>l) ) continue; // check bounds along the axis
        double cR = d.norm2() - cH*cH; //
        if( cR>r*r ) continue;         // check bounds in the plane
        selection.push_back(i);
        n++;
    }
    return n;
}

int Builder2::select_in_box( const Vec3d& p0, const Vec3d& fw, const Vec3d& up, const Vec3d& Lmin, const Vec3d& Lmax ){
    Vec3d lf = cross(fw,up);
    lf.normalize();
    int n=0;
    for(int i=0;i<verts.size();i++){
        Vec3d  d  =  verts[i].pos-p0;
        Vec3d  Td{ lf.dot(d), up.dot(d), fw.dot(d) }; 
        if( Td.isBetween( Lmin, Lmax ) ){ 
            selection.push_back(i);
            n++;
        }
    }
    return n;
}

void Builder2::makeSelectrionUnique(){
    std::sort( selection.begin(), selection.end() );
    auto it = std::unique( selection.begin(), selection.end() );
    selection.resize( std::distance(selection.begin(), it) );
}

int Builder2::findClosestVert(const Vec3d& p0,int i0,int n){
    if(n==-1) n=verts.size();
    double r2min=1e+300;
    int imin=-1;
    for(int j=0;j<n;j++){ double r2=(verts[i0+j].pos-p0).norm2(); if(r2<r2min){r2min=r2; imin=j;} }
    return imin;
}

int Builder2::findVert(const Vec3d& p0, double Rmax, int n, int* sel ){
    //printf( "Mesh::Builder2::findVert() p(%g,%g,%g) Rmax: %g  \n", p0.x,p0.y,p0.z, Rmax );
    if(n==-1){ n  =selection.size(); sel=selection.data(); }
    double r2min=Rmax;
    int imin=-1;
    for(int ii=0;ii<n;ii++){ 
        int i = sel[ii];
        double r2=(verts[i].pos-p0).norm2(); 
        if(r2<r2min){r2min=r2; imin=i;} 
    }
    // if( bExitError && (imin<0)){ 
    //     printf( "ERROR in findVert() imin=%i cannot find vert near p(%g,%g,%g) within Rmax: %g => Exit() \n", imin, p0.x,p0.y,p0.z, Rmax );
    //     printSelectedVerts(); 
    //     exit(0); 
    // };
    return imin;
}

    int Builder2::selectVertsAlongLine( Vec3d p0, Vec3d p1, double r, bool bSort ){
        std::vector<int> sel; // we need to make local selection so we can add to global selection
        sel.clear();
        Vec3d ax = p1-p0;
        double l = ax.normalize();
        double R2 = r*r;
        std::vector<double> xalong;
        for(int i=0; i<verts.size(); i++){
            Vec3d d = verts[i].pos-p0;
            double x = ax.dot(d); // axial component
            if( (x<0)|| (x>l) ) continue;
            d.add_mul( ax, -x ); // orthogonal component
            double r2 = d.norm2(); 
            if( r2>R2 ) continue;
            sel.push_back(i);
            if(bSort) xalong.push_back(x);
        }
        if(bSort){  // we need to sort selection by xalong
            sortArrayByAnother( sel.size(), sel.data(), xalong.data() );
        }
        selection.insert( selection.end(), sel.begin(), sel.end() );
        return sel.size();
    };

    int Builder2::selectVertsAlongPolyline( double r, bool bSort ){
        std::vector<int> sel = selection; // we need to backup selection so we can accumulate all found points into global selection
        for(int i=1; i<sel.size(); i++){
            Vec3d p0 = verts[sel[i-1]].pos;
            Vec3d p1 = verts[sel[i  ]].pos;
            selectVertsAlongLine( p0, p1, r, bSort );
        }
        return selection.size();
    }

    int Builder2::plateBetweenVertStrips( int n, int* ivs1, int* ivs2, int nsub ){
        if( n<2 ) return 0;
        for( int i=0; i<n; i++){
            edge( ivs1[i], ivs2[i] ); // this is just quick for debugig
            /// TODO: later we will do subdivision and faces / triangles creation
        }
        return 0;
    }

    int Builder2::plateBetweenEdges( int nsub, double r, bool bSort ){
        printf( "plateBetweenEdges() nsub=%i  r=%f  bSort=%i  sel.size=%i\n", nsub, r, bSort, selection.size() );
        std::vector<int> sel = selection;
        std::vector<int> strip1;
        std::vector<int> strip2;
        if( sel.size()==3 ){ // corner ( triangle )
            selection.clear(); Builder2::selectVertsAlongLine( verts[sel[1]].pos, verts[sel[0]].pos, r, bSort ); strip1 = selection;
            selection.clear(); Builder2::selectVertsAlongLine( verts[sel[1]].pos, verts[sel[2]].pos, r, bSort ); strip2 = selection;
            int n = _min( strip1.size(), strip2.size() );
            plateBetweenVertStrips( n-1, strip1.data()+1, strip2.data()+1, nsub ); 
            /// TODO: we ignore the first vert so we can use Quad plateBetweenVertStrips, later we need to make triangle version
            selection.insert( selection.end(), strip1.begin()+1, strip1.end() );
        } else if( sel.size()==4 ){ // quad
            selection.clear(); Builder2::selectVertsAlongLine( verts[sel[0]].pos, verts[sel[1]].pos, r, bSort ); strip1 = selection;
            selection.clear(); Builder2::selectVertsAlongLine( verts[sel[3]].pos, verts[sel[2]].pos, r, bSort ); strip2 = selection;
            int n = _min( strip1.size(), strip2.size() );
            plateBetweenVertStrips( n, strip1.data(), strip2.data(), nsub );
            selection.insert( selection.end(), strip1.begin(), strip1.end() );
        }else{
            printf( "ERROR in plateBetweenEdges() selection.size()=%i \n", selection.size() );
            printSelectedVerts();
            //exit(0);
        }
        return 0;
    }

    Vec2i Builder2::conect_vertex( int iv, int stickType, int n, int* iverts ){
        int ie = edges.size();
        for(int i=0;i<n;i++){ /*if (iv == iverts[i]) continue;*/ edge(iv,iverts[i],stickType); } // The check `if (iv == iverts[i]) continue;` is removed here as the caller (make_anchor_point) now ensures the anchor point itself is not in the selection list.
        return Vec2i{ie,ie+n-1};
    }

    int Builder2::conected_vertex( const Vec3d& p, int stickType, int n, int* iverts ){
        int iv = vert( p );
        conect_vertex( iv, stickType, n, iverts );
        return iv;
    };

    int Builder2::make_anchor_point( const Vec3d& p, int stickType, double Rcolapse, double r, const Vec3d* fw, double l ){
        selection.clear();
        if( fw ){ select_in_cylinder( p, *fw, r, l ); }
        else    { select_in_sphere(p, r);             }
        int iv=-1;
        if(Rcolapse>0){ iv = findVert(p, Rcolapse); }  // use existing vert if found
        if(iv<0      ){ iv = vert( p );             }  // create new vert if not found

    // Root cause fix: Remove the anchor point itself from the selection list
    // if it was included by select_in_sphere/cylinder (which it will be if it's an existing vertex)
    auto it = std::remove(selection.begin(), selection.end(), iv);
    selection.erase(it, selection.end());
    conect_vertex( iv, stickType, selection.size(), selection.data() );
        return iv;
    }

    int Builder2::make_anchor_points( int nv, Vec3d* vs, int* ivrts, int anchorType, double Rcolapse, double r, const Vec3d* fw, double l ){
        int nv0 = verts.size();
        for(int i=0; i<nv; i++){ ivrts[i]=make_anchor_point( vs[i], anchorType, Rcolapse, r, fw, l ); }
        return verts.size()-nv0;
    }



int Builder2::bondsBetweenVertRanges( Vec2i v1s, Vec2i v2s, double Rmax, int et ){
    double R2max = Rmax*Rmax;
    int n1 = v1s.b-v1s.a;
    int n2 = v2s.b-v2s.a;
    int nb =0;
    for(int i=0; i<n1; i++){
        Vec3d& p1 = verts[v1s.a+i].pos;
        for(int j=0; j<n2; j++){
            Vec3d& p2 = verts[v2s.a+j].pos;
            double r2 = (p1-p2).norm2();
            //printf( "bondsBetweenVertRanges[%i,%i] r=%g Rmax=%g \n", i, j, sqrt(r2), Rmax );
            if( r2<R2max ){ 
                edge( v1s.a+i, v2s.a+j, et ); 
                nb++;
            }
        }
    }
    return nb;
}

int Builder2::vstrip(Vec3d p0, Vec3d p1, int n, int et ){
    Vec3d d=p1-p0; 
    d.mul(1./n);
    Vec3d p = p0;
    int i0 = verts.size();
    //printf("vstri(%i){", i0 );
    for(int ii=0; ii<n+1; ii++){
        int i = vert(p);  // printf("%i ", i );
        if((et>-1)&&(ii>0))edge(i-1,i,et);  // long
        p.add(d);
    }
    //printf("}END\n" );
    return i0;
}

int Builder2::fstrip( int ip0, int ip1, int n, int ft, Vec2i et ){ 
    int i0 = tris.size();
    //printf("fstrip(%i,%i)\n", ip0, ip1 );
    for(int ii=0; ii<n+1; ii++){
        int it = ip0*1000+ii;
        //printf("it %i \n", it );
        int i = ip0+ii;
        int j = ip1+ii;
        if(et.x>-1)edge(j,i,et.x);         // perp
        //if(et.x>-1)edge(j,i,it);
        if(ii>0){
            if(et.y>-1)edge(i,j-1,et.y);   // diag
            //if(et.y>-1)edge(j-1,i,it);
            if(ft>-1){ 
                //tri(j-1,i,j,ft); 
                //printf("tri(%i,%i,%i)\n", i,j,j-1 );
                //printf("tri(%i,%i,%i)\n", i,j,j-1 );
                tri(j-1,i  ,j  ,ft);  //printf("tri(%i,%i,%i)\n", i,j,j-1 );
                //tri(j-1,i  ,i-1,ft);
                tri(i  ,j-1,i-1,ft);
                //tri(ip0,ip1,ip1-1,ft); printf("tri(%i,%i,%i)\n", ip0,ip1,ip1-1 );
                //tri(ip0,ip1,ip1-1, it ); printf("tri(%i,%i,%i)\n", ip0,ip1,ip1-1 );
            }
        }
    }
    return i0;
}

void Builder2::addPointCross( const Vec3d& p, double d ){
    stick( p + Vec3d{d,0,0}, p + Vec3d{-d, 0, 0} );
    stick( p + Vec3d{0,d,0}, p + Vec3d{ 0,-d, 0} );
    stick( p + Vec3d{0,0,d}, p + Vec3d{ 0, 0,-d} );
}

void Builder2::addArrow( const Vec3d& p1, const Vec3d& p2, double d ){
    stick( p1, p2 );
    Vec3d dir = p2 - p1;
    dir.normalize();
    Vec3d up, right;
    dir.getSomeOrtho(up, right);
    double head_len = d * 2.0;
    double head_wid = d;
    Vec3d head_base = p2 - dir * head_len;
    stick( p2, head_base + up * head_wid );
    stick( p2, head_base - up * head_wid );
    stick( p2, head_base + right * head_wid );
    stick( p2, head_base - right * head_wid );
}

/* Too Complicated logic inside
int vstrip(Vec3d p0, Vec3d p1, int n, Quat4i t=Quat4i{-1,-1,-1,-1}){
    Vec3d d=p1-p0; 
    d.mul(1./n);
    Vec3d p = p0;
    int oi=-1; 
    int i0 = verts.size();
    for(int ii=0; ii<n+1; ii++){
        int i = vert(p);
        int j = i-n;
        if(t.y>-1)edge(j,i,t.y);         // perp
        if(oi>0){
            if(t.x>-1)edge(i-1,i,t.x);   // long
            if(t.y>-1)edge(j-n,i,t.z);   // diag
            if(t.w>-1){ tri(j-1,j,i,t.w); tri(i,i-1,j,t.w); }
        p.add(d);
        oi=i;
    }
    return i0;
}
*/

void Builder2::box( Vec3d p, Vec3d ls, Mat3d rot ){
    int i000 = vert(p+rot.dotT({ ls.x, ls.y, ls.z}) );
    int i100 = vert(p+rot.dotT({-ls.x, ls.y, ls.z}) );
    int i010 = vert(p+rot.dotT({ ls.x,-ls.y, ls.z}) );
    int i110 = vert(p+rot.dotT({-ls.x,-ls.y, ls.z}) );
    int i001 = vert(p+rot.dotT({ ls.x, ls.y,-ls.z}) );
    int i101 = vert(p+rot.dotT({-ls.x, ls.y,-ls.z}) );
    int i011 = vert(p+rot.dotT({ ls.x,-ls.y,-ls.z}) );
    int i111 = vert(p+rot.dotT({-ls.x,-ls.y,-ls.z}) );
    // x
    edge(i000,i100);
    edge(i001,i101);
    edge(i010,i110);
    edge(i011,i111);
    // y
    edge(i000,i010);
    edge(i001,i011);
    edge(i100,i110);
    edge(i101,i111);
    // z
    edge(i000,i001);
    edge(i010,i011);
    edge(i100,i101);
    edge(i110,i111);

    /*
    double as[4]{-1.0,-1.0, 1.0,1.0};
    double bs[4]{-1.0, 1.0,-1.0,1.0};
    for(int ii=0; ii<4; ii++){
        Vec3d p;
        int i,j;
        // z
        rot.dot_to_T( {ls.a*as[ii],ls.b*bs[ii], 0.0}, p );
        i=vert( p-rot.c*ls.z ); j=vert( p+rot.c*ls.z );
        edge(i,j);
        // y
        rot.dot_to_T( {ls.a*as[ii],0.0,ls.c*bs[ii]}, p );
        i=vert( p-rot.b*ls.y ); j=vert( p+rot.b*ls.y );
        edge(i,j);
        // x
        rot.dot_to_T( {0.0,ls.b*as[ii],ls.c*bs[ii]}, p );
        i=vert( p-rot.a*ls.x ); j=vert( p+rot.a*ls.x );
        edge(i,j);
    }
    */
}

void Builder2::frustrumFace( const Vec3d& p0, const Mat3d& rot, double La, double Lb, double h, double Lbh, double Lah  ){
    printf("frustrumFace: La: %f Lb: %f h: %f Lbh: %f Lah: %f \n", La,Lb,h,Lbh);
    int i0 = verts.size();
    Vec3d p;
    p=p0+rot.a* La       + rot.c*0.1; vert( p+rot.b*Lb       ); vert( p-rot.b*Lb       );
    p=p0+rot.a* (La-Lah) + rot.c*h;   vert( p+rot.b*(Lb-Lbh) ); vert( p-rot.b*(Lb-Lbh) );
    p=p0+rot.a*-(La-Lah) + rot.c*h;   vert( p+rot.b*(Lb-Lbh) ); vert( p-rot.b*(Lb-Lbh) );
    p=p0+rot.a*-La       + rot.c*0.1; vert( p+rot.b*Lb       ); vert( p-rot.b*Lb       );

    int i,j;
    for(int ii=0; ii<4; ii++){
        int i1=i0+ii*2;
        edge(i1,i1+1);
        if(ii<3){
            edge(i1  ,i1+2);
            edge(i1+1,i1+3);
        }
    }
}

void Builder2::snapBoxFace( const Vec3d& p0, const Mat3d& rot, double La, double Lb ){
    //printf("snapBoxFace: La: %f Lb: %f h: %f Lbh: %f Lah: %f \n", La,Lb,h,Lbh);
    int i00=findVert( p0-rot.a*La-rot.b*Lb, R_snapVert );
    int i01=findVert( p0-rot.a*La+rot.b*Lb, R_snapVert );
    int i10=findVert( p0+rot.a*La-rot.b*Lb, R_snapVert );
    int i11=findVert( p0+rot.a*La+rot.b*Lb, R_snapVert );
    //printf("snapBoxFace: ivs: %i %i %i %i \n", i00,i01,i10,i11 );
    int ies[4];
    ies[0]=findOrAddEdges( {i00,i01} );
    ies[1]=findOrAddEdges( {i01,i11} );
    ies[2]=findOrAddEdges( {i11,i10} );
    ies[3]=findOrAddEdges( {i10,i00} );
    //printf("snapBoxFace: ies: %i %i %i %i \n", ies[0],ies[1],ies[2],ies[3] );
    polygon(4, ies); // top
}

void Builder2::snapFrustrumFace( const Vec3d& p0, const Mat3d& rot, double La, double Lb, double h, double Lbh, double Lah, bool bFace ){
    //printf("frustrumFace: La: %f Lb: %f h: %f Lbh: %f Lah: %f \n", La,Lb,h,Lbh);
    int iv[8];
    Vec3d p;
    p=p0+rot.a* La                ; iv[0]=findVert( p+rot.b*Lb, R_snapVert ); iv[1]= findVert( p-rot.b*Lb, R_snapVert );
    p=p0+rot.a* (La-Lah) + rot.c*h; iv[2]=    vert( p+rot.b*(Lb-Lbh) );       iv[3]=     vert( p-rot.b*(Lb-Lbh) );
    p=p0+rot.a*-(La-Lah) + rot.c*h; iv[4]=    vert( p+rot.b*(Lb-Lbh) );       iv[5]=     vert( p-rot.b*(Lb-Lbh) );
    p=p0+rot.a*-La                ; iv[6]=findVert( p+rot.b*Lb, R_snapVert ); iv[7]= findVert( p-rot.b*Lb, R_snapVert );
    int ie_v01 = edge(iv[0], iv[2]); int ie_v02 = edge(iv[1], iv[3]); //  |   |
    int ie_h1  = edge(iv[2], iv[3]);                                  //   ---
    int ie_v11 = edge(iv[2], iv[4]); int ie_v12 = edge(iv[3], iv[5]); //  |   |
    int ie_h2  = edge(iv[4], iv[5]);                                  //   ---
    int ie_v21 = edge(iv[4], iv[6]); int ie_v22 = edge(iv[5], iv[7]); //  |   |
    if(bFace){
        { int ies[4]{ie_h1, ie_v11, ie_h2, ie_v12};  polygon(4, ies); } // Front face (quad)
        { int ies[4]{ findEdgeByVerts({iv[0], iv[1]}),  ie_v01, ie_h1,  ie_v02  };  polygon(4, ies); } // top
        { int ies[4]{ findEdgeByVerts({iv[0], iv[6]}),  ie_v01, ie_v11, ie_v21  };  polygon(4, ies); } // left
        { int ies[4]{ findEdgeByVerts({iv[1], iv[7]}),  ie_v02, ie_v12, ie_v22  };  polygon(4, ies); } // right            
        { int ies[4]{ findEdgeByVerts({iv[6], iv[7]}),  ie_v21, ie_h2,  ie_v22  };  polygon(4, ies); } // botton
    }
}

void Builder2::prismFace( const Vec3d& p0, const Mat3d& rot, double La, double Lb, double h, double Lbh ){
    //printf("prismFace: La: %f Lb: %f h: %f Lbh: %f  \n", La,Lb,h,Lbh);
    // ToDo: We must find the points on cube where to attach the edges, not to create a new vertexes.
    int i0 = verts.size();
    Vec3d p;
    p=p0+rot.a* La + rot.c*0.1; vert( p+rot.b*Lb       ); vert( p-rot.b*Lb       );
    p=p0+rot.c* h;              vert( p+rot.b*(Lb-Lbh) ); vert( p-rot.b*(Lb-Lbh) );
    p=p0+rot.a*-La + rot.c*0.1; vert( p+rot.b*Lb       ); vert( p-rot.b*Lb       );

    int i,j;
    for(int ii=0; ii<3; ii++){
        int i1=i0+ii*2;
        edge(i1,i1+1);
        if(ii<2){
            edge(i1  ,i1+2);
            edge(i1+1,i1+3);
        }
    }
}

void Builder2::snapPrismFace( const Vec3d& p0, const Mat3d& rot, double La, double Lb, double h, double Lbh, bool bFace ){
    //prismFace( p0, rot, La, Lb, h, Lbh );
    printf("snapPrismFace: La: %f Lb: %f h: %f Lbh: %f  \n", La,Lb,h,Lbh);
    // ToDo: We must find the points on cube where to attach the edges, not to create a new vertexes.
    //int i0 = verts.size();
    int iv[6];
    Vec3d p;
    //p=p0+rot.a* La + rot.c*0.1; vert( p+rot.b*Lb       ); vert( p-rot.b*Lb       );
    p=p0+rot.a* La; iv[0]=findVert( p+rot.b*Lb,  R_snapVert ); iv[1]=findVert( p-rot.b*Lb , R_snapVert ); 
    p=p0+rot.c* h;  iv[2]=    vert( p+rot.b*(Lb-Lbh) );        iv[3]=    vert( p-rot.b*(Lb-Lbh) );
    //p=p0+rot.a*-La + rot.c*0.1; vert( p+rot.b*Lb       ); vert( p-rot.b*Lb       );
    p=p0+rot.a*-La; iv[4]=findVert( p+rot.b*Lb,  R_snapVert ); iv[5]=findVert( p-rot.b*Lb , R_snapVert ); 
    int ie_h1 =findEdgeByVerts({iv[0],iv[1]});
    int ie_v11=edge(iv[0],iv[2]); int ie_v12=edge(iv[1],iv[3]);  //   |   |
    int ie_h2 =edge(iv[2],iv[3]);                                //    ---
    int ie_v21=edge(iv[2],iv[4]); int ie_v22=edge(iv[3],iv[5]);  //   |   |
    int ie_h3 =findEdgeByVerts({iv[4],iv[5]});
    if(bFace){
        { int ies1[4]{ie_h1,ie_v11,ie_h2,ie_v12};  polygon( 4, ies1 ); }
        { int ies2[4]{ie_h2,ie_v21,ie_h3,ie_v22};  polygon( 4, ies2 ); }
    }
}



void Builder2::quad( Quat4i q, int face_type, int edge_type ){
    if( face_type>-2 ){ tri( q.x, q.y, q.z, face_type ); tri( q.x, q.w, q.y, face_type ); }
    if( edge_type>-3 ){ edge( q.y, q.z, edge_type );  }
    if( edge_type>-2 ){ edge( q.x, q.y, edge_type ); edge( q.y, q.w, edge_type ); edge( q.w, q.z, edge_type ); edge( q.z, q.x, edge_type ); }
}

int Builder2::rope( int ip0, int ip1, int typ, int nseg ){   
    //printf( "MeshBuilder2::rope(%i,%i,n=%i,t=%i)\n", ip0,ip1, n,typ );
    int i0 = verts.size();
    Vec3d p0,d;
    if( nseg!=1 ){
        p0 = verts[ip0].pos;
        d  = verts[ip1].pos-p0;
        if( nseg<0 ){ 
            double l = d.norm(); nseg=(int)(l/max_size)+1;   
            //printf( "rope() l=%g n=%i ]n", l, n ); 
        }         
        //printf( "rope() n=%i d(%g,%g,%g)  p0(%g,%g,%g) \n", n,  d.x,d.y,d.z,  p0.x,p0.y,p0.z );   
        d.mul(1./(double)nseg);  // we need this only it n>1
        //printf( "rope() n=%i d(%g,%g,%g)  p0(%g,%g,%g) \n", n,  d.x,d.y,d.z,  p0.x,p0.y,p0.z );   
    }
    int oi=ip0;
    for(int ii=1; ii<nseg; ii++){
        Vec3d p = p0+d*ii;
        int i=vert(p);
        edge(oi,i, typ);
        //Vec3d& v = verts.back().pos;
        //printf( "rope(%i,%i)[%i] p(%g,%g,%g)\n",  ip0,ip1,ii, v.x,v.y,v.z );
        oi=i;
    }
    edge( oi,ip1, typ );
    // ToDo: implementation by LINE_STRIP ?
    return i0;
};

void Builder2::ropes( int n, int nseg, const Vec2i* ends, int typ ){
    for(int ii=0; ii<n; ii++){ Vec2i ed = ends[ii]; rope( ed.x, ed.y, typ, nseg ); }
}

int Builder2::ring( Vec3d p, Vec3d a, Vec3d b, Vec2d cs, int n ){
    Vec2d rot = Vec2d{1.0,0.0};
    int i0  = verts.size();
    //int ip0 = ip0 + n;
    int ip0 = i0 + n;
    for(int i=0; i<n; i++){
        int ip = vert(p + a*rot.x + b*rot.y );
        edge( ip0, ip, edge_type.x );
        rot.mul_cmplx( cs );
        ip0 = ip;
    }
    return i0;
}

int Builder2::ring( Vec3d p, Vec3d ax, Vec3d up, double R, int n ){
    double dphi = 2*M_PI/n;
    if( n<0 ){ double dL   = R*dphi; int n2=(int)(dL/max_size)+1; }
    Vec2d cs; cs.fromAngle( dphi );
    Vec3d a,b; 
    b.set_cross(ax,up); b.normalize();
    a.set_cross(b ,ax); a.normalize();
    return ring(p,a*R,b*R,cs,n );
}

void Builder2::tube( Vec3d p0, Vec3d p1, Vec3d up, Vec2d R, Vec2i n ){
    Vec3d ax = p1-p0; double L = ax.normalize();
    if( n.x<0 ){ n.x=(int)(L/max_size)+1; }
    if( n.y<0 ){ double r = fmax(R.x,R.y); n.y=(int)(r/max_size)+1; }
    double dL = L/n.y;
    Vec2d cs; cs.fromAngle( 2*M_PI/n.x );
    Vec3d a,b; 
    b.set_cross(ax,up); b.normalize();
    a.set_cross(b ,ax); a.normalize();
    int oe;
    double dR = (R.y-R.x)/n.y;
    for(int iy=0; iy<n.y; iy++){
        int e = edges.size();
        double r = R.x + dR*iy;
        ring( p0+ax*(dL*iy), a*r, b*r, cs, n.x );
        if(iy>0){
            for(int ix=0; ix<n.x; ix++){
                Vec2i e1 = edges[e+ix].f.xy();
                Vec2i e2 = edges[oe+ix].f.xy(); 
                tri( e1.x, e1.y, e2.x, face_type ); 
                tri( e1.y, e2.y, e2.x, face_type );
                edge( e1.x, e2.x, edge_type.y );        
                //edge( e1.y, e2.y, edge_type.y );  // it would do it twice
                edge( e1.y, e2.x, edge_type.y );    // diagonal edge
            }
        }
        oe=e;
    }
}

//int strip( Vec3d p0, Vec3d p1, Vec3d d, int n,  ){
//    for(int i=0;i<n;i++){
//    }
//}


//int vstrip(Vec3d p0, Vec3d p1, int n, int et=-1 ){
//int fstrip( int ip1, int ip0, int n, int ft=-1, Vec2i et={-1,-1} ){ 


int Builder2::plate( Vec3d p00, Vec3d p01, Vec3d p10, Vec3d p11, Quat4i t, Vec2i n, int fillType ){
    Vec3d dx0=p01-p00; double lx0=dx0.norm();
    Vec3d dx1=p11-p10; double lx1=dx1.norm();
    Vec3d dy0=p10-p00; double ly0=dy0.norm();
    Vec3d dy1=p11-p01; double ly1=dy1.norm();
    if(n.x<0){ int n1=(int)(lx0/max_size)+1; int n2=(int)(lx1/max_size)+1;  n.x=_max(n1,n2); }
    if(n.y<0){ int n1=(int)(ly0/max_size)+1; int n2=(int)(ly1/max_size)+1;  n.y=_max(n1,n2); }
    //n.x=5;
    //n.y=5;
    //Quat4i t_;
    int i0=verts.size();
    int oiv;
    for(int iy=0;iy<n.y+1;iy++){
        double cy=iy/(double)n.y;  double my=1-cy; 
        Vec3d p0y = p00*my+p10*cy;
        Vec3d p1y = p01*my+p11*cy;
        //if(iy==0){ t_=Quat4i{t.x,-2,-2,-2}; }else{t=t;}
        //strip(p0,p1, n.x, t);
        int iv = vstrip( p0y, p1y, n.x, t.x );
        if(iy>0){
            fstrip( oiv, iv, n.x, t.w, {t.y,t.z} );
            //fstrip( oiv, iv, n.x, -2, {t.y,t.z} );

        }
        oiv=iv;
    }
    return i0;
}

// ToDo: This is too complicated, put we should remove it or move it elsewhere
int Builder2::plate_quad( int ip00, int ip01, int ip10, int ip11, Quat4i typs, Vec2i n, int fillType ){ 
    printf( "plate_quad(%i,%i;%i,%i) typs{%i,%i,%i,%i} n{%i,%i} fillType=%i \n", ip00,ip01,ip10,ip11, typs.x,typs.y,typs.z,typs.w, n.x, n.y, fillType ); 
    int i0 = verts.size();
    Vec3d p00,p01,p10,p11;
    //DEBUG
    p00 = verts[ip00].pos; p01 = verts[ip01].pos; p10 = verts[ip10].pos; p11 = verts[ip11].pos;
    if(n.x<0){ int n1=(int)((p00-p01).norm()/max_size)+1; int n2=(int)((p10-p11).norm()/max_size)+1;  n.x=_max(n1,n2); }
    if(n.y<0){ int n1=(int)((p10-p00).norm()/max_size)+1; int n2=(int)((p11-p01).norm()/max_size)+1;  n.y=_max(n1,n2); }
    printf( "p00(%g,%g,%g) p01(%g,%g,%g) p10(%g,%g,%g) p11(%g,%g,%g)\n", p00.x,p00.y,p00.z,  p01.x,p01.y,p01.z,  p10.x,p10.y,p10.z,  p11.x,p11.y,p11.z );
    n.x=3; n.y=3;
    printf( "plate_quad n(%i,%i)\n", n.x, n.y);
    int ex[n.x];
    int ey[n.y]; 
    //DEBUG
    // --- edges
    //typs.x = 234234;
    ex[0    ]=edges.size(); int i0x = rope( ip00, ip01, typs.x, n.y );
    ex[n.x-1]=edges.size(); int i1x = rope( ip10, ip11, typs.x, n.y );
    ey[0    ]=edges.size(); int i0y = rope( ip00, ip10, typs.x, n.x );
    ey[n.y-1]=edges.size(); int i1y = rope( ip01, ip11, typs.x, n.x );
    //DEBUG

    //rope( i0y+5, i1y+5, typs.y, n.x );
    // WARRNING: This will not work properly because the points generated here for x-strips and y-strips are generated twice, they are not colapsed => we should generate points first, only then add edges
    for(int iy=1; iy<n.y-1; iy++){ ex[iy]=rope( i0y+iy, i1y+iy, typs.y, n.x ); }
    for(int ix=1; ix<n.x-1; ix++){ ex[ix]=rope( i0x+ix, i1x+ix, typs.y, n.y ); }
    
    /*
    //DEBUG
    for(int iy=1; iy<n.y; iy++){
        for(int ix=1; ix<n.x; ix++){
            Vec2i e1= edges[ex[ix+1]].f.xy();
            Vec2i e2= edges[ex[ix  ]].f.xy();
            if(fillType>0){
                tri( e1.x, e1.y, e2.x, typs.z );
                tri( e1.y, e2.y, e2.x, typs.z );
                edge( e1.y, e2.x, typs.y );        // diagonal edge
            }else{
                tri( e2.x, e2.y, e1.x, typs.z );
                tri( e2.y, e1.y, e1.x, typs.z );
                edge( e2.y, e1.x, typs.y );
            }
        }
    }
    */
    //DEBUG
    // ToDo: implementation by TRIANGLE_STRIP ?
    return i0;
};

int Builder2::export_pos( Vec3d* ps, int i0, int i1 ){   if(i1<0){ i1=verts.size()-i1; };
    for(int i=i0; i<=i1; i++){ ps[i]=verts[i].pos; }
    return i1-i0+1;
}

int Builder2::export_pos( float4* ps, int i0, int i1 ){     if(i1<0){ i1=verts.size()-i1; };
    for(int i=i0; i<=i1; i++){ ps[i]=*(float4*)&verts[i].lo;}
    return i1-i0+1;
}

int Builder2::export_edges( Vec2i* eds, int i0, int i1 ){   if(i1<0){ i1=edges.size()-i1; };
    for(int i=i0; i<=i1; i++){ eds[i]=edges[i].lo; }
    return i1-i0+1;
}

int Builder2::export_tris( Quat4i* tri, int i0, int i1 ){  if(i1<0){ i1=edges.size()-i1; };
    for(int i=i0; i<=i1; i++){ tri[i]=tris[i]; }
    return i1-i0+1;
}

void Builder2::write_obj( const char* fname, uint8_t mask ){
    printf( "Mesh::Builder2::write_obj(%s)\n", fname );
    FILE * pFile;
    pFile = fopen(fname,"w");
    if (pFile == NULL) {
        printf("Error opening file %s\n", fname);
        return;
    }
    bool bVerts    = mask & ObjMask::Verts;
    bool bNormals  = mask & ObjMask::Normals;
    bool bUVs      = mask & ObjMask::UVs;
    bool bTris     = mask & ObjMask::Tris;
    bool bPolygons = mask & ObjMask::Polygons;
    bool bEdges    = mask & ObjMask::Edges;
    fprintf(pFile, "# SimpleSimulationEngine MeshBuilder2 OBJ export\n");
    if(bVerts){
        fprintf(pFile, "# Vertices: %zu\n", verts.size());
        for(const auto& v : verts){ fprintf(pFile, "v %f %f %f\n", v.pos.x, v.pos.y, v.pos.z); }
    }
    if(bNormals){
        fprintf(pFile, "# Vertex Normals: %zu\n", verts.size());
        for(const auto& v : verts){ fprintf(pFile, "vn %f %f %f\n", v.nor.x, v.nor.y, v.nor.z); }
    }
    if(bUVs){
        fprintf(pFile, "# Texture Coords: %zu\n", verts.size());
        for(const auto& v : verts){ fprintf(pFile, "vt %f %f\n", v.uv.x, v.uv.y); }
    }
    if(bTris){
        fprintf(pFile, "# Triangles: %zu\n", tris.size());
        for(const auto& t : tris){
            if(bUVs && bNormals){ fprintf(pFile, "f %i/%i/%i %i/%i/%i %i/%i/%i\n", t.x+1, t.x+1, t.x+1, t.y+1, t.y+1, t.y+1, t.z+1, t.z+1, t.z+1); } 
            else if (bNormals)  { fprintf(pFile, "f %i//%i %i//%i %i//%i\n", t.x+1, t.x+1, t.y+1, t.y+1, t.z+1, t.z+1);                            } 
            else if (bUVs)      { fprintf(pFile, "f %i/%i %i/%i %i/%i\n", t.x+1, t.x+1, t.y+1, t.y+1, t.z+1, t.z+1);                               } 
            else                { fprintf(pFile, "f %i %i %i\n", t.x+1, t.y+1, t.z+1);                                                             }
        }
    }
    if(bPolygons){
        fprintf(pFile, "# Polygons (from chunks)\n");
        for(const auto& ch : chunks){
            if(ch.w == (int)ChunkType::face){
                fprintf(pFile, "f");
                int* ivs = strips.data() + ch.x;
                for(int i=0; i<ch.z; i++){
                    int iv = ivs[i] + 1;
                    if(bUVs && bNormals){ fprintf(pFile, " %i/%i/%i", iv, iv, iv); } 
                    else if (bNormals  ){ fprintf(pFile, " %i//%i", iv, iv); } 
                    else if (bUVs      ){ fprintf(pFile, " %i/%i", iv, iv); } 
                    else                { fprintf(pFile, " %i", iv); }
                }
                fprintf(pFile, "\n");
            }
        }
    }
    if(bEdges){
        fprintf(pFile, "# Edges: %zu\n", edges.size());
        for(const auto& e : edges){ fprintf(pFile, "l %i %i\n", e.x+1, e.y+1); }
    }
    fclose(pFile);
    printf( "Mesh::Builder2::write_obj(%s) DONE\n", fname );
}

void Builder2::read_obj(const char* fname, uint8_t mask) {
    printf("Mesh::Builder2::read_obj(%s)\n", fname);
    FILE* pFile = fopen(fname, "r");
    if (pFile == NULL) { printf("ERROR in Builder2::read_obj(%s) Cannot open file\n", fname); exit(0); }

    bool bVerts    = mask & ObjMask::Verts;
    bool bNormals  = mask & ObjMask::Normals;
    bool bUVs      = mask & ObjMask::UVs;
    bool bTris     = mask & ObjMask::Tris;
    bool bPolygons = mask & ObjMask::Polygons;
    bool bEdges    = mask & ObjMask::Edges;

    int vert_offset = verts.size();
    std::vector<Vec3d> file_normals;
    std::vector<Vec2d> file_uvs;
    char line_buffer[2048];
    while (fgets(line_buffer, sizeof(line_buffer), pFile)) {
        char type[3] = {0};
        sscanf(line_buffer, "%2s", type);

        if ( type[0]=='v'){
            if ( type[1]=='n' ) {
                if(!bNormals) continue;
                Vec3d n;
                sscanf(line_buffer, "vn %lf %lf %lf", &n.x, &n.y, &n.z);
                file_normals.push_back(n);
            } else if ( type[1]=='t' ) {
                if(!bUVs) continue;
                Vec2d uv;
                sscanf(line_buffer, "vt %lf %lf", &uv.x, &uv.y);
                file_uvs.push_back(uv);
            }else{
                if(!bVerts) continue;
                Vec3d p;
                sscanf(line_buffer, "v %lf %lf %lf", &p.x, &p.y, &p.z);
                vert(p);
            }
        } else if ( type[0]=='f' ) {
            bool bFaces = bTris || bPolygons;
            if(!bFaces) continue;
            std::vector<int> face_v_indices;
            std::vector<int> face_vn_indices;
            std::vector<int> face_vt_indices;
            const char* p = line_buffer + 1;
            while (*p) {
                while (*p && isspace(*p)) p++;
                if (*p == '\0' || *p == '\r' || *p == '\n') break;
                int v_idx = 0, vt_idx = 0, vn_idx = 0;
                if (sscanf(p, "%d/%d/%d", &v_idx, &vt_idx, &vn_idx) == 3) {
                    face_v_indices.push_back(v_idx); face_vt_indices.push_back(vt_idx); face_vn_indices.push_back(vn_idx);
                } else if (sscanf(p, "%d//%d", &v_idx, &vn_idx) == 2) {
                    face_v_indices.push_back(v_idx); face_vn_indices.push_back(vn_idx);
                } else if (sscanf(p, "%d/%d", &v_idx, &vt_idx) == 2) {
                    face_v_indices.push_back(v_idx); face_vt_indices.push_back(vt_idx);
                } else if (sscanf(p, "%d", &v_idx) == 1) {
                    face_v_indices.push_back(v_idx);
                }
                while (*p && !isspace(*p)) p++;
            }
            if (face_v_indices.size() >= 3) {
                bool is_tri = (face_v_indices.size() == 3);
                if (is_tri && !bTris) continue;
                if (!is_tri && !bPolygons) continue;

                int n = face_v_indices.size();
                int iedges[n];
                int ivs[n];
                for (int i = 0; i < n; ++i) {
                    ivs[i] = face_v_indices[i] - 1 + vert_offset;
                }
                if (bNormals && face_vn_indices.size() == n && !file_normals.empty()) {
                    for (int i = 0; i < n; ++i) { if (ivs[i] < verts.size() && (face_vn_indices[i]-1) < file_normals.size()) verts[ivs[i]].nor = file_normals[face_vn_indices[i]-1]; }
                }
                if (bUVs && face_vt_indices.size() == n && !file_uvs.empty()) {
                    for (int i = 0; i < n; ++i) { if (ivs[i] < verts.size() && (face_vt_indices[i]-1) < file_uvs.size()) verts[ivs[i]].uv = file_uvs[face_vt_indices[i]-1]; }
                }
                for (int i = 0; i < n; ++i) {
                    iedges[i] = findOrAddEdges({ivs[i], ivs[(i + 1) % n]});
                }
                polygon(n, iedges);
            }
        } else if ( type[0]=='l' ) {
            if(!bEdges) continue;
            std::vector<int> line_v_indices;
            const char* p = line_buffer + 1;
            int v_idx;
            while(*p){
                while(*p && isspace(*p)) p++;
                if(*p == '\0' || *p == '\r' || *p == '\n') break;
                if(sscanf(p, "%d", &v_idx) == 1){
                    line_v_indices.push_back(v_idx - 1 + vert_offset);
                }
                while(*p && !isspace(*p)) p++;
            }
            if(line_v_indices.size() >= 2){
                for(size_t i=0; i < line_v_indices.size() - 1; ++i){
                    edge(line_v_indices[i], line_v_indices[i+1]);
                }
            }
        }
    }
    fclose(pFile);
    printf("Mesh::Builder2::read_obj(%s) DONE\n", fname);
}

void Builder2::printSelection(){
    printf( "Mesh::Builder2::printSelection() n=%i mode=%i :", selection.size(), selection_mode );
    for(int ii=0; ii<selection.size(); ii++){  printf( "%i ", selection[ii] ); }
    printf( "\n" );
}

void Builder2::printSelectedVerts(){
    printf( "Mesh::Builder2::printSelectedVerts() n=%i\n", selection.size() );
    for(int ii=0; ii<selection.size(); ii++){
        int i=selection[ii];
        printf( "%i -> %i pos: %16.10f %16.10f %16.10f\n", ii, i, verts[i].pos.x, verts[i].pos.y, verts[i].pos.z );
    }
}

void Builder2::printVerts(){
    printf( "Mesh::Builder2::printVerts() n=%i\n", verts.size() );
    for(int i=0; i<verts.size(); i++){
        printf( "%3i pos: %16.10f %16.10f %16.10f\n", i, verts[i].pos.x, verts[i].pos.y, verts[i].pos.z );
    }
}

void Builder2::printEdges(){
    printf( "Mesh::Builder2::printEdges() n=%i\n", edges.size() );
    for(int i=0; i<edges.size(); i++){
        printf( "%3i -> ivs: %3i %3i t: %i t2: %i \n", i, edges[i].x, edges[i].y, edges[i].z, edges[i].z  );
    }
}

void Builder2::printChunkRange( int ich, int ich2 ){
    printf( "Mesh::Builder2::printChunkRange() ich %i ich2 %i\n", ich, ich2 );
    if(ich2<0){ ich2=ich+1; }
    for(int i=ich; i<ich2; i++){
        Quat4i ch = chunks[i];
        printf( "chunk[%3i] -> is: %3i %3i t: %i t2: %i \n", i, ch.x, ch.y, ch.z, ch.w  );
    }
}

Vec3d Builder2::getCOG(int n, const int* ivs) const {
    //printf( "getCOG() n %i \n", n );
    Vec3d cog = Vec3dZero;
    for(int i=0; i<n; i++){
        //printf( "getCOG() %i %f %f %f\n", i, verts[ivs[i]].pos.x, verts[ivs[i]].pos.y, verts[ivs[i]].pos.z );
        cog.add( verts[ivs[i]].pos );
    }
    cog.mul(1./n);
    return cog;
}


void Builder2::move_verts( const std::vector<int>& indices, const Vec3d& shift ){
    for( int iv : indices ){
        if (iv >= 0 && iv < verts.size()) {
            verts[iv].pos.add(shift);
        }
    }
}

void Builder2::scale_verts( const std::vector<int>& indices, const Vec3d& p, const Vec3d& s ){
    Vec3d inv_s;
    if(s.x!=0) inv_s.x=1/s.x;
    if(s.y!=0) inv_s.y=1/s.y;
    if(s.z!=0) inv_s.z=1/s.z;
    for( int iv : indices ){
        if (iv >= 0 && iv < verts.size()) {
            Vert& v = verts[iv];
            v.pos.sub(p);
            v.pos.mul(s);
            v.pos.add(p);
            v.nor.mul(inv_s);
            v.nor.normalize();
        }
    }
}

void Builder2::rotate_verts( const std::vector<int>& indices, const Vec3d& p, const Mat3d& rot ){
    for( int iv : indices ){
        if (iv >= 0 && iv < verts.size()) {
            Vert& v = verts[iv];
            v.pos.sub(p);
            rot.dot_to(v.pos, v.pos);
            v.pos.add(p);
            rot.dot_to(v.nor, v.nor);
        }
    }
}

int Builder2::duplicateBlock( int iblock ){
    if( iblock < 0 || iblock >= blocks.size() ){ return -1; }
    Quat4i i0 = blocks[iblock];
    Quat4i i1;
    if( (iblock+1) < blocks.size() ){ i1 = blocks[iblock+1]; } else { i1 = latsBlock(); }
    int nVerts = i1.x - i0.x;
    int nEdges = i1.y - i0.y;
    int nTris  = i1.z - i0.z;
    int vert_offset = verts.size();
    for(int i=0; i<nVerts; i++){ verts.push_back( verts[i0.x + i] ); }
    for(int i=0; i<nEdges; i++){ Quat4i old_e = edges[i0.y + i]; edge( old_e.x + vert_offset, old_e.y + vert_offset, old_e.w, old_e.z ); }
    for(int i=0; i<nTris;  i++){ Quat4i old_t = tris [i0.z + i]; tri ( old_t.x + vert_offset, old_t.y + vert_offset, old_t.z + vert_offset, old_t.w ); }
    return block();
}








void Builder2::alling_polygons( int n, const int* ivs1, int* ivs2, int ipiv ){

    // compute centers of gravity
    Vec3d cog1 = getCOG( n, ivs1);
    Vec3d cog2 = getCOG( n, ivs2);
    // compute axis
    Vec3d ax = cog2 - cog1;
    ax.normalize();
    Vec3d u = verts[ ivs1[ipiv] ].pos - cog1;
    u.makeOrthoU(ax);
    u.normalize();
    Vec3d v; v.set_cross(ax, u);
    v.normalize();

    //printf( "alling_polygons() cog1 %f %f %f cog2 %f %f %f \n", cog1.x, cog1.y, cog1.z, cog2.x, cog2.y, cog2.z );
    //printf( "alling_polygons() ax %f %f %f u %f %f %f v %f %f %f \n", ax.x, ax.y, ax.z, u.x, u.y, u.z, v.x, v.y, v.z );

    //Mat3d fromDirUp( const VEC& dir, const VEC& up );
    
    // Project points to uv plane
    Vec2d uv1[n], uv2[n];
    for(int i=0; i<n; i++){
        Vec3d d1 = verts[ivs1[i]].pos - cog1;
        Vec3d d2 = verts[ivs2[i]].pos - cog2;
        uv1[i] = Vec2d{d1.dot(u), d1.dot(v)};
        uv2[i] = Vec2d{d2.dot(u), d2.dot(v)};
        //printf( "%i uv1 %f %f uv2 %f %f \n", i, uv1[i].x, uv1[i].y, uv2[i].x, uv2[i].y );
        // normalize to direction only
        uv1[i].normalize();
        uv2[i].normalize();
    }
    
    // 6) match points by nearest uv projection
    int map[n];
    for(int i=0; i<n; i++){
        double dmin = -1.0;
        int    imin = -1;
        Vec2d uvi = uv1[i];
        for(int j=0; j<n; j++){
            double d = uvi.dot(uv2[j]);
            //printf( "d=%f\n", d );
            if(d > dmin){
                dmin = d;
                imin = ivs2[j];
            }
        }
        map[i] = imin;
    }
    for(int i=0; i<n; i++){
        ivs2[i] = map[i];
    }
}



int Builder2::bridge_quads( Quat4i q1, Quat4i q2, int nseg, Quat4i stickTypes, Quat4i mask, bool bAlling ){
    printf( "====Mesh::Builder2::bridge_quads() nseg %i BEFORE q1: %i %i %i %i q2 = %i %i %i %i \n", nseg, q1.x, q1.y, q1.z, q1.w, q2.x, q2.y, q2.z, q2.w );
    if( bAlling ) alling_polygons(4, q1.array, q2.array, 0);
    //printf( "Mesh::Builder2::bridge_quads() AFTER q1: %i %i %i %i q2 = %i %i %i %i \n", q1.x, q1.y, q1.z, q1.w, q2.x, q2.y, q2.z, q2.w );
    Vec3d A1 = verts[q1.x].pos; Vec3d A2 = verts[q2.x].pos;
    Vec3d B1 = verts[q1.y].pos; Vec3d B2 = verts[q2.y].pos;
    Vec3d C1 = verts[q1.z].pos; Vec3d C2 = verts[q2.z].pos;
    Vec3d D1 = verts[q1.w].pos; Vec3d D2 = verts[q2.w].pos;
    int dnp = 4;
    int i00start = verts.size();
    int i00      = i00start;
    double dc = 1.0/((double)nseg);
    int oA=q1.x,oB=q1.y,oC=q1.z,oD=q1.w; // Start with q1 vertices
    for (int i=0; i<nseg; i++){ // Loop nseg times to create nseg segments
        printf( "Mesh::Builder2::bridge_quads() iseg %i\n", i );
        int iA,iB,iC,iD;
        if( i < (nseg-1) ){ // Create intermediate rings
            // vertices
            double c  = (i+1)*dc;
            double mc = 1-c;
            iA = vert( A1*mc + A2*c );
            iB = vert( B1*mc + B2*c );
            iC = vert( C1*mc + C2*c );
            iD = vert( D1*mc + D2*c );
            // ring edges
            printf( "Mesh::Builder2::bridge_quads()   ring-edges is(%i,%i,%i,%i)\n", iA,iB,iC,iD );
            edge( iA,iB, stickTypes.y );
            edge( iB,iC, stickTypes.y );
            edge( iC,iD, stickTypes.y );
            edge( iD,iA, stickTypes.y );
        }else{ // This is the last segment (i == nseg - 1), connect to q2
            iA=q2.x;
            iB=q2.y;
            iC=q2.z;
            iD=q2.w;
        }
        // longitudinal edges
        printf( "Mesh::Builder2::bridge_quads()   longitudinal-edges os(%i,%i,%i,%i) is(%i,%i,%i,%i)\n", oA,oB,oC,oD, iA,iB,iC,iD );
        edge( oA,iA,stickTypes.x );
        edge( oB,iB,stickTypes.x );
        edge( oC,iC,stickTypes.x );
        edge( oD,iD,stickTypes.x );
        // spiral edges
        if(mask.x){
            printf( "Mesh::Builder2::bridge_quads()   spiral-edges mask.x=%i os(%i,%i,%i,%i) is(%i,%i,%i,%i)\n", mask.x, oA,oB,oC,oD, iA,iB,iC,iD );
            edge( oA,iB, stickTypes.z );
            edge( oB,iC, stickTypes.z );
            edge( oC,iD, stickTypes.z );
            edge( oD,iA, stickTypes.z );
        }
        if(mask.y){
            printf( "Mesh::Builder2::bridge_quads()   spiral-edges mask.y=%i os(%i,%i,%i,%i) is(%i,%i,%i,%i)\n", mask.y, oA,oB,oC,oD, iA,iB,iC,iD );
            edge( oB,iA, stickTypes.z );
            edge( oC,iB, stickTypes.z );
            edge( oD,iC, stickTypes.z );
            edge( oA,iD, stickTypes.z );
        }
        // internal edges
        if(mask.z){
            printf( "Mesh::Builder2::bridge_quads()   internal-edges mask.z=%i os(%i,%i) is(%i,%i)\n", mask.z, oA,oB, iC,iD );
            edge( oA,iC, stickTypes.z );
            edge( oB,iD, stickTypes.z );
        }
        if(mask.w){
            printf( "Mesh::Builder2::bridge_quads()   internal-edges mask.w=%i os(%i,%i) is(%i,%i)\n", mask.w, oC,oD, iA,iB );
            edge( oC,iA, stickTypes.z );
            edge( oD,iB, stickTypes.z );
        }
        oA=iA;
        oB=iB;
        oC=iC;
        oD=iD;
        // i00+=dnp;
    }
    //return ibloc;
    return i00;
}

int Builder2::extrudeVertLoop( int n, int* iverts, Vec3d d, bool bEdges, bool bFace, bool bTris, bool bSort ){
    int iv0=verts.size();
    int ie0=edges.size();
    // if(bEdges || bFace){
    // }
    int   ivs[n];
    Vec2i es [n]; 
    int   ies[n];
    for( int i=0; i<n; i++ ){
        if( bSort ){ 
            es[i]=Vec2i{ iverts[i], iverts[(i+1)%n] }; // old verts ( already created )
        }else{
            ivs[i] = iverts[i];
            es[i]=Vec2i{iv0+i,iv0+((i+1)%n) };    // new verts ( not yet created )
        }
        //printf( "extrudeVertLoop() es[%i]=%i %i \n", i, es[i].x, es[i].y );
    }
    if(bSort){ sortPotentialEdgeLoop( n, es, ivs ); };
    for( int i=0; i<n; i++ ){
        //printf( "extrudeVertLoop() ivs[%i]=%i verts.size() %i \n", i, ivs[i], verts.size() );
        int iv = ivs[i];
        ivs[i] = vert( verts[iv].pos+d );
    }
    if(bEdges){
        for( int i=0; i<n; i++ ){
            ies[i] = edge( es[i].x, es[i].y );
        }
    }
    if(bFace){
        //polygonChunk( n, int* igs, ivs, bTris );
        //int ich = chunk( Quat4i{is0, is0+n, n, (int)ChunkType::face} );
        int i0 = strips.size();
        int ich = chunk( Quat4i{i0, i0+n, n, (int)ChunkType::face} );
        for(int i=0; i<n; i++){ strips.push_back( ivs[i] ); } // store verts
        for(int i=0; i<n; i++){ strips.push_back( ies[i] ); } // store edges
        if( bPolygonToTris){ polygonToTris( ich ); }
        return ich;
    }
    return iv0;
}

int Builder2::extrudeFace( int ich, double L, Quat4i stickTypes, Quat4i maks ){
    int ivs[4];
    int n    = loadChunk( ich, ivs );
    Vec3d nr = getChunkNormal( ich ); 
    int ich2 = extrudeVertLoop( n, ivs, nr*5.0, true, true, true, false );
    bridge_quads( *(Quat4i*)getChunkStrip( ich ), *(Quat4i*)getChunkStrip( ich2 ), n, stickTypes, maks );
    return ich2;
}

int Builder2::loadChunk( int ich, int* iedges, int* iverts ){
    Quat4i ch = chunks[ich];
    if( ch.w==(int)ChunkType::face ){
        for(int i=0; i<ch.z; i++){
            if(iedges){ iedges[i]=strips[ch.x+i     ]; }
            if(iverts){ iverts[i]=strips[ch.x+ch.z+i]; }
        }
        return ch.z;
    }
    return 0;
};

//Vec3d getChunkNormal( Quat4i ch ){
Vec3d Builder2::getChunkNormal( int ich ){
    Quat4i ch = chunks[ich];
    if( ch.w==(int)ChunkType::face ){
        const Vec3d& a = verts[strips[ch.x  ]].pos;
        const Vec3d& b = verts[strips[ch.x+1]].pos;
        const Vec3d& c = verts[strips[ch.x+2]].pos;
        Vec3d nr = cross(b-a, c-a); 
        nr.normalize(); 
        return nr;
    }
    return Vec3dZero;
}


// ============== From  SpaceCraft2Mesh2.h


/**
 * Creates a girder in the truss structure.
 * 
 * @param p0 starting point
 * @param p1 ending point
 * @param up up-vector 
 * @param n  number of segments along the girder 
 * @param width The width of the girder.
 * @return The index of the created block containing the starting indexes of the points and edges
 */
int Builder2::girder1( Vec3d p0, Vec3d p1, Vec3d up, int n, double width, Quat4i stickTypes, bool bCaps ){
    // ToDo: ad p0,p1 etc. - maybe we should rather specify indexes of existing verts rather than positions of new verts ?
    //                       No. this is done in girder1_caps
    //printf( "Mesh::girder1() n=%i p0(%g,%g,%g) p1(%g,%g,%g) up(%g,%g,%g) stickTypes(%i,%i,%i,%i)\n", n, p0.x,p0.y,p0.z, p1.x,p1.y,p1.z, up.x,up.y,up.z,   stickTypes.x,stickTypes.y,stickTypes.z,stickTypes.w );
    
    int ip0=-1,ip1=-1;
    if( bCaps ){
        ip0 = vert( p0 ); 
        ip1 = vert( p1 );
    }
    
    Vec3d dir = p1-p0;
    double length = dir.normalize();
    up.makeOrthoU(dir);
    Vec3d side; side.set_cross(dir,up); side.normalize();
    //print(dir); print(up); print(side); printf("dir up side \n");
    double dl = length/(2*n + 1);
    //int dnb = 2+4+4+4;
    int dnp = 4;
    int i00start = verts.size();
    int i00      = i00start;
    //int ibloc = block();   // this is better to call manually from outside
    //Quat4i istart;
    for (int i=0; i<n; i++){
        int i01=i00+1; int i10=i00+2; int i11=i00+3;
        vert( p0 + side*-width + dir*(dl*(1+2*i  )) );
        vert( p0 + side*+width + dir*(dl*(1+2*i  )) );
        vert( p0 + up  *-width + dir*(dl*(1+2*i+1)) );
        vert( p0 + up  *+width + dir*(dl*(1+2*i+1)) );
        edge( i00,i01,stickTypes.y );
        edge( i10,i11,stickTypes.y );
        edge( i00,i10,stickTypes.z );
        edge( i00,i11,stickTypes.z );
        edge( i01,i10,stickTypes.z );
        edge( i01,i11,stickTypes.z );
        if( i<(n-1) ){
            edge( i10,i00+dnp,stickTypes.w );
            edge( i10,i01+dnp,stickTypes.w );
            edge( i11,i00+dnp,stickTypes.w );
            edge( i11,i01+dnp,stickTypes.w );
            edge( i00,i00+dnp,stickTypes.x );
            edge( i01,i01+dnp,stickTypes.x );
            edge( i10,i10+dnp,stickTypes.x );
            edge( i11,i11+dnp,stickTypes.x );
        }
        i00+=dnp;
    }
    if( bCaps ){
        edge( i00start+0 ,ip0,stickTypes.x );
        edge( i00start+1 ,ip0,stickTypes.x );
        edge( i00start+2 ,ip0,stickTypes.x );
        edge( i00start+3 ,ip0,stickTypes.x );
        int i00end = i00-dnp;
        edge( i00end  +0 ,ip1,stickTypes.x );
        edge( i00end  +1 ,ip1,stickTypes.x );
        edge( i00end  +2 ,ip1,stickTypes.x );
        edge( i00end  +3 ,ip1,stickTypes.x );
    }
    //return ibloc;
    return i00;
}

/// this generate flat zig-zag triangle strip made of equal sized triangles. 
/// The upper and lower edges are parallel to each other 
/// The result vary if n is odd or even
///   1. for even n, the edges are tilted with respect to dir=p1-p0 because p1 is on the upper edge, and p0 is on the lower edge 
///   2. for odd n, both p1 and p0 are on the upper edge, and the edges are parallel to p1-p0
int Builder2::triangle_strip( Vec3d p0, Vec3d p1, Vec3d up, int n, double width, int stickType, bool bCaps ){
    printf( "Mesh::triangle_strip() n=%i p0(%g,%g,%g) p1(%g,%g,%g) up(%g,%g,%g)\n", n, p0.x,p0.y,p0.z, p1.x,p1.y,p1.z, up.x,up.y,up.z );
    Vec3d dir = p1-p0;
    if( (n&1)==0) dir.add_mul(up, -width);
    double length = dir.normalize();
    double dl     = length / (n+1);
    int i00 = verts.size();
    int ip0=-1,ip1=-1;
    if (bCaps) {
        ip0=vert(p0);
        ip1=vert(p1);
        i00 += 2;
    }
    int i00start = i00;
    int ip;
    for (int i = 0; i <n; i++) {
        Vec3d p = p0 + dir * dl*(i+1.0);
        if ( (i&1)==0 ){ p.add(up*width); }
        ip=vert(p);
        if (i > 0) { edge(ip, ip-1, stickType); } // zig-zag
        if (i > 1) { edge(ip, ip-2, stickType); } // beam mains
    }
    if (bCaps) {
        edge(i00start   , ip0, stickType);
        edge(i00start +1, ip0, stickType);
        edge(ip       -1, ip1, stickType);
        edge(ip         , ip1, stickType);
    }
    return i00;
}



/* @brief : plot planar fill between two strips of vertexes (e.g. between two girders or ropes)
*
*/
// int plateOnGriders( Builder2& mesh, Vec2i ns, Vec2i iblocks, Vec2i byN, Vec2i offs, Vec2d span1, Vec2d span2, Quat4i stickTypes ){
int Builder2::plateOnGriders( Vec2i ns, Vec2i prange1, Vec2i prange2, Vec2i byN, Vec2i offs, Vec2d span1, Vec2d span2, Quat4i stickTypes ){
    //printf( "plateOnGriders() ns(%i,%i) iblocks(%i,%i) byN(%i,%i) offs(%i,%i) span1(%6.3f,%6.3f) span2(%6.3f,%6.3f) \n", ns.x,ns.y, iblocks.x,iblocks.y, byN.x,byN.y,  offs.x,offs.y, span1.x,span1.y,  span2.x,span2.y );
    int n1   = (prange1.y - prange1.x)/byN.x;
    int n2   = (prange2.y - prange2.x)/byN.y;
    double step = 1.0/(ns.x-1);
    if(offs.x<0){
        int i0 = prange1.x + byN.x*(n1/2);
        offs.x=findClosestVert( verts[ prange2.x + byN.y*(n2/2) ].pos,i0,byN.x);
    }
    if(offs.y<0){
        int i0 = prange2.x + byN.y*(n2/2);
        offs.y=findClosestVert( verts[ prange1.x + byN.x*(n1/2) ].pos,i0,byN.y);
    }
    for(int i=0;i<ns.x; i++){
        double c = i*step;
        int i1 = prange1.x + offs.x + byN.x*(int)(  ( span1.x*(1-c)+span1.y*c)* n1 + 0.5 );
        int i2 = prange2.x + offs.y + byN.y*(int)(  ( span2.x*(1-c)+span2.y*c)* n2 + 0.5 );
        //printf( "plateOnGriders()[%i] (%i/%i) (%i/%i)\n", i, i1,n1, i2,n2 );
        edge( i1,i2,stickTypes.x );
    }
    //return ibloc;
    return 0;
}




int Builder2::girder1_caps( int ip0, int ip1, int kind ){
    //printf( "Truss::girder1_caps() ip0=%i ip1=%i ipbeg=%i ipend=%i \n", ip0, ip1, ipbeg, ipend );
    // it is expected that we call this immediately after girder1 - so we know the indexes of the first and last point of the girder block
    int ipbeg = blocks.back().x;  // index of the first point of the girder block
    int ipend = verts.size()-4;   // index of the last point  of the girder block
    int ie0 = edges.size();
    edge( ip0,ipbeg+0,kind );
    edge( ip0,ipbeg+1,kind );
    edge( ip0,ipbeg+2,kind );
    edge( ip0,ipbeg+3,kind );
    edge( ip1,ipend+0,kind );
    edge( ip1,ipend+1,kind );
    edge( ip1,ipend+2,kind );
    edge( ip1,ipend+3,kind );
    return ie0;
}

int Builder2::girder1( int ip0, int ip1, Vec3d up, int n, double width, Quat4i stickTypes ){
    int ib = girder1( verts[ip0].pos, verts[ip1].pos, up, n, width, stickTypes );
    girder1_caps( ip0, ip1, stickTypes.x );
    return ib;
}

/**
 * Creates a wheel-shaped truss structure.
 *
 * @param p0 starting point ( center of the wheel )
 * @param p1 ending point   ( on the perimeter of the wheel )
 * @param ax axis vector
 * @param n The number of segments along the perimeter of the wheel.
 * @param width The width of the wheel rim.
 * @return The index of the created block containing the starting indexes of the points and edges
 */
int Builder2::wheel( Vec3d p0, Vec3d p1, Vec3d ax, int n, Vec2d wh, Quat4i stickTypes ){
    printf( "Truss::wheel() n=%i p0(%g,%g,%g) p1(%g,%g,%g) ax(%g,%g,%g) stickTypes(%i,%i,%i,%i) \n", n, p0.x,p0.y,p0.z, p1.x,p1.y,p1.z, ax.x,ax.y,ax.z,   stickTypes.x,stickTypes.y,stickTypes.z,stickTypes.w );
    //int kind_long   = 0;
    //int kind_perp   = 1;
    //int kind_zigIn  = 2;
    //int kind_zigOut = 3;
    Vec3d dir = p1-p0;
    double r = dir.normalize();
    ax.makeOrthoU(dir);
    Vec3d side; side.set_cross(dir,ax);
    //double dl = length/(2*n + 1);
    //int dnb = 2+4+4+4;
    //print(dir); print(ax); print(side); printf("dir up side \n");
    int dnp = 4;
    int i00 = verts.size();
    int i000 = i00;

    Vec2d  rot = {1.0,0.0};
    Vec2d drot; drot.fromAngle( M_PI/n );
    //int ibloc = blocks.size();
    //blocks.push_back( {i00,edges.size()} );
    for (int i=0; i<n; i++){
        int i01=i00+1; int i10=i00+2; int i11=i00+3;

        Vec3d R = dir*rot.a + side*rot.b;
        //vert( p0 + R*r );
        vert( p0 + R*(r+wh.x) );
        vert( p0 + R*(r-wh.x) );
        // points.push_back( p0 +  R*(r-width) );
        // points.push_back( p0 +  R*(r+width) );
        rot.mul_cmplx(drot);
        R       = dir*rot.a + side*rot.b;
        vert( p0 + ax*+wh.y + R*r );
        vert( p0 + ax*-wh.y + R*r );
        // points.push_back( p0 + ax*-width + R*r );
        // points.push_back( p0 + ax*+width + R*r );
        rot.mul_cmplx(drot);

        edge( i00,i01,stickTypes.y );
        edge( i10,i11,stickTypes.y );
        edge( i00,i10,stickTypes.z );
        edge( i00,i11,stickTypes.z );
        edge( i01,i10,stickTypes.z );
        edge( i01,i11,stickTypes.z );
        // edges .push_back( (TrussEdge){i00,i01,kind_perp}  );
        // edges .push_back( (TrussEdge){i10,i11,kind_perp}  );
        // edges .push_back( (TrussEdge){i00,i10,kind_zigIn} );
        // edges .push_back( (TrussEdge){i00,i11,kind_zigIn} );
        // edges .push_back( (TrussEdge){i01,i10,kind_zigIn} );
        // edges .push_back( (TrussEdge){i01,i11,kind_zigIn} );
        if( i<(n-1) ){
            edge( i10,i00+dnp,stickTypes.w );
            edge( i10,i01+dnp,stickTypes.w );
            edge( i11,i00+dnp,stickTypes.w );
            edge( i11,i01+dnp,stickTypes.w );
            edge( i00,i00+dnp,stickTypes.x );
            edge( i01,i01+dnp,stickTypes.x );
            edge( i10,i10+dnp,stickTypes.x );
            edge( i11,i11+dnp,stickTypes.x );
            // edges.push_back( (TrussEdge){i10,i00+dnp,kind_zigOut} );
            // edges.push_back( (TrussEdge){i10,i01+dnp,kind_zigOut} );
            // edges.push_back( (TrussEdge){i11,i00+dnp,kind_zigOut} );
            // edges.push_back( (TrussEdge){i11,i01+dnp,kind_zigOut} );
            // edges.push_back( (TrussEdge){i00,i00+dnp,kind_long} );
            // edges.push_back( (TrussEdge){i01,i01+dnp,kind_long} );
            // edges.push_back( (TrussEdge){i10,i10+dnp,kind_long} );
            // edges.push_back( (TrussEdge){i11,i11+dnp,kind_long} );
        }else{
            edge( i10,i000+0,stickTypes.w );
            edge( i10,i000+1,stickTypes.w );
            edge( i11,i000+0,stickTypes.w );
            edge( i11,i000+1,stickTypes.w );
            edge( i00,i000+0,stickTypes.x );
            edge( i01,i000+1,stickTypes.x );
            edge( i10,i000+2,stickTypes.x );
            edge( i11,i000+3,stickTypes.x );
            // edges.push_back( (TrussEdge){i10,i000+0,kind_zigOut} );
            // edges.push_back( (TrussEdge){i10,i000+1,kind_zigOut} );
            // edges.push_back( (TrussEdge){i11,i000+0,kind_zigOut} );
            // edges.push_back( (TrussEdge){i11,i000+1,kind_zigOut} );
            // edges.push_back( (TrussEdge){i00,i000+0,kind_long} );
            // edges.push_back( (TrussEdge){i01,i000+1,kind_long} );
            // edges.push_back( (TrussEdge){i10,i000+2,kind_long} );
            // edges.push_back( (TrussEdge){i11,i000+3,kind_long} );
        }
        i00+=dnp;
    }
    return i000;
}

int Builder2::ngon( Vec3d p0, Vec3d p1, Vec3d ax, int n,  int stickType ){
    Vec3d dir = p1-p0;
    double r = dir.normalize();
    ax.makeOrthoU(dir);
    Vec3d side; side.set_cross(dir,ax);
    Vec2d  rot = {1.0,0.0};
    Vec2d drot; drot.fromAngle( 2*M_PI/n );
    int i0 = verts.size();
    for (int i=0; i<n; i++){
        Vec3d R = dir*rot.a + side*rot.b;
        vert( p0 + R*r );
        if (i > 0) {
            edge(i0+i-1, i0+i, stickType);
        }
        rot.mul_cmplx(drot);
    }
    edge(i0+n-1, i0, stickType);  // Close the loop
    return i0;
}


int Builder2::rope( Vec3d p0, Vec3d p1, int n, int ropeType, int anchorType, double Rcolapse, double r ){
    int i0 = make_anchor_point( p0, anchorType, Rcolapse, r );
    int i1 = make_anchor_point( p1, anchorType, Rcolapse, r );
    printf( "Builder2::rope() i0=%i i1=%i n=%i ropeType=%i anchorType=%i Rcolapse=%g r=%g \n", i0,i1,n,ropeType,anchorType,Rcolapse,r );
    rope( i0, i1, ropeType, n );
    return i0;
}

int Builder2::ropes( int nv, Vec3d* vs, int ne, int nseg, const Vec2i* ends, int ropeType, int anchorType, double Rcolapse, double r ){
    int ivrts[nv];
    make_anchor_points( nv, vs, ivrts, anchorType, Rcolapse, r );
    for( int i=0; i<ne; i++ ){
        Vec2i e = ends[i];
        rope( ivrts[e.x], ivrts[e.y], ropeType, nseg );
    }
    return nv;
};


/**
 * Creates a panel of truss elements between four corner points. Adds the points and edges to the Truss object.
 * // TODO: make also triangular panel
 * 
 * @param p00 The first corner point.
 * @param p01 The second corner point.
 * @param p10 The third corner point.
 * @param p11 The fourth corner point.
 * @param n The number of subdivisions along each side of the panel.
 * @param width The width of the truss elements.
 */
int Builder2::panel( Vec3d p00, Vec3d p01, Vec3d p10, Vec3d p11, Vec2i n, double width, Quat4i stickTypes ){
    // ToDo: ad p00,p01,p10,p11 etc. - maybe we should rather specify indexes of existing verts rather than positions of new verts ?
    //printf( "Mesh::panel() n(%i,%i) w=%g p00(%g,%g,%g) p01(%g,%g,%g) p10(%g,%g,%g) p11(%g,%g,%g) \n", n.x,n.y, p00.x,p00.y,p00.z, p01.x,p01.y,p01.z, p10.x,p10.y,p10.z, p11.x,p11.y,p11.z );
    //int kind_long   = 0;
    //int kind_perp   = 1;
    //int kind_zigIn  = 2;
    //int kind_zigOut = 3;
    Vec2d step = {1.0/n.a,1.0/n.b};
    int di = 2*n.a-1;
    //int ibloc = block();   // this is better to call manually from outside
    int i0  = verts.size();
    //int i00 = verts.size();
    for (int ib=0; ib<n.b; ib++){
        double db,mb;
        db = ib*step.b;     mb=1-db;
        Vec3d p0  = p00*mb + p10*db;
        Vec3d p1  = p01*mb + p11*db;
        db += 0.5*step.b; mb=1-db;
        Vec3d p0_ = p00*mb + p10*db;
        Vec3d p1_ = p01*mb + p11*db;
        for (int ia=0; ia<n.a; ia++){
            double da,ma;
            da = ia*step.a; ma = 1-da;
            Vec3d p   = p0 *ma + p1 *da;
            //points.push_back( p              );
            vert( p );
            int bi = i0+di; if( ib==n.b-2 )bi-=ia;
            int dia = 2;    if( ib==n.b-1 )dia=1;
            if (ia<(n.a-1)) edge( i0,i0+dia,stickTypes.y );
            if (ib<(n.b-1)) edge( i0,bi    ,stickTypes.y );
            if( (ia<(n.a-1))&&(ib<(n.b-1)) ){ // diagonal
                Vec3d p_  = p0_*ma + p1_*da;
                da += 0.5*step.a; ma=1-da;
                Vec3d p__ = p0_*ma + p1_*da;
                Vec3d up; up.set_cross( p_-p, p__-p ); up.normalize();
                //points.push_back( p__ + up*width );
                vert( p__ + up*width );
                if( ia<(n.a-2) ) edge( i0+1,i0+1+dia,stickTypes.z );
                if( ib<(n.b-2) ) edge( i0+1,bi+1    ,stickTypes.z );
                edge( i0+1,i0     ,stickTypes.w );
                edge( i0+1,i0+dia ,stickTypes.w );
                edge( i0+1,bi     ,stickTypes.w );
                if( ib==n.b-2 )dia=1;
                edge( i0+1,bi+dia ,stickTypes.w );
                i0++;
            }
            i0++;
        }
    }
    return i0;
}

void Builder2::facingNodes( const CMesh& cmesh, int nnod, const Vec3d* points, Vec2i* out_chs, int nplane, const int* planes, const int* planeVs ){
    //Vec2i chs[nnod];
    bool bSpecialNodes = (planes!=nullptr);
    //bSpecialNodes=false;
    for(int i=0; i<nnod; i++){
        Quat4i i0s = addCMesh( cmesh, !bSpecialNodes, points[i] ); // bFaces=false, we only want the wireframe initially
        //printf("  addCMesh() i0s(%i,%i,%i)\n", i0s.a, i0s.b, i0s.c);
        if(bSpecialNodes){ out_chs[i] = addFaces( nplane, planes, planeVs, true, i0s.x ); }
        else             { out_chs[i] = Vec2i{ i0s.w, (int)chunks.size()-1};                }
        printf( "facingNodes() %i chs(%i,%i)\n", i, out_chs[i].a, out_chs[i].b );
        //printf("  addFaces() chs(%i,%i) range\n", chs[i].a, chs[i].b);
    }
}


void Builder2::bridgeFacingPolygons( Vec3d p1, Vec3d p2, const Vec2i chr1, const Vec2i chr2, int nseg, Quat4i stickTypes, Quat4i maks ){
    printf(" --- findMostFacingNormal() chr1(%i,%i) chr2(%i,%i)   p1(%g,%g,%g) p2(%g,%g,%g)\n", chr1.x,chr1.y, chr2.x,chr2.y, p1.x,p1.y,p1.z, p2.x,p2.y,p2.z);
    Vec3d hray = p2-p1;
    hray.normalize();
    int ich1 = findMostFacingNormal(hray, chr1, 0.0, true );
    int ich2 = findMostFacingNormal(hray, chr2, 0.0, true );
    // if((ich1<0)||(ich2<0)){  
    //     printf("ERROR in Builder2::bridgeFacingPolygons(): ich1,2 %3i %3i \n", ich1, ich2 ); 
    //         for(int i=0; i<nch; i++){
    //             int current_chunk_id = chs[i];
    //             //printf("  [DEBUG] findMostFacingNormal: Processing chunk ID %d (of %d total)\n", current_chunk_id, nch);
    //         }
    //     Vec3d nr = polygonNormal(current_chunk_id); // This will now print its own debug info
    //     exit(0); 
    // }
    //printf(" ich1,2 %3i %3i  @iq1,2 %p %p \n", ich1, ich2, iq1, iq2 );
    int* iq1 = getChunkStrip( ich1 );
    int* iq2 = getChunkStrip( ich2 );
    //printf(" ich1,2 %3i %3i  @iq1,2 %p %p \n", ich1, ich2, iq1, iq2 );
    bridge_quads( *(Quat4i*)iq1, *(Quat4i*)iq2, nseg, stickTypes, maks );
}

void Builder2::bridgeFacingPolygons( int nrod, const Vec2i* edges, const Vec3d* points, int nseg, const Vec2i* chs,  Quat4i stickTypes, Quat4i maks ){
    //printf(" === Builder2::bridgeFacingPolygons() nrod=%i\n", nrod );
    for(int i=0; i<nrod; i++){
        Vec2i e = edges[i];
        //printf("Builder2::findMostFacingNormal() e(%i,%i)");
        bridgeFacingPolygons( points[e.x], points[e.y], chs[e.x], chs[e.y], nseg, stickTypes, maks );
        // Vec2i e = edges[i];
        // //printf(" ===== findMostFacingNormal() e(%i,%i) chs.a(%i,%i) chs.b(%i,%i)\n", e.x, e.y, chs[e.x].a, chs[e.x].b, chs[e.y].a, chs[e.y].b);
        // Vec3d hray = points[e.x] - points[e.y];
        // hray.normalize();
        // int ich1 = findMostFacingNormal(hray, chs[e.x], 0.0, true );
        // int ich2 = findMostFacingNormal(hray, chs[e.y], 0.0, true );
        // if((ich1<0)||(ich2<0)){ 
        //     printf("ERROR in oct_nodes: ich1,2 %3i %3i \n", ich1, ich2 ); exit(0); 
        // }
        // //printf(" ich1,2 %3i %3i  @iq1,2 %p %p \n", ich1, ich2, iq1, iq2 );
        // int* iq1 = getChunkStrip( ich1 );
        // int* iq2 = getChunkStrip( ich2 );
        // //printf(" ich1,2 %3i %3i  @iq1,2 %p %p \n", ich1, ich2, iq1, iq2 );
        // bridge_quads( *(Quat4i*)iq1, *(Quat4i*)iq2, nseg, stickTypes, maks );
    }
}

int Builder2::checkAllPointsConnected(bool bExit,bool bPrint) const{
    if (verts.empty()) {
        if(bPrint) printf("Mesh::Builder2::checkAllPointsConnected(): No vertices to check.\n");
        return 0;
    }
    std::vector<bool> connected(verts.size(), false);
    for (const auto& edge : edges) {
        if (edge.x < verts.size()) connected[edge.x] = true;
        if (edge.y < verts.size()) connected[edge.y] = true;
    }
    int unconnected_count = 0;
    for (size_t i = 0; i < verts.size(); ++i) {
        if (!connected[i]) {
            if (bPrint) {
                printf("WARNING: Mesh::Builder2::checkAllPointsConnected() - Vertex %zu at (%f, %f, %f) is not connected to any edge.\n",  i, verts[i].pos.x, verts[i].pos.y, verts[i].pos.z);
            }
            unconnected_count++;
        }
    }
    if (unconnected_count > 0) {  
        if(bPrint) printf("ERROR: Mesh::Builder2::checkAllPointsConnected() found %d unconnected vertices.\n", unconnected_count);
        if(bExit) exit(0);
    }
    return unconnected_count;
}

} // namespace Mesh