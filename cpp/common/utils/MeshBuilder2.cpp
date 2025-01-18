#include "MeshBuilder2.h"

namespace Mesh {

void Builder2::clear(){ blocks.clear(); verts.clear(); edges.clear(); tris.clear(); chunks.clear(); strips.clear(); }
void Builder2::printSizes(){printf( "MeshBuilder::printSizes() blocks=%i verts=%i edges=%i tris=%i chunks=%i strips=%i \n", blocks.size(), verts.size(), edges.size(), tris.size(), chunks.size(), strips.size() );}

// int Builder2::getOtherEdgeVert(int ie, int iv){
//     Quat4i& e = edges[ie];
//     return e.x==iv ? e.y : e.x;
// }

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
    int ie0    = iedges[0];
    Quat4i& e0 = edges[ie0];
    int iv    = e0.y;
    int oe    = ie0; 
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

void Builder2::buildVerts2Edge(){
    vert2edge.clear();
    for( int i=0; i<edges.size(); i++ ){
        Vec2i e = edges[i].lo;
        uint64_t key = symetric_id( e );
        vert2edge[key] = i;
    }
}

int Builder2::polygon( int n, int* iedges ){
    int ivs[n];
    if( !sortEdgeLoop( n, iedges, ivs ) ) return false;
    int i0 = strips.size();
    int ich = chunk( Quat4i{i0, i0+n, n, (int)ChunkType::face} );
    for(int i=0; i<n; i++){ strips.push_back( ivs[i]       ); } // store verts
    for(int i=0; i<n; i++){ strips.push_back( selection[i] ); } // store edges
    if( bPolygonToTris){ polygonToTris( ich ); }
    return ich;
}

int Builder2::polygonToTris( int i ){
    Quat4i ch = chunks[i];
    int n = ch.z;
    int i0 = ch.x;
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

int Builder2::pickVertex( const Vec3d& ray0, const Vec3d& hRay, double R ){
    double tmin =  1e+300;
    int imin    = -1;
    for(int i=0; i<verts.size(); i++){
        //if(ignore)if(ignore[i])continue;
        double ti = raySphere( ray0, hRay, R, verts[i].pos );
        if(ti<tmin){ imin=i; tmin=ti; }
    }
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
    switch( (SelectionMode)selection_mode ){
        case SelectionMode::vert: return pickVertex    ( ro, rh, Rmax ); break;
        case SelectionMode::edge: return pickEdgeSelect( ro, rh, Rmax ); break;
        case SelectionMode::face: return pickTriangle  ( ro, rh, true ); break;
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

void Builder2::selectRect( const Vec3d& p0, const Vec3d& p1, const Mat3d& rot  ){
    selset.clear();
    //printf( "pickSelect() selection_mode %i  \n", selection_mode);
    switch( (SelectionMode)selection_mode ){
        case SelectionMode::edge: selectRectEdge( p0, p1, rot ); break;
    }
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
    if(n==-1){
        n  =selection.size();
        sel=selection.data();
    }
    double r2min=Rmax;
    int imin=-1;
    for(int ii=0;ii<n;ii++){ 
        int i = sel[ii];
        double r2=(verts[i].pos-p0).norm2(); if(r2<r2min){r2min=r2; imin=i;} 
    }
    if( bExitError && (imin<0)){ 
        printf( "ERROR in findVert() imin=%i cannot find vert near p(%g,%g,%g) within Rmax: %g => Exit() \n", imin, p0.x,p0.y,p0.z, Rmax );
        printSelectedVerts(); 
        exit(0); 
    };
    return imin;
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

// void snapFrustrumFace( const Vec3d& p0, const Mat3d& rot, double La, double Lb, double h, double Lbh, double Lah, bool bFace=true ){
//     //printf("frustrumFace: La: %f Lb: %f h: %f Lbh: %f Lah: %f \n", La,Lb,h,Lbh);
//     int iv[8];
//     Vec3d p;
//     p=p0+rot.a* La                ; iv[0]=findVert( p+rot.b*Lb, R_snapVert       ); iv[1]= findVert( p-rot.b*Lb, R_snapVert       );
//     p=p0+rot.a* (La-Lah) + rot.c*h; iv[2]=   vert( p+rot.b*(Lb-Lbh) );              iv[3]=     vert( p-rot.b*(Lb-Lbh) );
//     p=p0+rot.a*-(La-Lah) + rot.c*h; iv[4]=   vert( p+rot.b*(Lb-Lbh) );              iv[5]=     vert( p-rot.b*(Lb-Lbh) );
//     p=p0+rot.a*-La                ; iv[6]=findVert( p+rot.b*Lb, R_snapVert       ); iv[7]= findVert( p-rot.b*Lb, R_snapVert       );
//     edge(iv[0],iv[2]); edge(iv[1],iv[3]); //  |   |
//     edge(iv[2],iv[3]);                    //   ---
//     edge(iv[2],iv[4]); edge(iv[3],iv[5]); //  |   |
//     edge(iv[4],iv[5]);                    //   ---
//     edge(iv[4],iv[6]); edge(iv[5],iv[7]); //  |   |
//     if(bFace){
//         // AI: implement this 
//     }
// }

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
        { int ies[4]{ findEdgeByVerts({iv[6], iv[7]}),  ie_v21, ie_h2,  ie_v22  };  polygon(4, ies); }  // botton
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

int Builder2::rope( int ip0, int ip1, int typ, int n ){   
    //printf( "MeshBuilder2::rope(%i,%i,n=%i,t=%i)\n", ip0,ip1, n,typ );
    int i0 = verts.size();
    Vec3d p0,d;
    if( n!=1 ){
        p0 = verts[ip0].pos;
        d  = verts[ip1].pos-p0;
        if( n<0 ){ 
            double l = d.norm(); n=(int)(l/max_size)+1;   
            //printf( "rope() l=%g n=%i ]n", l, n ); 
        }         
        //printf( "rope() n=%i d(%g,%g,%g)  p0(%g,%g,%g) \n", n,  d.x,d.y,d.z,  p0.x,p0.y,p0.z );   
        d.mul(1./(double)n);  // we need this only it n>1
        //printf( "rope() n=%i d(%g,%g,%g)  p0(%g,%g,%g) \n", n,  d.x,d.y,d.z,  p0.x,p0.y,p0.z );   
    }
    int oi=ip0;
    for(int ii=1; ii<n; ii++){
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

void Builder2::printSelectedVerts(){
    printf( "Mesh::Builder2::printSelectedVerts() n=%i\n", selection.size() );
    for(int ii=0; ii<selection.size(); ii++){
        int i=selection[ii];
        printf( "%i -> %i pos: %16.10f %16.10f %16.10f\n", ii, i, verts[i].pos.x, verts[i].pos.y, verts[i].pos.z );
    }
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


int Builder2::rope( Vec3d p0, Vec3d p1, int n,  int stickType ){
    int i0 = vert( p0 );
    int i1 = vert( p1 );
    rope( i0, i1, stickType, n );
    return i0;
}


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


} // namespace Mesh