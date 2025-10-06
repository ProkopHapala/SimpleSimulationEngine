
#ifndef MechMesh2D_h
#define MechMesh2D_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "geom2D.h"
#include "globals.h"




/*

Idea
-----

visco-elestic simulation with Multi-Materials on using arbitrary triangular mesh
 * mass is stored in triangle COG
 * each triangle has different material


ToDo:
-----
    1) Swap edges
    2) Equation of states
    3) dcog/dvert   -   https://www.mathopenref.com/coordcentroid.html





  Proper way how to calculate forces on vertexes :
---------------------------------------------------
F_vert = dE/d(r_vert) = Sum_i (dE/dV_i)*( dV_i/d(r_vert) )
    (dE/dV_i)           ...    obtain from Equations of State (  EOS )
    ( dV_i/d(r_vert) )  ...    obtain form geometry


add 3) dcog/dvert  Derivation:

*/



#include <vector>
#include <algorithm>
#include "arrayAlgs.h"


inline bool linkEdges( const Vec2i& a, const Vec2i& b, Vec3i t ){
    if( a.a == b.a ){ t.c=a.a; t.a=a.b; t.b=b.b; return true; }
    if( a.a == b.b ){ t.c=a.a; t.a=a.b; t.b=b.a; return true; }
    if( a.b == b.a ){ t.c=a.b; t.a=a.a; t.b=b.b; return true; }
    if( a.b == b.b ){ t.c=a.b; t.a=a.a; t.b=b.a; return true; }
    return false;
}



namespace Mesh2D{

struct Edge{
    Vec2i verts;
    Vec2i tris;
    inline bool addTri(int it){
        if(tris.a<0){ tris.a=it; return true; }
        if(tris.b<0){ tris.b=it; return true; }
        return false;
    }

};

struct Tri{
    Vec3i verts;
    Vec3i edges;
};

class Builder{ public:
    std::vector<Vec2d> verts;
    std::vector<Edge>  edges;
    std::vector<Tri>   tris;

    void fromArrayVerts( int n, Vec2d* vs ){
        for(int i=0; i<n; i++){ verts.push_back( vs[i] ); }
    }
    void fromArrayEdges( int n, Vec2i* es ){
        for(int i=0; i<n; i++){ edges.push_back( {es[i], (Vec2i){-1,-1}} ); }
    }
    //void fromArrayTris( int n, Vec3i* vs ){
    //    for(int i=0; i<n; i++){ tris.push_back( ts[i] ); }
    //}

    void sortEdges(){
        for(Edge& e:edges){ e.verts.order(); } // make sure vert with lower index is first
        //struct { bool operator()(const Edge& a, const Edge& b)const{  return ((a.verts.a<<16) + a.verts.b) < ((b.verts.a<<16) + b.verts.b); } } edgeComparator;
        struct { bool operator()(const Edge& a, const Edge& b)const{  return a.verts.a < b.verts.a; } } edgeComparator;
        std::sort( edges.begin(), edges.end(), edgeComparator);
    }

    bool insertTri2edverts( const Vec2i& a, const Vec2i& b ){
        Vec3i t;
        if ( linkEdges( a, b, t ) ){ tris.push_back( (Tri){t,(Vec3i){-1,-1,-1}} ); return true; }
        return false;
    }

    bool insertTri2edge( int ia, int ib ){
        Edge& ea = edges[ia];
        Edge& eb = edges[ib];
        Vec3i t;
        if ( linkEdges( ea.verts, ea.verts, t ) ){
            int itri = tris.size();
            ea.addTri( itri );
            eb.addTri( itri );
            tris.push_back( (Tri){t,(Vec3i){ia,ib,-1}} );
            return true;
        }
        return false;
    }

    inline bool unfoldEdgeVertIndex(int ie, Vec2i& vs){
        //int ie = v2eds[i0+permut[k]];
        bool swaped = ie&1; //  first or second endge
        ie>>=1;
        vs = edges[ie].verts;
        if(swaped) _swap(vs.a,vs.b);
        return swaped;
    }

    void makeTris( ){

        int ne = edges.size();
        int nv = verts.size();

        Vec2i* ed2vs = new Vec2i[ne];
        int*   v2eds = new int  [ne*2];
        int*   vNs   = new int  [nv];
        int*   v0s   = new int  [nv];

        DEBUG
        for(int i=0; i<ne; i++){ ed2vs[i] = edges[i].verts; };
        DEBUG
        objects2cells( ne*2, nv, (int*)ed2vs, vNs, v0s, v2eds );
        DEBUG
        const int perVertMax = 16;
        float ranks [perVertMax];
        int   permut[perVertMax];
        Vec2d dirs  [perVertMax];
        DEBUG
        for(int iv=0; iv<nv; iv++){
            int i0  = v0s[iv];
            int nve = vNs[iv];
            if(perVertMax<=nve){ printf("ERROR in Mesh2D::Builder::makeTris(): number of edges per vertex(%i) > perVertMax(%i) \n", nve,perVertMax); exit(0); };
            printf("-----------\n");
            for(int k=0; k<nve; k++){
                /*
                int ie = v2eds[i0+k];
                int swaped = ie&1; //  first or second endge
                ie>>=1;
                Vec2i vs = ed2vs[ie];
                if(swaped) _swap(vs.a,vs.b);
                */
                Vec2i vs;
                int ie = v2eds[i0+k];
                unfoldEdgeVertIndex(ie, vs);
                //printf( "[%i|%i] %i (%i,%i)\n", k, i0,   ie, vs.a, vs.b );
                Vec2d va   = verts[vs.b] - verts[vs.a];
                va.normalize();
                dirs[k] = va;
                double phi = va.y;
                if(va.x>0) phi = -2-phi;
                ranks[k] = phi;
                printf( "[%i] %g |  %i(%i,%i) (%g,%g)\n", k, phi,   ie,  vs.a, vs.b, va.x,va.y );
            }
            for(int k=0; k<nve; k++){ permut[k]=k; }
            quickSort( ranks, permut, 0, nve);
            //insertTri2edge( permut[m-1], permut[0] );
            int ok = permut[nve-1];
            for(int kk=0; kk<nve; kk++){
                /*
                printf( "--[%i] %g \n", k, ranks[permut[k]] );
                int ie = v2eds[i0+permut[k]];
                int swaped = ie&1; //  first or second endge
                ie>>=1;
                Vec2i vs = ed2vs[ie];
                if(swaped) _swap(vs.a,vs.b);
                Vec2d va   = verts[vs.b] - verts[vs.a];
                //insertTri2edge( permut[k-1], permut[k] );
                */
                int k = permut[kk];
                Vec2d v1 ,v2;
                Vec2i vs1,vs2;
                int ie1,ie2;
                v1=dirs[k ];
                v2=dirs[ok];
                double s = v1.cross(v2);
                printf("s %g \n", s);
                if( s > 0 ){
                    ie1=v2eds[i0+k ]; unfoldEdgeVertIndex(ie1, vs1);
                    ie2=v2eds[i0+ok]; unfoldEdgeVertIndex(ie2, vs2);
                }
            }
        }
        DEBUG
        delete [] ed2vs;
        delete [] v2eds;
        delete [] vNs;
        delete [] v0s;
    }


    /*
    void makeTris( ){
        const int perVertMax = 16;
        float ranks [perVertMax];
        int   permut[perVertMax];
        // sortEdges();
        // we assume edges are sorted
        int i  = 0;
        int iv = edges[0].verts.a;
        int ne = edges.size();
        for(int j=1; j<ne+1; j++){
            int jv;
            if(j==ne){ jv=-1; }else{ jv= edges[j].verts.a; };
             // handle last edge case
            if(jv!=iv){
                // -- rank edges connected to this vertex by angle round circle
                int m = j-i;
                printf( " === %i %i  %i \n", i,j, m );
                //for(int k=0; k<m; k++){  printf(" e[%i](%i,%i) \n", k, edges[i+k].verts.a, edges[i+k].verts.b ); };
                if(perVertMax<m){ printf("ERROR in Mesh2D::Builder::makeTris(): number of edges per vertex(%i) > perVertMax(%i) \n",m,perVertMax); exit(0); };
                for(int k=0; k<m; k++){
                    Vec2i ivs  = edges[i+k].verts;
                    Vec2d va   = verts[ivs.b] - verts[ivs.a];
                    va.normalize();
                    double phi = va.y;
                    if(va.x>0) phi = -2-phi;
                    ranks[k] = phi;
                    printf( "[%i] %g \n", k, phi );
                }
                quickSort( ranks, permut, 0, m);
                insertTri2edge( permut[m-1], permut[0] );
                for(int k=1; k<m; k++){
                    insertTri2edge( permut[k-1], permut[k] );
                }
                i=j;
            }
        }
        delete [] ed2vs;
        delete [] v2eds;
        delete [] vNs;
        delete [] v0s;
    }
    */
};


};
















class MechMesh2D{ public:

    int nverts     =0;  // number of vertexes
    Vec2d* pos     =0;  // position of vertexes
    Vec2d* vels    =0;  // velocity of vertexes
    Vec2d* fvs     =0;  // acceleration
    double* invMass=0;  // invMass divided by 3 (because of centroid coordinates)

    int nedges   =0;
    Vec2i*  e2v  =0;   // edge -> vertex
    Vec2i*  e2t  =0;   // edge -> triangle
    Vec2d*  efs  =0;   // force ond edges ?
    double* els  =0;   // lengths

    int ntris    =0;   //
    Vec3i* t2v   =0;   // triangle -> vertex
    Vec3i* t2e   =0;   // triangle -> edge
    Vec2d*  ps   =0;

    //Vec2d*  dcog =0;  // dcog/dvert ... chnage of centre of mass (COG) by movement of vertex
    int*    material =0;
    double* pressure =0;
    double* volume   =0;
    double* mass     =0;
    double* Us       =0; // internal energy

void realloc(int nverts_, int nedges_, int ntris_ ){
    nedges=nverts_;
    nedges=nedges_;
    ntris =ntris_;

    _realloc( pos  ,nverts); // position of vertexes
    _realloc( vels ,nverts); // velocity of vertexes
    _realloc( fvs  ,nverts); // velocity of vertexes

    _realloc( e2v  ,nedges);   // edge -> vertex
    _realloc( e2t  ,nedges);   // edge -> triangle
    _realloc( efs  ,nedges);   // force ond edges ?
    _realloc( els  ,nedges);   // lengths

    _realloc( t2v      ,ntris);   // triangle -> vertex
    _realloc( t2e      ,ntris);   // triangle -> edge
    _realloc( ps       ,ntris);
    //_realloc( dcog     ,ntris);      // dcog/dvert ... chnage of centre of mass (COG) by movement of vertex
    _realloc( material ,ntris);
    _realloc( mass     ,ntris);
    _realloc( pressure ,ntris);
    _realloc( volume   ,ntris);
    _realloc( Us       ,ntris); // internal energy
}

void updateVerts(double dt){
    for(int i=0; i<nverts; i++){
        Vec2d& vi = vels[i];
        vi.add_mul(fvs[i],dt*invMass[i]);
        pos[i].add_mul(vi,dt);
    }
}


/*
double momentum2verts(){
    for(int i=0; i<ntris; i++){
        const  Vec3i& vi = t2v[i];
        const  Vec2d  pi = ps [i];
        //double invm      = invMass3[i];
        //pi.mul( invMass3[i] );
        fvs[vi.a].add( pi );
        fvs[vi.b].add( pi );
        fvs[vi.c].add( pi );

        //const Vec2d& dcog = dcogs[i];
        //vs[vi.a] += pi*dcog.a; // maybe we need inverse?
        //vs[vi.b] += pi*dcog.b;
        //vs[vi.c] += pi*dcog.c;

        //const Vec2d& invmi = invm[i]; // inverse mass of vertex with respect to volume ( ToDo - may be 2x2 tensor ?)
        //vs[vi.a] += pi*invmi.a;
        //vs[vi.b] += pi*invmi.b;
        //vs[vi.c] += pi*invmi.c;
    }
}

double updateMomentum(double dt){
    for(int i=0; i<ntris; i++){
        const Vec3i& e = t2e[i];
        ps[i].add_mul( efs[e.a] + efs[e.b] + efs[e.c], dt ); // accelerate volume by forces
    }
}

double pressure2edges(){
    for(int i=0; i<nedges; i++){
        const Vec3i& e  = e2v[i];
        const Vec3i& t  = e2t[i];
        double dp = pressure[t.a] - pressure[t.b]; // ToDo: check orientation ???
        Vec2d f;
        f.set_perp(pos[e.b]-pos[e.a]);
        f.mul(dp); // dp*l/|dr| = dp
        efs[i] = f;
    }
}

double updateVolume(){
    // https://en.wikipedia.org/wiki/Adiabatic_process
    for(int i=0; i<ntris; i++){
        const Vec3i& t = t2v[i];
        double V   = triangleArea(pos[t.a],pos[t.b],pos[t.c]);
        volume  [i] = V;
        pressure[i] = EOS_u2p( material[i], Us[i], V );
    }
};

*/

/*
double getPressure(int i){

    const Vec3i& t = t2v[i];
    double V   = triangleArea(pos[t.a],pos[t.b],pos[t.c]);
    volume  [i] = V;

    return Us[i]/Vs[i];
}
*/

void getForces(){
    // F_vert = dE/d(r_vert) = Sum_i (dE/dV_i)*( dV_i/d(r_vert) )
    for(int i=0; i<ntris; i++){
        const  Vec3i& vi = t2v[i];
        //double dEdV = get_pressure( i ); // pressure
        double dEdV;
        { // get pressure
            const  Vec3i& t = t2v[i];
            double V        = triangleArea(pos[t.a],pos[t.b],pos[t.c]);
            // U = TS - pV     (better use Helmholtz energy?)
            dEdV           = Us[i]/V; // pressure from EOS
        }

        Vec2d  dVdr;
        const Vec2d& pa = pos[vi.a];
        const Vec2d& pb = pos[vi.b];
        const Vec2d& pc = pos[vi.c];

        dtriangleArea(  pb, pc, dVdr );   fvs[vi.a].add_mul( dVdr, dEdV );
        dtriangleArea(  pc, pa, dVdr );   fvs[vi.a].add_mul( dVdr, dEdV );
        dtriangleArea(  pa, pb, dVdr );   fvs[vi.a].add_mul( dVdr, dEdV );
        // NOTE : WARRNING - be ware of chirality of triangles !!!!
    }
}


void step(double dt){
/*
    updateVolume  ();
    pressure2edges();
    updateMomentum(dt);
    momentum2verts();
    updateVerts   (dt);
*/

    //updateVolume  ();
    getForces   ();
    updateVerts (dt);

}

//double EOS(int imat, double mass, double V, double U, double& p, double& T){
//double EOS(int imat, double mass, double V, double& p, double& T){
//    // https://en.wikipedia.org/wiki/Internal_energy
//    // U = pV = nRT
//
//    T = U/(C*mass);
//   p = U/V;
//};

/*
//double EOS(int imat, double mass, double V, double U, double& p, double& T){
double EOS_U2p(int imat, double U, double V){
    // https://en.wikipedia.org/wiki/Internal_energy
    // https://en.wikipedia.org/wiki/Adiabatic_process
    // U = pV = nRT
    p = U/V;
};

bool mergeCond(int i, int j){ // check if merging does not cause too much error
    double dE =
    pressure[i]*volume[i] + pressure[j]*volume[j]  -  (pressure[i]+pressure[j])/ + *volume[j]  // pressure energy loss
    ps[i].norm2(i)/mass   + ps[j].norm2()/mass[j]  -  ;                 // kinetic energy loss
}

*/

bool trySwapEdge(int i, int j, int ie){
    if ( material[i] != material[j] ) return false;

    const Vec2i& vi = e2v[ie];
    const Vec2d dp  = pos[vi.b]- pos[vi.a];
    double l2_new    = dp.norm2();

    double minLfactor = 1.5;
    if( sq(els[ie]*minLfactor) > l2_new ){ // check if new edge sufficiently shorter

        // physical properties of new volume
        Vec2d  p     = ps[i] + ps[j];

        double press,temp;
        //double press = EOS( material[i], mass[i], volume[i], Us[i], press, temp );
        //double press = EOS( material[i], mass[i], volume[i], press, temp );
        //double press = EOS_u2p( material[i], volume[i], press );
        press = Us[i]/volume[i];
        double m     = mass[i] + mass[j];

        double maxEerr2 = 0.01;
        double dEp=0,dEk=0;
        //double dEp = press*V     - ( pressure[i]*volume[i] + pressure[j]*volume[j] );
        //double dEk = p.norm2()*m - ( ps[i].norm2()*mass[i] + ps[j].norm2()*mass[j] );

        if( dEp*dEp + dEk*dEk < maxEerr2 ){  // check energy error due to merging
            // swap edge

        }

    }
}


void remesh(){
/*
reorganize mesh such that for two triangles sharing an edge, the longer edges are switched to shorter edge
     --------       -------
    /     . /      /\     /
   /    .  /      /  \    /
  /   .   /  --> /    \   /
 /  .    /      /      \  /
/ .     /      /        \ /
 --------      --------

 * this is done only for same materials
 * momentum and pressure must be updated

*/
}


};

#endif

