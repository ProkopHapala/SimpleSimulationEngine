
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "testUtils.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "SDL_utils.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"

#include "raytrace.h"
#include "Forces.h"

#include "Molecule.h"
#include "MMFFmini.h"
#include "NBFF.h"
#include "MMFFparams.h"
#include "MMFFBuilder.h"
#include "DynamicOpt.h"
#include "QEq.h"

#include "Draw3D_Molecular.h"  // it needs to know MMFFparams

//#include "NBSRFF.h"
#include "IO_utils.h"

#include "AppSDL2OGL_3D.h"









// ============= Graph analysis



//int splitGraphs( int nb, Vec2i* bonds, int a0, int b0 ){
int splitGraphs( int nb, Vec2i* bonds, int b0, std::unordered_set<int>& innodes ){
    //printf( "splitGraphs \n" );
    std::unordered_set<int> exbonds; // excluded bonds
    //std::unordered_set<int> innodes;
    exbonds.insert(b0);
    //innodes.insert(a0);
    int n0;
    do{
        n0=innodes.size();
        for( int ib=0; ib<nb; ib++ ){
            //printf( "ib %i n0 %i \n", ib, n0 );
            if( exbonds.find(ib) != innodes.end() ) continue;
            const Vec2i& b = bonds[ib];
            int ia=-1;
            if     ( innodes.find(b.a) != innodes.end() ){ ia=b.b; }
            else if( innodes.find(b.b) != innodes.end() ){ ia=b.a; }
            if(ia>=0){
                innodes.insert(ia);
                exbonds.insert(ib);
            }
        }
    }while( innodes.size()>n0 );
    /*
    int i=0;
    for( int ia : innodes){
        printf( "splitGraphs[%i] %i \n", i, ia );
        i++;
    }
    exit(0);
    */
    return innodes.size();
}
















// ==========================
// TestAppMMFFmini
// ==========================

void plotSurfPlane( Vec3d normal, double c0, Vec2d d, Vec2i n ){
    Vec3d da,db;
    normal.getSomeOrtho( da,db );
    da.mul( d.a/da.norm() );
    db.mul( d.b/db.norm() );
    //glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos(normal, {0.0,0.0,0.0} );
    //glColor3f(0.0f,1.0f,0.0f); Draw3D::drawVecInPos(da*10, {0.0,0.0,0.0} );
    //glColor3f(0.0f,0.0f,1.0f); Draw3D::drawVecInPos(db*10, {0.0,0.0,0.0} );
    Draw3D::drawRectGridLines( n*2, (da*-n.a)+(db*-n.b) + normal*c0, da, db );
}

void drawAngle(int i, MMFFmini& ff){
    //Vec2i b=ff.ang2bond[i];
    float sz=0.1;
    Vec3i a=ff.ang2atom[i];
    Draw3D::drawArrow(ff.apos[a.x],ff.apos[a.y],sz);
    Draw3D::drawArrow(ff.apos[a.y],ff.apos[a.z],sz);
    //if(b.i&SIGN_MASK){ Draw3D::drawArrow(apos[a.x],apos[a.y]); }else{ Draw3D::drawArrow(apos[a.x],apos[a.y]); };
    //if(b.j&SIGN_MASK){   };
}

void drawTors(int i, MMFFmini& ff){
    float sz=0.1;
    //Vec2i b=ff.ang2bond[i];
    Quat4i t=ff.tors2atom[i];
    Draw3D::drawArrow(ff.apos[t.x],ff.apos[t.y],sz);
    Draw3D::drawArrow(ff.apos[t.y],ff.apos[t.z],sz);
    Draw3D::drawArrow(ff.apos[t.z],ff.apos[t.w],sz);
    //if(b.i&SIGN_MASK){ Draw3D::drawArrow(apos[a.x],apos[a.y]); }else{ Draw3D::drawArrow(apos[a.x],apos[a.y]); };
    //if(b.j&SIGN_MASK){   };
}



/*
void rotateAtoms( Vec3d p0, Vec3d ax, int n, int* selection, Vec3d* atoms, double angle ){
    Vec2d rot; rot.fromAngle(angle);
    for(int i=0; i<n; i++){
        int ia  = selection[i];
        Vec3d d = selection - p0;
        d.rotate()

    }
}
*/






bool bPrint = true;


inline double cos_func(double ca, double cb, double cab ){  return ( cab - ca*cb )/sqrt( (1-ca*ca)*(1-cb*cb) ); }

inline void ddot_num(Vec3d h, Vec3d h2, Vec3d& f, double d){
    h.normalize(); h2.normalize();
    h.x+=d; h.normalize(); f.x=h.dot(h2); h.x-=2*d; h.normalize(); f.x-=h.dot(h2); h.x+=d;
    h.y+=d; h.normalize(); f.y=h.dot(h2); h.y-=2*d; h.normalize(); f.y-=h.dot(h2); h.y+=d;
    h.z+=d; h.normalize(); f.z=h.dot(h2); h.z-=2*d; h.normalize(); f.z-=h.dot(h2); h.z+=d;
    f.mul(1/(d*2));
}

//inline double evalTorq(Vec3d& ha,Vec3d& hb,Vec3d& hab,   Vec3d& f1, Vec3d& f2, Vec3d& f2, Vec3d& f4 ){
inline double evalTorq(Vec3d& ha,Vec3d& hb,Vec3d& hab,   Vec3d& fa, Vec3d& fb, Vec3d& fab ){

    double invra  = 1.0d/ha .normalize();
    double invrb  = 1.0d/hb .normalize();
    double invrab = 1.0d/hab.normalize();

    //double invra3 = invra*invra*invra;
    //double invrb3 = invrb*invrb*invrb;
    //double invrc3 = invrc*invrc*invrc;

    // check
    Vec3d haT = ha - hab*hab.dot(ha); haT.normalize();
    Vec3d hbT = hb - hab*hab.dot(hb); hbT.normalize();
    double cT = haT.dot(hbT);


    double ca   = hab.dot(ha);
    double cb   = hab.dot(hb);
    double cab  = ha .dot(hb);
    double sa2  = (1-ca*ca);
    double sb2  = (1-cb*cb);
    double invs = 1/sqrt( sa2*sb2 );
    double c    = ( cab - ca*cb )*invs;  //  c = <  ha - <ha|hab>hab   | hb - <hb|hab>hab    >

    //double c   = cos_func(ca,cb,cab);

    //printf( "<fa|fb> %g cT %g cS %g \n", cs.x, cT, cS );

    // derivatives to get forces

    double invs3 = invs*invs*invs;
    double dcab  = invs;                           // dc/dcab = dc/d<ha|hb>
    double dca   = (1-cb*cb)*(ca*cab - cb)*invs3;  // dc/dca  = dc/d<ha|hab>
    double dcb   = (1-ca*ca)*(cb*cab - ca)*invs3;  // dc/dca  = dc/d<hb|hab>

    fa =Vec3dZero;
    fb =Vec3dZero;
    fab=Vec3dZero;

    Mat3Sd J;

    J.from_dhat(ha);    // -- by ha
    J.mad_ddot(hab,fa, dca ); // dca /dha = d<ha|hab>/dha
    J.mad_ddot(hb ,fa, dcab); // dcab/dha = d<ha|hb> /dha

    J.from_dhat(hb);    // -- by hb
    J.mad_ddot(hab,fb, dcb ); // dcb /dhb = d<hb|hab>/dha
    J.mad_ddot(ha ,fb, dcab); // dcab/dhb = d<hb|ha> /dha

    J.from_dhat(hab);         // -- by hab
    J.mad_ddot(ha,fab, dca);  // dca/dhab = d<ha|hab>/dhab
    J.mad_ddot(hb,fab, dcb);  // dcb/dhab = d<hb|hab>/dhab
    // derivative cab = <ha|hb>

    fa .mul(invra *invra *invra );
    fb .mul(invrb *invrb *invrb );
    fab.mul(invrab*invrab*invrab);

    if(bPrint){

        double d=0.001d;
        double dcaE  =( cos_func(ca+d,cb  ,cab  ) - cos_func(ca-d,cb  ,cab  ) )/(2*d);
        double dcbE  =( cos_func(ca  ,cb+d,cab  ) - cos_func(ca  ,cb-d,cab  ) )/(2*d);
        double dcabE =( cos_func(ca  ,cb  ,cab+d) - cos_func(ca  ,cb  ,cab-d) )/(2*d);


        Vec3d fE=Vec3dZero,f=Vec3dZero,ferr;
        Vec3d h;
        h=ha;

        ddot_num(ha,hab,fE,d); f=Vec3dZero; J.from_dhat(ha   );   J.dhat_dot(hab,f);
        printf( "d<a|ab>/da err %g : anal  (%g,%g,%g)  | num  (%g,%g,%g) \n", (f-fE).norm(), f.x,f.y,f.z,   fE.x,fE.y,fE.z );
        ddot_num(ha,hb,fE,d);  f=Vec3dZero; J.from_dhat(ha   );   J.dhat_dot(hb,f);
        printf( "d<a|b>/dha  err %g : anal  (%g,%g,%g)  | num  (%g,%g,%g) \n",(f-fE).norm(), f.x,f.y,f.z,   fE.x,fE.y,fE.z );

        ddot_num(hb,hab,fE,d); f=Vec3dZero; J.from_dhat(hb   );   J.dhat_dot(hab,f);
        printf( "d<b|ab>/da err %g : anal  (%g,%g,%g)  | num  (%g,%g,%g) \n",(f-fE).norm(), f.x,f.y,f.z,   fE.x,fE.y,fE.z );
        ddot_num(hb,ha,fE,d);  f=Vec3dZero; J.from_dhat(hb   );   J.dhat_dot(ha,f);
        printf( "d<a|b>/db  err %g : anal  (%g,%g,%g)  | num  (%g,%g,%g) \n",(f-fE).norm(), f.x,f.y,f.z,   fE.x,fE.y,fE.z );

        ddot_num(hab,ha,fE,d); f=Vec3dZero; J.from_dhat(hab   );   J.dhat_dot(ha,f);
        printf( "d<a|ab>/dab err %g : anal  (%g,%g,%g)  | num  (%g,%g,%g) \n",(f-fE).norm(), f.x,f.y,f.z,   fE.x,fE.y,fE.z );
        ddot_num(hab,hb,fE,d);  f=Vec3dZero; J.from_dhat(hab   );   J.dhat_dot(hb,f);
        printf( "d<b|ab>/dab err %g : anal  (%g,%g,%g)  | num  (%g,%g,%g) \n",(f-fE).norm(), f.x,f.y,f.z,   fE.x,fE.y,fE.z );

        printf( "dca,dcb,dcab : anal  (%g,%g,%g)  | num  (%g,%g,%g) \n", dca,dcb,dcab,   dcaE,dcbE,dcabE );
        printf( "c,cT   %g %g       ca,cb,cab (%g, %g, %g) \n",  c,cT,  ca,cb, cab );
        printf( " <ha|fa> %g <hb|fb> %g <hab|fa> %g <hab|fb> %g \n", ha.dot(fa), hb.dot(fb), hab.dot(fa), hab.dot(fb) );

    }
    //if(bPrint) printf( "n", fa.dot(fb), fa.dot(fb), fa.dot(fb), fa.dot(fb), );

    /*
    aforce[ia.x].add(fa*-1);
    aforce[ia.y].add(fa-fab);
    aforce[ia.z].add(fab-fb);
    aforce[ia.w].add(fb);
    */
    return c;
}


void doChecks(){
    Vec3d ha      =(Vec3d){1.0, 1.0,-0.2};
    Vec3d hb      =(Vec3d){1.0,-1.0,+0.3};
    Vec3d hab     =(Vec3d){0.1, 0.2, 1.0};
    //ha.normalize();
    //hb.normalize();
    //hab.normalize();
    Vec3d fa,fb,fab;
    evalTorq( ha,hb,hab,  fa,fb,fab );
    bPrint=false;

    auto getEF_a  = [&](Vec3d p,Vec3d& f)->double{ return evalTorq( p,hb,hab,  f,fb,fab ); };
    auto getEF_b  = [&](Vec3d p,Vec3d& f)->double{ return evalTorq( ha,p,hab,  fa,f,fab ); };
    auto getEF_ab = [&](Vec3d p,Vec3d& f)->double{ return evalTorq( ha,hb,p,   fa,fb,f  ); };

    Vec3d f,fE;
    printf("ha  "); checkDeriv3d(getEF_a , ha, 0.001d,  fE, f );
    printf("hb  "); checkDeriv3d(getEF_b , hb, 0.001d,  fE, f );
    printf("hab "); checkDeriv3d(getEF_ab, hab,0.001d,  fE, f );

    //exit(0);
    /*
    double ds = 0.02;
    Vec3d& h  = ha0;
    Vec3d& f  =
    h   = ha;
    for(int i=0;i<3;i++){
        double E0,E1;
        h.array[i]-=ds*0.5;  E0= torqForce( ha,hb,hab,   fa,fb,fab );
        h.array[i]+=ds;      E1= torqForce( ha,hb,hab,   fa,fb,fab );
        double fE.array[i]=(E1-E0)/ds;
    }
    torqForce( ha,hb,hab, fa,fb,fab );

    printf( "fE(%g,%g,%g) fE(%g,%g,%g)\n", fE.x, fE.y, fE.z, f.x,f.y,f.z );
    */
    /*
    for(int i=0; i<n; i++){
        ha.x+=dx;
        double E  = torqForce( ha,hb,hab,   f1,f2,f2,f4 );
        double fx = (E-oE)/dx;
        printf( f1 );
    }
    */
}

// ===========================================
// ================= MAIN CLASS ==============
// ===========================================

class TestAppMMFFmini : public AppSDL2OGL_3D { public:
	Molecule    mol;
	MMFFparams  params;
    MMFFmini    ff;
    NBFF       nff;
    MM::Builder builder;
    DynamicOpt  opt;

    int* atypes = 0;

    bool bNonBonded = true;


    std::vector<int> selection;


    //std::unordered_map<std::string,int> atomTypeDict;

    //Mat3d lvec;

    bool bConverged = false;
    bool bRunRelax  = false;

    int     fontTex;
    int     ogl_sph=0;
    int     ogl_mol=0;

    char str[256];

    Vec3d ray0;
    int ipicked  = -1, ibpicked = -1, iangPicked = -1;
    int perFrame =  50;


    double drndv =  10.0;
    double drndp =  0.5;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppMMFFmini( int& id, int WIDTH_, int HEIGHT_ );

	int makeMoleculeInline();
	int makeMoleculeInlineBuilder( bool bPBC );
	int loadMoleculeMol( const char* fname, bool bAutoH, bool loadTypes );
	int loadMoleculeXYZ( const char* fname, bool bAutoH );

	void drawSystem( );

};


int TestAppMMFFmini::loadMoleculeXYZ( const char* fname, bool bAutoH ){
    params.loadAtomTypes( "common_resources/AtomTypes.dat" );
    params.loadBondTypes( "common_resources/BondTypes.dat" );
    builder.bindParams(&params);
    int nheavy = builder.load_xyz( fname, bAutoH, true, true );
    readMatrix( "common_resources/polymer.lvs", 3, 3, (double*)&builder.lvec );



    //builder.printAtoms ();
    //builder.printConfs ();
    builder.printAtomConfs();
    //exit(0);
    builder.export_atypes(atypes);

    builder.bDEBUG = true;
    //builder.autoBonds ();             builder.printBonds ();
    builder.autoBondsPBC();             builder.printBonds ();  // exit(0);
    //builder.autoBondsPBC(-0.5, 0, -1, {0,0,0});             builder.printBonds ();  // exit(0);
    builder.autoAngles( 0.5, 0.5 );     builder.printAngles();


    //bNonBonded = false;
    //exit(0);
    builder.toMMFFmini( ff, &params );

    builder.saveMol( "data/polymer.mol" );

    return nheavy;
}


int TestAppMMFFmini::loadMoleculeMol( const char* fname, bool bAutoH, bool bLoadTypes ){

    /// should look in   test_SoftMolecularDynamics.cpp

    if(bLoadTypes){
        printf( "bLoadTypes==True : load atom and bond types from file \n" );
        params.loadAtomTypes( "common_resources/AtomTypes.dat" );
        params.loadBondTypes( "common_resources/BondTypes.dat");
        //builder.params = &params;
    }else{
        printf( "bLoadTypes==False : make default Atom names dictionary \n" );
        params.initDefaultAtomTypeDict( );
    }
    //mol.atomTypeDict  = &params.atomTypeDict;
    //mol.atomTypeNames = &params.atomTypeNames;
    mol.bindParams(&params);
    mol.loadMol( fname );
    //mol.loadMol_const( "common_resources/propylacid.mol");
    //mol.loadMol_const( "/home/prokop/Dropbox/TEMP/ERC2021/Molecules/chain--frag-4---N-----.mol" );
    //exit(0);
    //int iH = 1;
    int iH = params.atomTypeDict["H"];
    int nh     = mol.countAtomType(iH); printf( "nh %i\n", nh );
    int nheavy = mol.natoms - nh;
    if(bAutoH){
        printf( "bAutoH==True : MMFFBuilder will re-calculate hydrogens, pi-orbitals and free electron pairs \n" );
        builder.insertFlexibleMolecule_ignorH( &mol, Vec3dZero, Mat3dIdentity, iH );
        //builder.setConfs( 0, 0, 0, mol.natoms-nh );
        //for(int i=0;i<(mol.natoms-nh);i++){ builder.makeSPConf(i,0,0); }
        //for(int i=0;i<mol.natoms;i++)     { builder.makeSPConf(i,0,0); }
    }else{
        printf( "bAutoH==False : Angles assigned by simple algorithm Molecule::autoAngles \n" );
        //mol.bondsOfAtoms();
        params.assignREs( mol.natoms, mol.atomType, mol.REQs );
        mol.autoAngles(true);
        Vec3d cog = mol.getCOG_av();
        mol.addToPos( cog*-1.0d );
        builder.insertMolecule(&mol, Vec3dZero, Mat3dIdentity, false );
        builder.toMMFFmini( ff, &params );
    }

    //builder.sortAtomsOfBonds();
    builder.tryAddConfsToAtoms(0, nh);
    builder.tryAddBondsToConfs();
    //for(int i=0; i<nh; i++){ builder.addConfToAtom(i); }
    //builder.tryAddBondsToConfs();

    //mol.printAtomInfo();
    //mol.printAtom2Bond();
    //mol.printAngleInfo();
    builder.printAtoms();
    //builder.printBonds();
    //builder.printAngles();
    //builder.printConfs();

    //bNonBonded = false;      // ToDo : WARRNING : this is just hack, because builder.sortBonds() does not seem to work, we have to switch off Non-Bonding interactions
    builder.trySortBonds();
    //builder.sortBonds();
    builder.printBonds();
    builder.printAngles();
    builder.printConfs();
    builder.toMMFFmini( ff, &params );

    //Draw3D::shapeInPoss( ogl_sph, ff.natoms, ff.apos, 0 );

    /*
    ogl_mol = glGenLists(1);
    glNewList( ogl_mol, GL_COMPILE );
        Draw3D::drawLines( mol.nbonds, (int*)mol.bond2atom, mol.pos );
    glEndList();
    */

    return nheavy;
}

int TestAppMMFFmini::makeMoleculeInlineBuilder( bool bPBC ){
    //const int natom=4,nbond=3,nang=2,ntors=1;
    //const int natom=4,nbond=3,nang=0,ntors=0;
    const int natom=4,nbond=4;
    Vec3d apos0[] = {
        {-2.0,0.0,0.0},  // 0
        {-1.0,2.0,0.0},  // 1
        {+1.0,2.0,0.0},  // 2
        {+2.0,0.0,1.0},  // 3
    };
    Vec2i bond2atom[] = {
        {0,1},  // 0
        {1,2},  // 1
        {2,3},  // 2
        {3,0},  // 3  // PBC
    };
    // ============ Build molecule
    MM::Atom brushAtom        { -1, -1, -1 , Vec3dZero, MM::Atom::defaultREQ };
    MM::Bond brushBond        { -1, {-1,-1}, 1.5,  25.0 };
    builder.capBond = MM::Bond{ -1, {-1,-1}, 1.07, 15.0 };

    builder.insertAtoms( natom, brushAtom,  apos0    );
    builder.insertBonds( nbond, brushBond, bond2atom );
    builder.setConfs   ( 0, 0 );
    builder.autoAngles ( 2.5, 1.25 );

    // instert aditional dihedral
    MM::Dihedral brushDihedral{ -1,   {-1,-1,-1},    3, 0.5 };  println(brushDihedral);
    builder.insertDihedralByAtom( {0,1,2,3}, brushDihedral );
    builder.trySortBonds();

    builder.toMMFFmini( ff, &params );

    builder.lvec.a = (Vec3d){  5.0,0.0,0.0 };
    builder.lvec.b = (Vec3d){  0.0,5.0,0.0 };
    builder.lvec.c = (Vec3d){  0.0,0.0,5.0 };
    if(bPBC){    // --- Periodic Boundary Conditions
        ff.initPBC();                 // as far as good, pbc-shifts are curenlty zero, so no change
        ff.pbcShifts[1] = builder.lvec.a*-1.; // make bond 3 from nighboring cell
        ff.printBondParams();
    }

    return natom;
}

int TestAppMMFFmini::makeMoleculeInline(){
    printf( "----- makeMoleculeInline \n" );
    const int natom=5+2,nbond=4+3,nang=6, ntors=0;
    Vec3d apos0[] = {
        { 0.5, 0.5, 0.5},  // 0
        {-1.0,+1.0,+1.0},  // 1
        {+1.0,-1.0,+1.0},  // 2
        {-1.0,-1.0,-1.0},  // 3
        {+1.0,+1.0,-1.0},  // 4
        {-1.0,-1.0,-2.0},  // 5
        {+1.0,+1.0,-2.0}   // 6
    };
    Vec2i bong2atom[] = {
        {0,1},  // 0
        {0,2},  // 1
        {0,3},  // 2
        {0,4},  // 3
        {5,6},  // 4
        {3,5},  // 5
        {4,6}   // 6
    };
    Vec2i ang2bond[] = {
        {0,1},
        {0,2},
        {0,3},
        {1,2},
        {1,3},
        {2,3}
    };
    double a0s[] ={
        2.0,
        2.0,
        2.0,
        2.0,
        2.0,
        2.0
    };
    Vec3i tors2bond[] = { };
    double l0    = 1.5;
    double Kbond = 10.0;
    double Kang  = 3.0;
    double Ktors = 1.0;
    int tors_n = 3;
    ff.realloc(natom,nbond,nang,ntors);
    printf( "----- Atoms \n" );
    for(int i=0; i<ff.natoms; i++){
        ff.apos[i] = apos0[i];
    }
    printf( "----- Bonds \n" );
    for(int i=0; i<ff.nbonds; i++){
        ff.bond2atom[i] = bong2atom[i];
        ff.bond_k [i] = Kbond;
        ff.bond_l0[i] = l0;
    }
    printf( "----- Angles \n" );
    for(int i=0; i<ff.nang; i++){
        ff.ang2bond[i] = ang2bond[i];
        double a0 = -a0s[i]/2.0; // NOTE: we use half-angle
        ff.ang_cs0[i] = { cos(a0), sin(a0) };
        ff.ang_k  [i] = Kang;
    }
    printf( "----- Dihedrals \n" );
    for(int i=0; i<ff.ntors; i++){
        ff.tors2bond[i] = tors2bond[i];
        ff.tors_k   [i] = Ktors;
        ff.tors_n   [i] = tors_n;
    }
    ff.angles_bond2atom  ();
    ff.torsions_bond2atom();
    return natom;
};

//void TestAppMMFFmini::makeAtoms(){}
//template<typename T> std::function<T(const T&,const T&         )> F2;

TestAppMMFFmini::TestAppMMFFmini( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    doChecks();

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    // ------ using molecule from mol-file does not seem to work now - There is some problem with AtomConfs
    // >> This algorithm assumes all atoms with conf precede atoms without confs in the array
    // >>   ERROR in builder.sortBonds() => exit

    int nheavy = 0;
    nheavy = loadMoleculeXYZ( "common_resources/polymer.xyz", false );
    //nheavy = loadMoleculeXYZ( "common_resources/polymer-noH.xyz", false );
    //nheavy = loadMoleculeMol(  "/home/prokop/Dropbox/TEMP/ERC2021/Molecules/chain--frag-4---N-----.mol", false, true);
    //nheavy = loadMoleculeMol( "common_resources/propylacid-q.mol", false, true);   // use old method loading whole .mol file with hydrogens // currently distorted molecule :-(
    //nheavy = loadMoleculeMol( "common_resources/propylacid.mol", true, false);   // use new method  // currently makes NaNs :-(
    //nheavy = makeMoleculeInlineBuilder();     // piece of polyethylene
    //nheavy = makeMoleculeInline();

    nff.bindOrRealloc( ff.natoms, ff.nbonds, ff.apos, ff.aforce, 0, ff.bond2atom );
    builder.export_REQs( nff.REQs );

    if(bNonBonded){
        if( !checkPairsSorted( nff.nmask, nff.pairMask ) ){
            printf( "ERROR: nff.pairMask is not sorted => exit \n" );
            exit(0);
        };
    }else{
        printf( "WARRNING : we ignore non-bonded interactions !!!! \n" );
    }

    opt.bindOrAlloc( 3*ff.natoms, (double*)ff.apos, 0, (double*)ff.aforce, 0 );
    //opt.setInvMass( 1.0 );
    opt.cleanVel( );

    // ======== Test before we run
    nff.printAtomParams();
    ff.printAtomPos();
    ff.printBondParams();
    ff.printAngleParams();
    double E = ff.eval(true);
    printf( "iter0 E = %g \n", E );
    printf("TestAppMMFFmini.init() DONE \n");
    //exit(0);

    Draw3D::makeSphereOgl( ogl_sph, 3, 1.0 );

    //selection.insert( selection.end(), {12, 16, 14, 6, 2, 3,   20,18,31,25,26} );
    //selection.insert( selection.end(), {13,29,30} );
    //splitGraphs( ff.nbonds, ff.bond2atom, 12, 22 );
    std::unordered_set<int> innodes;
    //innodes.insert(12);
    innodes.insert(11);
    splitGraphs( ff.nbonds, ff.bond2atom, 22, innodes );
    for( int i : innodes ){ selection.push_back(i); }

}

void TestAppMMFFmini::drawSystem( ){
    //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLines ( ff.nbonds, (int*)ff.bond2atom, ff.apos );
    glColor3f(0.0f,0.0f,0.0f); Draw3D::bondsPBC(  ff.nbonds, ff.bond2atom, ff.apos, &builder.bondPBC[0], builder.lvec );
    glColor3f(0.5f,0.0f,0.0f); Draw3D::atomLabels(  ff.natoms,  ff.apos, fontTex );
    glColor3f(0.0f,0.0f,1.0f); Draw3D::bondLabels( ff.nbonds,       ff.bond2atom, ff.apos, fontTex, 0.02 );
    //glColor3f(1.0f,0.0f,0.0f); Draw3D::vecsInPoss( ff.natoms, ff.aforce, ff.apos, 300.0              );
    //Draw3D::atomsREQ  ( ff.natoms, ff.apos,   nff.REQs, ogl_sph, 1.0, 0.25, 1.0 );
    Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, 0.25, 1.0 );
}

void TestAppMMFFmini::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //printf( "====== Frame # %i \n", frameCount );

    Draw3D::drawAxis(  10. );

	if( ogl_mol ){
        glCallList( ogl_mol );
        return;
        //exit(0);
    }

	//ibpicked = world.pickBond( ray0, camMat.c , 0.5 );
    ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
    Draw3D::drawPointCross( ray0, 0.1 );
    //Draw3D::drawVecInPos( camMat.c, ray0 );
    if(ipicked>=0) Draw3D::drawLine( ff.apos[ipicked], ray0);

    double Ftol = 1e-4;

    //ff.apos[0].set(-2.0,0.0,0.0);
    perFrame = 1;
    //perFrame = 50;
    //bRunRelax = false;

    //Vec3d::rotate( selection.size(), &selection[0], ff.apos, ff.apos[12], (ff.apos[13]-ff.apos[12]).normalized(), 0.1 );
    Vec3d::rotate( selection.size(), &selection[0], ff.apos, ff.apos[13], (ff.apos[11]-ff.apos[13]).normalized(), 0.1 );

    if(bRunRelax){

        //builder.lvec.a.mul(1.001);
        builder.lvec.a.mul(0.999);
        for(int itr=0; itr<perFrame; itr++){
            //if(bConverged) break;
            //printf( "======= frame %i \n", frameCount );

            //printf( "DEBUG run 1 \n" );

            // rotate arom[0]
            //ff.apos[0] = ff.apos[1] + (ff.apos[0]-ff.apos[1]).rotate( 2*M_PI/perFrame, ff.apos[2]-ff.apos[1] );

            double E=0;
            ff.cleanAtomForce();
            E += ff.eval(false);
            if(bNonBonded){
                //E += nff.evalLJQ_sortedMask();   // This is fast but does not work in PBC
                E += nff.evalLJQ_pbc( builder.lvec, {1,1,1} );
            }
            //Vec3d cog,fsum,torq;
            //checkForceInvariatns( ff.natoms, ff.aforce, ff.apos, cog, fsum, torq );
            //printf( "fsum %g torq %g   cog (%g,%g,%g) \n", fsum.norm(), torq.norm(), cog.x, cog.y, cog.z );


            //for(int i=0; i<world.natoms; i++){ world.aforce[i].set(0.0d); }
            //printf( "DEBUG x.1 \n" );
            //world.eval_bonds(true);
            //world.eval_angles();
            //printf( "DEBUG x.2 \n" );
            //world.eval_angles();
            //printf( "DEBUG x.3 \n" );
            //world.eval_LJq_On2();

            //exit(0);
            if(ipicked>=0){
                Vec3d f = getForceSpringRay( ff.apos[ipicked], (Vec3d)cam.rot.c, ray0, -1.0 );
                //printf( "f (%g,%g,%g)\n", f.x, f.y, f.z );
                ff.aforce[ipicked].add( f );
            };

            for(int i=0; i<ff.natoms; i++){
                //ff.aforce[i].add( getForceHamakerPlane( ff.apos[i], {0.0,0.0,1.0}, -3.0, 0.3, 2.0 ) );
                //printf( "%g %g %g\n",  world.aforce[i].x, world.aforce[i].y, world.aforce[i].z );
            }


            //exit(0);

            //for(int i=0; i<world.natoms; i++){ world.aforce[i].add({0.0,-0.01,0.0}); }
            //int ipivot = 0;
            //world.aforce[ipivot].set(0.0);
            //opt.move_LeapFrog(0.01);
            //opt.move_MDquench();

            //opt.move_GD(0.001);
            double f2 = opt.move_FIRE();
            //exit(0);

            /*
            printf( "E %g |F| %g |Ftol %g \n", E, sqrt(f2), Ftol );
            if(f2<sq(Ftol)){
                bConverged=true;
                printf( "CONVERGED \n" );
            }
            */
        }
    }

    Draw3D::drawTriclinicBox(builder.lvec, Vec3dZero, Vec3dOne );
    glColor3f(0.6f,0.6f,0.6f); Draw3D::plotSurfPlane( (Vec3d){0.0,0.0,1.0}, -3.0, {3.0,3.0}, {20,20} );

    if(builder.bPBC){
        //printf( "draw PBC \n" );
        Draw3D::drawPBC( (Vec3i){1,1,0}, builder.lvec, [&](){drawSystem();} );
    }else{
        drawSystem();
    }

    for(int i=0; i<selection.size(); i++){
        int ia = selection[i];
        glColor3f( 0.f,1.f,0.f );
        Draw3D::drawSphereOctLines( 8, 0.3, ff.apos[ia] );
    }


    if(iangPicked>=0){
        glColor3f(0.,1.,0.);
        Draw3D::angle( ff.ang2atom[iangPicked], ff.ang_cs0[iangPicked], ff.apos, fontTex );
    }

    /*
    printf("==========\n");
    for(int i=0; i<world.natoms; i++){
        printf("iatom %i (%g,%g,%g) (%g,%g,%g) \n", i, world.apos[i].x,world.apos[i].y,world.apos[i].z, world.aforce[i].x,world.aforce[i].y,world.aforce[i].z  );
    }
    if(frameCount>=10){STOP = true;}
    */

};


void TestAppMMFFmini::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_p:  first_person = !first_person; break;
                //case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;

                case SDLK_LEFTBRACKET:  ff.i_DEBUG=(ff.i_DEBUG+1)%ff.nang; break;
                case SDLK_RIGHTBRACKET: ff.i_DEBUG=(ff.i_DEBUG+1)%ff.nang; break;

                //case SDLK_v: for(int i=0; i<world.natoms; i++){ ((Vec3d*)opt.vel)[i].add(randf(-drndv,drndv),randf(-drndv,drndv),randf(-drndv,drndv)); } break;
                //case SDLK_p: for(int i=0; i<world.natoms; i++){ world.apos[i].add(randf(-drndp,drndp),randf(-drndp,drndp),randf(-drndp,drndp)); } break;

                //case SDLK_LEFTBRACKET:  if(ibpicked>=0) world.bond_0[ibpicked] += 0.1; break;
                //case SDLK_RIGHTBRACKET: if(ibpicked>=0) world.bond_0[ibpicked] -= 0.1; break;

                //case SDLK_a: world.apos[1].rotate(  0.1, {0.0,0.0,1.0} ); break;
                //case SDLK_d: world.apos[1].rotate( -0.1, {0.0,0.0,1.0} ); break;

                //case SDLK_w: world.apos[1].mul( 1.1 ); break;
                //case SDLK_s: world.apos[1].mul( 0.9 ); break;


                case SDLK_SPACE: bRunRelax=!bRunRelax; break;

                case SDLK_g: iangPicked=(iangPicked+1)%ff.nang;
                    printf( "ang[%i] cs(%g,%g) k %g (%i,%i,%i)\n", iangPicked, ff.ang_cs0[iangPicked].x, ff.ang_cs0[iangPicked].y, ff.ang_k[iangPicked],
                        ff.ang2atom[iangPicked].a,ff.ang2atom[iangPicked].b,ff.ang2atom[iangPicked].c );
                    break;

            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    ipicked = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.natoms, ff.apos );
                    printf( "picked atom %i \n", ipicked );
                    break;
                case SDL_BUTTON_RIGHT:
                    //ibpicked = ff.pickBond( ray0, (Vec3d)cam.rot.c , 0.5 );
                    //printf("ibpicked %i \n", ibpicked);
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    ipicked = -1;
                    break;
                case SDL_BUTTON_RIGHT:
                    //ibpicked = -1;
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
    //STOP = false;
}

void TestAppMMFFmini::drawHUD(){
    glDisable ( GL_LIGHTING );

}

// ===================== MAIN

TestAppMMFFmini * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppMMFFmini( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















