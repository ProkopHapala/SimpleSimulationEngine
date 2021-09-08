#ifndef MMFFBuilder_h
#define MMFFBuilder_h

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>

#include "macroUtils.h"
#include "testUtils.h"

//#include "Molecule.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"

#include "molecular_utils.h"

//#include "Molecule.h"
//#include "MMFF.h"
//#include "MMFFmini.h"
#include "MMFFparams.h"

// =============== Structs for Atoms, Bonds etc...

namespace MM{

static const double const_eVA2_Nm = 16.0217662;

struct Atom{
    constexpr const static Vec3d HcapREQ    = (Vec3d){ 1.4870, sqrt(0.000681    ), 0 };
    constexpr const static Vec3d defaultREQ = (Vec3d){ 1.7,    sqrt(0.0037292524), 0 };

    // this breaks {<brace-enclosed initializer list>} in C++11
    //int type  = -1;
    //int frag  = -1;
    //int iconf = -1;
    //Vec3d pos;
    //Vec3d REQ = defaultREQ;   // constexpr Vec3d{1.7,sqrt(0.0037292524),0}

    int type;
    int frag;
    int iconf;
    Vec3d pos;
    Vec3d REQ;   // constexpr Vec3d{1.7,sqrt(0.0037292524),0}

    //Atom() = default;

    void print()const{ printf( " Atom{ t %i c %i f %i REQ(%g,%g,%g) pos(%g,%g,%g)}", type, iconf, frag, REQ.x, REQ.y, REQ.z, pos.x,pos.y,pos.z ); }

    Atom() = default;
    Atom(const Vec3d& pos_):type(0),frag(-1),iconf(-1),REQ(defaultREQ),pos(pos_){};
    Atom(const Vec3d& pos_,const Vec3d& REQ_):type(0),frag(-1),iconf(-1),REQ(REQ_),pos(pos_){};
    Atom(int type_,int frag_,int iconf_,const Vec3d& pos_,const Vec3d& REQ_):type(type_),frag(frag_),iconf(iconf_),REQ(REQ_),pos(pos_){};
};

#define N_NEIGH_MAX 4
enum class NeighType: int {
    pi    = -2,
    epair = -3,
    H     = -4
};


struct AtomConf{

    int iatom=-1;
    uint8_t n     =0;
    uint8_t nbond =0;
    uint8_t npi   =0; // pi bonds
    uint8_t ne    =0; // electron pairs
    uint8_t nH    =0; //
    int neighs[N_NEIGH_MAX]; // neighs  - NOTE: bonds not atoms !!!!

    //AtomConf() = default;

    inline int findNeigh(int ib){
        for(int i=0; i<N_NEIGH_MAX; i++){
            //printf( "MM:AtomConf.findNeigh()[%i] ng %i ja %i \n", i, neighs[i], ja );
            if(neighs[i]==ib) return i;
        }
        return -1;
    }

    inline bool addNeigh(int ia, uint8_t& ninc ){
        if(n>=N_NEIGH_MAX)return false;
        if(ia>=0){ neighs[nbond]=ia; }else{ neighs[N_NEIGH_MAX-(n-nbond)-1]=ia; };
        ninc++;
        n++;
        //printf( "bond.addNeigh n==%i ninc==%i\n", n, ninc );
        return true;
    };

    inline bool addBond (int i){ return addNeigh(i,nbond); };
    inline bool addH    (     ){ return addNeigh((int)NeighType::H    ,nH ); };
    inline bool addPi   (     ){ return addNeigh((int)NeighType::pi   ,npi); };
    inline bool addEpair(     ){ return addNeigh((int)NeighType::epair,ne ); };

    inline void clearNonBond(){ n=nbond; npi=0;ne=0;nH=0; };
    inline void clearBond   (){ nbond=0; n=npi+ne+nH;     };
    inline void setNonBond(int npi_,int ne_){ npi=npi_; ne=ne_; n=nbond+npi+ne+nH;  }
    inline void init0(){ for(int i=0; i<N_NEIGH_MAX; i++)neighs[i]=-1; nbond=0; clearNonBond(); }

    void print()const{ printf( " AtomConf{ ia %i, n %i nb %i np %i ne %i nH %i (%i,%i,%i,%i) }", iatom, n, nbond, npi, ne, nH , neighs[0],neighs[1],neighs[2],neighs[3] ); }

    AtomConf() = default;
    //AtomConf(int iatom_,int npi_)       :iatom(iatom_),npi(npi_),ne(0  ),nH(0),nbond(0),n(npi_    ){};
    AtomConf(int iatom_,int npi_,int ne_):iatom{iatom_},npi{npi_},ne{ne_},nH{0},nbond{0},n{npi_+ne_}{ for(int i=0;i<N_NEIGH_MAX;i++)neighs[i]=-1; };
    //AtomConf(const MMFFAtomConf&) = default;
    //AtomConf(std::initializer_list<MMFFAtomConf>) {};
};

struct Bond{
    // --- this breaks {<brace-enclosed initializer list>} in C++11
    int    type  = -1;
    Vec2i  atoms = (Vec2i){-1,-1};
    double l0=1,k=0;
    //Vec3<int8_t> ipbc; // for periodic boundary conditions

    //int    type;
    //Vec2i  atoms;
    //double l0,k;

    inline int getNeighborAtom(int ia)const{
        if     (ia==atoms.i){ return atoms.j; }
        else if(ia==atoms.j){ return atoms.i; }
        return -1;
    }

    void print()const{ printf( " Bond{t %i a(%i,%i) l0 %g k %g}", type, atoms.i, atoms.j, l0, k ); };

    Bond()=default;
    Bond(int type_, Vec2i atoms_, double l0_, double k_):type(type_),atoms(atoms_),l0(l0_),k(k_){};
};

struct Angle{

    // --- this breaks {<brace-enclosed initializer list>} in C++11
    //int type     = -1;
    //Vec2i  bonds = (Vec2i){-1,-1};
    //double a0    = 0;
    //double k     = 0;
    //Angle()=default;

    int type;
    Vec2i  bonds;
    double a0;
    double k;

    void print()const{ printf( " Angle{t %i b(%i,%i) a0 %g k %g}", type, bonds.i, bonds.j, a0, k ); }

    Angle()=default;
    Angle( int type_, Vec2i bonds_, double a0_, double k_):type(type_), bonds(bonds_),a0(a0_),k(k_){ };
};


struct Dihedral{

    // --- this breaks {<brace-enclosed initializer list>} in C++11
    //int type     = -1;
    //Vec3i  bonds = (Vec3i){-1,-1,-1};
    //int    n=0;
    //double k=0;

    int    type;
    Vec3i  bonds;
    int    n;
    double k;

    //Dihedral()=default;

    void print()const{ printf( " Dihedral{t %i b(%i,%i,%i) n %i k %g}", type, bonds.a, bonds.b,bonds.c, n, k ); }

    Dihedral()=default;
    Dihedral( int type_, Vec3i  bonds_, int n_, double k_ ):type(type_), bonds(bonds_), n(n_), k(k_){};
};


//struct Molecule{
//    Molecule * mol ;
//    Vec3i      i0  ;
//};
#ifdef Molecule_h
struct Fragment{
    //int imolType;
    Vec2i atomRange;
    Vec2i bondRange;
    Vec2i angRange;
    Vec2i dihRange;
    Vec3d  pos;
    Quat4d rot;
    Molecule * mol;     // ToDo : we can replace this by MolID to leave dependence on Molecule_h
    Vec3d    * pos0s;

    Fragment()=default;
    Fragment(Molecule* mol_, Vec3d pos_, Quat4d rot_, Vec2i atomRange_ ):
        pos(pos_),rot(rot_),
        atomRange{atomRange_},bondRange{0,0},angRange{0,0},dihRange{0,0},
        mol{mol_},pos0s{mol_->pos}{};
};
#endif // Molecule_h

class Builder{  public:

    bool bDEBUG = false;

    //static int iDebug = 0;
    std::vector<Atom>       atoms;
    std::vector<Bond>       bonds;
    std::vector<Angle>      angles;
    std::vector<Dihedral>   dihedrals;

    std::vector<AtomConf>  confs;
    //std::vector<int>  atom_neighs;

    bool bPBC=false;
    Mat3d lvec = Mat3dIdentity;  // lattice vectors for PBC (periodic boundary conditions)
    std::vector<Vec3i> bondPBC;

    MMFFparams* params = 0;  // Each function which needs this can take it as parameter
    std::vector       <std::string>*     atomTypeNames = 0;
    std::unordered_map<std::string,int>* atomTypeDict  = 0;

#ifdef Molecule_h
    //std::vector<Molecule> mols;
    std::vector<Fragment>   frags;
    std::vector<Molecule*>              molTypes;
    std::unordered_map<std::string,int> molTypeDict;
    std::unordered_map<size_t,size_t> fragTypes;
    std::unordered_map<size_t,size_t> mol2molType;
#endif // Molecule_h

    Vec3d defaultREQ  {  1.5, 0.0, 0.0 };
    Bond  defaultBond { -1, {-1,-1}, 1.5, 1.0 };
    //Angle defaultAngle{ -1, {-1,-1}, 0.0, 0.5 };
    Angle defaultAngle{ -1, {-1,-1}, M_PI, 0.5 };

    Bond bondBrush = defaultBond;

    Atom capAtom      = (Atom){ (int)NeighType::H,     -1,-1, {0,0,0}, Atom::HcapREQ };
    Atom capAtomEpair = (Atom){ (int)NeighType::epair, -1,-1, {0,0,0}, {0,0,0} };
    Atom capAtomPi    = (Atom){ (int)NeighType::pi,    -1,-1, {0,0,0}, {0,0,0} };
    Bond capBond      = (Bond){ -1,  {-1,-1},  1.07, 100/const_eVA2_Nm };
    Vec3d    capUp   = (Vec3d){0.0d,0.0d,1.0d};
    bool bDummyPi    = false;
    bool bDummyEpair = false;


    // =================== Functions =====================

    void clear(){
        atoms.clear(); //printf("DEBUG a.1 \n");
        bonds.clear(); //printf("DEBUG a.2 \n");
        angles.clear();
        dihedrals.clear();
#ifdef Molecule_h
        //mols .clear(); //printf("DEBUG a.3 \n");
        frags.clear(); //printf("DEBUG a.4 \n");
        fragTypes.clear();
#endif // Molecule_h
    }

    void bindParams( MMFFparams* params_ ){
        params = params_;
        atomTypeNames = &params->atomTypeNames;
        atomTypeDict  = &params->atomTypeDict;
    }

    void initDefaultAtomTypeDict(){
        makeDefaultAtomTypeDict( atomTypeNames, atomTypeDict );
    }

    // ============= Add Capping Hydrogens

    const AtomConf* getAtomConf(int ia)const{
        int ic=atoms[ia].iconf;
        if(ic>=0){ return &confs[ic]; }
        return 0;
    }

    int getBondToNeighbor( int ia, int ja )const {
        const AtomConf* conf = getAtomConf(ia);
        if(conf){
            for(int i=0; i<conf->nbond; i++){
                int ib  = conf->neighs[i];
                if(ib<0) continue;
                int jai = bonds[ib].getNeighborAtom(ia);
                if(jai==ja){ return ib; }
            }
        }
        return -1;
    }

    inline int getBondByAtoms(int i, int j)const{
        int ib;
        ib = getBondToNeighbor( i, j ); if( ib>=0 ) return ib;
        ib = getBondToNeighbor( j, i ); if( ib>=0 ) return ib;
        return -1;
    }

    AtomConf* addConfToAtom( int ia, const AtomConf* conf=0 ){
        //printf( "MM::Builder.addConfToAtom ia %i \n", ia );
        int ic = confs.size();
        atoms[ia].iconf = ic;
        if   (conf){ confs.push_back( *conf      );                       }
        else       { confs.push_back(  AtomConf()); confs.back().init0(); }
        confs.back().iatom=ia;
        return &confs.back();
    }
    int tryAddConfsToAtoms( int i0=0, int imax=-1 ){
        if(imax<0){ imax=atoms.size(); }
        int n=0;
        for(int ia=0;ia<imax;ia++){
            int ic=atoms[ia].iconf;
            if(ic<0){
                addConfToAtom( ia, 0 );
                n++;
            }
        }
        return n;
    }

    AtomConf* insertAtom(const Atom& atom, AtomConf* conf ){
        atoms.push_back(atom);
        return addConfToAtom( atoms.size()-1, conf );
    }
    void insertAtom(const Atom& atom ){ atoms.push_back(atom); }

    void insertAtom( int ityp, const Vec3d& pos, Vec3d* REQ=0, int npi=-1, int ne=0 ){
        if(REQ==0)REQ=&defaultREQ;
        int iconf=-1;
        if(npi>=0){
            //printf( "insertAtom npi>0 => make Conf \n" );
            iconf=confs.size();
            confs.push_back( AtomConf(atoms.size(), npi, ne ) );
        }
        atoms.push_back( Atom    ( ityp,-1,iconf, pos, *REQ )  );
    }

    void tryAddBondToAtomConf( int ib, int ia, bool bCheck ){
        int ic = atoms[ia].iconf;
        //printf( "MM::Builder.addBondToAtomConf ia %i ib %i ic %i \n", ia, ib, ic );
        if(ic>=0){
            if(bCheck){
                int ing = confs[ic].findNeigh(ib);
                //printf( "ing %i \n", ing );
                if( 0<ing ){
                    //printf( " ia %i ib %i ing %i RETURN \n", ia, ib, ing );
                    return; // neighbor already present in conf
                }
            }
            //printf( "MM::Builder.addBondToAtomConf ia %i ib %i ADDED \n", ia, ib );
            confs[ic].addBond(ib);
        }
    }

    void addBondToConfs( int ib, bool bCheck ){
        const Bond& bond = bonds[ib];
        tryAddBondToAtomConf( ib, bond.atoms.i, bCheck );
        tryAddBondToAtomConf( ib, bond.atoms.j, bCheck );
    }
    void tryAddBondsToConfs( int i0=0, int imax=-1 ){
        if(imax<0) imax=bonds.size();
        for(int ib=i0; ib<imax; ib++){
            addBondToConfs( ib, true );
        }
    }

    void insertBond(const Bond& bond ){
        int ib = bonds.size();
        bonds.push_back(bond);
        tryAddBondToAtomConf( ib, bond.atoms.i, false );
        tryAddBondToAtomConf( ib, bond.atoms.j, false );
        /*
        int ic = atoms[bond.atoms.i].iconf;
        int jc = atoms[bond.atoms.j].iconf;
        //if(ic>=0){ confs[ic].addBond(bond.atoms.j); }
        //if(jc>=0){ confs[jc].addBond(bond.atoms.i); }
        //printf( "insertBond %i(%i,%i) to c(%i,%i) l0 %g k %g\n", ib, bond.atoms.i,bond.atoms.j, ic, jc, bond.l0, bond.k );
        if(ic>=0){
            confs[ic].addBond(ib);
            //printf( "   i.conf " ); println(confs[ic]);
        }
        if(jc>=0){
            confs[jc].addBond(ib);
            //printf( "   j.conf " ); println(confs[jc]);
        }
        */
    }

    //void addCap(int ia,Vec3d& hdir, Atom* atomj, int btype){
    void addCap(int ia,Vec3d& hdir, Atom* atomj ){
        int ja=atoms.size();
        //capAtom;
        Atom atom_tmp;
        if(atomj==0){
            atom_tmp=capAtom;
            atomj=&atom_tmp;
        }
        //if(btype<0) btype=capBond.type;
        atomj->pos = atoms[ia].pos + hdir;
        //atoms.push_back( *atomj );
        insertAtom(*atomj);
        //bonds.push_back( (Bond){btype,{ia,ja}} );
        capBond.atoms.set(ia,ja);
        insertBond( capBond );
        //printf("addCap %i \n", ia );
        //int ic = atoms[ia].iconf;
        //confs[ic].addBond(ja);
        //confs[ic].addBond(bonds.size());
        //println(confs[ic]);
    }

    void makeConfGeom(int nb, int npi, Vec3d* hs){
        Mat3d m;
        if(nb==3){ // defined by 3 sigma bonds
            m.b.set_cross( hs[1]-hs[0], hs[2]-hs[0] );
            m.b.mul( -1/m.b.norm() );
            if(npi==0){ // sp3 no-pi
                if( 0 < m.b.dot( hs[0]+hs[1]+hs[2] ) ){ m.b.mul(-1.); }
                hs[3]=m.b;
            }else{
                hs[3]=m.b;
            }
        }else if(nb==2){ // defined by 2 sigma bonds
            m.fromCrossSafe( hs[0], hs[1] );
            if      (npi==0){ // -CH2- like sp3 no-pi
                const double cb = 0.81649658092; // sqrt(2/3)
                const double cc = 0.57735026919;  // sqrt(1/3)
                hs[nb  ] = m.c*cc+m.b*cb;
                hs[nb+1] = m.c*cc-m.b*cb;
            }else if(npi==1){ // =CH- like  sp 1-pi
                hs[nb  ] = m.c;
                hs[nb+1] = m.b;
                //printf("like =CH- H(%g,%g,%g) pi(%g,%g,%g,) \n", hs[nb].x, hs[nb].y, hs[nb].z, hs[nb+1].x, hs[nb+1].y, hs[nb+1].z );
            }else{            // #C- like sp 2-pi
                hs[nb  ] = m.c;
                hs[nb+1] = m.b;
            }
        }else if(nb==1){
            m.c = hs[0]; m.c.normalize();
            m.c.getSomeOrtho(m.b,m.a);
            if      (npi==0){ // -CH3 like sp3 no-pi
                const double ca = 0.81649658092;  // sqrt(2/3)
                const double cb = 0.47140452079;  // sqrt(2/9)
                const double cc =-0.33333333333;  // 1/3
                hs[nb  ] = m.c*cc + m.b*(cb*2) ;
                hs[nb+1] = m.c*cc - m.b* cb    + m.a*ca;
                hs[nb+2] = m.c*cc - m.b* cb    - m.a*ca;
            }else if(npi==1){ // =CH2 like sp2 1-pi
                const double ca = 0.87758256189;  // 1/2
                const double cc =-0.5;            // sqrt(1/8)
                hs[nb  ] = m.c*cc + m.a*ca;
                hs[nb+1] = m.c*cc - m.a*ca;
                hs[nb+2] = m.b;
            }else{            // #CH sp  2-pi
                hs[nb  ] = m.c*-1;
                hs[nb+1] = m.b;
                hs[nb+2] = m.a;
            }
        }else if(nb==0){
            m.c = hs[0]; m.c.normalize();
            m.c.getSomeOrtho(m.b,m.a);
            if      (npi==0){ //  CH4 like sp3 no-pi
                const double ca = 0.81649658092;  // sqrt(2/3)
                const double cb = 0.47140452079;  // sqrt(2/9)
                const double cc =-0.33333333333;  // 1/3
                hs[nb  ] = m.c*cc + m.b*(cb*2) ;
                hs[nb+1] = m.c*cc - m.b* cb    + m.a*ca;
                hs[nb+2] = m.c*cc - m.b* cb    - m.a*ca;
                hs[nb+3] = m.c;
            }
        }
    }

    void makeSPConf(int ia,int npi,int ne){
        //if(nH==0){ // ToDo : check reasonable limits npi, nh
        int ic = atoms[ia].iconf;
        AtomConf& conf = confs[ic];
        conf.clearNonBond();
        int nb   = conf.nbond;
        int ncap = 4-nb-npi;   // number of possible caps
        int nH   = ncap-ne;
        //printf("-- "); println(conf);
        //printf( "ia %i nb,npi %i,%i   n,nH,ne %i,%i,%i \n", ia,   nb,npi,  n,nH,ne );
        //Mat3d m;
        Vec3d hs[4];
        for(int i=0;i<nb;i++){
            int ib = conf.neighs[i];
            int ja = bonds[ib].getNeighborAtom(ia);
            hs[i]  = atoms[ja].pos - atoms[ia].pos;
            hs[i].normalize();
        }
        makeConfGeom(conf.nbond, npi, hs);
        bool Hmask[]{1,1,1,1};
        //if(nH!=ncap) Hmask[rand()%ncap]=0;
        //bool breverse = (nH==2)&&(ncap==3);
        bool breverse;
        if(ncap<4){
            if(ne>0) Hmask[rand()%ncap]=0;
            breverse = (ne>1);
        }else{
            for(int i=0;i<ne;i++)Hmask[3-i]=0;
            breverse = 0;
        }
        //printf( "makeSPConf: atom[%i] ncap %i nH %i nb %i npi %i ne %i Hmask{%i,%i,%i,%i}  \n", ia, ncap, nH, nb,npi,ne,  (int)Hmask[0],(int)Hmask[1],(int)Hmask[2],(int)Hmask[3] );
        for(int i=0; i<ncap; i++){
            if     (Hmask[i]!=breverse){ addCap(ia,hs[i+nb],&capAtom       ); }
            else if(bDummyEpair       ){ addCap(ia,hs[i+nb],&capAtomEpair  ); }
        }
        if(bDummyPi){ for(int i=0; i<npi; i++){ addCap(ia,hs[i+ncap+nb],&capAtomPi); } }
        conf.npi=npi;
        conf.ne =ne;
        //printf("-> "); println(conf);
    }

    bool tryMakeSPConf(int ia){
        const AtomConf* conf = getAtomConf(ia);
        //printf("tryMakeSPConf %i conf %li\n", ia, (long)conf  );
        if(conf){
            //printf("tryMakeSPConf: proceed !!! \n"  );
            makeSPConf(ia,conf->npi,conf->ne);
            return true;
        }
        return false;
    }

    int makeAllConfsSP(){
        int n=0,na=atoms.size();
        for(int i=0;i<na;i++){
            if(tryMakeSPConf(i)){n++;}
        }
        return n;
    }

    // ============= Angles

    void addAnglesToBond( int ib, int n,const int* neighs, double a0, double k ){
        for(int j=0; j<n; j++){
            angles.push_back( (Angle){-1,  (Vec2i){ neighs[ib], neighs[j]}, a0,k} );
            //printf("[%li]",angles.size()); println(angles.back());
        }
    }

    void addAnglesUpToN( int n, const int* neighs, double a0, double k ){
        for(int i=0; i<n; i++){ addAnglesToBond( i, i, neighs, a0, k ); }
    }

    bool addAnglesToAtom( int ia, double ksigma, double kpi ){
        const AtomConf* conf = getAtomConf(ia);
        if(conf==0) return false;
        int nsigma = conf->nbond;
        //printf("addAnglesToAtom[%i] nsigma %i npi %i \n", ia, nsigma, conf->npi  );
        // ------ Pi Bonds
        if( bDummyPi && (conf->npi>0) ){
            //printf( "addAnglesToAtom[%i] angles to dummy Pi-orbital \n", ia );
            //nsigma -= conf->npi; // ToDo : Check this
            for(int i=0;i<conf->npi;i++){ addAnglesToBond( i+nsigma, i+nsigma, conf->neighs, M_PI_2, kpi ); }
        }
        // ------- Sigma bonds
        static const double a0s[]{ 0.0, 0.0, M_PI, 120*M_PI/180, 109.5*M_PI/180 };
        static const double Ks []{ 0.0, 0.0, 1.4,           1.2,            1.0 };
        if(nsigma>=2){
            int iangType = nsigma+conf->ne;
            double a0 = a0s[iangType];
            ksigma   *=  Ks[iangType];
            //printf( "addAnglesToAtom[%i] ns %i npi %i a0,ks %g %g   {%g,%g,%g,%g} %g \n", ia, nsigma, conf->npi, a0, ksigma, a0s[0],a0s[1],a0s[2],a0s[3] , a0s[nsigma] );
            addAnglesUpToN( nsigma, conf->neighs, a0, ksigma );
            if(bDEBUG){
                printf("addAnglesToAtom[%i] ", ia); printAtomConf(ia);
                printf( " Sigma(%i,a0[%i] %g, k %g)", ia,  nsigma,conf->npi, iangType,  a0*(180/M_PI), ksigma );
                if(bDummyPi && (conf->npi>0) ){ printf( " Pi(n%i,a0 %g, k %g)", conf->npi,  M_PI_2*(180/M_PI), kpi ); }else{ puts(""); };
            }
        }
        return true;
    }

    void autoAngles(double ksigma, double kpi){
        for(int i=0; i<atoms.size(); i++){
            if(atoms[i].iconf>=0){
                addAnglesToAtom( i, ksigma, kpi );
            }
        }
    }

    // =============== Dihedrals

    bool insertDihedralByAtom(const Quat4i& ias, Dihedral& tors ){
        int ib1 = getBondByAtoms(ias.x,ias.y); if(ib1<0) return false;
        int ib2 = getBondByAtoms(ias.y,ias.z); if(ib2<0) return false;
        int ib3 = getBondByAtoms(ias.z,ias.w); if(ib3<0) return false;
        tors.bonds.set(ib1,ib2,ib3);
        dihedrals.push_back(tors);
        return true;
    }

    /*
    void setAtoms( const Atom* brushAtom, const Vec3d* ps=0, int imin=0, int imax=-1,  ){
        if(imax<0){ imax=atoms.size(); }
        for(int i=imin;i<imax;i++){
            if(ps       )brushAtom.pos = ps[i];
            if(brushAtom)builder.insertAtom( *brushAtom, true);
        }
    }
    void setBondTypes( Bond brushBond, const Vec2i* bond2atom=0, int imin=0, int imax=-1 ){
        if(imax<0){ imax=bonds.size(); }
        for(int i=imin;i<imax;i++){
            if(bond2atom)brushBond.atoms=bond2atom[i];
            builder.insertBond(brushBond);
        }
    }
    */
    void insertAtoms( int n, Atom brushAtom, const Vec3d* ps, bool withConf=true ){
        for(int i=0;i<n;i++){
            brushAtom.pos = ps[i];
            if(withConf){ insertAtom( brushAtom, 0 ); }
            else        { insertAtom( brushAtom    ); }
        }
    }
    void insertBonds( int n, Bond brushBond, const Vec2i* bond2atom ){
        for(int i=0;i<n;i++){
            brushBond.atoms=bond2atom[i];
            insertBond( brushBond );
        }
    }
    void setConfs( int npi, int ne, int imin=0, int imax=-1 ){
        if(imax<0){ imax=atoms.size(); }
        for(int i=imin;i<imax;i++){
            makeSPConf( i, npi, ne );
        }
    }

    void touchingAtoms( int i0, int imax, const Vec3d& p, double R0, double Rfac, std::vector<int>& found ){
        for(int i=i0; i<imax; i++){  // for pbc we need all atom pairs
            const Atom& A = atoms[i];
            Vec3d dp = A.pos - p; // pbc here
            double R = (R0 + A.REQ.x)*Rfac;
            //printf( "[%i,%i] R %g \n", i0-1, i, dp.norm() );
            if(  dp.norm2() < (R*R) ){
                printf( "[%i,%i] R %g \n", i0-1, i, dp.norm() );
                found.push_back(i);
            }
        }
    }

    void autoBonds( double R=-0.5, int i0=0, int imax=-1 ){
        // ToDo : periodic boundary conditions
        if(imax<0)imax=atoms.size();
        bool byParams = (R<0);
        double Rfac=-R;
        //if( byParams && (params==0) ){ printf("ERROR in MM::Builder.autoBonds() byParams(R<0) but params==NULL \n"); exit(0); }
        for(int i=i0; i<imax; i++){
            const Atom& A = atoms[i];
            for(int j=i+1; j<imax; j++){  // for pbc we need all atom pairs
                const Atom& B = atoms[j];
                Vec3d dp = B.pos - A.pos; // pbc here
                if(byParams){ R = (B.REQ.x + A.REQ.x)*Rfac; }
                if(  dp.norm2() < (R*R) ){
                    bondBrush.atoms={i,j};
                    insertBond( bondBrush );
                }
            }
        }
    }

    inline Vec3d pbcShift( Vec3i G ){ return lvec.a*G.a + lvec.b*G.b + lvec.c*G.c; }

    void autoBondsPBC( double R=-0.5, int i0=0, int imax=-1, Vec3i npbc=Vec3iOne ){
        bPBC = true;
        // ToDo : periodic boundary conditions
        if(imax<0)imax=atoms.size();
        bool byParams = (R<0);
        double Rfac=-R;
        //if( byParams && (params==0) ){ printf("ERROR in MM::Builder.autoBonds() byParams(R<0) but params==NULL \n"); exit(0); }
        std::vector<int> found;
        for(int i=i0; i<imax; i++){
            const Atom& A = atoms[i];
            R = A.REQ.x;
            int ipbc=0;
            //printf( "#==== Atom[%i] \n", i );
            for(int ix=-npbc.x;ix<=npbc.x;ix++){
                for(int iy=-npbc.y;iy<=npbc.y;iy++){
                    for(int iz=-npbc.z;iz<=npbc.z;iz++){
                        int   j0=i+1;
                        //Vec3d vpbc = lvec.a*ix + lvec.b*iy + lvec.c*iz;
                        //Vec3d vpbc; lvec.dot_to_T( {(double)ix,(double)iy,(double)iz} );
                        //Vec3d p = A.pos - pbcShift( {ix,iy,iz} );
                        Vec3d p = A.pos - lvec.lincomb( ix, iy, iz );
                        //printf( "# pbc[%i,%i,%i][%i] v(%g,%g,%g)\n", ix,iy,iz, ipbc,  vpbc.x, vpbc.y, vpbc.z );
                        found.clear();
                        touchingAtoms( j0, imax, p, R, Rfac, found );
                        for(int j:found){
                            bondBrush.atoms={i,j};
                            insertBond( bondBrush );
                            bondPBC.push_back( {ix,iy,iz} );
                        }
                        ipbc++;
                    }
                }
            }
        }
    }

    // =============== Configurations

    int checkBond2Conf(bool bPrint)const{
        for(int i=0;i<bonds.size(); i++){
            //printf("checkBond2Conf b[%i]\n", i );
            const Bond& b = bonds[i];
            int i_ = getBondByAtoms(b.atoms.i,b.atoms.j);
            if(i_!= i){
                if(bPrint){
                    printf( "MMFFbuilder.checkBond2Conf: getBondByAtoms(bond[%i/%li]) returned %i \n", i,bonds.size(), i_ );
                }
                return i;
            }
        }
        return -1;
    }

    int checkConf2Bond(bool bPrint)const{
        int nb=0;
        std::vector<int> ng(atoms.size(), 0);
        for(const Bond& b: bonds){ ng[b.atoms.i]++; ng[b.atoms.j]++; };
        for(int ia=0;ia<atoms.size(); ia++){
            //printf("checkConf2Bond[%i] \n", ia );
            const AtomConf* conf = getAtomConf(ia); // we need to modify it
            if(conf==0){
                if( nb<bonds.size() ){
                    if(bPrint){
                        printf( "MMFFbuilder.checkConf2Bond: atom[%i/%li].conf==null nb(%i)<bonds.size(%li) \n", ia, atoms.size(), nb,bonds.size()  );
                    }
                    return ia;
                } else continue;
            }
            int nbconf = conf->nbond;
            if(nbconf != ng[ia] ){
                    if(bPrint){
                        printf( "MMFFbuilder.checkConf2Bond: atom[%i/%li].conf.nbond==%i but counted %i bonds \n", ia, atoms.size(), nbconf, ng[ia] );
                        println( (*conf) );
                    }
                    return ia;
            }
            for(int j=0; j<nbconf; j++){
                int ib = conf->neighs[j];
                int ja = bonds[ib].getNeighborAtom(ia);
                if(ja<0){
                    if(bPrint){
                        printf( "MMFFbuilder.checkConf2Bond: atom[%i/%li].neighs[%i/%i]->bonds[%i/%li].getNeighborAtom(%i) returned %i \n", ia,atoms.size(), j,nbconf, ib,bonds.size(), ia, ja );
                        println( (*conf)   );
                        println( bonds[ib] );
                    }
                    return ia;
                }
            }
            //printf("checkConf2Bond[%i] nb %i \n", ia, nb );
            nb+=nbconf;
        }
        return -1;
    }

    bool checkBondsSorted( int iPrint=0 )const{
        int ia=-1,ja=-1;
        if(iPrint>1)printf("checkBondsSorted %li \n", bonds.size() );
        for(int i=0;i<bonds.size(); i++){
            const Vec2i& b = bonds[i].atoms;
            if(iPrint>1)printf( "pair[%i] %i,%i | %i %i  | %i %i %i \n", i, b.i, b.j,   ia,ja ,   b.i>=b.j,  b.i<ia, b.j<=ja );
            if(b.i>=b.j){ if(iPrint>0){ printf("b.i>=b.j b[%i](%i,%i) ia,ja(%i,%i)\n", i,b.i,b.j,ia,ja); }; return false; }
            if(b.i<ia)  { if(iPrint>0){ printf("b.i<ia   b[%i](%i,%i) ia,ja(%i,%i)\n", i,b.i,b.j,ia,ja); }; return false; }
            else if (b.i>ia){ia=b.i; ja=-1; };
            if(b.j<=ja){  if(iPrint>0){ printf("b.j<=ja  b[%i](%i,%i) ia,ja(%i,%i)\n", i,b.i,b.j,ia,ja); }; return false; }
            ja=b.j;
        }
        if(iPrint>1)printf("checkBondsSorted DONE !\n");
        return true;
    }

    void sortAtomsOfBonds(){
        for(int i=0; i<bonds.size(); i++){ bonds[i].atoms.order(); }
    }


    bool sortBonds(){
        //printf( "sortBonds \n" );
        // sort bonds so that
        //   1) (b.i<b.j)
        //   1) if(bk.i) (b.i<b.j)
        //   1) (b.i<b.j)

        //int bsort    = new[bonds.size()];
        //Bond * bback = new Bond[bonds.size()];
        //int *   invBsort = new int     [bonds.size()];

        // use smart pointer to solve problems with delete[] when return on fail
        std::unique_ptr<Bond[]> bback   (new Bond[bonds.size()]);
        std::unique_ptr<int []> invBsort(new int [bonds.size()]);

        int nga[N_NEIGH_MAX];
        int ngb[N_NEIGH_MAX];

        sortAtomsOfBonds();
        //printBonds();

        int nb=0;
        for(int ia=0; ia<atoms.size(); ia++ ){
            // assume atoms with conf are first, capping are later
            //const AtomConf* conf = getAtomConf(ia);
            if(nb>=bonds.size())break;
            AtomConf* conf = (AtomConf*)getAtomConf(ia); // we need to modify it
            if(!conf){
                printf( "ERROR in MMFF.sortBonds(): atom[%i/%li] without conf (confs.size(%li)) met before all bonds enumerated nb(%i)<bonds.size(%li) \n", ia, atoms.size(), confs.size(), nb, bonds.size() );
                printf( " This algorithm assumes all atoms with conf precede atoms without confs in the array \n" );
                printf( " => return \n" );
                return false;
            }
            int nbconf=conf->nbond;
            int * neighs = conf->neighs;
            //printf( "ia %i nb %i conf.nb %i\n", ia, nb, nbconf );
            for(int i=0;i<nbconf;i++ ){
                int ib=neighs[i];
                if(ib<0){ printf("ERROR in MMFF.sortBonds(): atom[%i].condf inconsistent nbond=%i neigh[%i]<0 \n", ia, conf->nbond, i ); return false; }
                int ja = bonds[ib].getNeighborAtom(ia);
                //if(ja<ia)continue; // this bond was processed before
                nga[i]=ja;
                ngb[i]=ib;
            }
            int ja=-1;
            for(int i=0;i<nbconf;i++ ){      // take bonds on atom in order
                int ipick = selectMinHigher(ja, nbconf, nga );
                ja=nga[ipick];
                //neighs[i] = ngb[ipick]; // make conf sorted
                //printf( " atom[%i].neigh[%i] %i \n", ia, i, ja  );
                //printf( " atom[i %i -> j %i] ng %i \n", ia, ja, i  );
                if(ja<ia)continue;      // this bond was processed before (Hopefully)
                int ib = ngb[ipick];

                //bsort   [nb]=ib;
                bback[nb]   = bonds[ib];
                invBsort[ib]=nb;
                //printf( " bond[%i] -> bond[%i] \n", ib, nb );
                nb++;
            }
            // clean conf so it can be re-inserted
            conf->nbond=0;
            conf->n-=nbconf;
        }
        bonds.clear();
        for(int i=0; i<nb;i++){
            bback[i].atoms.order();
            //bonds[i]=bback[i];
            insertBond( bback[i] );
            //printf( " bond[%i] (%i,%i) \n", i, bback[i].atoms.i, bback[i].atoms.j );
        }
        for(int i=0; i<angles.size();i++){
            Vec2i& bs = angles[i].bonds;
            bs.a = invBsort[bs.a];
            bs.b = invBsort[bs.b];
        }
        for(int i=0; i<dihedrals.size();i++){
            Vec3i& bs = dihedrals[i].bonds;
            bs.a = invBsort[bs.a];
            bs.b = invBsort[bs.b];
            bs.c = invBsort[bs.c];
        }
        return true;
    }

    bool trySortBonds( int iPrint=0 ){
        bool sorted = checkBondsSorted( iPrint );
        if( !sorted ){
            if( !sortBonds() ){
                printf( " ERROR in builder.sortBonds() => exit \n" );
                exit(0);
            }
        }
        return sorted;
    }

    void printAtoms(){
        printf(" # MM::Builder.printAtoms() \n");
        for(int i=0; i<atoms.size(); i++){
            printf("atom[%i]",i); atoms[i].print(); puts("");
        }
    }
    void printBonds(){
        printf(" # MM::Builder.printBonds() \n");
        for(int i=0; i<bonds.size(); i++){
            printf("bond[%i]",i); bonds[i].print(); if(bPBC)printf(" pbc(%i,%i,%i)",bondPBC[i].x,bondPBC[i].y,bondPBC[i].z); puts("");
        }
    }
    void printAngles(){
        printf(" # MM::Builder.printAngles() \n");
        for(int i=0; i<angles.size(); i++){
            printf("angle[%i]", i); angles[i].print(); puts("");
        }
    }
    void printConfs(){
        printf(" # MM::Builder.printConfs() \n");
        for(int i=0; i<confs.size(); i++){
            printf("conf[%i]", i); confs[i].print(); puts("");
        }
    }
    void printAtomConf(int i){
        const Atom& A = atoms[i];
        printf("atom[%i] T %i ic %i ", i, A.type, A.iconf);
        if(A.iconf>=0){
            const AtomConf& c = confs[A.iconf];
            printf(" Conf[%i] n %i nb %i npi %i ne %i nH %i ", A.iconf, c.n, c.nbond, c.npi, c.ne, c.nH );
        }
    }
    void printAtomConfs(){
        printf(" # MM::Builder.printAtomConfs() \n");
        for(int i=0; i<atoms.size(); i++){ printAtomConf(i); puts(""); }
    }

    int write2xyz( FILE* pfile, const char* comment="#comment" ){
        //write2xyz( pfile, atoms.size(), int* atypes, Vec3d* pos, const std::unordered_map<std::string,int>& atomTypeDict, const char* comment="#comment" ){
        fprintf(pfile, "%i\n",  atoms.size() );
        fprintf(pfile, "%s \n", comment      );
        for(int i=0; i<atoms.size(); i++){
            int ityp         = atoms[i].type;
            const Vec3d&  pi = atoms[i].pos;
            //printf( "write2xyz %i %i (%g,%g,%g) %s \n", i, ityp, pi.x,pi.y,pi.z, params->atypes[ityp].name );
            fprintf( pfile, "%s   %15.10f   %15.10f   %15.10f \n", (*atomTypeNames)[ityp].c_str(), pi.x,pi.y,pi.z );
        };
        return atoms.size();
    }

    int save2xyz( const char * fname, const char* comment="#comment"  ){
        FILE* pfile = fopen(fname, "w");
        if( pfile == NULL ) return -1;
        int n = write2xyz(pfile, comment );
        fclose(pfile);
        return n;
    }

    int saveMol( const char* fname ){
        FILE* pfile = fopen(fname, "w");
        if( pfile == NULL ) return -1;
        fprintf(pfile, "\n" );
        fprintf(pfile, "SimpleSimulationEngine::MMFFBuilder:: saveMol()\n" );
        fprintf(pfile, "\n" );
        fprintf(pfile, "%3i%3i  0  0  0  0  0  0  0  0999 V2000\n", atoms.size(), bonds.size() );
        for(int i=0; i<atoms.size(); i++){
            int ityp         = atoms[i].type;
            const Vec3d&  pi = atoms[i].pos;
            fprintf( pfile, "%10.4f%10.4f%10.4f %-3s 0  0  0  0  0  0  0  0  0  0  0  0\n",  pi.x,pi.y,pi.z, (*atomTypeNames)[ityp].c_str() );
        }
        for(int i=0; i<bonds.size(); i++){
            int ityp          = bonds[i].type;   if(ityp==-1) ityp=1;
            const Vec2i&  ats = bonds[i].atoms;
            fprintf( pfile, "%3i%3i%3i  0  0  0  0\n",  ats.a+1, ats.b+1, ityp );
        }
        fclose(pfile);
        return atoms.size() + bonds.size();
    }

    int load_xyz( const char * fname, bool noH=false, bool bConf=true, bool bDebug=false ){
        if(bDebug)printf( "MM::Builder.load_xyz(%s)\n", fname );
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        int natoms; Vec3d pos,REQ=defaultREQ; char at_name[8]; int npi,ne=0;
        const int nbuf=1024;
        char buff[nbuf]; char* line;
        line = fgets( buff, nbuf, pFile ); // number of atoms
        sscanf( line, "%i", &natoms );
        if(bDebug)printf( "natoms %i \n", natoms );
        line = fgets( buff, nbuf, pFile ); // comment, ignore
        int n0 = atoms.size();
        for(int i=0; i<natoms; i++){
            line     = fgets( buff, nbuf, pFile ); // comment, ignore
            int nret = sscanf( line,       "%s %lf %lf %lf %lf %i  ",     at_name, &pos.x, &pos.y, &pos.z, &REQ.z, &npi );
            if(bDebug)printf   (  ".xyz[%i] %s %lf %lf %lf %lf %i\n", i, at_name,  pos.x,  pos.y,  pos.z,  REQ.z,  npi  );
            if( nret < 5 ){ REQ.z=0;  };
            if( nret < 6 ){ npi  =-1; };
            auto it = atomTypeDict->find( at_name );
            if( it != atomTypeDict->end() ){
                int ityp=it->second;
                if(params){
                    params->assignRE( ityp, REQ );
                    ne = params->atypes[ityp].nepair();
                }
                if( noH && (at_name[0]='H') && (at_name[1]='\0')  ) continue;
                insertAtom( it->second, pos, &REQ, npi, ne );
            }
        }
        return atoms.size() - n0;
    }

    inline void natom_def(int& n,int i0)const{ if(n<0){ n=atoms .size()-i0; }; }
    inline void nbond_def(int& n,int i0)const{ if(n<0){ n=bonds .size()-i0; }; }
    inline void nang_def (int& n,int i0)const{ if(n<0){ n=angles.size()-i0; }; }
    inline void ndih_def (int& n,int i0)const{ if(n<0){ n=dihedrals.size()-i0; }; }

    void export_REQs(Vec3d* REQs, int i0=0, int n=-1)const{
        natom_def(n,i0);
        //_allocIfNull(REQs,n);
        for(int i=0; i<n; i++){ REQs[i]= atoms[i0+i].REQ; }
    }

    void export_apos(Vec3d* apos, int i0=0, int n=-1)const{
        natom_def(n,i0);
        //_allocIfNull(apos,n);
        for(int i=0; i<n; i++){ apos[i]= atoms[i0+i].pos; }
    }

    void export_atypes(int*& atypes, int i0=0, int n=-1)const{
        natom_def(n,i0);
        _allocIfNull(atypes,n);
        for(int i=0; i<n; i++){ atypes[i]= atoms[i0+i].type; }
    }

    void export_bonds(Vec2i* b2a, double* l0s=0, double* ks=0, int i0=0, int n=-1)const{
        nbond_def(n,i0);
        for(int i=0; i<n; i++){
            const Bond& b  = bonds[i0+i];
            b2a[i] = b.atoms;
            if(ks )ks [i] = b.k;
            if(l0s)l0s[i] = b.l0;
        }
    }

    void export_angles(Vec2i* a2b, double* a0s=0, Vec2d* cs0s=0, double* ks=0, int i0=0, int n=-1)const{
        nang_def(n,i0);
        for(int i=0; i<n; i++){
            const Angle& a = angles[i0+i];
            a2b[i] = a.bonds;
            if(ks )ks [i] = a.k;
            if(a0s)a0s[i] = a.a0;
            if(cs0s)cs0s[i].fromAngle( a.a0 * 0.5 ); // NOTE: we divide angle by 2
        }
    }

    void export_dihedrals(Vec3i* d2b, int* tn=0, double* ks=0, int i0=0, int n=-1)const{
        ndih_def(n,i0);
        for(int i=0; i<n; i++){
            const Dihedral& d = dihedrals[i0+i];
            d2b[i] = d.bonds;
            if(ks)ks[i] = d.k;
            if(tn)tn[i] = d.n;
        }
    }

#ifdef Molecule_h
    void clearMolTypes( bool deep ){
        if(deep){ for(Molecule* mol : molTypes ){ mol->dealloc(); delete mol; } }
        molTypeDict.clear();
        molTypes.clear();
    }

    void assignAtomREQs( const MMFFparams* params ){
        for(int i=0; i<atoms.size(); i++){
            //mmff->aLJq [i]  = atoms[i].type;
            int ityp = atoms[i].type;
            atoms[i].REQ.x = params->atypes[ityp].RvdW;
            atoms[i].REQ.y = params->atypes[ityp].EvdW;
            atoms[i].REQ.z = 0;
            //atomTypes[i]  = atoms[i].type;
        }
    }

    int loadMolTypeXYZ(const char* fname, const MMFFparams* params ){
        Molecule* mol = new Molecule();      //printf( "DEBUG 1.1.1 \n" );
        mol->atomTypeDict = &params->atomTypeDict; //printf( "DEBUG 1.1.2 \n" );
        //printf("mol->atypNames %i %i \n", mol->atypNames, &params->atypNames );
        mol->loadXYZ( fname );             //printf( "DEBUG 1.1.3 \n" );
        if(params) params->assignREs( mol->natoms, mol->atomType, mol->REQs ); //printf( "DEBUG 1.1.4 \n" );
        int ityp = molTypes.size();
        mol2molType[(size_t)mol]=ityp;
        molTypes.push_back(mol);  //printf( "DEBUG 1.1.5 \n" );
        return molTypes.size()-1; //printf( "DEBUG 1.1.6 \n" );
    }

    int registerRigidMolType( int natoms, Vec3d* pos, Vec3d* REQs, int* atomType ){
        Molecule* mol = new Molecule();
        mol->allocate( natoms, 0 );
        for(int i=0; i<mol->natoms; i++){ mol->pos[i]=pos[i]; mol->REQs[i]=REQs[i]; mol->atomType[i]=atomType[i]; }
        int ityp = molTypes.size();
        mol2molType[(size_t)mol]=ityp;
        molTypes.push_back(mol);
        return molTypes.size()-1;
    }

    int loadMolType(const std::string& fname, const std::string& label, const MMFFparams* params ){
        //printf( "fname:`%s` label:`%s` \n", fname.c_str(), label.c_str()  );
        int itype = loadMolTypeXYZ( fname.c_str(), params );
        //printf( "fname:`%s` label:`%s` itype %i \n", fname.c_str(), label.c_str(), itype  );
        molTypeDict[label] = itype;
        return itype;
    };

    void insertFlexibleMolecule_ignorH( Molecule * mol, const Vec3d& pos, const Mat3d& rot, int iH = 1 ){
        int natom0  = atoms.size();
        int nbond0  = bonds.size();
        std::vector<int> atomInds(mol->natoms);
        std::vector<int> bondInds(mol->nbonds);
        for(int i=0; i<mol->natoms; i++){
            if( mol->atomType[i]==iH ) continue;
            atomInds[i] = atoms.size();
            Vec3d  REQi = mol->REQs[i];   REQi.y = sqrt(REQi.y);
            Vec3d p; rot.dot_to(mol->pos[i],p); p.add( pos );
            atoms.push_back( (Atom){mol->atomType[i], -1, -1, p, REQi } );
        }
        for(int i=0; i<mol->nbonds; i++){
            //bonds.push_back( (Bond){mol->bondType[i], mol->bond2atom[i] + ((Vec2i){natom0,natom0}), defaultBond.l0, defaultBond.k } );
            bondInds[i] = bonds.size();
            const Vec2i& b = mol->bond2atom[i];
            bonds.push_back( Bond(mol->bondType[i], { atomInds[b.a],atomInds[b.b] }, defaultBond.l0, defaultBond.k ) );
        }
        for(int i=0; i<mol->nang; i++){
            const Vec2i& ang = mol->ang2bond[i];
            angles.push_back( (Angle){ 1, { bondInds[ang.a],bondInds[ang.b] }, defaultAngle.a0, defaultAngle.k } );
        }
    }

    void insertFlexibleMolecule( Molecule * mol, const Vec3d& pos, const Mat3d& rot ){
        //printf( "# MM::Builder::insertFlexibleMolecule \n" );
        int natom0  = atoms.size();
        int nbond0  = bonds.size();
        for(int i=0; i<mol->natoms; i++){
            //Vec3d LJq = (Vec3d){0.0,0.0,0.0};  // TO DO : LJq can be set by type
            //Vec3d LJq = (Vec3d){1.0,0.03,0.0}; // TO DO : LJq can be set by type
            Vec3d  REQi = mol->REQs[i];   REQi.y = sqrt(REQi.y);
            Vec3d p; rot.dot_to(mol->pos[i],p); p.add( pos );
            //printf( "insertAtom[%i] pos(%g,%g,%g) -> p(%g,%g,%g)\n", i, mol->pos[i].x, mol->pos[i].y, mol->pos[i].z, p.x,p.y,p.z );
            atoms.push_back( (Atom){mol->atomType[i], -1, -1, p, REQi } );
        }
        for(int i=0; i<mol->nbonds; i++){
            //bonds.push_back( (Bond){mol->bondType[i], mol->bond2atom[i] + ((Vec2i){natom0,natom0}), defaultBond.l0, defaultBond.k } );
            bonds.push_back( Bond(mol->bondType[i], mol->bond2atom[i] + ((Vec2i){natom0,natom0}), defaultBond.l0, defaultBond.k ) );
        }
        for(int i=0; i<mol->nang; i++){
            double alfa0 = defaultAngle.a0;
            if( mol->ang0s ) alfa0 = mol->ang0s[i];
            angles.push_back( (Angle){ 1, mol->ang2bond[i] + ((Vec2i){nbond0,nbond0}), alfa0, defaultAngle.k } );
            //printf( "angle[%i|%i,%i] %g|%g %g \n", i, angles.back().bonds.a, angles.back().bonds.b, angles.back().a0, alfa0, angles.back().k );
        }
    }

    int insertRigidMolecule( Molecule * mol, const Vec3d& pos, const Mat3d& rot ){
        int natoms0 = atoms.size();
        Quat4d qrot; qrot.fromMatrix(rot);
        int ifrag = frags.size();
        //printf( "insertMolecule mol->natoms %i \n", mol->natoms );
        for(int i=0; i<mol->natoms; i++){
            //Vec3d REQi = (Vec3d){1.0,0.03,mol->}; // TO DO : LJq can be set by type
            //atoms.push_back( (Atom){mol->atomType[i],mol->pos[i], LJq } );
            Vec3d  REQi = mol->REQs[i];   REQi.y = sqrt(REQi.y); // REQi.z = 0.0;
            Vec3d  p; rot.dot_to(mol->pos[i],p); p.add( pos );
            atoms.push_back( (Atom){mol->atomType[i], ifrag, -1, p, REQi } );
        }
        //frags.push_back( (Fragment){natoms0, atoms.size()-natoms0, pos, qrot, mol}  );
        frags.push_back( Fragment{ mol, pos, qrot,  {natoms0, atoms.size()-natoms0} } );
        //size_t mol_id = static_cast<size_t>(mol);
        size_t mol_id = (size_t)(mol);
        auto got = fragTypes.find(mol_id);
        if ( got == fragTypes.end() ) {
            fragTypes[ mol_id ] = frags.size()-1; // WTF ?
        }else{}
        return ifrag;
    }

    int insertMolecule( Molecule * mol, const Vec3d& pos, const Mat3d& rot, bool rigid, int noH=-1 ){
        //int natom0  = atoms .size();
        //int nbond0  = bonds .size();
        //int nangle0 = angles.size();
        //mols.push_back( (MMFFmol){mol, (Vec3i){natom0,nbond0,nangle0} } );
        if( rigid ){
            return insertRigidMolecule( mol, pos, rot );
        }else{
            if(noH>0){ insertFlexibleMolecule_ignorH( mol, pos, rot, noH ); }
            else     { insertFlexibleMolecule       ( mol, pos, rot      ); }
            return -1;
        }
    }

    int insertMolecule( int itype, const Vec3d& pos, const Mat3d& rot, bool rigid ){
        return insertMolecule( molTypes[itype], pos, rot, rigid );
    };

    int insertMolecule( const std::string& molName, const Vec3d& pos, const Mat3d& rot, bool rigid ){
        //printf( "insertMolecule molName %s itype %i \n", molName.c_str(), molTypeDict[molName] );
        return insertMolecule( molTypes[ molTypeDict[molName] ], pos, rot, rigid );
    };

#endif // Molecule_h


#ifdef ForceField_h
    void toForceField( ForceField& ff ){
        if(iDebug>0) printf( " MMFFbuilder.toForceField na %li nb %li nA %li nd %li \n", atoms.size(), bonds.size(), angles.size(), dihedrals.size() );
        //mmff->deallocate();
        ff.realloc( atoms.size(), bonds.size(), angles.size(), dihedrals.size() );
        for(int i=0; i<atoms.size(); i++){
            ff.apos [i]  = atoms[i].pos;
            if(iDebug>0){ printf("[%i]", i); atoms[i].print(); if( atoms[i].iconf>=0){confs[atoms[i].iconf].print();} puts(""); }
        }
        for(int i=0; i<bonds.size(); i++){
            const Bond& b  = bonds[i];
            const Vec2i& ib    = b.atoms;
            ff.bond2atom[i]    = ib;
            //if(params){
            //    params->getBondParams( atoms[ib.x].type, atoms[ib.y].type, bonds[i].type, ff.bond_l0[i], ff.bond_k[i] );
            //}else{
            //    //printf( "no params \n" );
            //    ff.setBondParam(i, b.l0, b.k );
            //}
            ff.setBondParam(i, b.l0, b.k );
            if(iDebug>0){  printf( "bond[%i] (%i,%i) %g %g | %g %g\n", i, ff.bond2atom[i].i, ff.bond2atom[i].j, ff.bond_l0[i], ff.bond_k[i], b.l0, b.k ); }
            //bondTypes[i]       = bonds[i].type;
        }
        for(int i=0; i<angles.size(); i++){
            const Angle& a  = angles[i];
            ff.ang2bond[i] = a.bonds;
            ff.setAngleParam(i, a.a0, a.k );
            if(iDebug>0){  printf( "angle[%i] (%i,%i) (%g,%g) %g\n", i, ff.ang2bond[i].i, ff.ang2bond[i].j, ff.ang_cs0[i].x, ff.ang_cs0[i].y, ff.ang_k[i] ); }
        }
        for(int i=0; i<dihedrals.size(); i++){
            const Dihedral& d  = dihedrals[i];
            ff.tors2bond[i] = d.bonds;
            ff.setTorsParam( i, d.n, d.k );
            if(iDebug>0){ printf( "dihedrals[%i] (%i,%i,%i) %i %g\n", i, ff.tors2bond[i].a, ff.tors2bond[i].b, ff.tors2bond[i].c, ff.tors_n[i], ff.tors_k[i] ); }
        }
        ff.angles_bond2atom();
        ff.torsions_bond2atom();
        //exit(0);
    }
#endif // ForceField_h


void updatePBC( Vec3d* pbcShifts ){
    for(int i=0; i<bonds.size(); i++){
        pbcShifts[i] = pbcShift( bondPBC[i] );
    }
}



#ifdef MMFFmini_h
    void toMMFFmini( MMFFmini& ff, const MMFFparams* params ){
        ff.realloc( atoms.size(), bonds.size(), angles.size(), dihedrals.size() );
        export_apos     ( ff.apos );
        export_bonds    ( ff.bond2atom,   ff.bond_l0, ff.bond_k );
        export_angles   ( ff.ang2bond, 0, ff.ang_cs0, ff.ang_k  );
        export_dihedrals( ff.tors2bond,   ff.tors_n,  ff.tors_k );

        if( bPBC ){ ff.initPBC(); updatePBC( ff.pbcShifts ); }

        ff.angles_bond2atom();
        ff.torsions_bond2atom();
    }

    // ----  OLD version ----
    /*
    void toMMFFmini( MMFFmini& ff, const MMFFparams* params ){
        //printf( "na %i nb %i nA %i \n", atoms.size(), bonds.size(), dihedrals.size() );
        //mmff->deallocate();
        ff.realloc( atoms.size(), bonds.size(), angles.size(), dihedrals.size() );
        for(int i=0; i<atoms.size(); i++){
            ff.apos [i]  = atoms[i].pos;
            //println(atoms[i]);
            //if( atoms[i].iconf>=0 ) println(confs[atoms[i].iconf]);
        }
        for(int i=0; i<bonds.size(); i++){
            const Bond& b  = bonds[i];
            const Vec2i& ib    = b.atoms;
            ff.bond2atom[i]    = ib;
            if(params){
                params->getBondParams( atoms[ib.x].type, atoms[ib.y].type, bonds[i].type, ff.bond_l0[i], ff.bond_k[i] );
            }else{
                //printf( "no params \n" );
                ff.setBondParam(i, b.l0, b.k );
            }
            //printf( "bond[%i] (%i,%i) %g %g | %g %g\n", i, ff.bond2atom[i].i, ff.bond2atom[i].j, ff.bond_l0[i], ff.bond_k[i], b.l0, b.k );
            //bondTypes[i]       = bonds[i].type;
        }
        //printf( "toMMFFmini . Angles.size() %i \n", angles.size() );
        for(int i=0; i<angles.size(); i++){
            const Angle& a  = angles[i];
            ff.ang2bond[i] = a.bonds;
            ff.setAngleParam(i, a.a0, a.k );
            //printf( "angle[%i] (%i,%i) (%g,%g) %g\n", i, ff.ang2bond[i].i, ff.ang2bond[i].j, ff.ang_cs0[i].x, ff.ang_cs0[i].y, ff.ang_k[i] );
        }
        for(int i=0; i<dihedrals.size(); i++){
            const Dihedral& d  = dihedrals[i];
            ff.tors2bond[i] = d.bonds;
            ff.setTorsParam( i, d.n, d.k );
            //printf( "dihedrals[%i] (%i,%i,%i) %i %g\n", i, ff.tors2bond[i].a, ff.tors2bond[i].b, ff.tors2bond[i].c, ff.tors_n[i], ff.tors_k[i] );
        }
        ff.angles_bond2atom();
        ff.torsions_bond2atom();
        //exit(0);
    }
    */
#endif // MMFFmini_h

#ifdef MMFF_h
    void toMMFF( MMFF * mmff, MMFFparams* params ){
        //mmff->deallocate();
        mmff->allocate( atoms.size(), bonds.size(), angles.size(), 0 );
        //int * atomTypes = new int[atoms.size()];
        //int * bondTypes = new int[bonds.size()];
        for(int i=0; i<atoms.size(); i++){
            mmff->atypes[i] = atoms[i].type;
            mmff->atom2frag[i] = atoms[i].frag;
            mmff->apos [i]  = atoms[i].pos;
            mmff->aREQ [i]  = atoms[i].REQ;
            //atomTypes[i]  = atoms[i].type;
            //printf( "iatom %i atype %i ifrag %i pos (%g,%g,%g) REQ (%g,%g,%g) \n", i, atoms[i].type, atoms[i].frag, atoms[i].pos.x,atoms[i].pos.y,atoms[i].pos.z, atoms[i].REQ.x,atoms[i].REQ.y,atoms[i].REQ.z );
        }
        for(int i=0; i<bonds.size(); i++){
            mmff->bond2atom[i] = bonds[i].atoms;
            Vec2i ib           = bonds[i].atoms;
            params->getBondParams( atoms[ib.x].type, atoms[ib.y].type, bonds[i].type, mmff->bond_0[i], mmff->bond_k[i] );
            //bondTypes[i]       = bonds[i].type;
        }
        for(int i=0; i<angles.size(); i++){
            mmff->ang2bond[i] = angles[i].bonds;
            mmff->ang_0[i] = {1.0,0.0}; // TODO FIXME
            mmff->ang_k[i] = 0.5;       // TODO FIXME
        }
        if( frags.size()>0 ){
            mmff->allocFragment( frags.size() );
            for(int i=0; i<frags.size(); i++){
                MM::Fragment& fragi = frags[i];
                mmff->frag2a  [i] = fragi.atomRange.x;
                mmff->fragNa  [i] = fragi.atomRange.y;
                mmff->fapos0s [i] = fragi.mol->pos;
                double * posi= (mmff->poses + i*8);
                *(Vec3d *)(posi  )= fragi.pos;
                *(Quat4d*)(posi+4)= fragi.rot;
            }
        }
        //params.fillBondParams( mmff->nbonds, mmff->bond2atom, bondTypes, atomTypes, mmff->bond_0, mmff->bond_k );
        //delete [] atomTypes;
        //delete [] bondTypes;
    }
#endif // MMFF_h

#ifdef EFF_h
//#if defined(EFF_h) || defined(EFF_old_h)
    void toEFF( EFF& ff, const EFFAtomType* params, double esize, double dpair ){
        //int ne = bonds.size() * 2; // ToDo
        int ne = 0;
        int na = 0;
        DEBUG
        for(int i=0; i<atoms.size(); i++){
            int ityp = atoms[i].type;
            if( ityp==capAtomEpair.type || ityp==capAtomPi.type ) continue;
            na++;
            ne += params[ ityp ].ne;
            printf( "[%i] ityp %i ne %i  %i \n", i, ityp, params[ityp].ne, ne );
        }
        DEBUG
        printf( "na %i ne %i | %i \n", atoms.size(), ne, bonds.size()*2 );
        ff.realloc( na, ne );
        DEBUG
        for(int i=0; i<na; i++){
            int ityp = atoms[i].type;
            if( ityp==capAtomEpair.type || ityp==capAtomPi.type ) continue;
            //printf( "[%i] ityp %i \n" );
            ff.apos [i]  = atoms[i].pos;
            //ff.aQ [i]  = params[ ityp ].ne; // ToDo
            ff.aPars[i]  = EFF::default_AtomParams[ityp];
        }
        DEBUG
        for(int i=0; i<bonds.size(); i++){
            const MM::Bond& b  = bonds[i];
            const Vec2i& ib    = b.atoms;
            double c1=0.5-dpair;
            double c2=1-c1;
            int i2 = i*2;
            ff.epos [i2] = atoms[ib.a].pos*c1 + atoms[ib.b].pos*c2;
            ff.esize[i2] = esize;
            ff.espin[i2] = 1;
            i2++;
            ff.epos [i2] = atoms[ib.a].pos*c2 + atoms[ib.b].pos*c1;
            ff.esize[i2] = esize;
            ff.espin[i2] = -1;
            //break;
        }
        DEBUG
        // ToDo:  pi-bonds & e-pairs

    }

    /*
    int countValenceElectrons(){
        int ne=0;
        for(int i=0; i<atoms.size(); i++){

        }
    }
    */
#endif  // EFF_h


#ifdef EFF_old_h
    void toEFF_old( EFF& ff, const EFFAtomType* params, double esize, double dpair ){
        //int ne = bonds.size() * 2; // ToDo
        int ne = 0;
        int na = 0;
        DEBUG
        for(int i=0; i<atoms.size(); i++){
            int ityp = atoms[i].type;
            if( ityp==capAtomEpair.type || ityp==capAtomPi.type ) continue;
            na++;
            ne += params[ ityp ].ne;
            printf( "[%i] ityp %i ne %i  %i \n", i, ityp, params[ityp].ne, ne );
        }
        DEBUG
        printf( "na %i ne %i | %i \n", atoms.size(), ne, bonds.size()*2 );
        ff.realloc( na, ne );
        DEBUG
        for(int i=0; i<na; i++){
            int ityp = atoms[i].type;
            if( ityp==capAtomEpair.type || ityp==capAtomPi.type ) continue;
            //printf( "[%i] ityp %i \n" );
            ff.apos [i]  = atoms[i].pos;
            ff.aQ   [i]  = params[ ityp ].ne; // ToDo
        }
        DEBUG
        for(int i=0; i<bonds.size(); i++){
            const MM::Bond& b  = bonds[i];
            const Vec2i& ib    = b.atoms;
            double c1=0.5-dpair;
            double c2=1-c1;
            int i2 = i*2;
            ff.epos [i2] = atoms[ib.a].pos*c1 + atoms[ib.b].pos*c2;
            ff.esize[i2] = esize;
            ff.espin[i2] = 1;
            i2++;
            ff.epos [i2] = atoms[ib.a].pos*c2 + atoms[ib.b].pos*c1;
            ff.esize[i2] = esize;
            ff.espin[i2] = -1;
            //break;
        }
        DEBUG
        // ToDo:  pi-bonds & e-pairs

    }

    /*
    int countValenceElectrons(){
        int ne=0;
        for(int i=0; i<atoms.size(); i++){

        }
    }
    */
#endif  // EFF_old_h

}; // MMFFBuilder


} // namespace MMFF

#endif // MMFFBuilder_h
