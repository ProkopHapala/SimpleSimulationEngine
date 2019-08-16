#ifndef MMFFBuilder_h
#define MMFFBuilder_h

#include <string>
#include <vector>
#include <unordered_map>

#include "Molecule.h"
#include "MMFF.h"
#include "MMFFmini.h"
#include "MMFFparams.h"

//namespace MMFF{

class MMFFAtom{ public:
    int type;
    int frag;  //=-1;
    int iconf; //
    Vec3d pos;
    Vec3d REQ;

    //void copyTo(MMFFAtom& c){ c.type=type; }
    //MMFFAtom(){};
    //MMFFAtom(int type_,int ifrag_,Vec3d pos_,Vec3d REQ_): type(type_),ifrag(ifrag_),pos(pos_),REQ(REQ_){};

};

#define N_NEIGH_MAX 8
enum class NeighType: int {
    pi    = -1,
    epair = -2,
    H     = -3
};


class MMFFAtomConf{ public:
    //int npi;       // number of pi bonds
    //int nepair;    // number of electron pairs
    //int nH;        // number of capping hydrogens
    int iatom=-1;
    int n     =0;
    int nbond =0;
    int npi   =0; // pi bonds
    int ne    =0; // electron pairs
    int nH    =0; //
    int neighs[N_NEIGH_MAX]; // neighs  - NOTE: bonds not atoms !!!!

    MMFFAtomConf() = default;
    MMFFAtomConf(const MMFFAtomConf&) = default;
    MMFFAtomConf(std::initializer_list<MMFFAtomConf>) {};

    inline bool addNeigh(int ia, int& ninc ){
        if(n>=N_NEIGH_MAX)return false;
        neighs[n]=ia;
        ninc++;
        n++;
        return true;
    };

    //inline bool addBond( Vec2i b ){
    //    if     (b.i==iatom){ addNeigh(b.j,nbond); }
    //    else if(b.j==iatom){ addNeigh(b.i,nbond); }
    //    return false;
    //}
    inline bool addBond (int i){ return addNeigh(i,nbond); };
    inline bool addH    (     ){ return addNeigh((int)NeighType::H    ,nH ); };
    inline bool addPi   (     ){ return addNeigh((int)NeighType::pi   ,npi); };
    inline bool addEpair(     ){ return addNeigh((int)NeighType::epair,ne ); };
    inline int  clearNonBond(){ n=nbond; npi=0;ne=0;nH=0; };

};

class MMFFBond{ public:
    int    type;
    Vec2i  atoms;
    double l0,k;
    inline int getNeighborAtom(int ia){ if(ia==atoms.i){return atoms.j;}else if(ia==atoms.j){return atoms.i;}else{ return -1; } }
};

class MMFFAngle{ public:
    int type;
    Vec2i  bonds;
    double a0,k;
};

class MMFFmol{ public:
    Molecule * mol;
    Vec3i    i0;
};

class MMFFfrag{ public:
    //int imolType;
    int atom0, natom;
    Vec3d  pos;
    Quat4d rot;
    Molecule * mol;
    Vec3d    * pos0s;
};

class MMFFBuilder{  public:

    MMFFparams*  params = NULL;
    std::unordered_map<std::string,int> molTypeDict;
    std::vector<Molecule*> molTypes;

    std::vector<MMFFAtom>  atoms;
    std::vector<MMFFBond>  bonds;
    std::vector<MMFFAngle> angles;
    std::vector<MMFFmol>   mols;
    std::vector<MMFFfrag>  frags;

    std::vector<MMFFAtomConf>  confs;
    //std::vector<int>  atom_neighs;

    std::unordered_map<size_t,size_t> fragTypes;
    std::unordered_map<size_t,size_t> mol2molType;

    MMFFBond  defaultBond;
    MMFFAngle defaultAngle;

    MMFFAtom capAtom; // = (MMFFAtom){,,};
    MMFFAtom capAtomEpair;
    MMFFAtom capAtomPi;
    MMFFBond capBond;
    Vec3d    capUp   = (Vec3d){0.0d,0.0d,1.0d};
    bool bDummyPi    = false;
    bool bDummyEpair = false;

    void clearMolTypes( bool deep ){
        if(deep){ for(Molecule* mol : molTypes ){ mol->dealloc(); delete mol; } }
        molTypeDict.clear();
        molTypes.clear();
    }

    void clear(){
        //printf( "DEBUG MMFFBuilder.clear \n");
        atoms.clear(); //printf("DEBUG a.1 \n");
        bonds.clear(); //printf("DEBUG a.2 \n");
        angles.clear();
        mols .clear(); //printf("DEBUG a.3 \n");
        frags.clear(); //printf("DEBUG a.4 \n");
        fragTypes.clear();
        //printf( "DEBUG MMFFBuilder.clear DONE \n");
    }

    int loadMolType(const char* fname ){
        Molecule* mol = new Molecule();      //printf( "DEBUG 1.1.1 \n" );
        mol->atypNames = &params->atypNames; //printf( "DEBUG 1.1.2 \n" );
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

    int loadMolType(const std::string& fname, const std::string& label ){
        //printf( "fname:`%s` label:`%s` \n", fname.c_str(), label.c_str()  );
        int itype = loadMolType( fname.c_str() );
        //printf( "fname:`%s` label:`%s` itype %i \n", fname.c_str(), label.c_str(), itype  );
        molTypeDict[label] = itype;
        return itype;
    };

    int insertMolecule( Molecule * mol, const Vec3d& pos, const Mat3d& rot, bool rigid ){
        int natom0  = atoms.size();
        int nbond0  = bonds.size();
        int nangle0 = angles.size();
        mols.push_back( (MMFFmol){mol, (Vec3i){natom0,nbond0,nangle0} } );

        int natoms0 = atoms.size();
        if( rigid ){
            Quat4d qrot; qrot.fromMatrix(rot);
            int ifrag = frags.size();
            printf( "insertMolecule mol->natoms %i \n", mol->natoms );
            for(int i=0; i<mol->natoms; i++){
                //Vec3d REQi = (Vec3d){1.0,0.03,mol->}; // TO DO : LJq can be set by type
                //atoms.push_back( (MMFFAtom){mol->atomType[i],mol->pos[i], LJq } );
                Vec3d  REQi = mol->REQs[i];   REQi.y = sqrt(REQi.y); // REQi.z = 0.0;
                Vec3d  p; rot.dot_to(mol->pos[i],p); p.add( pos );
                atoms.push_back( (MMFFAtom){mol->atomType[i], ifrag, -1, p, REQi } );
            }
            frags.push_back( (MMFFfrag){natoms0, atoms.size()-natoms0, pos, qrot, mol}  );
            //size_t mol_id = static_cast<size_t>(mol);
            size_t mol_id = (size_t)(mol);
            auto got = fragTypes.find(mol_id);
            if ( got == fragTypes.end() ) {
                fragTypes[ mol_id ] = frags.size()-1; // WTF ?
            }else{}
            return ifrag;
        }else{
            for(int i=0; i<mol->natoms; i++){
                //Vec3d LJq = (Vec3d){0.0,0.0,0.0}; // TO DO : LJq can be set by type
                //Vec3d LJq = (Vec3d){1.0,0.03,0.0}; // TO DO : LJq can be set by type
                Vec3d  REQi = mol->REQs[i];   REQi.y = sqrt(REQi.y);
                Vec3d p; rot.dot_to(mol->pos[i],p); p.add( pos );
                atoms.push_back( (MMFFAtom){mol->atomType[i], -1, -1, p, REQi } );
            }
            for(int i=0; i<mol->nbonds; i++){
                bonds.push_back( (MMFFBond){mol->bondType[i], mol->bond2atom[i] + ((Vec2i){natom0,natom0}), defaultBond.l0, defaultBond.k } );
            }
            for(int i=0; i<mol->nang; i++){
                angles.push_back( (MMFFAngle){ 1, mol->ang2bond[i] + ((Vec2i){nbond0,nbond0}), defaultAngle.a0, defaultAngle.k } );
            }
            return -1;
        }
    }

    int insertMolecule( int itype, const Vec3d& pos, const Mat3d& rot, bool rigid ){
        return insertMolecule( molTypes[itype], pos, rot, rigid );
    };

    int insertMolecule( const std::string& molName, const Vec3d& pos, const Mat3d& rot, bool rigid ){
        printf( "insertMolecule molName %s itype %i \n", molName.c_str(), molTypeDict[molName] );
        return insertMolecule( molTypes[ molTypeDict[molName] ], pos, rot, rigid );
    };

    void assignAtomTypes(){
        for(int i=0; i<atoms.size(); i++){
            //mmff->aLJq [i]  = atoms[i].type;
            int ityp = atoms[i].type;
            atoms[i].REQ.x = params->atypes[ityp].RvdW;
            atoms[i].REQ.y = params->atypes[ityp].EvdW;
            atoms[i].REQ.z = 0;
            //atomTypes[i]  = atoms[i].type;
        }
    }

    void toMMFF( MMFF * mmff ){
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
                MMFFfrag& fragi = frags[i];
                mmff->frag2a  [i] = fragi.atom0;
                mmff->fragNa  [i] = fragi.natom;
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

    void toMMFFmini( MMFFmini * mmff ){
        //mmff->deallocate();
        mmff->realloc( atoms.size(), bonds.size(), angles.size(), 0 );
        for(int i=0; i<atoms.size(); i++){
            mmff->apos [i]  = atoms[i].pos;
        }
        for(int i=0; i<bonds.size(); i++){
            mmff->bond2atom[i] = bonds[i].atoms;
            Vec2i ib           = bonds[i].atoms;
            params->getBondParams( atoms[ib.x].type, atoms[ib.y].type, bonds[i].type, mmff->bond_l0[i], mmff->bond_k[i] );
            //bondTypes[i]       = bonds[i].type;
        }
        for(int i=0; i<angles.size(); i++){
            mmff->ang2bond[i] = angles[i].bonds;
            mmff->ang_cs0[i] = {1.0,0.0}; // TODO FIXME
            mmff->ang_k[i] = 0.5;       // TODO FIXME
        }
    }

    // ============= Add Capping Hydrogens

    bool getAtomConf(int ia,MMFFAtomConf* conf){ int ic=atoms[ia].iconf; if(ic>=0){ conf=&confs[ic];return true;}return false; };

    MMFFAtomConf* insertAtom(const MMFFAtom& atom, bool bConf ){
        atoms.push_back(atom);
        if(bConf){
            int ic = confs.size();
            int ia = atoms.size()-1;
            printf( "insertAtom ia %i ic %i \n", ia, ic );
            atoms.back().iconf = ic;
            confs.push_back(MMFFAtomConf());
            MMFFAtomConf& c = confs.back();
            c.iatom = ia;
            return &c;
        }
        return 0;
    }

    void insertBond(const MMFFBond& bond ){
        int ib = bonds.size();
        bonds.push_back(bond);
        int ic = atoms[bond.atoms.i].iconf;
        int jc = atoms[bond.atoms.j].iconf;
        //if(ic>=0){ confs[ic].addBond(bond.atoms.j); }
        //if(jc>=0){ confs[jc].addBond(bond.atoms.i); }
        printf( "insertBond %i(%i,%i) to %i,%i\n", ib, bond.atoms.i,bond.atoms.j, ic, jc );
        if(ic>=0){ confs[ic].addBond(ib); }
        if(jc>=0){ confs[jc].addBond(ib); }
    }

    void addCap(int ia,Vec3d& hdir, MMFFAtom* atomj, int btype){
        int ja=atoms.size();
        //capAtom;
        MMFFAtom atom_tmp;
        if(atomj==0){
            atom_tmp=capAtom;
            atomj=&atom_tmp;
        }
        if(btype<0) btype=capBond.type;
        atomj->pos = atoms[ia].pos + hdir;
        atoms.push_back( *atomj );
        bonds.push_back( (MMFFBond){btype,{ia,ja}} );
        int ic = atoms[ia].iconf;
        //confs[ic].addBond(ja);
        confs[ic].addBond(bonds.size());
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
        }
    }

    void makeSPConf(int ia,int npi,int ne){
        //if(nH==0){ // ToDo : check reasonable limits npi, nh
        int ic = atoms[ia].iconf;
        MMFFAtomConf& conf = confs[ic];
        conf.clearNonBond();
        int nb = conf.nbond;
        int n  = 4-nb-npi;   // number
        int nH = n-ne;
        printf( "ia %i nb %i npi %i ne %i n %i nH %i ne %i \n", ia, nb,npi,ne,n,nH,ne );
        //Mat3d m;
        Vec3d hs[4];
        for(int i=0;i<nb;i++){
            int ib = conf.neighs[i];
            int ja = bonds[ib].getNeighborAtom(ia);
            hs[i]  = atoms[ja].pos - atoms[ia].pos;
            hs[i].normalize();
        }
        makeConfGeom(conf.nbond, npi, hs);
        bool Hmask[]{1,1,1};
        if(nH!=n) Hmask[rand()%n]=0;
        bool breverse = (nH==2)&&(n==3);
        for(int i=0; i<n; i++){
            if     (Hmask[i])   { addCap(ia,hs[i+nb],&capAtom       ,0); }
            else if(bDummyEpair){ addCap(ia,hs[i+nb],&capAtomEpair  ,0); }
        }
        if(bDummyPi){
            for(int i=0; i<npi; i++){ addCap(ia,hs[i+n+nb],&capAtomPi,0); }
        }
    }

    void makeSPConf(int ia){
        MMFFAtomConf conf;
        getAtomConf(ia,&conf);
        makeSPConf(ia,conf.npi,conf.ne);
    }

    // ============= Angles

    void addAnglesToBond( int ib, int n, int* neighs, double a0, double k ){
        for(int j=0; j<n; j++){
            angles.push_back( (MMFFAngle){defaultAngle.type,(Vec2i){ ib, neighs[j]}, a0,k} );
        }
    }
    void addAnglesUpToN( int n, int* neighs, double a0, double k ){
        for(int i=0; i<n; i++){ addAnglesToBond( i, i, neighs, a0, k ); }
    }
    void addAnglesToAtom( int ia, double ksigma, double kpi ){
        MMFFAtomConf& conf = confs[ atoms[ia].iconf ];
        int nsigma = conf.nbond;
        if(bDummyPi && (conf.npi>0)){
            nsigma -= conf.npi;
            for(int i=0;i<conf.npi;i++){ addAnglesToBond( i+nsigma, i+nsigma, conf.neighs, M_PI_2, kpi ); }
        }
        constexpr double a0s[]{ 0.0, M_PI, 120*M_PI/180, 109.5*M_PI/180 };
        double a0=a0s[nsigma];
        addAnglesUpToN( nsigma, conf.neighs, a0, ksigma );
    }

    //void flipCap(){}
};

int write2xyz( FILE* pfile, MMFF * mmff, MMFFparams * params ){
    fprintf(pfile, "%i\n", mmff->natoms );
    fprintf(pfile, "#comment \n");
    for(int i=0; i<mmff->natoms; i++){
        int ityp   = mmff->atypes[i];
        Vec3d&  pi = mmff->apos[i];
        //printf( "write2xyz %i %i (%g,%g,%g) %s \n", i, ityp, pi.x,pi.y,pi.z, params->atypes[ityp].name );
        fprintf( pfile, "%s   %15.10f   %15.10f   %15.10f \n", params->atypes[ityp].name, pi.x,pi.y,pi.z );
    };
    return mmff->natoms;
}

int save2xyz( char * fname, MMFF * mmff, MMFFparams * params ){
    FILE* pfile = fopen(fname, "w");
    if( pfile == NULL ) return -1;
    int n = write2xyz(pfile, mmff, params );
    fclose(pfile);
    return n;
}

//} namespace MMFF

#endif
