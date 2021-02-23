
#ifndef Molecule_h
#define Molecule_h

#include <unordered_map>
#include <vector>
#include <string>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
//#include "IOutils.h"

#include "MMFFparams.h"

inline int atomName2int(char ch ){
    int i = -1;
    switch(ch){
        case 'H': i=1;  break;
        case 'B': i=5;  break;
        case 'C': i=6;  break;
        case 'N': i=7;  break;
        case 'O': i=8;  break;
        case 'F': i=9;  break;
        case 'P': i=15; break;
        case 'S': i=16; break;
    }
    return i;
}

void cpstr( const char* str, char* tmp, int i0, int n ){
    tmp[n]='\0';
    for(int i=0; i<n; i++){
        tmp[i]=str[i0+i];
    }
}

double getDouble(const char* str, int i0, int i1 ){
    const int n=i1-i0;
    char  tmp[n+1];
    cpstr( str, tmp, i0, n );
    //double f; //sscanf( "%lf", &f );
    return atof(tmp);
}

double getInt(const char* str, int i0, int i1 ){
    const int n=i1-i0;
    char  tmp[n+1];
    cpstr( str, tmp, i0, n );
    //printf("str: %s",   str );
    //printf("tmp: %s\n", tmp );
    //int i; sscanf( "%i", &i );
    return atoi(tmp);
}


inline int otherAtom( Vec2i ib, int ia ){ return (ia==ib.x)?ib.y:ib.x; }

#define NBMAX 4

class Molecule{ public:
    int  natoms=0, nbonds=0;
    Vec3d  * pos       = NULL;
    Vec2i  * bond2atom = NULL;
    //double * charges   = NULL;
    Vec3d  * REQs      = NULL;
    int    * atomType  = NULL;
    int    * bondType  = NULL;

    int   * atom_nb    = NULL;
    int   * atom2bond  = NULL;

    int nang = 0;
    Vec2i  * ang2bond   = NULL;
    double * ang0s      = NULL;
    //int * ang2atom = NULL;

    const MMFFparams* params = 0;
    // ToDo: this is redudant can be changed to params
    const std::vector       <std::string    >* atomTypeNames=0;
    const std::unordered_map<std::string,int>* atomTypeDict =0;

    void allocate( int natoms_, int nbonds_ ){
        natoms=natoms_; nbonds=nbonds_;
        _realloc( pos       , natoms );
        //_realloc( charges   , natoms );
        _realloc( REQs      , natoms );
        _realloc( atomType  , natoms );
        _realloc( bondType  , nbonds );
        _realloc( bond2atom , nbonds );
    }

    void bondsOfAtoms(){
        atom_nb   = new int[natoms];
        atom2bond = new int[natoms*NBMAX];
        for(int i=0; i<natoms; i++){ atom_nb[i]=0; }
        for(int i=0; i<nbonds; i++){
            Vec2i ib = bond2atom[i];
            int nb1 = atom_nb[ib.x];
            int nb2 = atom_nb[ib.y];
            if( (nb1>=4)||(nb2>=4) ){
                printf( "ERROR in Molecule.bondsOfAtoms() [ib %i] bond count(%i,%i) > NBMAX(%i) \n", i, nb1, nb2, NBMAX );
                exit(0);
            }
            atom2bond[ib.x*NBMAX+nb1] = i;
            atom2bond[ib.y*NBMAX+nb2] = i; // -i to denote directionality ?
            atom_nb[ib.x]=nb1+1;
            atom_nb[ib.y]=nb2+1;
        }
    }

    /*

    void sortAtomsOfBonds(){

    }

    void sortAtomNeighs(int ia){
        int nb = atom_nb[ia];
        int* ibs = atom2bond + i*NBMAX;
        int nghs[NBMAX];
        int bs  [NBMAX];
        for(int i=0; i<nb i++){

        }
        for(int i=0; i<nb i++){
            int ib = selectMinHigher(ja, nbconf, ibs );
        }
    }

    void remakeBondsSorted(){
        if( (atom_nb==0)||(atom2bond==0) )bondsOfAtoms();
        int ib=0;

        for(int i=0; i<natoms; i++){
            int nb = atom_nb[i];
            for(int ij=0; ij<nb; ij++){
                int j = atom2bond[i*NBMAX+ij];
                bond2atom[ib] = (Vec2i){i,j};
                ib++;
            }
        }
    }
    */

    void findBonds_brute(double Rmax){
        double R2max = sq(Rmax);
        std::vector<Vec2i> pairs;
        for(int i=0;i<natoms;i++){
            const Vec3d& pi = pos[i];
            for(int j=0;j<i;j++){
                Vec3d d; d.set_sub(pos[j],pi);
                double r2 = d.norm2();
                //printf( "%i,%i %g %g  (%g,%g,%g)  (%g,%g,%g) \n", i, j, r2, R2max, pi.x,pi.y,pi.z,    pos[j].x,pos[j].y,pos[j].z );
                if(r2<R2max){
                    pairs.push_back({i,j});
                }
            }
        }
        nbonds=pairs.size();
        _realloc(bond2atom,nbonds);
        for(int i=0;i<nbonds;i++){  bond2atom[i]=pairs[i]; }
    }

    void bindParams( const MMFFparams* params_ ){
        params = params_;
        atomTypeNames = &params->atomTypeNames;
        atomTypeDict  = &params->atomTypeDict;
    }

    int countAtomType( int ityp )const{
        int n=0;
        for(int i=0; i<natoms; i++){ if(atomType[i]==ityp)n++; }
        return n;
    }

    inline double getAng0( int ia ){
        int ityp = atomType[ia];
        if(params) ityp = params->atypes[ityp].iZ; // ToDo : later this can be improved
        int nb   = atom_nb [ia];
        printf( "getAng0 iT %i iZ %i nb %i \n", atomType[ia], ityp, nb );
        switch(ityp){
            case 1: return M_PI; break; // H
            case 6: // C
                if     (nb==2){ return M_PI;               }
                else if(nb==3){ return M_PI*0.66666666666; }
                else if(nb==4){ return M_PI*0.60833333333; }
                else          { return M_PI;               }
                break;
            case 7: // N
                if     (nb==2){ return M_PI*0.6666666666;  }
                else if(nb==3){ return M_PI*0.60833333333; }
                else          { return M_PI;               }
                break;
            case 5: return M_PI; break; // B
            case 8:  // O
            case 14: // Si
            case 16: // S
                return M_PI*0.60833333333;
            default: return M_PI;
        }
        return M_PI;
    }

    int autoAngles( bool bAng0=false ){
        // create angle between each pair of bonds of an atom
        nang = 0;
        if( (atom_nb==0)||(atom2bond==0) )bondsOfAtoms();
        for(int i=0; i<natoms; i++){ int nb = atom_nb[i]; nang+=(nb*(nb-1))>>1; }
        ang2bond       = new Vec2i [nang];
        if(bAng0)ang0s = new double[nang];
        int * a2b = atom2bond;
        int iang = 0;
        for(int ia=0; ia<natoms; ia++){
            int nb = atom_nb[ia];
            double a0;
            if((bAng0)&&(nb>1)){
                a0 = getAng0( ia );
                printf( "atom[%i] a0 %g \n ", ia, a0 );
            }
            for(int i=0; i<nb; i++){
                int ib1 = a2b[i];
                for(int j=0; j<i; j++){
                    ang2bond[iang]={ib1,a2b[j]};
                    if(bAng0){
                        ang0s[iang] = a0;
                        printf( "autoAngles[%i][ia%i|%i,%i] %g[rad] | nZ %i nb %i \n", iang, ia,i,j, ang0s[iang], atomType[ia], nb );
                    }
                    iang++;
                }
            }
            a2b+= NBMAX;
        }

        return iang;
    }

    Vec3d getCOG_av()const{
        Vec3d cog = (Vec3d){0.0,0.0,0.0};
        for(int i=0; i<natoms; i++){
            cog.add(pos[i]);
        }
        cog.mul(1.0d/natoms);
        return cog;
    }

    void addToPos( Vec3d dp ){ for(int i=0; i<natoms; i++){ pos[i].add(dp); } }

    /*
    // DEPRECATED
    int loadMol_old( const char* fname ){
        // xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
        // http://www.daylight.com/meetings/mug05/Kappler/ctfile.pdf
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        char buff[1024];
        char * line;
        int nl;
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        sscanf( line, "%i %i\n", &natoms, &nbonds );
        //printf("%i %i \n", natoms, nbonds );

        allocate(natoms,nbonds);
        for(int i=0; i<natoms; i++){
            //char ch;
            char at_name[8];
            double junk;
            line = fgets( buff, 1024, pFile );  //printf("%s",line);
            sscanf( line, "%lf %lf %lf %s %lf %lf\n", &pos[i].x, &pos[i].y, &pos[i].z,  at_name, &junk, &REQs[i].z );
            //printf(       "%lf %lf %lf %s %lf %lf\n",  pos[i].x,  pos[i].y,  pos[i].z,  at_name,  junk,  REQs[i].z );
            // atomType[i] = atomChar2int( ch );
            auto it = atomTypeDict->find( at_name );
            if( it != atomTypeDict->end() ){
                atomType[i] = it->second;
            }else{
                //atomType[i] = atomChar2int( at_name[0] );
                atomType[i] = -1;
            }
        }
        for(int i=0; i<nbonds; i++){
            line = fgets( buff, 1024, pFile );  //printf("%s",line);
            sscanf( line, "%i %i %i\n", &bond2atom[i].x, &bond2atom[i].y, &bondType[i] );
            printf(       "%i %i %i\n",  bond2atom[i].x,  bond2atom[i].y,  bondType[i] );
            bond2atom[i].x--;
            bond2atom[i].y--;
        }
        return natoms + nbonds;
    }
    */

    void assignAtomType(int i, const char * at_name ){
        if(atomTypeDict){
            auto it = atomTypeDict->find( at_name );
            if( it != atomTypeDict->end() ){
                atomType[i] = it->second;
                if(1==it->second)REQs[i].x=1.0; // Hydrogen is smaller
            }else{
                //atomType[i] = atomChar2int( at_name[0] );
                atomType[i] = -1;
            }
            //printf( "at_name[%i] %s -> %i \n", i, at_name, atomType[i] );
        }else{
            atomType[i] = atomName2int( at_name[0] );
        }
    }


    int loadMol( const char* fname ){
        // 0        10         20
        //   -13.0110  -15.2500   -0.0030 N   0  0  0  0  0  0  0  0  0  0  0  0
        // xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
        // http://www.daylight.com/meetings/mug05/Kappler/ctfile.pdf
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        char buff[1024];
        char * line;
        int nl;
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        //sscanf( line, "%i %i\n", &natoms, &nbonds );
        //printf("natoms, nbonds %i %i \n", natoms, nbonds );
        //line = fgets( buff, 1024, pFile );
        natoms = getInt( line, 0, 3 );
        nbonds = getInt( line, 3, 6 );
        printf("loadMol(%s):  natoms, nbonds %i %i \n", fname, natoms, nbonds );
        //exit(0);
        //return 0;
        allocate(natoms,nbonds);
        int istr_x      =0;
        int istr_y      =10;
        int istr_z      =20;
        int istr_name   =25;
        int istr_charge =30;

        for(int i=0; i<natoms; i++){
            //char ch;
            double junk;
            char at_name[8];
            line = fgets( buff, 1024, pFile );  //printf("%s",line);
            sscanf( line, "%lf %lf %lf %s %lf %lf\n", &pos[i].x, &pos[i].y, &pos[i].z,  at_name, &junk, &(REQs[i].z) );
            printf(       "%lf %lf %lf %s %lf %lf\n",  pos[i].x,  pos[i].y,  pos[i].z,  at_name,  junk,   REQs[i].z   );
            REQs[i].x = 1.5; // [A]  van der Waals radius default
            REQs[i].y = 0;   // [eV] van der Waals binding energy default
            assignAtomType( i, at_name );
        }
        for(int i=0; i<nbonds; i++){
            line = fgets( buff, 1024, pFile );  //printf("%s",line);
            bond2atom[i].x = getInt( line, 0, 3 );
            bond2atom[i].y = getInt( line, 3, 6 );
            bondType[i]    = getInt( line, 6, 9 );
            //sscanf( line, "%i %i %i\n", &bond2atom[i].x, &bond2atom[i].y, &bondType[i] );
            //printf(       "%i %i %i\n",  bond2atom[i].x,  bond2atom[i].y,  bondType[i] );
            bond2atom[i].x--;
            bond2atom[i].y--;
        }
        return natoms + nbonds;
    }

    int loadXYZ(const char* fname ){
        // xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
        // http://www.daylight.com/meetings/mug05/Kappler/ctfile.pdf
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        char buff[1024];
        char * line;
        int nl;
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        sscanf( line, "%i \n", &natoms );
        //printf("natoms %i \n", natoms );
        //allocate(natoms,0);
        allocate(natoms,nbonds);
        line = fgets( buff, 1024, pFile ); // comment
        for(int i=0; i<natoms; i++){
            //char ch;
            char at_name[8];
            double junk;
            line = fgets( buff, 1024, pFile );  //printf("%s",line);
            double Q;
            int nret = sscanf( line, "%s %lf %lf %lf %lf\n", at_name, &pos[i].x, &pos[i].y, &pos[i].z, &Q );
            if( nret >= 5 ){  REQs[i].z=Q; }else{ REQs[i].z=0; };
            printf(       "mol[%i] %s %lf %lf %lf  %lf    *atomTypeDict %i\n", i,  at_name, pos[i].x,  pos[i].y,  pos[i].z,   REQs[i].z, atomTypeDict );
            // atomType[i] = atomChar2int( ch );
            auto it = atomTypeDict->find( at_name );
            if( it != atomTypeDict->end() ){
                atomType[i] = it->second;
            }else{
                //atomType[i] = atomChar2int( at_name[0] );
                atomType[i] = -1;
            }
            //printf( " i %i name %s ityp %i \n", i, at_name, atomType[i] );
        }

        //printf( "atypNames.size() %i \n", atypNames->size() );
        //for ( auto element : *atypNames ){
	    //    printf(  "atypNames[%s]-> %i \n", element.first.c_str(), element.second );
        //}
        //printf("loadXYZ DONE \n");
        return natoms;
    }

    int loadXYZ_bas(const char* fname ){
        // xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
        // http://www.daylight.com/meetings/mug05/Kappler/ctfile.pdf
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        char buff[1024];
        char * line;
        int nl;
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        sscanf( line, "%i \n", &natoms );
        //printf("natoms %i \n", natoms );
        //allocate(natoms,0);
        nbonds=0;
        allocate(natoms,nbonds);
        //line = fgets( buff, 1024, pFile ); // comment
        for(int i=0; i<natoms; i++){
            //char ch;
            //char at_name[8];
            //double junk;
            line = fgets( buff, 1024, pFile );  //printf("%s",line);
            double Q;
            int nret = sscanf( line, "%i %lf %lf %lf %lf\n", &atomType[i], &pos[i].x, &pos[i].y, &pos[i].z, &Q );
            if( nret >= 5 ){  REQs[i].z=Q; }else{ REQs[i].z=0; };
            printf(       "mol[%i] %i %lf %lf %lf  %lf    n", i,  atomType[i], pos[i].x,  pos[i].y,  pos[i].z,   REQs[i].z  );
            // atomType[i] = atomChar2int( ch );
        }
        return natoms;
    }

    void printAtom2Bond()const{
        if(atom2bond==0){ printf("# Molecule.printAtom2Bond() atom2bond == NULL ; return \n"); return; };
        printf("# Molecule.printAtom2Bond()\n");
        int * a2b = atom2bond;
        for(int ia=0; ia<natoms; ia++){
            int nb = atom_nb[ia];
            printf("%i :", ia);
            for(int i=0; i<nb; i++){
                //printf(" %i", a2b[i] );
                printf(" %i", otherAtom( bond2atom[a2b[i]], ia) );
            }
            printf("\n");
            a2b+= NBMAX;
        }
    }

    void printAtomInfo()const{
        printf(" # Molecule.printAtomInfo() : \n" );
        for(int i=0; i<natoms; i++){
            printf( "atom[%i] t %i pos(%g,%g,%g) REQs(%g,%g,%g) \n", i, atomType[i], pos[i].x,pos[i].y,pos[i].z, REQs[i].x, REQs[i].y, REQs[i].z );
        }
    }

    void printAngleInfo()const{
        if(ang2bond==0){ printf("# Molecule.printAngleInfo() ang2bond == NULL ; return \n"); return; };
        printf(" # Molecule.printAngleInfo()\n");
        for(int i=0; i<nang; i++){
            double a0=0; if(ang0s)a0=ang0s[i];
            printf( "angle[%i|%i,%i] a0 %g[rad] \n", i, ang2bond[i].a, ang2bond[i].b, a0 );
        }
    }

    void dealloc(){
        _dealloc( pos       );
        _dealloc( bond2atom );
        _dealloc( REQs      );
        _dealloc( atomType  );
        _dealloc( bondType  );
        _dealloc( atom_nb   );
        _dealloc( atom2bond );
        _dealloc( ang2bond  );
    }

};

#endif






