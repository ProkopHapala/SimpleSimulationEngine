
#ifndef Molecule_h
#define Molecule_h

#include <unordered_map>
#include <vector>
#include <string>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

inline int atomChar2int(char ch ){
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
    Vec2i * ang2bond   = NULL;
    //int * ang2atom = NULL;

    std::unordered_map<std::string,int>* atypNames = NULL;

    void allocate( int natoms_, int nbonds_ ){
        natoms=natoms_; nbonds=nbonds_;

        if(pos      ==NULL) pos       = new Vec3d [natoms];
        //if(charges  ==NULL) charges   = new double[natoms];
        if(REQs     ==NULL) REQs      = new Vec3d [natoms];
        if(bond2atom==NULL) bond2atom = new Vec2i [nbonds];
        if(atomType ==NULL) atomType  = new int   [natoms];
        if(bondType ==NULL) bondType  = new int   [nbonds];

        /*
        if(pos      )delete [] pos;       pos       = new Vec3d [natoms];
        //if(charges  ==NULL) charges   = new double[natoms];
        if(REQs     )delete [] REQs;      REQs      = new Vec3d [natoms];
        if(bond2atom)delete [] bond2atom; bond2atom = new Vec2i [nbonds];
        if(atomType )delete [] atomType;  atomType  = new int   [natoms];
        if(bondType )delete [] bondType;  bondType  = new int   [nbonds];
        */
    }

    int bondsOfAtoms(){
        atom_nb   = new int[natoms];
        atom2bond = new int[natoms*NBMAX];
        for(int i=0; i<natoms; i++){ atom_nb[i]=0; }
        for(int i=0; i<nbonds; i++){
            Vec2i ib = bond2atom[i];
            int nb1 = atom_nb[ib.x];
            int nb2 = atom_nb[ib.y];
            atom2bond[ib.x*NBMAX+nb1] = i;
            atom2bond[ib.y*NBMAX+nb2] = i; // -i to denote directionality ?
            atom_nb[ib.x]=nb1+1;
            atom_nb[ib.y]=nb2+1;
        }
    }

    void printAtom2Bond(){
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

    void printAtomInfo(){
        printf("Molecule::printAtomInfo : \n" );
        for(int i=0; i<natoms; i++){
            printf( "%i %f %f %f \n", i, REQs[i].x, REQs[i].y, REQs[i].z );
        }
    }

    int autoAngles(){
        // create angle between each pair of bonds of an atom
        nang = 0;
        for(int i=0; i<natoms; i++){ int nb = atom_nb[i]; nang+=(nb*(nb-1))>>1; }
        ang2bond = new Vec2i[nang];
        int * a2b = atom2bond;
        int iang = 0;
        for(int ia=0; ia<natoms; ia++){
            int nb = atom_nb[ia];
            for(int i=0; i<nb; i++){
                int ib1 = a2b[i];
                for(int j=0; j<i; j++){
                    ang2bond[iang]={ib1,a2b[j]};
                    iang++;
                }
            }
            a2b+= NBMAX;
        }

        return iang;
    }

    Vec3d getCOG_av(){
        Vec3d cog = (Vec3d){0.0,0.0,0.0};
        for(int i=0; i<natoms; i++){
            cog.add(pos[i]);
        }
        cog.mul(1.0d/natoms);
        return cog;
    }

    void addToPos( Vec3d dp ){ for(int i=0; i<natoms; i++){ pos[i].add(dp); } }

    int loadMol( char* fname ){
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
        printf("%i %i \n", natoms, nbonds );
        allocate(natoms,nbonds);
        for(int i=0; i<natoms; i++){
            //char ch;
            char at_name[8];
            double junk;
            line = fgets( buff, 1024, pFile );  //printf("%s",line);
            sscanf( line, "%lf %lf %lf %s %lf %lf\n", &pos[i].x, &pos[i].y, &pos[i].z,  at_name, &junk, &REQs[i].z );
            printf(       "%lf %lf %lf %s %lf %lf\n",  pos[i].x,  pos[i].y,  pos[i].z,  at_name,  junk,  REQs[i].z );
            // atomType[i] = atomChar2int( ch );
            auto it = atypNames->find( at_name );
            if( it != atypNames->end() ){
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


    int loadMol_const( char* fname ){
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
        sscanf( line, "%i %i\n", &natoms, &nbonds );
        printf("%i %i \n", natoms, nbonds );
        allocate(natoms,nbonds);

        int istr_x      =0;
        int istr_y      =10;
        int istr_z      =20;
        int istr_name   =25;
        int istr_charge =30;

        for(int i=0; i<natoms; i++){
            //char ch;
            char at_name[8];
            line = fgets( buff, 1024, pFile );  //printf("%s",line);
            sscanf( line, "%lf %lf %lf %s\n", &pos[i].x, &pos[i].y, &pos[i].z,  at_name );
            printf(       "%lf %lf %lf %s\n",  pos[i].x,  pos[i].y,  pos[i].z,  at_name );
            // atomType[i] = atomChar2int( ch );
            auto it = atypNames->find( at_name );
            if( it != atypNames->end() ){
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


    int loadXYZ( char* fname ){
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
        printf("%i \n", natoms );
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
            if( nret == 5 ){  REQs[i].z=Q; };
            printf(       "%s %lf %lf %lf  %lf\n",  at_name, pos[i].x,  pos[i].y,  pos[i].z,   REQs[i].z );
            // atomType[i] = atomChar2int( ch );
            auto it = atypNames->find( at_name );
            if( it != atypNames->end() ){
                atomType[i] = it->second;
            }else{
                //atomType[i] = atomChar2int( at_name[0] );
                atomType[i] = -1;
            }
        }
        return natoms;
    }


};

#endif






