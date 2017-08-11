#ifndef MMFFparams_h
#define MMFFparams_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

#include <string>
#include <unordered_map>

inline uint64_t getBondTypeId( uint16_t at1, uint16_t at2, uint8_t order ){
    if (at1>at2){ SWAP(at1,at2,uint16_t); }
    return pack64( at1, at2, order, 0 );
}

class BondType{ public:
    double length;
    double stiffness;
    uint16_t at1,at2;
    uint8_t  order;
};

class AtomType{ public:
    char      name[8];
    uint8_t   iZ;         // proton number
    uint8_t   neval;      // number of valence electrons
    uint8_t   valence;    // sum of bond orders of all bonds
    uint8_t   sym;        // sp, sp2, sp3 ... tetrahedral, triangle, linear, kink, square, octahedral
    uint32_t  color;
    double    RvdW;
    double    EvdW;

    char* toString( char * str ){
        sprintf( str, "printf: %s %i %i %i %i %lf %lf %x", name,  iZ,   neval,  valence,   sym,    RvdW, EvdW,   color );
        return str;
    }

    void fromString( char * str ){
        int iZ_, neval_, valence_, sym_;
        //char sclr[6];
        sscanf( str, " %s %i %i %i %i %lf %lf %x \n", name, &iZ_, &neval_, &valence_, &sym_,  &RvdW, &EvdW, &color );
        iZ=iZ_; neval=neval_; valence=valence_; sym=sym_;
        //printf( "printf: %s %i %i %i %i %lf %lf %x \n", name,  iZ,   neval_,  valence,   sym,    RvdW, EvdW,   color );
        //char ss[256]; printf("%s\n", toString(ss) );
    }

};

// ======================
// ====   MMFFparams
// ======================

class MMFFparams{ public:

    // http://www.science.uwaterloo.ca/~cchieh/cact/c120/bondel.html

    std::unordered_map<std::string,int>   atypNames;
    std::vector<AtomType>                    atypes;
    std::unordered_map<uint64_t,BondType>     bonds;

    double default_bond_length      = 2.0;
    double default_bond_stiffness   = 1.0;

    int loadAtomTypes(char * fname){

        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        char buff[1024];
        char * line;
        int nl;

        AtomType atyp;
        int i;
        for(int i; i<0xFFFF; i++){
            line = fgets( buff, 1024, pFile );
            if(line==NULL) break;
            atyp.fromString( line );
            atypes.push_back(atyp);
            atypNames[atyp.name] = atypes.size()-1;
        }
        return i;
    }

    int loadBondTypes(char * fname){

        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        char buff[1024];
        char * line;
        BondType bt;
        //line = fgets( buff, 1024, pFile ); //printf("%s",line);
        //sscanf( line, "%i %i\n", &n );
        int i;
        for(int i; i<0xFFFF; i++){
            line = fgets( buff, 1024, pFile );
            if(line==NULL) break;
            //printf("%s",line);
            sscanf(  line, "%i %i %i %lf %lf\n", &bt.at1, &bt.at2, &bt.order, &bt.length, &bt.stiffness );
            printf(        "%i %i %i %lf %lf\n",  bt.at1,  bt.at2,  bt.order,  bt.length,  bt.stiffness );
            uint64_t id = getBondTypeId( bt.at1, bt.at2, bt.order );
            //printf( ":: (%i,%i,%i) -> %i \n", bt.at1, bt.at2, bt.order, id );
            //bt.at1--; bt.at2--;
            bonds[id]=bt;
        }
        return i;
    }

    bool getBondParams( int atyp1, int atyp2, int btyp, double& l0, double& k ){
        uint64_t id  = getBondTypeId( atypes[atyp1].iZ, atypes[atyp2].iZ, btyp );
        //printf( "(%i,%i,%i) -> %i \n", atypes[atyp1].iZ, atypes[atyp2].iZ, btyp, id );
        //uint64_t id  = getBondTypeId( atyp1, atyp2, btyp );
        auto it = bonds.find(id);
        if   ( it == bonds.end() ){ l0=default_bond_length; k=default_bond_stiffness; return false;}
        else                      { l0=it->second.length;   k=it->second.stiffness;   return true; }
    }

    void fillBondParams( int nbonds, Vec2i * bond2atom, int * bondOrder, int * atomType, double * bond_0, double * bond_k ){
        printf("fillBondParams: %i\n", nbonds);
        for(int i=0; i<nbonds; i++){
            Vec2i ib = bond2atom[i];
            getBondParams( atomType[ib.x], atomType[ib.y], bondOrder[i], bond_0[i], bond_k[i] );
            //printf( "%i (%i %i) %i %g %g \n", i, atomType[ib.x], atomType[ib.y], bondOrder[i], bond_0[i], bond_k[i] );
        }
    }

};

#endif﻿
