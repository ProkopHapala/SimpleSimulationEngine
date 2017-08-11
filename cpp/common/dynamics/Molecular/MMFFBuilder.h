#ifndef MMFFBuilder_h
#define MMFFBuilder_h

#include <vector>

#include "Molecule.h"
#include "MMFF.h"


class MMFFAtom{ public:
    int type;
    Vec3d pos;
    Vec3d LJq;
};

class MMFFBond{ public:
    int    type;
    Vec2i  atoms;
    //double l0,k;
};

class MMFFAngle{ public:
    int type;
    Vec2i  bonds;
    //double a0,k;
};

class MMFFmol{ public:
    Molecule * mol;
    Vec3i    i0;
};

class MMFFBuilder{  public:
    std::vector<MMFFAtom>  atoms;
    std::vector<MMFFBond>  bonds;
    std::vector<MMFFAngle> angles;
    std::vector<MMFFmol>   mols;

    void insertMolecule( Molecule * mol, Vec3d pos, Mat3d rot ){
        int natom0  = atoms.size();
        int nbond0  = bonds.size();
        int nangle0 = angles.size();
        mols.push_back( (MMFFmol){mol, (Vec3i){natom0,nbond0,nangle0} } );
        for(int i=0; i<mol->natoms; i++){
            //Vec3d LJq = (Vec3d){0.0,0.0,0.0}; // TO DO : LJq can be set by type
            Vec3d LJq = (Vec3d){1.0,0.03,0.0}; // TO DO : LJq can be set by type
            Vec3d p; rot.dot_to(mol->pos[i],p); p.add( pos );
            atoms.push_back( (MMFFAtom){mol->atomType[i], p, LJq } );
        }
        for(int i=0; i<mol->nbonds; i++){
            bonds.push_back( (MMFFBond){mol->bondType[i], mol->bond2atom[i] + ((Vec2i){natom0,natom0}) } );
        }
        for(int i=0; i<mol->nang; i++){
            angles.push_back( (MMFFAngle){ 1, mol->ang2bond[i] + ((Vec2i){nbond0,nbond0}) } );
        }
    }

    void assignAtomTypes( MMFFparams * params ){
        for(int i=0; i<atoms.size(); i++){
            //mmff->aLJq [i]  = atoms[i].type;
            int ityp = atoms[i].type;
            atoms[i].LJq.x = params->atypes[ityp].RvdW;
            atoms[i].LJq.y = params->atypes[ityp].EvdW;
            atoms[i].LJq.z = 0;
            //atomTypes[i]  = atoms[i].type;
        }
    }

    void toMMFF( MMFF * mmff, MMFFparams * params ){
        mmff->deallocate();
        mmff->allocate( atoms.size(), bonds.size(), angles.size(), 0 );
        //int * atomTypes = new int[atoms.size()];
        //int * bondTypes = new int[bonds.size()];
        for(int i=0; i<atoms.size(); i++){
            mmff->atypes[i] = atoms[i].type;
            mmff->apos [i]  = atoms[i].pos;
            mmff->aLJq [i]  = atoms[i].LJq;
            //atomTypes[i]  = atoms[i].type;
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
        //params.fillBondParams( mmff->nbonds, mmff->bond2atom, bondTypes, atomTypes, mmff->bond_0, mmff->bond_k );
        //delete [] atomTypes;
        //delete [] bondTypes;
    }
};

int write2xyz( FILE* pfile, MMFF * mmff, MMFFparams * params ){
    fprintf(pfile, "%i\n", mmff->natoms );
    fprintf(pfile, "#comment \n");
    for(int i=0; i<mmff->natoms; i++){
        int ityp   = mmff->atypes[i];
        Vec3d&  pi = mmff->apos[i];
        printf( "%i %i (%g,%g,%g) %s \n", i, ityp, pi.x,pi.y,pi.z, params->atypes[ityp].name );
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

#endif
