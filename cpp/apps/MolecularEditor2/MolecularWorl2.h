#ifndef MolecularWorld_h
#define MolecularWorld_h

/*
 - molecules (both rigid and soft) are converted to atoms and moved

*/


﻿class MolecularWorld2 : public MMFF{
    public:
    
    int        nx,ny,nz;
    CubicRuler ruler;
    std::vector<int> grid[nx*ny*nz];
    
    void atoms2grid(){
    
    }
    
    void shortRangeForce(){
    
    }

}

#endif
