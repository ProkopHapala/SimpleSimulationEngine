
#ifndef MolecularConfiguration_h
#define MolecularConfiguration_h

#include "datatypes.h"
#include "Vec3.h"
#include "quaternion.h"
#include "BatchBuff.h"

template<int N>
struct MolecularConfiguration{
    float poses[N*8]; //  pos.xyz,E,
};

float confDistance( int ntypes, int2* types, float* confI_, float* confJ_, float cpos, float crot ){
    double dist2 = 0;

    for(int ityp=0; ityp<ntypes; ityp++){
        int2 typ=types[ityp];
        Quat4f* confI=(Quat4f*)(confI_+typ.x*8);
        for(int i=0; i<typ.y; i++){
            const Quat4f& mposi = confI[0];
            const Quat4f& qroti = confI[1];
            float r2,r2min = 1e+30;
            Quat4f* confJ=(Quat4f*)(confJ_+typ.x*8);
            for(int j=0; j<typ.y; j++){
                r2  = cpos * mposi.f.dist2 ( confJ[0].f );
                r2 += crot * qroti.dist_cos( confJ[1]   );
                if(r2<r2min) r2min=r2;
                confJ+=2;
            }
            dist2+=r2;
            confI+=2;
        }
    }
    return dist2;
}

class MolecularConfigurationMemory { public:
    int    nTypes = 0;
    int    nMol   = 0;
    int2*  types  = 0;
    //std::vector<AtomicConfiguration> confs;
    BatchBuff<float> confs;

    void init(int nTypes_, int* types_){
        nTypes=nTypes_;
        //types=types_;
        nMol=0;
        _realloc( types, nTypes );
        for(int i=0; i<nTypes; i++){
            types[i].x  = nMol;
            types[i].y  = types_[i];
            nMol       += types_[i];
        }
        confs.setBatchSize( 8*nMol );
    }

};

#endif
