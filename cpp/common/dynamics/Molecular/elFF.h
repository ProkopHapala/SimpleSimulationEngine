
#ifndef elFF_h
#define elFF_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Forces.h"

class ElFF { public:

    int maxA2W = 8;
    int maxW2A = 8;

    Vec3d* apos = 0; // atomi position
    Vec3d* wpos = 0; // position of wavepackets

    double*  acsums = 0; // summ of coeffs from  Sum_k{C_ki}
    //double*  rhos   = 0; // density matrix Sum_k Cki*Ckj
//    double*     = 0; //

    double* wcs = 0; // wave coefs
    int*    wfNs  = 0;
    int*    a2w   = 0; // all waves on atom[i]

    int*    ao2wc = 0; // atomic orbital -> wave function coeficient
    int*    wc2ac = 0; // wave function coeficient -> atomic orbital



    int*    w2a   = 0; // wave  -> atom 
    int*    w2b   = 0; // wave  -> basis
    int*    w2c   = 0; // wave  -> coeff   w2c = w2a*4 + w2b
    int*    wna   = 0; // number of atoms per
    int*    anw   = 0; /

// minimize spread of wave pockets
/*

- wave functions are spatially constrained. This would mean either cut-off distance, or limited number of coefficeitns
- limited number of coefficents is more efficient, and leads to simpler datastructures (const size), but it is questionable how form a bond in such case
- It may spowned stochastically in neighboring atoms during the relaxation

*/

void getH( Vec3d d, int i, int j ){
/*
    xy
    xx
    yy
    zz
    switch i:
*/
}


/*
int eval_dS(int i){
    // it is better implemented in atomic space
    for(int ia=0; ia++; ia<0){
        for(int ii=0; ii++; ii<0){
            int i = a2w[ii];
            
            for(int jj=0; jj<ii; jj<0){
                int j = a2w[jj];
                
            }
        }
    }
}
*/


int evalDerivs(int i){

// Sij
double Sij = 0;
int nk   = wfNs[j];
int ioff = i*maxW2A;
found
for(int k=0;k<nk; k++ ){ // loop over coefs
    int kk = ioff+k;
    int ak = w2a[kk]; // atom index
    double cik = wcs[kk];
    for(jj=0;jj<anw[ak];j++){
        j       = a2w[kk];
        Sij[j] +=  cik * wcs[joff+jk];
    }
};


// dS
/*
    int nk   = wfNs[j];
    int ioff = i*maxW2A;
    for(int k=0;k<nk; k++ ){ // loop over coefs
        int kk = ioff+k;
        double cik = wcs[kk];
        //int ja   = w2a[kk]; // atom index
        int ibas = w2b[kk];
        // Sij
        // 
        double sumCli = acsums[ibas]; 
    };
*/

}


}

#endif
