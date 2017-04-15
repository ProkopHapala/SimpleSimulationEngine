#ifndef  RigidMolecule_h
#define  RigidMolecule_h

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#define R2SAFE  1.0e-2
#define F2MAX   10.0

//  ============== Globals

constexpr float R_MAX = 1.8;
constexpr float R2MAX = R_MAX*R_MAX;

constexpr int nAtoms = 4;
constexpr int nMols  = 32;
float pos   [nMols*8];
float vel   [nMols*8];
float force [nMols*8];

Vec3f atomsT[nAtoms*nMols];

//Vec3f molAtoms   [nAtoms*3] = { 0.0f,0.0f,0.0f,   1.0f,0.0f,0.0f,   0.0f,1.0f,0.0f,   0.0f,0.0f,1.0f  };
Vec3f molAtoms     [nAtoms] = { 1.0f,1.0f,1.0f,   -1.0f,-1.0f,1.0f,   -1.0f,1.0f,-1.0f,   1.0f,-1.0f,-1.0f  };

void initParticles( int n, float * pos, double step, double drnd ){
    float root_n = pow(n,0.33333);
    int nz    = (int)(floor(root_n));
    int ny    = (int)(floor(sqrt((n/nz)+1)));
    int nx    = ((int)(n/(nz*ny)))+1;
    int nrest = n - nz*ny*nx;
    printf( "n %i %g   (%i,%i,%i) \n", n, root_n, nx, ny, nz );
    int i = 0;
    for(int iz=0; iz<nz; iz++){
        for(int iy=0; iy<ny; iy++){
            for(int ix=0; ix<nx; ix++){
                if(i<n){
                    int i8 = i<<3;
                    ((Vec3f*)(pos+i8))-> set(
                        (ix-0.5*nx)*step + randf(-drnd,drnd),
                        (iy-0.5*ny)*step + randf(-drnd,drnd),
                        (iz-0.5*nz)*step + randf(-drnd,drnd)
                    );
                    //((Quat4f*)(pos+i8+4))->fromUniformS3( {randf(), randf(), randf()} );
                    ((Quat4f*)(pos+i8+4))->setOne();
                }
                i++;
            }
        }
    }
}

inline void acum_force( const Vec3f& p1, const Vec3f& p2, Vec3f& f ){
    Vec3f dp; dp.set_sub( p2, p1 );
    float r2 = dp.norm2();
    if (r2 > R2MAX) return;
    //Draw2D::drawLine(p1,p2);
    float ir2 = 1/( r2 + R2SAFE );
    //float fr  = (0.2-ir2)*(R2MAX-r2);
    float fr  = (0.7-ir2)*(R2MAX-r2);
    f.add_mul( dp, fr );
}

void addAtomForce( int n, Vec3f* pos, Vec3f* force ){
    for(int i=0; i<n; i++){
        Vec3f f; f.set(0.0);
        Vec3f& pi = pos[i];
        for(int j=0; j<n; j++){ acum_force( pi, pos[j], f ); }
        force[i] = f;
    }
}

void transformAtoms( const Vec3f& pos, const Quat4f& rot, int npoints, Vec3f * points, Vec3f * Tpoints ){
    Mat3f T;
    rot.toMatrix( T);
    for( int i=0; i<npoints; i++ ){
        Vec3f Tp;
        T.dot_to_T(   points[i],  Tp );
        //T.dot_to(   points[i],  Tp );
        Tpoints[i].set_add( pos, Tp  );
    }
}

void transformAllAtoms( int nMols, int nAtoms, float * pos, Vec3f * molAtoms, Vec3f * atomsT ){
    int iatom0 = 0;
    for(int i=0; i<nMols; i++){
        int i8 = i<<3;
        transformAtoms( *(Vec3f*)(pos+i8), *(Quat4f*)(pos+i8+4), nAtoms, molAtoms, atomsT+iatom0 );
        iatom0  += nAtoms;
    }
}




/*
void forceFromPoints( int npoints, Vec3f * points, Vec3f * forces,  const Quat4f& q,  Vec3f& fp, Quat4f& fq ){
    //printf( "forceFromPoints ----------------\n" );
    for( int i=0; i<npoints; i++ ){
        //printf( " %i   %f %f %f   %f %f %f \n", i,   points[i].x, points[i].y, points[i].z,   forces[i].x, forces[i].y, forces[i].z  );
        q .addForceFromPoint( points[i], forces[i], fq );
        fp.add( forces[i] );
        //printf( " %i   %f %f %f   %f %f %f %f \n", i,   forces[i].x, forces[i].y, forces[i].z,    fp.x, fp.y, fp.z,  fq.x, fq.y, fq.z, fq.w  );
    }
}

void cleanPointForce( int npoints, Vec3d * forces ){	for( int i=0; i<npoints; i++ ){  forces[i].set(0.0);  }	}

float move_leap_frog( int n, Vec2f * pos, Vec2f * vel, Vec2f * force, float dt, float damp ){
    float cdamp = 1 - damp*dt;
    double f2max = 0;
    for(int i=0; i<n; i++){
        float f2 = force[i].norm2();
        if(f2>f2max) f2max = f2;
    }
    if(f2max > F2MAX ) dt *= sqrt( sqrt(F2MAX/f2max) );
    for(int i=0; i<n; i++){
        vel[i].mul(cdamp);
        vel[i].add_mul( force[i], dt );
        pos[i].add_mul( vel[i]  , dt );
        //wrap( pos[i] );
    }
    return f2max;
}
*/

#endif

