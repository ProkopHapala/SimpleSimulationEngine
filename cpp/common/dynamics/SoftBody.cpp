
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include <SoftBody.h>
#include <Cholesky.h>

// ==== Dynamics

void SoftBody::cleanForces( ){
	for( int i=0; i<npoints; i++ ){ forces[i]=Vec3dZero; }
}

void SoftBody::evalPointForces( ){
	for( int i=0; i<npoints; i++ ){
		//forces[i].set(0.0);
		evalPointForce( i, gravity, airFlow );
		//printf( " %i  %f %f %f   %f %f %f \n", i, forces[i].x, forces[i].y, forces[i].z, mass[i], invMass[i], drag[i] );
	}
}

void SoftBody::evalBondForces( ){
	if(fmax>1e+12){ for( int i=0; i<nbonds; i++ ){	addBondForce  ( bonds[i]       ); } }
	else          { for( int i=0; i<nbonds; i++ ){	addBondForceGS( bonds[i], fmax ); } }
}

void SoftBody::evalKinkForces( ){
    if(kinks){
        for( int i=0; i<nkink; i++ ){
            //addKinkForce( kinks[i] );
            dampKink( kinks[i] );
        }
	}
}

/*
void SoftBody::evalForces( ){
	//printf( "==============\n" );
	//for( int i=0; i<nbonds;  i++ ){	addBondForce( bonds[i] ); }
	evalPointForces( );
	evalStickForces( );
	//for( int i=0; i<npoints; i++ ){ printf( " point force %i  %f %f %f   %f %f %f \n", i, forces[i].x, forces[i].y, forces[i].z, mass[i], invMass[i], drag[i] ); }
}
*/

void SoftBody::evalForceLinearizedBonds( ){
	for( int i=0; i<nbonds; i++ ){
        const Bond& b  = bonds[i];
        Vec3d hat = dirs[i];
        double dl = hat.dot( disps[b.j] - disps[b.i] );
        if( dl>0 ){ dl*=b.type->kTens; }else{ dl*=b.type->kPress; };
        hat.mul( dl+f0s[i] );
        forces[b.i].sub(hat);
        forces[b.j].add(hat);
	}
}

void SoftBody::linearizedBonds( ){
    //for( int i=0; i<; i++){ disps[i].set(0.0); };
	for( int i=0; i<nbonds; i++ ){
        const Bond& b = bonds[i];
        Vec3d d  = points[b.j] - points[b.i];
        double l = d.normalize();
        dirs[i]  = d;
        double dl = l-b.l0;
        if( dl>0 ){ dl*=b.type->kTens; }else{ dl*=b.type->kPress; };
        f0s[i]   = dl;
	}
}

void SoftBody::disp2pos( ){
    for( int i=0; i<npoints; i++){
        points[i].add( disps[i] );
        disps[i].set(0.0);
    };
}

void SoftBody::applyConstrains(){
    if(fix){
	for( int i=0; i<nfix; i++ ){
		int ip = fix[i];
		if( (ip>0)&&(ip<npoints) ){
            forces    [ip].set(0.0);
            velocities[ip].set(0.0);
		}
	}
	}
}

void SoftBody::move_LeapFrog( ){
	for (int i=0; i<npoints; i++){
        Vec3d f=forces[i];
        Vec3d v=velocities[i];
        // damping
		v.mul( 1-damp*dt );
		//double cfv = f.dot(v);
		//if(cfv<0){ f.mul( 1-damp_fvdot*cfv ); } // this would break action-reaction law
		v        .add_mul( f, invMass[i] * dt );
		points[i].add_mul( v, dt );
		velocities[i] = v;
/*
		printf( " (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) %3.3f %3.3f %3.3f \n",
            points[i].x, points[i].y, points[i].z,
            velocities[i].x, velocities[i].y, velocities[i].z,
            forces[i].x, forces[i].y, forces[i].z,
            invMass[i] ,  dt,  damp
        );
*/
	}
}

int SoftBody::relaxStepGS( double errMax ){
    // Gauss-Seidel relaxation step
	int ndone=0;
	for( int i=0; i<nbonds; i++ ){
        const Bond& bond = bonds[i];
        Vec3d d; d.set_sub( points[bond.i], points[bond.j] );
        double l  = d.norm();
        double dl = l - bond.l0; // todo - we shoud ignore ropes in pressure
        if(dl>errMax){
            double f  = bond.getForce( dl );
            d.mul(dl/l);
            points[bond.j].add( d );
            points[bond.i].sub( d );
        }
    }
    return ndone;
}

void SoftBody::step( ){
    //for (int i=0; i<npoints; i++){ forces[i].set(0.0); }
    //evalForces( );
    cleanForces    ();
    evalPointForces();
    evalBondForces ();
    applyConstrains();
    move_LeapFrog( );
}



void SoftBody::rhs_ProjectiveDynamics(
    double dt, Vec3d* ps, Vec3d* b
    //int*    neighBs, int* neighBs_sizes, int np,
    //Vec2i*  bonds, int nbonds,
    //double* masses, double dt, double* ks,
    //Vec3d*  points, double* l0s,
    //Vec3d*  pnew, Vec3d* b
) {
    double idt2 = 1.0 / (dt * dt);
    for (int i = 0; i < npoints; i++) {
        Vec3d bi;
        bi.set_mul( ps[i], mass[i] * idt2 );
        //int neighB_start = (i == 0) ? 0 : neighBs_sizes[i-1];
        //int neighB_end   = neighBs_sizes[i];
        int nb0 = i*N_MAX_NEIGH;
        int nb1 = nb0 + N_MAX_NEIGH;
        for (int nb  = nb0; nb < nb1; nb++) {
            int ib   = neighBs[nb];
            const Bond& b = bonds[ib]; 
            double k = b.k;
            int i_ = bonds[ib].i;
            int j_ = bonds[ib].j;
            int j = (i_ == i) ? j_ : i_;
            Vec3d  d = ps[i] - ps[j];
            bi.add_mul(d,k*b.l0/d.norm());  
        }
        b[i] = bi;
    }
}

void SoftBody::run_cholesky(
    int niter
    //int np, int nbonds,
    //Vec3d* points, Vec3d* velocity, int* fixed, int nfixed,
    //int* neighBs, int* neighBs_sizes,
    //Vec2i* bonds, double* masses, double* ks, double* l0s,
    //double* LDLT_L, double* LDLT_D, int* neighsLDLT,
    //void (*eval_forces)(Vec3d* points, Vec3d* velocity, Vec3d* force, int np),
    //double dt, int niter
){
    Vec3d* ps_cor  = new Vec3d[npoints];
    Vec3d* ps_pred = new Vec3d[npoints];
    //Vec3d* force   = new Vec3d[npoints];
    Vec3d* b       = new Vec3d[npoints];

    //double* zz = new double[npoints*3];
    double* yy = new double[npoints*3];

    memcpy(ps_cor, points, npoints * sizeof(Vec3d));

    const int m=3; 

    for (int iter = 0; iter < niter; iter++) {
        // Evaluate forces
        //eval_forces(points, velocity, force, npoints);

        // Predict step
        double dt2 = dt*dt;
        for (int i = 0; i < npoints; i++) {    ps_pred[i] = points[i] + velocities[i]*dt + forces[i]*dt2; }

        //for (int i = 0; i < nfixed; i++) {  ps_pred[fixed[i]] = points[fixed[i]];}  // Apply fixed constraints

        rhs_ProjectiveDynamics( dt, ps_pred, b);   // Compute right-hand side

        // Solve using LDLT decomposition (assuming solve_LDLT_sparse is adapted for Vec3d)
        //Lingebra::solve_LDLT_sparse( npoints ,3, LDLT_L, LDLT_D, (double*)b, (double*)ps_cor, neighsLDLT );

        Lingebra::forward_substitution_sparse( npoints,m,  LDLT_L, (double*)b, yy, neighsLDLT );
        for (int i=0; i< npoints; i++){  double c=(1/LDLT_D[i]);   for(int s=0;s<m;s++){ yy[i*m+s]*=c; }; } // Diagonal 
        Lingebra::forward_substitution_transposed_sparse( npoints,m, LDLT_L, yy, (double*)ps_cor, neighsLDLT );

        // Compute residual
        // double res = 0.0;
        // for (int i = 0; i < npoints; i++) {
        //     Vec3d d = ps_cor[i] - points[i];
        //     double l = d.norm();
        //     if (l > res) res = l;
        // }
        // printf("residual[%d] : %f\n", iter, res);

        // Update velocity and points
        for (int i = 0; i < npoints; i++) {
            Vec3d dv = (ps_cor[i] - ps_pred[i])*(1/dt);
            velocities[i].add(dv);
            points[i] = ps_cor[i];
        }
    }

    // free(ps_cor);
    // free(ps_pred);
    // free(b);
    delete []ps_cor; 
    delete []ps_pred; 
    delete []b;
    delete []yy;
}


// ==== Setup

void SoftBody::deallocateAll( ){
    if( points    ) delete points;       //printf( "DEBUG delete points\n" );
    if( velocities) delete velocities;   //printf( "DEBUG delete velocities\n" );
    if( forces    ) delete forces ;      //printf( "DEBUG delete forces\n" );
    if( mass      ) delete mass;         //printf( "DEBUG delete mass\n" );
    if( drag      ) delete drag;         //printf( "DEBUG delete drag\n" );
    if( invMass   ) delete invMass;      //printf( "DEBUG delete invMass\n" );
    if( bonds     ) delete bonds;        //printf( "DEBUG delete bonds\n" );
    if( fix       ) delete fix;          //printf( "DEBUG delete fix\n" );
}

void SoftBody::allocate( int npoints_, int nbonds_, int nfix_ ){
    npoints=npoints_; nbonds=nbonds_; nfix=nfix_;
	points     = new Vec3d [npoints];
	mass       = new double[npoints]; for(int i=0; i<npoints; i++ ){ mass[i]=0.0d; }
	drag       = new double[npoints]; for(int i=0; i<npoints; i++ ){ drag[i]=0.0d; }
	velocities = new Vec3d [npoints ];
	forces     = new Vec3d [npoints ];
	invMass    = new double[npoints ];
    gravity.set(0.0,-9.81,0.0);
    airFlow.set(0.0,0.0,0.0);

	bonds  = new Bond[ nbonds ];

    if(nfix>0)fix = new int[nfix]; for(int i=0; i<nfix; i++ ){ fix[i]=-1; }
}

void SoftBody::setPoints( int npoints_,  Vec3d  * points_, double * mass_, double * drag_ ){
	npoints = npoints_;
	if( points_ ){ points = points_; }else{  points=new Vec3d [npoints]; }
	if( mass_   ){ mass   = mass_;   }else{  mass  =new double[npoints]; for(int i=0; i<npoints; i++ ){ mass[i]=0.0d; } }
	if( drag_   ){ drag   = drag_;   }else{  drag  =new double[npoints]; for(int i=0; i<npoints; i++ ){ drag[i]=0.0d; } }
	velocities = new Vec3d [npoints ];
	forces     = new Vec3d [npoints ];
	invMass    = new double[npoints ];
    gravity.set(0.0,-9.81,0.0);
    airFlow.set(0.0,0.0,0.0);
}

void SoftBody::setConstrains( int nfix_, int  * fix_  ){
    nfix = nfix_;
    if( fix_ != NULL ){ fix = fix_;  }else{ fix = new int[nfix]; for(int i=0; i<nfix; i++ ){ fix[i]=-1; } };
}

void SoftBody::setBonds( int nbonds_, int * ips, int * its, BondType * bts ){
    nbonds = nbonds_;
    bonds  = new Bond[ nbonds ];
    if(ips&&its&&bts){
        for(int i=0; i<nbonds; i++){
            int i2=i<<1;
            bonds[i].id     = i;
            bonds[i].i      = ips[i2  ];
            bonds[i].j      = ips[i2+1];
            bonds[i].type   = &bts[its[i]];
            //printf("%i %i (%i,%i) (%3.3g,%3.3g) (%3.3g,%3.3g) %3.3g \n", i, bonds[i].id, bonds[i].i, bonds[i].j, bonds[i].type->kTens, bonds[i].type->kPress, bonds[i].type->sPress, bonds[i].type->sPress, bonds[i].type->linearDensity );
            bonds[i].broken = false;
        }
    };
}

int SoftBody::findBonds( double lmax, BondType * bt ){
    int nbmax = npoints*(npoints-1)/2;
    int * bs  = new int[nbmax*2];
    int   n2  = 0;
    for(int i=0;i<npoints; i++){
        for(int j=0;j<i; j++){
            double r = getBondLength( i, j ); // TODO : this can be optimized
            if(r<lmax){
                bs[n2  ]=i;
                bs[n2+1]=j;
                n2+=2;
            }
        };
    };
    nbonds = n2>>1;
    bonds = new Bond[nbonds];
    for(int i=0; i<nbonds; i++){
        int i2=i<<1;
        bonds[i].id     = i;
        bonds[i].i      = bs[i2  ];
        bonds[i].j      = bs[i2+1];
        bonds[i].type   = bt;
        bonds[i].broken = false;
    }
    delete [] bs;
    return nbonds;
}

#include "testUtils.h"
int SoftBody::findKinks( double damp, double kstiff ){
    int*   npls = new int  [npoints];
    Vec2i* ijs  = new Vec2i[npoints];
    DEBUG
    for(int i=0; i<npoints; i++){ npls[i]=0; }
    auto func = [&](int ip, int ib){ int ni=npls[ip]; if(ni<2){ ijs[ip].array[ni]=ib; }; npls[ip]=ni+1; };
    for(int i=0; i<nbonds; i++){
        const Bond& b = bonds[i];
        //ni=npls[b.i];
        //npls[b.i]++;
        //npls[b.j]++;
        func(b.i,i);
        func(b.j,i);
    }
    nkink=0;
    for(int i=0; i<npoints; i++){ if(npls[i]==2) nkink++; }
    //kinks = new Kink[nkink];
    //printf( "nkink %i \n", nkink );
    _realloc( kinks, nkink );
    int ik=0;
    for(int i=0; i<npoints; i++){
        if(npls[i]==2){
            int a,b,ib;
            ib=ijs[i].a; if(bonds[ib].i==i){a=bonds[ib].j;}else{a=bonds[ib].i;}
            ib=ijs[i].b; if(bonds[ib].i==i){b=bonds[ib].j;}else{b=bonds[ib].i;}
            //printf( "kink[%i] (%i,%i, %i) \n", ik, a,b, i );
            kinks[ik]=(Kink){ a,b,i,  damp, kstiff  };
            ik++;
        }
    }
    delete [] ijs;
    delete [] npls;
    //exit(0);
    return nkink;
};


void SoftBody::prepareBonds( bool l0_fromPos ){
    for( int i=0; i<nbonds; i++ ){
        Bond& bond = bonds[i];
        if( l0_fromPos ){  bond.l0 = getBondLength( bond.i, bond.j ); }
        double dmass = bond.getMass() * 0.5;
        double ddrag = bond.getDrag() * 0.5;
        mass[bond.i] += dmass;
        mass[bond.j] += dmass;
        drag[bond.i] += ddrag;
        drag[bond.j] += ddrag;
        //printf( " bond %i %i %e %e %e \n", i, bond.type.id, bond.type.linearDensity, bond.type.kPress, bond.type.kTens );
        //printf( " bond %i %f %f  %f \n", i, bond.l0, dmass, ddrag );
    }
}

void SoftBody::preparePoints(  bool clearVelocity, double constDrag, double constMass ){
    for (int i=0; i<npoints; i++){
        if( clearVelocity ) velocities[i].set(0.0);
        if( constDrag >= 0 ) drag[i] = constDrag;
        if( constMass >  0 ) mass[i] = constMass;
		invMass[i] = 1/mass[i];
		//printf( " poins %i %e %e   (%e,%e,%e) \n", i, mass[i], invMass[i], velocities[i].x, velocities[i].y, velocities[i].z  );
	}
}



/*
void allocate_RVF          ( ){  }
void deallocate_RVF        ( ){ delete pos; delete vel; delete force; }
void allocate_pointParams  ( ){ mass = new double[npoints]; drag = new double[npoints]; invMass = new double[npoints]; }
void deallocate_pointParams( ){ delete mass; delete drag; delete invMass; }
void init_Axuliary         ( ){ for( int i=0; i<npoints; i++ ){ invMass = 1/mass; } }

void allocate_Bonds        ( ){ bonds = new Bonds[nbonds]; }
void deallocate_Bonds      ( ){ delete bonds; }

void init( int npoints_, int nbonds_, int nfix_, bool own_points_, bool own_mass_ ){
    own_points = own_points_;
    own_mass   = own_mass_;
    own_fix    = own_fix_;
    npoints = npoints_;
    nbonds  = nbonds_;
    nfix    = nfix_;
    if(own_points){
        pos = new Vec3d[npoints]; vel = new Vec3d[npoints]; force = new Vec3d[npoints];
    }


}
*/


/*
void SoftBody::bondFromType( int ib, int * bondTypes, const BondTypes& bondTypesBooks ){
	int itype = bondTypes[ib];
	int ib2 = ib<<1;
	int i   = bonds[ib2  ];
	int j   = bonds[ib2+1];
	Vec3d d; d.set_sub( points[i], points[j] );
	l0s[ib]			= d.norm();
	double area     = bondTypesBooks.area   [itype];
	double bondMass = bondTypesBooks.density[itype] * area * l0s[ib];
	kTens [ ib ]    = bondTypesBooks.kTens  [itype] * area;
	kPres [ ib ]    = bondTypesBooks.kPres  [itype] * area;
	//double tensStrength = bondTypes_tensStrength[i] * area;
	//double presStrength = bondTypes_presStrength[i] * area;
	mass[i] += bondMass*0.5;
	mass[j] += bondMass*0.5;
	//printf( " %i %i   %i   %i %i   %f %f %f %f  \n", ib, ib2, itype, i, j,   bondMass, l0s[ib], kTens[ib], kPres[ib]   );
}

*/

/*
void SoftBody::draw( float forceScale ){
	glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
	for( int ib=0; ib<nbonds;  ib++ ){
		int ib2 = ib<<1;
		Vec3f pi,pj, d;
		convert( points[ bonds[ib2   ] ], pi );
		convert( points[ bonds[ib2+1 ] ], pj );
		d.set_sub( pi, pj );
		float l  = d.norm(); // this should be optimized
		float dl = ( l - l0s[ib] ) / l;
		if( dl > 0 ){
			float f = (float)(kTens[ib]*dl)*forceScale;
			glColor3f( f,0,0 );
		}else{
			float f = (float)(kPres[ib]*dl)*forceScale;
			glColor3f( 0,0,-f );
		}
		glVertex3d( pi.x, pi.y, pi.z );
		glVertex3d( pj.x, pj.y, pj.z );
		//printf( " %f %f %f   %f %f %f \n", pi.x, pi.y, pi.z,   pj.x, pj.y, pj.z  );
	}
	glEnd();
}
*/


/*

void SoftBody::init(
	int npoints_, int nbonds_, int nfix_,
	Vec3d  * points_, double * mass_, double * drag_,
	int    * bonds_, int * bondTypes, const BondTypes& bondTypeBooks,
	int    * fix_
){
	npoints = npoints_;  nbonds = nbonds_; nfix = nfix_;
	points  = points_; bonds   = bonds_; fix = fix_;
	if( mass_ != NULL ){ mass = mass_; }else{  mass=new double[npoints]; for(int i=0; i<npoints; i++ ){ mass[i]=0; } }
	if( drag_ != NULL ){ drag = drag_; }else{  drag=new double[npoints]; for(int i=0; i<npoints; i++ ){ drag[i]=0; } }
	velocities = new Vec3d[npoints];
	forces     = new Vec3d[npoints];
	invMass    = new double[ npoints ];

	l0s   = new double[nbonds];
	kTens = new double[nbonds];
	kPres = new double[nbonds];
	for (int i=0; i<nbonds; i++ ){ bondFromType( i, bondTypes, bondTypeBooks ); }
	for (int i=0; i<npoints; i++){
		invMass[i] = 1/mass[i]; velocities[i].set(0.0);
		// printf( " %f   %f %f %f \n", invMass[i], velocities[i].x, velocities[i].y, velocities[i].z  );
	}
}

void SoftBody::init(
	int npoints_, int nbonds_, int nfix_,
	Vec3d  * points_, double * mass_,    double * drag_,
	int    * bonds_,  double * kTens_,   double * kPres_,   double * l0s_,
	int    * fix_
){
	npoints = npoints_;  nbonds = nbonds_; nfix = nfix_;
	points = points_;    mass    = mass_;
	bonds   = bonds_;    kTens      = kTens_;      kPres = kPres_; l0s = l0s_;
	if( drag_ != NULL ){ drag = drag_; }else{  drag=new double[npoints]; for(int i=0; i<npoints; i++ ){ drag[i]=0; } }
	fix     = fix_;
	velocities = new Vec3d[npoints];
	forces     = new Vec3d[npoints];
	invMass    = new double[ npoints ];
	for (int i=0; i<npoints; i++){ invMass[i] = 1/mass[i]; velocities[i].set(0.0); }
}

*/

