
#ifndef  DrawVegetation_h
#define  DrawVegetation_h

int pushCylinderTris( int n, float r1, float r2, Vec3f base, Vec3f tip, Vec3f up, std::vector<Vec3f>& verts, std::vector<Vec3f>* normals ){
	int nvert=0;

	Vec3f dir,left;
	dir.set_sub( tip, base );
	dir.normalize();
    left = dir.getOrtho(up);

    float alfa = 2*M_PI/n;
    Vec2f rot,drot;
    rot .set(1.0f,0.0f);
    drot.set( cos( alfa ), sin( alfa ) );

	Vec3f q; q.set(dir); q.add_mul( up, -(r1-r2) );
	float pnab =  dir.dot( q )/q.norm();
	float pnc  =  sqrt( 1 - pnab*pnab );

	Vec3f op,opn;
	for(int i=0; i<=n; i++ ){
		Vec3f p,pn;
		p .set( rot.x*up.x + rot.y*left.x, rot.x*up.y + rot.y*left.y, rot.x*up.z + rot.y*left.z );
		pn.set( pnab*p.x   + pnc*dir.x   , pnab*p.y   + pnc*dir.y   , pnab*p.z   + pnc*dir.z    );
		if( i>0 ){
            verts  .push_back(base+op*r1); verts  .push_back(tip+p*r2); verts  .push_back(tip+op*r2);
            verts  .push_back(base+op*r1); verts  .push_back(tip+p*r2); verts  .push_back(base+p*r1);
            if(normals){
                normals->push_back(opn); normals->push_back(pn); normals->push_back(opn);
                normals->push_back(opn); normals->push_back(pn); normals->push_back(pn);
            }
		}
		op=p; opn=pn;
        rot.mul_cmplx( drot );
	}
	return nvert;
};

void pushTris_Treestep( int level, Vec3f pos, Vec3f dir, Vec3f up, std::vector<Vec3f>& verts, std::vector<Vec3f>* normals ){
    float l = dir.norm();
    static const float drnd = 0.15;
    dir.add( randf(-drnd,drnd)*l, randf(-drnd,drnd)*l, randf(-drnd,drnd)*l );
    float f = randf(0.5,0.9);
    dir.mul( f );
    up=cross(dir,up); up.normalize();
    Vec3f pos_ = pos + dir;
    pushCylinderTris( 6, 0.1*l, 0.1*l*f, pos, pos_, (Vec3f){0.0f,1.0f,0.0f},verts,normals);
    //pushCylinderTris(verts,normals);
    if( level>0 ){
        level--;
        pushTris_Treestep( level, pos_, dir+up*0.3*l, up, verts, normals );
        pushTris_Treestep( level, pos_, dir-up*0.3*l, up, verts, normals );
    }
}

void tree_step( int level, Vec3f pos, Vec3f dir, std::vector<Vec3f>& branches, std::vector<Vec3f>& leafs ){
    static const float drnd = 0.6;
    //dir.x *= randf(1.0-drnd,1.0);
    //dir.y *= randf(1.0-drnd,1.0);
    //dir.z *= randf(1.0-drnd,1.0);
    float l = dir.norm();
    dir.add( randf(-drnd,drnd)*l, randf(-drnd,drnd)*l, randf(-drnd,drnd)*l );
    dir.mul( randf(0.5,0.9) );
    Vec3f pos_ = pos + dir;
    branches.push_back(pos );
    branches.push_back(pos_);
    if( level==0 ){
        leafs.push_back(pos_);
    }else{
        level--;
        tree_step( level, pos_, dir,  branches, leafs );
        tree_step( level, pos_, dir,  branches, leafs );
    }
}

#endif
