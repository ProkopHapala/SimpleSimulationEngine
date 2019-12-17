
#include "Truss.h"

//class Truss{ public:

void Truss::clear(){
    points.clear();
    edges.clear();
    blocks.clear();
    removed_points.clear();
    removed_edges.clear();
}

void Truss::sticksFormString( char * str ){
    //puts( str );
    char * pch = strtok (str,";\n");
    while (pch != NULL){
        //printf ("%s\n",pch);
        TrussEdge edge;
        edge.fromString( pch );
        edges.push_back(edge);
        //char str_tmp[32];
        //edges[edges.size()-1].print();
        //edges[edges.size()-1].toString(str_tmp);
        //printf ("%s\n", str_tmp );
        pch = strtok (NULL, ";\n");
    }
}

int Truss::loadXYZ( char* fname ){
    FILE * pFile = fopen(fname,"r");
    if( pFile == NULL ){
        printf("cannot find %s\n", fname );
        return -1;
    }
    char buff[1024];
    char * line;
    int nl;
    line = fgets( buff, 1024, pFile );
    int nvert, nedges;
    sscanf( line, "%i %i\n",  &nvert, &nedges );
    printf(       "%i %i \n",  nvert,  nedges );

    line = fgets( buff, 1024, pFile );
    sticksFormString(line); // load stricks // like comment in .xyz file

    for(int i=0; i<nvert; i++){
        char at_name[8];
        line = fgets( buff, 1024, pFile );  //printf("%s",line);
        Vec3d p;
        sscanf( line, "%s %lf %lf %lf \n", at_name, &p.x, &p.y, &p.z );
        printf(       "%s %lf %lf %lf \n", at_name,  p.x,  p.y,  p.z );
        points.push_back(p);
    }
    printf( "loadXYZ DONE\n");
    return nvert;
}

void Truss::affineTransform( Mat3d mat, bool T ){
    if( T ) { for(int i=0; i<points.size(); i++){ points[i] = mat.dotT(points[i]); }; }
    else    { for(int i=0; i<points.size(); i++){ points[i] = mat.dot (points[i]); }; }
}

int Truss::pickVertex( const Vec3d &ray0, const Vec3d &hRay ) const {
    double r2min=1e+300;
    int imin=0;
    for(int i=0; i<points.size(); i++){
        double t;
        double r2 = rayPointDistance2( ray0, hRay, points[i], t );
        if(r2<r2min){ imin=i; r2min=r2; }
    }
    return imin;
};

int Truss::pickVertex( const Vec3d& ray0, const Vec3d& hRay, double R ) const {
    //double tmin =  1e+300;
    double r2min =  R*R;
    int imin     = -1;
    for(int i=0; i<points.size(); i++){
        double t;
        double r2 = rayPointDistance2( ray0, hRay, points[i], t );
        //printf( "%i %i %g %g \n", i, imin, r2, r2min );
        if(r2<r2min){ imin=i; r2min=r2; }
        //double ti = raySphere( ray0, hRay, R, ps[i] );
        //if(ti<tmin){ imin=i; tmin=ti; }
    }
    return imin;
}

int Truss::pickEdge( const Vec3d& ray0, const Vec3d& hRay, double R ) const {
    double dist_min =  R;
    int    imin     = -1;
    for(int ib=0; ib<edges.size(); ib++){
        TrussEdge iat = edges[ib];
        double t1,t2;
        Vec3d p0 = points[iat.a];
        Vec3d d  = points[iat.b] - p0;
        double l = d.normalize();
        double dist = rayLine( ray0, hRay, p0, d, t1, t2 );
        if( (dist<dist_min) && (t2>0) && (t2<l) ){
            imin=ib; dist_min=dist;
        }
    }
    return imin;
};

void Truss::panel( Vec3d p00, Vec3d p01, Vec3d p10, Vec3d p11, Vec2i n, double width ){
    int kind_long   = 0;
    int kind_perp   = 1;
    int kind_zigIn  = 2;
    int kind_zigOut = 3;

    //int dnb = 2+4+4+4;
    int i0 = points.size();

    blocks.push_back( {i0,edges.size()} );

    Vec2d step = {1.0/n.a,1.0/n.b};

    int di = 2*n.a-1;

    int i00 = points.size();
    for (int ib=0; ib<n.b; ib++){
        double db,mb;
        db = ib*step.b;     mb=1-db;
        Vec3d p0  = p00*mb + p10*db;
        Vec3d p1  = p01*mb + p11*db;
        db += 0.5*step.b; mb=1-db;
        Vec3d p0_ = p00*mb + p10*db;
        Vec3d p1_ = p01*mb + p11*db;
        for (int ia=0; ia<n.a; ia++){
            double da,ma;
            da = ia*step.a; ma = 1-da;
            Vec3d p   = p0 *ma + p1 *da;

            points.push_back( p              );


            int bi = i0+di; if( ib==n.b-2 )bi-=ia;
            int dia = 2;    if( ib==n.b-1 )dia=1;

            if( ia<(n.a-1) ) edges.push_back( (TrussEdge){i0,i0+dia,kind_perp} );
            if( ib<(n.b-1) ) edges.push_back( (TrussEdge){i0,bi    ,kind_perp} );

            if( (ia<(n.a-1))&&(ib<(n.b-1)) ){ // diagonal
                Vec3d p_  = p0_*ma + p1_*da;
                da += 0.5*step.a; ma=1-da;
                Vec3d p__ = p0_*ma + p1_*da;
                Vec3d up; up.set_cross( p_-p, p__-p ); up.normalize();
                points.push_back( p__ + up*width );

                if( ia<(n.a-2) ) edges.push_back( (TrussEdge){i0+1,i0+1+dia,kind_perp} );
                if( ib<(n.b-2) ) edges.push_back( (TrussEdge){i0+1,bi+1    ,kind_perp} );

                edges .push_back( (TrussEdge){i0+1,i0     ,kind_zigIn} );
                edges .push_back( (TrussEdge){i0+1,i0+dia ,kind_zigIn} );
                edges .push_back( (TrussEdge){i0+1,bi     ,kind_zigIn} );
                if( ib==n.b-2 )dia=1;
                edges .push_back( (TrussEdge){i0+1,bi+dia ,kind_zigIn} );

                i0++;
            }

            i0++;
        }
    }
}

void Truss::girder1( Vec3d p0, Vec3d p1, Vec3d up, int n, double width ){
    int kind_long   = 0;
    int kind_perp   = 1;
    int kind_zigIn  = 2;
    int kind_zigOut = 3;
    Vec3d dir = p1-p0;
    double length = dir.normalize();
    up.makeOrthoU(dir);
    Vec3d side; side.set_cross(dir,up);

    //print(dir); print(up); print(side); printf("dir up side \n");
    double dl = length/(2*n + 1);
    //int dnb = 2+4+4+4;
    int dnp = 4;
    int i00 = points.size();

    blocks.push_back( {i00,edges.size()} );

    for (int i=0; i<n; i++){
        int i01=i00+1; int i10=i00+2; int i11=i00+3;
        points.push_back( p0 + side*-width + dir*(dl*(1+2*i  )) );
        points.push_back( p0 + side*+width + dir*(dl*(1+2*i  )) );
        points.push_back( p0 + up  *-width + dir*(dl*(1+2*i+1)) );
        points.push_back( p0 + up  *+width + dir*(dl*(1+2*i+1)) );
        edges .push_back( (TrussEdge){i00,i01,kind_perp}  );
        edges .push_back( (TrussEdge){i10,i11,kind_perp}  );
        edges .push_back( (TrussEdge){i00,i10,kind_zigIn} );
        edges .push_back( (TrussEdge){i00,i11,kind_zigIn} );
        edges .push_back( (TrussEdge){i01,i10,kind_zigIn} );
        edges .push_back( (TrussEdge){i01,i11,kind_zigIn} );
        if( i<(n-1) ){
            edges.push_back( (TrussEdge){i10,i00+dnp,kind_zigOut} );
            edges.push_back( (TrussEdge){i10,i01+dnp,kind_zigOut} );
            edges.push_back( (TrussEdge){i11,i00+dnp,kind_zigOut} );
            edges.push_back( (TrussEdge){i11,i01+dnp,kind_zigOut} );
            edges.push_back( (TrussEdge){i00,i00+dnp,kind_long} );
            edges.push_back( (TrussEdge){i01,i01+dnp,kind_long} );
            edges.push_back( (TrussEdge){i10,i10+dnp,kind_long} );
            edges.push_back( (TrussEdge){i11,i11+dnp,kind_long} );
        }
        i00+=dnp;
    }
}

void Truss::girder1_caps( int ip0, int ip1, int kind ){
    int ipbeg = blocks.back().x;
    int ipend = points.size()-4;
    edges.push_back( (TrussEdge){ip0,ipbeg+0,kind} );
    edges.push_back( (TrussEdge){ip0,ipbeg+1,kind} );
    edges.push_back( (TrussEdge){ip0,ipbeg+2,kind} );
    edges.push_back( (TrussEdge){ip0,ipbeg+3,kind} );
    edges.push_back( (TrussEdge){ip1,ipend+0,kind} );
    edges.push_back( (TrussEdge){ip1,ipend+1,kind} );
    edges.push_back( (TrussEdge){ip1,ipend+2,kind} );
    edges.push_back( (TrussEdge){ip1,ipend+3,kind} );
}

void Truss::wheel( Vec3d p0, Vec3d p1, Vec3d ax, int n, double width ){
    int kind_long   = 0;
    int kind_perp   = 1;
    int kind_zigIn  = 2;
    int kind_zigOut = 3;
    Vec3d dir = p1-p0;
    double r = dir.normalize();
    ax.makeOrthoU(dir);
    Vec3d side; side.set_cross(dir,ax);
    //double dl = length/(2*n + 1);
    //int dnb = 2+4+4+4;
    print(dir); print(ax); print(side); printf("dir up side \n");
    int dnp = 4;
    int i00 = points.size();
    int i000 = i00;

    Vec2d  rot = {1.0,0.0};
    Vec2d drot; drot.fromAngle( M_PI/n );

    blocks.push_back( {i00,edges.size()} );
    for (int i=0; i<n; i++){
        int i01=i00+1; int i10=i00+2; int i11=i00+3;

        Vec3d R = dir*rot.a + side*rot.b;
        points.push_back( p0 +  R*(r-width) );
        points.push_back( p0 +  R*(r+width) );
        rot.mul_cmplx(drot);
        R       = dir*rot.a + side*rot.b;
        points.push_back( p0 + ax*-width + R*r );
        points.push_back( p0 + ax*+width + R*r );
        rot.mul_cmplx(drot);

        edges .push_back( (TrussEdge){i00,i01,kind_perp}  );
        edges .push_back( (TrussEdge){i10,i11,kind_perp}  );
        edges .push_back( (TrussEdge){i00,i10,kind_zigIn} );
        edges .push_back( (TrussEdge){i00,i11,kind_zigIn} );
        edges .push_back( (TrussEdge){i01,i10,kind_zigIn} );
        edges .push_back( (TrussEdge){i01,i11,kind_zigIn} );
        if( i<(n-1) ){
            edges.push_back( (TrussEdge){i10,i00+dnp,kind_zigOut} );
            edges.push_back( (TrussEdge){i10,i01+dnp,kind_zigOut} );
            edges.push_back( (TrussEdge){i11,i00+dnp,kind_zigOut} );
            edges.push_back( (TrussEdge){i11,i01+dnp,kind_zigOut} );
            edges.push_back( (TrussEdge){i00,i00+dnp,kind_long} );
            edges.push_back( (TrussEdge){i01,i01+dnp,kind_long} );
            edges.push_back( (TrussEdge){i10,i10+dnp,kind_long} );
            edges.push_back( (TrussEdge){i11,i11+dnp,kind_long} );
        }else{
            edges.push_back( (TrussEdge){i10,i000+0,kind_zigOut} );
            edges.push_back( (TrussEdge){i10,i000+1,kind_zigOut} );
            edges.push_back( (TrussEdge){i11,i000+0,kind_zigOut} );
            edges.push_back( (TrussEdge){i11,i000+1,kind_zigOut} );
            edges.push_back( (TrussEdge){i00,i000+0,kind_long} );
            edges.push_back( (TrussEdge){i01,i000+1,kind_long} );
            edges.push_back( (TrussEdge){i10,i000+2,kind_long} );
            edges.push_back( (TrussEdge){i11,i000+3,kind_long} );
        }
        i00+=dnp;
    }
}

void Truss::makeGriders( int nEdges, TrussEdge* edges, Vec3d* points, GirderParams* params, Vec3d * ups ){
    for( int i=0; i<nEdges; i++ ){
        TrussEdge& ed = edges[i];
        girder1( points[ed.a], points[ed.b], ups[i], params[i].n, params[i].widthX );
    }
};

void Truss::makeGriders( Truss plan, GirderParams* params, Vec3d * ups, std::vector<Vec2i>* ends ){
    int np0 = points.size();
    blocks.push_back( {np0,edges.size()} );
    for( int i=0; i<plan.points.size(); i++ ){
        points.push_back( plan.points[i] );
    }
    int kind_end = 0;
    for( int i=0; i<plan.edges.size(); i++ ){
        TrussEdge& ed = plan.edges[i];
        int np0i = points.size();
        girder1( plan.points[ed.a], plan.points[ed.b], ups[i], params[i].n, params[i].widthX );
        int np0j = points.size();
        edges.push_back( (TrussEdge){np0+ed.a,np0i+0,kind_end} );
        edges.push_back( (TrussEdge){np0+ed.a,np0i+1,kind_end} );
        edges.push_back( (TrussEdge){np0+ed.a,np0i+2,kind_end} );
        edges.push_back( (TrussEdge){np0+ed.a,np0i+3,kind_end} );
        edges.push_back( (TrussEdge){np0+ed.b,np0j-1,kind_end} );
        edges.push_back( (TrussEdge){np0+ed.b,np0j-2,kind_end} );
        edges.push_back( (TrussEdge){np0+ed.b,np0j-3,kind_end} );
        edges.push_back( (TrussEdge){np0+ed.b,np0j-4,kind_end} );
        if( ends ){
            ends->push_back( (Vec2i){np0i+0,i} );
            ends->push_back( (Vec2i){np0i+1,i} );
            ends->push_back( (Vec2i){np0i+2,i} );
            ends->push_back( (Vec2i){np0i+3,i} );
            ends->push_back( (Vec2i){np0j-1,i} );
            ends->push_back( (Vec2i){np0j-2,i} );
            ends->push_back( (Vec2i){np0j-3,i} );
            ends->push_back( (Vec2i){np0j-4,i} );
        }
    }
};

void Truss::autoBridge(int n, Vec2i * ips, double rmax, int kind ){
    //printf( "autoBridge beg %i %i n=%i\n", points.size(),edges.size(), n );
    blocks.push_back( (Vec2i){points.size(),edges.size()} );
    double R2 = rmax*rmax;

    for(int i=0; i<n; i++){
        Vec2i& ipi = ips[i];
        Vec3d   pi = points[ipi.a];
        //printf( "ipi (%i,%i) pi (%g,%g,%g) \n", ipi.a, ipi.b, pi.x, pi.y, pi.z );
        for(int j=0; j<i; j++){
            Vec2i& ipj = ips[j];
            if( ipi.b != ipj.b ){
                Vec3d  d   = points[ipj.a] - pi;
                double l2  = d.norm2();
                if(l2<R2){ edges.push_back( (TrussEdge){ipi.a,ipj.a,kind} ); }
            }
        }
    }
    //printf( "autoBridge end %i %i\n", points.size(),edges.size() );
    //exit(0);
};

Vec2i* Truss::getIJs(){
    Vec2i* ijs = new Vec2i[edges.size()];
    for(int i=0; i<edges.size(); i++){
        ijs[i].set(edges[i].a,edges[i].b);
    }
    return ijs;
};





