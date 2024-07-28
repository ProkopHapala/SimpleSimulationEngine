
#include "Truss.h"

//class Truss{ public:

void Truss::clear(){
    points.clear();
    edges.clear();
    blocks.clear();
    removed_points.clear();
    removed_edges.clear();
}

/**
 * Parses a string and adds TrussEdges to Truss
 * 
 * @param str The string containing TrussEdge information separated by ';' or '\n'.
 */
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

/**
 * Loads the data of a truss from an XYZ file. The nodes are loaded form the body. The sticks are loaded from the comment line if format of three integers separated by ';' like "i1 j1 typ1; i2 j2 typ2; i3 j3 typ3;"
 * 
 * @param fname The file name of the XYZ file.
 * @return The number of vertices loaded from the file, or -1 if the file cannot be found.
 */
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
    sticksFormString(line); // load sticks from comment line in .xyz file

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

/**
 * Applies an affine transformation to the points of the truss.
 * 
 * @param mat The transformation matrix to apply.
 * @param T   If true, the transpose of the matrix is used for the transformation.
 *            If false, the matrix is used as is.
 */
void Truss::affineTransform( Mat3d mat, bool T, Vec3d p0, Vec3d p ){
    if( T ) { for(int i=0; i<points.size(); i++){ points[i] = mat.dotT(points[i]-p0)+p; }; }
    else    { for(int i=0; i<points.size(); i++){ points[i] = mat.dot (points[i]-p0)-p; }; }
}

/**
 * Picks the vertex in the truss that is closest to the given ray. Usefull e.g. for picking a vertex with the mouse.
 *
 * @param ray0 The starting point of the ray.
 * @param hRay The direction vector of the ray.
 * @return The index of the closest vertex in the truss.
 */
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

/**
 * Finds the index of the vertex in the Truss object that is closest to the given ray. Usefull e.g. for picking a vertex with the mouse.
 * 
 * @param ray0 The starting point of the ray.
 * @param hRay The direction vector of the ray.
 * @param R The radius of the sphere used for distance comparison.
 * @return The index of the closest vertex.
 */
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

/**
 * Finds the index of the closest edge to a given ray within a specified radius. Usefull e.g. for picking a edges with the mouse.
 *
 * @param ray0 The starting point of the ray.
 * @param hRay The direction vector of the ray.
 * @param R The radius within which to search for the closest edge.
 * @return The index of the closest edge, or -1 if no edge is found within the radius.
 */
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


void Truss::updateEdgesLengths(){
    for(TrussEdge& e : edges){ e.l0 = (points[e.a] - points[e.b]).norm(); }
};

void Truss::updateFacesAreas(){
    for(TrussFace& f : faces){
        Vec3d n; 
        //f.area  = normalAreaTriangle( n, points[f.a], points[f.b], points[f.c] ); // we want to be independent geom3d.h
        Vec3d a = points[f.a];
        Vec3d b = points[f.b];
        Vec3d c = points[f.c];
        n.set_cross( b-a, c-a );
        f.area = n.norm()*0.5;
    }
};

/*
int Truss::getMaxNeighs(){
    int nmax=0;
    for(TrussEdge& e : edges){ if(e.type>0) nmax++; }
    return nmax;
};
*/


/**
 * Creates a panel of truss elements between four corner points. Adds the points and edges to the Truss object.
 * // TODO: make also triangular panel
 * 
 * @param p00 The first corner point.
 * @param p01 The second corner point.
 * @param p10 The third corner point.
 * @param p11 The fourth corner point.
 * @param n The number of subdivisions along each side of the panel.
 * @param width The width of the truss elements.
 */
void Truss::panel( Vec3d p00, Vec3d p01, Vec3d p10, Vec3d p11, Vec2i n, double width ){
    printf( "Truss::panel() n(%i,%i) w=%g p00(%g,%g,%g) p01(%g,%g,%g) p10(%g,%g,%g) p11(%g,%g,%g) \n", n.x,n.y, p00.x,p00.y,p00.z, p01.x,p01.y,p01.z, p10.x,p10.y,p10.z, p11.x,p11.y,p11.z );
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

/**
 * Creates a girder in the truss structure.
 * 
 * @param p0 starting point
 * @param p1 ending point
 * @param up up-vector 
 * @param n  number of segments along the girder 
 * @param width The width of the girder.
 * @return The index of the created block containing the starting indexes of the points and edges
 */
int Truss::girder1( Vec3d p0, Vec3d p1, Vec3d up, int n, double width ){
    //printf( "Truss::girder1() n=%i p0(%g,%g,%g) p1(%g,%g,%g) up(%g,%g,%g) \n", n, p0.x,p0.y,p0.z, p1.x,p1.y,p1.z, up.x,up.y,up.z );
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

    int ibloc = blocks.size();
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
    return ibloc;
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

/*
function ngonTruss(p0::Vector{Float64}, p1::Vector{Float64}, ax::Vector{Float64}; n::Int=8, k::Float64=1.00 )
    #println("Truss::ngonTruss() n=$n p0=($p0[1],$p0[2],$p0[3]) p1=($p1[1],$p1[2],$p1[3]) ax=($ax[1],$ax[2],$ax[3])")

    dir  = p1 - p0
    r    = norm(dir)
    dir  = dir / norm(dir)
    ax   = make_ortho_u!(ax, dir)
    side = cross(dir, ax)

    rot  = 1.0 + 0.0im
    drot = cis(2*π / n)  # Using complex number representation for rotation

    points = Vector{Vector{Float64}}()
    edges  = Vector{TrussEdge}()

    io = n
    for i in 1:n
        R = dir * real(rot) + side * imag(rot)
        push!(points, p0 + R )
        push!(edges, TrussEdge(i, io, 1))
        rot *= drot  # Complex multiplication
        io=i
    end

    points = Matrix(hcat(points...)')

    bonds = Array{Tuple{Int, Int}, 1}(undef, 0)
    ks = Array{Float64, 1}(undef, 0)
    for e in edges
        push!(bonds, (e.i, e.j))
        push!(ks, k )
    end
    return points, bonds, ks
end
*/

int Truss::nGon( Vec3d p0, Vec3d p1, Vec3d ax, int n, int kind ){
    printf( "Truss::nGon() n=%i kind=%i p0(%g,%g,%g) p1(%g,%g,%g) ax(%g,%g,%g) \n", n, kind, p0.x,p0.y,p0.z, p1.x,p1.y,p1.z, ax.x,ax.y,ax.z );
    Vec3d dir = p1-p0;
    double r = dir.normalize();
    ax.makeOrthoU(dir);
    Vec3d side; side.set_cross(dir,ax);
    
    Vec2d  rot = {1.0,0.0};
    Vec2d drot; drot.fromAngle( 2.*M_PI/n );
    int i1   = points.size();
    int ibloc = blocks.size();
    blocks.push_back( {i1,edges.size()} );
    int i0 = i1+n-1;
    for (int i=0; i<n; i++){
        Vec3d R = dir*rot.a + side*rot.b;
        points.push_back( p0 +  R );
        edges .push_back( (TrussEdge){i0,i1,kind}  );
        rot.mul_cmplx(drot);
        i0=i1; i1++;
    }
    return ibloc;
}


/**
 * Creates a wheel-shaped truss structure.
 *
 * @param p0 starting point ( center of the wheel )
 * @param p1 ending point   ( on the perimeter of the wheel )
 * @param ax axis vector
 * @param n The number of segments along the perimeter of the wheel.
 * @param width The width of the wheel rim.
 * @return The index of the created block containing the starting indexes of the points and edges
 */
int Truss::wheel( Vec3d p0, Vec3d p1, Vec3d ax, int n, double width ){
    printf( "Truss::wheel() n=%i p0(%g,%g,%g) p1(%g,%g,%g) ax(%g,%g,%g) \n", n, p0.x,p0.y,p0.z, p1.x,p1.y,p1.z, ax.x,ax.y,ax.z );
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
    //print(dir); print(ax); print(side); printf("dir up side \n");
    int dnp = 4;
    int i00 = points.size();
    int i000 = i00;

    Vec2d  rot = {1.0,0.0};
    Vec2d drot; drot.fromAngle( M_PI/n );
    int ibloc = blocks.size();
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
    return ibloc;
}

/**
 * Creates griders from lists of points and edges.
 * 
 * @param nEdges number of edges
 * @param edges  array of TrussEdge objects  ToDo: we should rathe use Vec2i ( it is more general )
 * @param points array of node points
 * @param params array of GirderParams objects representing the parameters for each girder.
 * @param ups An array of Vec3d objects representing the up vectors for each girder.
 */
int Truss::makeGriders( int nEdges, TrussEdge* edges, Vec3d* points, GirderParams* params, Vec3d * ups ){
    int ib0 = blocks.size();
    for( int i=0; i<nEdges; i++ ){
        TrussEdge& ed = edges[i];
        girder1( points[ed.a], points[ed.b], ups[i], params[i].n, params[i].widthX );
    }
    return ib0;
};

/**
 * @brief Creates griders based on plan stored in a Truss object.
 *
 * This function adds griders to the truss based on the provided plan, girder parameters, ups vector, and ends vector.
 *
 * @param plan   truss plan to use for creating griders
 * @param params array of parameters for each girder
 * @param ups    array of up-vectors for each girder
 * @param ends   pointer to a vector of Vec2i to store the ends of each girder.
 * @return The index of the created block containing the starting indexes of the points and edges
 */
int Truss::makeGriders( Truss plan, GirderParams* params, Vec3d * ups, std::vector<Vec2i>* ends ){ 
    int np0 = points.size();
    int ib0 = blocks.size();
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
    return ib0;
};





// ToDo: makeCilinder     - for center of tanks and tubes
// ToDo: makeSphericalCap - for caps of tanks and exhaust nozzles
// ToDo: makeParaboloid   - for exhaust nozzle

int Truss::makeCylinder( Vec3d p0, Vec3d p1, double r0, double r1, int nPhi, int nL, double dStep, int edgeTyp, int faceTyp ){
    Vec3d ax   = p1-p0;  double L = ax.normalize();
    if( nPhi<0 ) nPhi = (int)(2*M_PI*fmax(r0,r1)/dStep); if(nPhi<6){ nPhi==6; } // make the elements are not longer than dStep
    if( nL  <0 ) nL   = (int)(L/dStep); if(nL<1)nL=1;
    Vec3d up,left;       ax.getSomeOrtho(up,left);
    Vec2d cph=Vec2dX, dph;
    double dr = (r1-r0)/nL;
    int ia0=0, ib0=0;
    for(int il=0; il<=nL; il++){ // longitudinal loop
        int ia_,ib_;
        double r = r0 + dr*il;   
        // one circle
        for(int iph=0; iph<=nPhi; iph++){ // loop over phi
            if(iph<nPhi){
                int ia = addPoint( p0 + left*(cph.x*r) + up*(cph.y*r) ); 
                addEdge( ia,(ia+1)%nPhi, edgeTyp );
                cph.mul_cmplx(dph);
            }
        } // end of loop over phi
        // between two circles
        if(il<nL){
            int ia0 = points.size()-1;
            int ia1 = ia0-nPhi;
            for(int iph=0; iph<nPhi; iph++){
                addEdge( iph+ia0,   iph         +ia1, edgeTyp );
                addEdge( iph+ia0, ((iph+1)%nPhi)+ia1, edgeTyp );
            }
        }
        // ToDo: add faces
        // ToDo: we can add edgest toward the center of the cylinder to make it more rigid
    }
    return 0;
}


/*
int Truss::makeCapsula( Vec3d p0, Vec3d p1, double r1, double r2, double theta1, double theta2, double dTheta, int nPhi, int edgeTyp, int faceTyp, bool capped ){
    // modified from Draw3D::drawCapsula()
    int nvert=0;
    Vec3d ax   = p1-p0;  double L = ax.normalize();
    Vec3d up,left;       ax.getSomeOrtho(up,left);
    Vec2d cph=Vec2dX, dph;
    dph.fromAngle( 2*M_PI/nPhi );
    // Cylinder
    Vec2d cth,dth;
    double dr = (r2-r1);
    double cv = sqrt(L*L+dr*dr);
    cth.set( L/cv, -dr/cv );
    
    // Central Cylinder part
    int ia0=0, ib0=0;
    for(int iph=0; iph<=nPhi; iph++){
        int ia_,ib_;
        if(iph<nPhi){
            int ia = addPoint( p0 + left*(cph.x*r1) + up*(cph.y*r1) ); 
            int ib = addPoint( p1 + left*(cph.x*r2) + up*(cph.y*r2) );
            if(iph==0){ ia0=ia; ib0=ib; }
            cph.mul_cmplx(dph);
        }
        if(iph>0){
            addEdge( ia_,ia, edgeTyp );
            addEdge( ib_,ib, edgeTyp );
            addEdge( ia, ib, edgeTyp );
            addEdge( ia, ib_,edgeTyp );
            ia_=ia; ib_=ib;
        }
        // ToDo: add faces
        nvert+=2;
    }
    
    double DTh,h;
    int nTheta;

    // Spherical Cap 1
    cph=Vec2dX;
    cth.set( L/cv, -dr/cv );
    DTh = (-theta1 - asin(cth.y));
    nTheta = (int)(fabs(DTh)/dTheta);
    dth.fromAngle( DTh/nTheta );
    //printf( " cth (%f,%f)  dth (%f,%f) \n", cth.x, cth.y,  dth.x, dth.y );
    r1/=cth.x;
    h  =-cth.y*r1;
    for(int ith=0; ith<(nTheta+1); ith++){
        Vec2d cth_ = Vec2d::mul_cmplx(cth,dth);
        for(int iph=0; iph<(nPhi+1); iph++){            
            int ia = addPoint( p0 + (left*(cph.x*r1) + up*(cph.y*r1))*cth.x  + ax*(h+cth.y*r1) ); 
            int ib = addPoint( p0 + (left*(cph.x*r1) + up*(cph.y*r1))*cth_.x + ax*(h+cth_.y*r1) );
            nvert+=2;
            cph.mul_cmplx(dph);

            // ToDo: add edges
            // ToDo: add faces

        }
        //printf( "%i cth (%f,%f)  cth_ (%f,%f) \n", ith, cth.x, cth.y,  cth_.x, cth_.y );
        cth=cth_;
    }

    // Spherical Cap 2
    cph=Vec2dX;
    cth.set( L/cv, -dr/cv );
    DTh    = (theta2-asin(cth.y));
    nTheta = (int)(fabs(DTh)/dTheta);
    dth.fromAngle(DTh/nTheta );
    r2/= cth.x;
    h  =-cth.y*r2;
    for(int ith=0; ith<(nTheta+1); ith++){
        Vec2d cth_ = Vec2d::mul_cmplx(cth,dth);
        for(int iph=0; iph<(nPhi+1); iph++){
            int ia = addPoint( p1 + (left*(cph.x*r2) + up*(cph.y*r2))*cth.x  + ax*(h+cth.y*r2) );
            int ib = addPoint(p1 + (left*(cph.x*r2) + up*(cph.y*r2))*cth_.x + ax*(h+cth_.y*r2) );
            nvert+=2;
            cph.mul_cmplx(dph);
        }
        cth=cth_;
    }
    return nvert;
}
*/

int Truss::addRope( int ip1, int ip2, int type, int nsub ){
    Vec3d& p1 =  points[ip1];
    Vec3d  d  = (points[ip2]-p1)*(1./nsub);
    int jp,ip=ip1;
    int np=0;
    for(int k=0; k<nsub; k++){
        if(k<(nsub-1)){
            points.push_back( p1+d*(k+1) ); np++;
            jp=points.size()-1;
        }else{
            jp=ip2;
        }
        edges.push_back( (TrussEdge){ip,jp,type} );
        ip=jp;
    }
    return np;
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

void Truss::massesToPoints( double* masses ){
    for( TrussEdge& e: edges ){
        masses[ e.a ] += e.mass*0.5;
        masses[ e.b ] += e.mass*0.5;
    }
    double c = 1./3;
    for( TrussFace& f: faces ){
        double m = f.mass*c;
        masses[ f.a ] += m;
        masses[ f.b ] += m;
        masses[ f.c ] += m;
    }
}

Vec2i* Truss::getIJs(){
    Vec2i* ijs = new Vec2i[edges.size()];
    for(int i=0; i<edges.size(); i++){
        ijs[i].set(edges[i].a,edges[i].b);
    }
    return ijs;
};





