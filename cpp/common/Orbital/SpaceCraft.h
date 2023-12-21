
#ifndef  SpaceCraft_h
#define  SpaceCraft_h

#include <vector>

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "geom3D.h"

#include "TriangleRayTracer.h"
#include "Radiosity.h"

#include "Truss.h"

#include "SpaceCraftComponents.h"


namespace SpaceCrafting{

class SpaceCraft : public CatalogItem { public:
    bool bPrint = true;

    SpaceCraftWorkshop* workshop = 0;

    std::vector<int>       LODs;  // levels of detail for OpenGL rendering
    Truss truss;
    std::vector<Node>      nodes;
    
    // Node-Connectors ( linear structures )
    std::vector<Rope>      ropes;
    std::vector<Girder>    girders;
    std::vector<Ring>      rings;
	std::vector<Gun>       guns;

    // Plate components  ( planar structures )     ToDo: maybe we should make just array 'plates' and store all plate-like components there (e.g. radiators, shields, collectors, etc.)
	std::vector<Radiator>  radiators;
	std::vector<Shield>    shields;

    // Volumetric structures
	std::vector<Tank>      tanks;    // Also Capacitors ?
	std::vector<Pipe>      pipes;    // Also Cables, should be attached to girders

    // Shell structures
	std::vector<Thruster>  thrusters;
	std::vector<Balloon>  balloons;
	std::vector<Rock>     rocks;
	// Truss * coarse = NULL;
	// Truss * fine   = NULL;

	// aux
	int pickedTyp = -1;

    // ==== functions

	void clear(){
        nodes.clear(); ropes.clear(); girders.clear(); rings.clear(); thrusters.clear(); guns.clear(); radiators.clear(); shields.clear(); tanks.clear(); pipes.clear();
        rocks.clear(); balloons.clear();
        truss.clear();
	};

	inline void linker2line( const NodeLinker& o, Line3d& l ){ l.a=nodes[o.p0].pos; l.b=nodes[o.p1].pos; }
	inline void plate2quad ( const Plate& o, Quad3d& qd )const{
        //qd.l1.a=nodes[ girders[o.g1].p0 ].pos;
        //qd.l1.b=nodes[ girders[o.g1].p1 ].pos;
        //qd.l2.a=nodes[ girders[o.g2].p0 ].pos;
        //qd.l2.b=nodes[ girders[o.g2].p1 ].pos;
        //Vec3d d;
        //d = p01-p00; qd.p01=qd.p00+d*o.g1span.x;  qd.p00.add_mul( d,o.g1span.y);
        //d = p11-p10; qd.p11=qd.p10+d*o.g2span.x;  qd.p10.add_mul( d,o.g2span.y);
        qd.l1.fromSubLine( nodes[ girders[o.g1].p0 ].pos, nodes[ girders[o.g1].p1 ].pos, o.g1span.x, o.g1span.y );
        qd.l2.fromSubLine( nodes[ girders[o.g2].p0 ].pos, nodes[ girders[o.g2].p1 ].pos, o.g2span.x, o.g2span.y );
	}

    inline double rayPlate( const Plate& plate, const Vec3d& ro, const Vec3d& rd, Vec3d& normal, const Vec3d& hX, const Vec3d& hY )const{
        Quad3d qd;
        //getQuad(qd);
        plate2quad(plate, qd);
        return qd.ray(ro, rd, normal, hX, hY );
	}

    inline double rayLinkLine( const NodeLinker& o, const Vec3d& ro, const Vec3d& rd, double rmax )const{
        Vec3d lp0=nodes[o.p0].pos;
        Vec3d lpd=nodes[o.p1].pos-lp0;
        double lmax = lpd.normalize();
        double t,l;
        double r = rayLine(ro, rd, lp0, lpd, t, l );
        //printf( "rayLinkLine t,l,r %g %g %g \n", t,l,r );
        if( (r>rmax) || (l>lmax) || (l<0) ) return t_inf;
        return t;
	}

int add_Node( const Vec3d& pos ){ 
    Node o;
    o.pos = pos;
    o.id  = nodes.size();
    nodes.push_back( o );
    //if(bPrint)printf( "Node (%g,%g,%g)  ->  %i\n",  pos.x, pos.y, pos.z, id );
    if(bPrint) o.print();
    return o.id;
};

int add_Rope( int p0, int p1, double thick, int matId=-1, const char* matn=0 ){
    Rope o;
    o.p0   = p0;
    o.p1   = p1;
    o.thick= thick;
    if( (matId<0) && matn){ o.face_mat = workshop->materials.getId(matn); }else{ o.face_mat = matId; }
    //if( (matId<0) && matn){ o.material = workshop->materials.get(matn); }else{ o.material = workshop->materials.vec[matId]; }
    o.id   = ropes.size();
    //if(bPrint)printf( "Rope (%i,%i) %g %s ->  %i\n", o.p0, o.p1, o.thick, matn, o.id );
    if(bPrint) o.print();
    ropes.push_back( o );
    return o.id;
};

int add_Girder  ( int p0, int p1, const Vec3d& up, int nseg, int mseg, const Vec2d& wh, Quat4i stickTypes=Quat4i{1,2,3,4} ){
    Girder o;
    o.p0   = p0;
    o.p1   = p1;
    o.up   = up;
    o.nseg = nseg;
    o.mseg = mseg;
    o.wh   = wh;
    o.st   = stickTypes;
    o.id   = girders.size();
    //if( (matId<0) && matn){ o.material = workshop->materials.get(matn); }else{ o.material = workshop->materials.vec[matId]; }
    //if(bPrint)printf( "Girder (%i,%i) (%g,%g,%g) (%i,%i), (%g,%g) st(%i,%i,%i,%i) ->  %i\n", o.p0, o.p1, o.up.x, o.up.y, o.up.z, o.nseg, o.mseg, o.wh.x, o.wh.y, o.st.x,o.st.y,o.st.z,o.st.w,   o.id );
    if(bPrint) o.print();
    girders.push_back( o );
    return o.id;
};

int add_Ring( Vec3d pos, Vec3d ax, Vec3d up, double R, int nseg, const Vec2d& wh, Quat4i stickTypes=Quat4i{1,2,3,4} ){
    Ring o;
    o.pose.pos = pos;
    o.pose.rot.fromDirUp(ax,up);
    o.R     = R;
    o.nseg  = nseg;
    o.wh    = wh;
    o.st    = stickTypes;
    o.id     = rings.size();
    //if( (matId<0) && matn){ o.material = workshop->materials.get(matn); }else{ o.material = workshop->materials.vec[matId]; }
    /*
    if(bPrint)printf( "Ring (%i,%i) (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) (%i), %g, (%g,%g) st(%i,%i,%i,%i) ->  %i\n", o.nseg,
        o.pose.pos.x, o.pose.pos.y, o.pose.pos.z,
        o.pose.rot.a.x, o.pose.rot.a.y, o.pose.rot.a.z,
        o.pose.rot.b.x, o.pose.rot.b.y, o.pose.rot.b.z,
        o.pose.rot.c.x, o.pose.rot.c.y, o.pose.rot.c.z,
        o.R, o.wh.x, o.wh.y, o.st.x,o.st.y,o.st.z,o.st.w, o.id );
    */
    if(bPrint) o.print();
    rings.push_back( o );
    return o.id;
};

int add_Radiator ( int g1, const Vec2d& g1span, int g2, const Vec2d& g2span, int face_mat=-1, double temperature=0 ){
    Radiator o;
    o.g1          = g1;
    o.g2          = g2;
    o.g1span      = g1span;
    o.g2span      = g2span;
    o.temperature = temperature;
    //if( (matId<0) && matn){ o.material = workshop->materials.get(matn); }else{ o.material = workshop->materials.vec[matId]; }
    o.face_mat = face_mat;
    o.id   = radiators.size();
    //if(bPrint)printf( "Radiator %i(%.2f,%.2f) %i(%.2f,%.2f) %f %i -> %i\n", o.g1, o.g1span.x, o.g1span.y,   o.g2, o.g2span.x, o.g2span.y,  o.temperature, o.face_mat, o.id );
    if(bPrint) o.print();
    radiators.push_back( o );
    return o.id;
};

int add_Shield ( int g1, const Vec2d& g1span, int g2, const Vec2d& g2span, double face_mat=-1 ){
    Shield o;
    o.g1          = g1;
    o.g2          = g2;
    o.g1span      = g1span;
    o.g2span      = g2span;
    o.face_mat = face_mat;
    o.id   = shields.size();
    //printf( "Radiator %i(%.2f,%.2f) %i(%.2f,%.2f) %f -> %i\n", o.g1, o.g1span.x, o.g1span.y,   o.g2, o.g2span.x, o.g2span.y, o.face_mat, o.id );
    if(bPrint) o.print();
    shields.push_back( o );
    //lua_pushnumber(L, o.id);
    return o.id;
};

int add_Tank ( Vec3d pos, Vec3d dir, Vec3d span, int face_mat=-1, int fill=-1, double full=0.0 ){
    Tank o;
    o.pose.pos   = pos;
    o.pose.rot.c.normalize();
    o.pose.rot.c.getSomeOrtho( o.pose.rot.b, o.pose.rot.a );
    o.span       = span;
    o.face_mat     = face_mat;
    o.commodityId  = fill;
    o.filled       = full;
    o.id   = tanks.size();
    if(bPrint)o.print();
    //printf( "Tank pos (%f,%f,%f) dir (%f,%f,%f) l %f r %f -> %i\n",  o.pose.pos.x,o.pose.pos.x,o.pose.pos.x,  o.pose.rot.c.x, o.pose.rot.c.y, o.pose.rot.c.z,   o.span.b, o.span.c, o.id );
    tanks.push_back( o );
    //lua_pushnumber(L, o.id);
    return o.id;
};

int add_Thruster ( Vec3d pos, Vec3d dir, Vec3d span, int type=-1 ){
    Thruster o;
    o. pose.pos   = pos;
    o. pose.rot.c = dir;
    o. pose.rot.c.normalize();
    o. pose.rot.c.getSomeOrtho( o.pose.rot.b, o.pose.rot.a );
    o. span       = span;
    o.type        = type;
    o.id   = thrusters.size();
    //printf( "Thruster pos (%f,%f,%f) dir (%f,%f,%f) l %f r %f -> %i\n",  o.pose.pos.x,o.pose.pos.x,o.pose.pos.x,  o.pose.rot.c.x, o.pose.rot.c.y, o.pose.rot.c.z,   o.span.b, o.span.c, o.id );
    if(bPrint) o.print();
    thrusters.push_back( o );
    //lua_pushnumber(L, o.id);
    return o.id;
};

int add_Gun     ( int suppId, const Vec2d& suppSpan, int type=-1 ){
    Gun o;
    o.suppId   = suppId;
    o.suppSpan = suppSpan;
    o.gunId    = type;
    o.id   = guns.size();
    //printf( "Gun %i(%.2f,%.2f) %s -> %i\n", o.suppId, o.suppSpan.x, o.suppSpan.y,   skind, o.id );
    if(bPrint) o.print();
    guns.push_back( o );
    //lua_pushnumber(L, o.id);
    return 1;
};
    

/**
 * Picks the closest component intersected by a ray in the scene. Used e.g. for mouse picking.
 * 
 * @param ro The origin of the ray.
 * @param rd The direction of the ray.
 * @param rmax The maximum distance to consider for picking (default: 1.0).
 * @return The index of the picked component, or -1 if no component is picked.
 */
	int pick(Vec3d ro, Vec3d rd, double rmax = 1.0 ){
        Quad3d qd;
        Vec3d  nr;
        Vec3d hX,hY;
        rd.getSomeOrtho(hX,hY);
        int imin=-1;
        pickedTyp=-1;
        double tmin = t_inf;
        for(int i=0; i<radiators.size(); i++){
            //plate2quad( plate, qd );
            double t = rayPlate( radiators[i], ro, rd, nr, hX, hY );
            //printf( "pick.radiators[%i] t %g tmin %g \n", i, t, tmin );
            if(t<tmin){ tmin=t; imin=i; pickedTyp=(int)ComponetKind::Radiator; }
        };
        for(int i=0; i<shields.size(); i++){
            double t = rayPlate( shields[i], ro, rd, nr, hX, hY );
            //printf( "pick.shields[%i] t %g tmin %g \n", i, t, tmin );
            if(t<tmin){ tmin=t; imin=i; pickedTyp=(int)ComponetKind::Shield; }
        }
        //double rmax = 5.0;
        for(int i=0; i<girders.size(); i++){
            double t = rayLinkLine( girders[i], ro, rd, rmax );
            //printf( "pick.girders[%i] t %g tmin %g \n", i, t, tmin );
            if(t<tmin){ tmin=t; imin=i; pickedTyp=(int)ComponetKind::Girder; }
        }
        for(int i=0; i<ropes.size(); i++){
            double t = rayLinkLine( ropes[i], ro, rd, rmax );
            //printf( "pick.girders[%i] t %g tmin %g \n", i, t, tmin );
            if(t<tmin){ tmin=t; imin=i; pickedTyp=(int)ComponetKind::Rope; }
        }
        //return tmin<0.9*t_inf;
        return imin;
	}

/**
 * Returns picked ship component
 * 
 * @param picked The index of the picked component.
 * @return The picked component. The type of the component is determined by the internal variable pickedTyp .
 */
	ShipComponent* getPicked( int picked ){
	    switch( (ComponetKind)pickedTyp){
            case ComponetKind::Radiator : return &radiators[picked]; break;
            case ComponetKind::Shield   : return &shields  [picked]; break;
            case ComponetKind::Girder   : return &girders  [picked]; break;
            case ComponetKind::Rope     : return &ropes    [picked]; break;
        }
        return 0;
    }

/**
 * Converts the SpaceCraft object to a Truss object. Nodes, ropes, girders and rings are converted to points and edges of the truss. Each ship component is assigned a range of points and edges in the truss.
 * 
 * @param truss The Truss object to populate with the converted data.
 */
	void toTruss(Truss& truss, bool bTanks=false )const {
        printf( "SpaceCraft::toTruss() \n" );
        int i=0;
        truss.newBlock();
        int ip0 = truss.points.size();
        for(Node o: nodes){
            truss.points.push_back( o.pos );
        }
        for(Rope o: ropes){
            truss.edges.push_back( (TrussEdge){o.p0,o.p1,0} );
        }
        for(Girder o: girders){
            //printf("DEBUG toTruss : girder #%i \n", i);
            truss.girder1( nodes[o.p0].pos, nodes[o.p1].pos, o.up, o.nseg, o.wh.a );
            truss.girder1_caps( o.p0, o.p1, 1 );

            Vec2i& bak = truss.blocks.back();
            o.poitRange  = {bak.x,truss.points.size()};
            o.stickRange = {bak.y,truss.edges .size()};
            i++;
        }
        i=0;
        for(Ring o: rings){
            //printf("DEBUG toTruss : ring #%i  %f   %f \n", i, o.nseg, o.wh.a );
            truss.wheel( o.pose.pos, o.pose.pos+o.pose.rot.b*o.R, o.pose.rot.c, o.nseg, o.wh.a );
            Vec2i& bak = truss.blocks.back();
            o.poitRange  = {bak.x,truss.points.size()};
            o.stickRange = {bak.y,truss.edges .size()};
            i++;
        }
        // ToDo: Shields
        // ToDo: Radiators
        // ToDo: Thrusters
        // ToDo: Tanks
        // ToDo: maybe Tanks would be better simulated as rigid-body objects, but than we need to couple the Truss (SoftBody) and RigidBody ?
        if( bTanks ){
            for(Tank o: tanks){
                //printf("DEBUG toTruss : tank #%i \n", i);
                Vec3d p0 = o.pose.pos + o.pose.rot.a*o.span.a*0.5;
                Vec3d p1 = o.pose.pos - o.pose.rot.a*o.span.a*0.5;
                int edgeMat = workshop->panelMaterials.vec[ o.face_mat ]->stickMaterialId;
                truss.makeCylinder( p0, p1, o.span.b, o.span.b, -1, -1, 1.0, o.face_mat, o.face_mat );
                Vec2i& bak = truss.blocks.back();
                o.poitRange  = {bak.x,truss.points.size()};
                o.stickRange = {bak.y,truss.edges .size()};
                i++;
            }
        }
        printf( "npoint %i nstick %i nblocks %i \n", truss.points.size(), truss.edges.size(), truss.blocks.size()  );
	}
	void toTruss(){ truss.clear(); toTruss(truss); };



    inline void plate2raytracer( const Plate& o, TriangleRayTracer& raytracer, double elemMaxSize, bool active )const{
        Quad3d qd;
        plate2quad(o, qd);
        raytracer.addTriangle( {qd.l1.a,qd.l1.b,qd.l2.b}, elemMaxSize, active );
        raytracer.addTriangle( {qd.l1.a,qd.l2.b,qd.l2.a}, elemMaxSize, active );
        //return fmin(t1,t2);
	}

	void toRayTracer( TriangleRayTracer& raytracer, double elemMaxSize )const{
        for(int i=0; i<radiators.size(); i++){ plate2raytracer( radiators[i], raytracer, elemMaxSize, true ); }
        for(int i=0; i<shields.size();   i++){ plate2raytracer( shields  [i], raytracer, elemMaxSize, true ); }
	}

    void updatePanelMaterials(){
        for (PanelMaterial* pm : workshop->panelMaterials.vec ){     
            double areaDensity = 0.0;
            for( PanelLayer& l : pm->layers ){ 
                areaDensity += l.materiallId = workshop->materials.vec[l.materiallId]->density * l.thickness; // [kg/m2]
            }   
            pm->areaDensity = areaDensity;
        }
    }

    void updateTrussEdges( Truss& truss ){
        for( TrussEdge& e : truss.edges ){
            Material* mat = workshop->materials.vec[ e.type ];
            double c = e.crossection / e.l0;
            e.kT = mat->Kpull * c;
            e.kP = mat->Kpush * c;
            e.mass = mat->density * e.crossection * e.l0;
        }
    }

    void updateTrussFaces( Truss& truss ){
        for( TrussFace& f : truss.faces ){
            PanelMaterial* mat = workshop->panelMaterials.vec[ f.type ];
            f.mass = mat->areaDensity * f.area;
        }
    }

};

} // namespace SpaceCrafting

#endif
