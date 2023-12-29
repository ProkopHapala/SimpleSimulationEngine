
#ifndef  SpaceCraft2Truss_h
#define  SpaceCraft2Truss_h

#include <vector>
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

//#include "Noise.h"
//#include "SphereSampling.h"
//#include "DrawSphereMap.h"

#include "SpaceCraft.h"
#include "Truss.h"
#include "TriangleRayTracer.h"
#include "Radiosity.h"

namespace SpaceCrafting{

    void plate2raytracer( const SpaceCraft& craft, const Plate& o, TriangleRayTracer& raytracer, double elemMaxSize, bool active ){
        Quad3d qd;
        craft.plate2quad(o, qd);
        raytracer.addTriangle( {qd.l1.a,qd.l1.b,qd.l2.b}, elemMaxSize, active );
        raytracer.addTriangle( {qd.l1.a,qd.l2.b,qd.l2.a}, elemMaxSize, active );
        //return fmin(t1,t2);
	}

	void toRayTracer( const SpaceCraft& craft, TriangleRayTracer& raytracer, double elemMaxSize ){
        for(int i=0; i<craft.radiators.size(); i++){ plate2raytracer( craft, *craft.radiators[i], raytracer, elemMaxSize, true ); }
        for(int i=0; i<craft.shields.size();   i++){ plate2raytracer( craft, *craft.shields  [i], raytracer, elemMaxSize, true ); }
	}

    void updateTrussEdges( const SpaceCraft& craft, Truss& truss ){
        for( TrussEdge& e : truss.edges ){
            Material* mat = craft.workshop->materials.vec[ e.type ];
            double c = e.crossection / e.l0;
            e.kT = mat->Kpull * c;
            e.kP = mat->Kpush * c;
            e.mass = mat->density * e.crossection * e.l0;
        }
    }

    void updateTrussFaces( const SpaceCraft& craft, Truss& truss ){
        for( TrussFace& f : truss.faces ){
            PanelMaterial* mat = craft.workshop->panelMaterials.vec[ f.type ];
            f.mass = mat->areaDensity * f.area;
        }
    }



/**
 * Converts the SpaceCraft object to a Truss object. Nodes, ropes, girders and rings are converted to points and edges of the truss. Each ship component is assigned a range of points and edges in the truss.
 * 
 * @param truss The Truss object to populate with the converted data.
 */
	void toTruss(SpaceCraft& craft, Truss& truss, bool bTanks=false ){
        printf( "SpaceCraft::toTruss() \n" );
        int i=0;
        truss.newBlock();
        int ip0 = truss.points.size();
        for(Node* o: craft.nodes){
            truss.points.push_back( o->pos );
        }
        for(Rope* o: craft.ropes){
            truss.edges.push_back( (TrussEdge){o->nodes.x->id,o->nodes.y->id,0} );
        }
        for(Girder* o: craft.girders){
            //printf("DEBUG toTruss : girder #%i \n", i);
            truss.girder1( o->nodes.x->pos, o->nodes.y->pos, o->up, o->nseg, o->wh.a );
            truss.girder1_caps( o->nodes.x->id, o->nodes.y->id, 1 );
            Vec2i& bak = truss.blocks.back();
            o->pointRange  = {bak.x,truss.points.size()};
            o->stickRange = {bak.y,truss.edges .size()};
            i++;
        }
        i=0;
        for(Ring* o: craft.rings){
            //printf("DEBUG toTruss : ring #%i  %f   %f \n", i, o->nseg, o->wh.a );
            truss.wheel( o->pose.pos, o->pose.pos+o->pose.rot.b*o->R, o->pose.rot.c, o->nseg, o->wh.a );
            Vec2i& bak = truss.blocks.back();
            o->pointRange  = {bak.x,truss.points.size()};
            o->stickRange = {bak.y,truss.edges .size()};
            i++;
        }
        // ToDo: Shields
        // ToDo: Radiators
        // ToDo: Thrusters
        // ToDo: Tanks
        // ToDo: maybe Tanks would be better simulated as rigid-body objects, but than we need to couple the Truss (SoftBody) and RigidBody ?
        if( bTanks ){
            for(Tank* o: craft.tanks){
                //printf("DEBUG toTruss : tank #%i \n", i);
                Vec3d p0 = o->pose.pos + o->pose.rot.a*o->span.a*0.5;
                Vec3d p1 = o->pose.pos - o->pose.rot.a*o->span.a*0.5;
                int edgeMat = craft.workshop->panelMaterials.vec[ o->face_mat ]->stickMaterialId;
                truss.makeCylinder( p0, p1, o->span.b, o->span.b, -1, -1, 1.0, o->face_mat, o->face_mat );
                Vec2i& bak = truss.blocks.back();
                o->pointRange  = {bak.x,truss.points.size()};
                o->stickRange = {bak.y,truss.edges .size()};
                i++;
            }
        }
        printf( "npoint %i nstick %i nblocks %i \n", truss.points.size(), truss.edges.size(), truss.blocks.size()  );
	}
	//void toTruss(  ){ truss.clear(); toTruss(truss); };

} // namespace SpaceCrafting

#endif
