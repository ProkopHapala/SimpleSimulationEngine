
#ifndef  clustering2D_h
#define  clustering2D_h

#include <vector>

#include "Vec3.h"
#include "geom3D.h"


class Cluster{ public:
    Vec3d center;
    std::vector<int> leafs;
    Cluster(Vec3d center_){ center=center_; }
};

class ClusterMap{ public:
    //std::vector<Rect2d> clusters;
    std::vector<Cluster> clusters;

    bool insert_dist( int o, const Vec3d& p, double R ){
        double R2 = R*R;
        bool isNew = true;
        int j=0;
        for( Cluster& c : clusters ){
            double r2 = p.dist2( c.center );
            //printf( " %i %i %f %f \n", o, j, r2, R2 ); j++;
            if(r2>R2) continue;
            //printf( " HIT \n", o, j, r2, R2 );
            c.leafs.push_back(o);
            return false;
            // insert p somewhere;
        };
        clusters.push_back( Cluster(p) );
        clusters.back().leafs.push_back(o);
        return true;
    }

};



class ClusterBox{ public:
    Box bbox;
    std::vector<int> leafs;
    ClusterBox(Box bbox_){ bbox=bbox_; }

    bool tryInsert( int o, const Vec3d& p, double R ){
        Box newbox=bbox;
        if     ( p.x<bbox.a.x ){ if( (bbox.b.x-p.x)<R ){ newbox.a.x=p.x; }else{ return false; } }
        else if( p.x>bbox.b.x ){ if( (p.x-bbox.a.x)<R ){ newbox.b.x=p.x; }else{ return false; } }
        if     ( p.y<bbox.a.y ){ if( (bbox.b.y-p.y)<R ){ newbox.a.y=p.y; }else{ return false; } }
        else if( p.y>bbox.b.y ){ if( (p.y-bbox.a.y)<R ){ newbox.b.y=p.y; }else{ return false; } }
        if     ( p.z<bbox.a.z ){ if( (bbox.b.z-p.z)<R ){ newbox.a.z=p.z; }else{ return false; } }
        else if( p.z>bbox.b.z ){ if( (p.z-bbox.a.z)<R ){ newbox.b.z=p.z; }else{ return false; } }
        //printf( " (%f,%f) ((%f,%f)(%f,%f)) -> ((%f,%f)(%f,%f))   ", p.x, p.y, bbox.a.x,bbox.a.y,bbox.b.x,bbox.b.y,  newbox.a.x,newbox.a.y,newbox.b.x,newbox.b.y );
        bbox = newbox;
        leafs.push_back(o);
        return true;
    }

};

class ClusterMapBox{ public:
    //std::vector<Rect2d> clusters;
    std::vector<ClusterBox> clusters;

    bool insert( int o, const Vec3d& p, double R, double span0 ){
        double R2 = R*R;
        int j=0;
        for( ClusterBox& c : clusters ){
            if( c.tryInsert(o,p,R) ){ return false; };
        };
        clusters.push_back( ClusterBox( (Box){ p-(Vec3d){span0,span0,span0}, p+(Vec3d){span0,span0,span0} } ) );
        clusters.back().leafs.push_back(o);
        return true;
    }

};


#endif

