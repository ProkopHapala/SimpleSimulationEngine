#ifndef PathFinder_h
#define PathFinder_h

#include <vector>
#include "Noise.h"
#include "arrayAlgs.h"
#include "Grid2DAlgs.h"

class Way{ public:
    int a,b; // centers
    std::vector<int>  path;
};

class PathFinder :public Grid2DAlg { public:
    double * height       = NULL;
    double * terrain_cost = NULL;

	//double * water   = NULL;
	double ch2=1.0,chminus=0.0,chplus=0.0;

	int    nContour=0,nContour_old=0;
	int    * contour1 = NULL;
	int    * contour2 = NULL;
	//double * water_    = NULL;
	//std::vector<int>    sinks;
	//std::vector<River*> rivers;

	std::vector<Vec2i>                 centers;
	std::unordered_map<uint64_t,Vec2i> pass;
	std::vector<Way*>                  paths;


	// PathFinder
    //std::vector<int> rivers;
    double* moveCosts = NULL;
    int* toBasin      = NULL;
    int* toTile       = NULL;

    void bind(double* height_, double* terrain_cost_){
        height=height_;
        terrain_cost=terrain_cost_;
    }

    // TODO: Inflow/outflow  +  Rivers
	void allocate(int n=-1){
        if( n<0 ) n=ntot;
		_realloc(contour1,n);
		_realloc(contour2,n);

		_realloc(toBasin  ,ntot);
		_realloc(toTile   ,ntot);
		_realloc(moveCosts,ntot);
    }

    void deallocate(){
        delete [] contour1; delete [] contour2;
        delete [] toBasin; delete toTile; delete moveCosts;
    }

    inline double heightCost(double dh){
        double val = ch2 * dh*dh;
        if( dh>0 ){ val+=chplus*dh; }else{ val-=chminus*dh; }
        return val;
    }

    void pepare(){
        //nContour=0,nContour_old=0;
        for(int i=0; i<ntot; i++){
            moveCosts[i] = 1e+300;
            toTile [i]   = -1;
            toBasin[i]   = -1;
        }
        nContour = centers.size();
        for(int i=0; i<nContour; i++){
            int ii = ip2i( centers[i] );
            contour2[i] = ii;
            //printf( "pepare %i: %i \n", i, ii );
            toBasin  [ii] = i;
            toTile   [ii] = ii;
            moveCosts[ii] = 0.0;
        }
    }

    void path_step(){
        printf( "========== path_step ====== nContour: %i nneigh %i\n", nContour, nneigh );
        nContour_old = nContour;  nContour = 0;
        int * tmp = contour1; contour1 = contour2; contour2 = tmp;
        for ( int ii=0; ii<nContour_old; ii++ ){
            int    i     = contour1[ii];
            Vec2i ip0    = i2ip(i);
            double cost0 = moveCosts[i];
            int ibas     = toBasin  [i];

            //printf( "%i: %i (%i,%i) %i %g\n", ii, i, ip0.x,ip0.y, ibas, cost0 );

            for(int ing=0; ing<nneigh; ing++){
                Vec2i ip = wrap_index( ip0 + neighs[ing] );
                if( validIndex( ip ) ){
                    //printf( "neigh %i: (%i,%i) %g \n", ing, ip0.x, ip0.y, neigh_dist[ing] );
                    extend_path( i, ip2i(ip), ibas, cost0+neigh_dist[ing] );
                };
            }
        }
    }

    void extend_path( int oi, int i, int ibas, double cost0 ){
        double cost  = cost0;
        if( height       ) cost += heightCost( height[i] - height[oi] );
        if( terrain_cost ) cost += terrain_cost[i];
        if( cost < moveCosts[i] ){
            moveCosts[i] = cost;
            toBasin  [i] = ibas;
            toTile   [i] = oi;
            contour2 [nContour] = i;
            nContour++;
        }
    }

    void findConnections(){
        Vec2i ip0;
        for(ip0.y=0; ip0.y<n.y; ip0.y++){
            for(ip0.x=0; ip0.x<n.x; ip0.x++){
                int i0   = ip2i(ip0);
                int ibas = toBasin[i0];
                for(int ing=0; ing<nneigh; ing++){
                    Vec2i dip = neighs[ing];
                    if( ( (dip.x>0)||(dip.y>0) ) ){
                        Vec2i jp = ip0 + dip;
                        if( !validIndex( jp ) ) continue;
                        int j    = ip2i( jp );
                        int jbas = toBasin[j];
                        if(jbas==ibas) continue;
                        //pass.insert( {pass.size(),(Vec2i){i0,j}} );
                        uint64_t id = symetric_id( (Vec2i){ ibas, jbas } );
                        auto found = pass.find(id);
                        if( found != pass.end() ){
                            double cost = moveCosts[i0] + moveCosts[j];
                            cost += heightCost( height[i0] - height[j] );
                            Vec2i ops = found->second;
                            double ocost = moveCosts[ops.x] + moveCosts[ops.y];
                            ocost += heightCost( height[ops.x] - height[ops.y] );
                            if(cost<ocost){
                                //printf( "replace bas(%i,%i)==(%i,%i) %i %i %g %g \n", ibas, jbas, id&0xFFFF, id>>32, i0, j, cost, ocost  );
                                found->second = Vec2i{i0,j};
                            };
                        }else{
                            //printf( "instert  bas(%i,%i)==(%i,%i) %i %i \n", ibas, jbas, id&0xFFFF, id>>32, i0, j  );
                            pass.insert({id,(Vec2i){i0,j}});
                        }
                    }
                }
            }
        }
    }

    void track( Way& way, int i1, int i2){
        std::vector<int> w1;
        int i;
        i = i1;
        w1.push_back(i);
        while(true){
            int i_ = toTile[i];
            if(i==i_) break;
            w1.push_back(i_);
            i=i_;
        };
        for(int i=w1.size()-1; i>=0; i--){
            way.path.push_back( w1[i] );
        }
        i = i2;
        way.path.push_back(i);
        while(true){
            int i_ = toTile[i];
            if(i==i_) break;
            way.path.push_back(i_);
            i=i_;
        };
    }

    void makePaths(){
        for( auto& item : pass ){
            Way* w = new Way();
            track( *w, item.second.x, item.second.y);
            paths.push_back( w );
        }
    }

};

#endif
