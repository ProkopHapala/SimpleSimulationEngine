
#ifndef  SphereTreeND_2_h
#define  SphereTreeND_2_h

#include <vector>

//#include "VecN.h"

inline double dist2limited( int n, double * xs, double * x0s, double dist2max ){
    double d2sum = 0;
    for(int i=0; i<n; i++){ double d = d-x0s[i]; d*=d; d2sum+=d; }
    //for(int i=0; i<n; i++){ double d=center[i]-d; d*=d; d2sum+=d;  if(d2sum>dist2max) return 1e+300; }
    return d2sum;
}

class SphereNodeND{
    public:
    //int      id;
    int              parrent;
    double *         center;
    std::vector<int> branches;

    SphereNodeND(){ };
    SphereNodeND(double* center_, int parrent_){  center=center_; parrent=parrent_; };

    inline void searchOp( int i, int n, double* xs, double r2max,     double& rmin, int& imin, int& nfound, int* found ){
        double  r2  = dist2limited( n, xs, center, r2max );
        if(r2>r2max) return;
        if(found)found[nfound]=i;
        nfound++;
        if(r2<rmin ){ rmin=r2; imin = i; }
    }


/*
    int search( double r2max, int n, double xs, double** data, int& closest, int * results ){
        int nfound = 0;
        double rmin = 1e+300;
        //closest = -1;
        for( int i : branches ){
            double   r2  = dist2limited( n, xs, data[i], r2max );
            if(r2>r2max) continue;
            if(results )results[nfound]=i;
            nfound++;
            if(r2<rmin ){ rmin=r2; closest = i; }
        }
        return nfound;
    };

    int search( double r2max, int n, double xs, SphereNodeND* data, int& imin, double& rmin, int * found ){
        int nfound = 0;
        double rmin = 1e+300;
        for( int i : branches ){
            double   r2  = dist2limited( n, xs, data[i].center, r2max );
            if(r2>r2max) continue;
            if(found)found[nfound]=i;
            nfound++;
            if(r2<rmin ){ rmin=r2; closest = i; }
        }
        return nfound;
    };

    */

};

/*
int search( double r2max, int n, double xs, int nfrom, int* from, SphereNodeND* data, int& imin, double& rmin, int * found ){
    int nfound = 0;
    double rmin = 1e+300;
    for( int i : branches ){
    for( int i = 0; i< nfrom; i++){

        double r2  = dist2limited( n, xs, data[i].center, r2max );
        if(r2>r2max) continue;
        if(found)found[nfound]=i;
        nfound++;
        if(r2<rmin ){ rmin=r2; closest = i; }
    }
    return nfound;
};
*/




class SphereTreeND{
    public:
    int nDim=0;
    int nLevel=0;
    double * RIs = NULL; // inner radius - child can be attached under this node ( but is not enclosed by it )
    double * ROs = NULL; // outer radius - point can interact with some child of this node ( but not necessarily inside the node RI )
    std::vector<SphereNodeND>* nodes = NULL;
    std::vector<double*>       leafs;

    int nMaxFound=0;
    //int nsup=0;
    int *sups=NULL,*subs=NULL;



    /*
    int findNeighAtLevel( int n, double * xs, int& level, int& iclosest ){
        for(int ilevel=0;ilevel<level;ilevel++){
            int nsub=0;
            int closest,
            for(int isup=0; isup<nsup; isup++ ){
                nsub += nodes[ilevel][isup]->findChildsCloserThan( ROs[ilevel], n, xs, &nodes[ilevel+1][0], closest, sups+nnew );
            }
            if( nsub == 0 ){
                return;
            }
            int* tmp = sups; sups = subs; sups = tmp;
        }
    }
    */

    /*
    int findNeighAtLevel( int n, double * xs, int& level, int& closest_ ){
        int ilevel=0,nsub=0;
        double rmin  = 1.0e+300d;
        double r2max = ROs[ilevel+1];
        std::vector<SphereNodeND>& nodes_sup = nodes[0];
        //std::vector<SphereNodeND>& nodes_sub = nodes[ilevel+1];
        int nsup = 0;
        for(int i=0; i<nodes_sup.size(); i++){
            double   r2  = dist2limited( n, xs, nodes_sub[i].center, r2max );
            if(r2>r2max) continue;
            if(results )results[nfound]=i;
            nsub++;
            if(r2<rmin ){ rmin=r2; closest = i; }
        }
        for(ilevel=0;ilevel<level;ilevel++){
            nsub        = 0;
            int closest = -1;
            std::vector<SphereNodeND>& nodes_sup = nodes[ilevel];
            std::vector<SphereNodeND>& nodes_sub = nodes[ilevel+1];
            for(int j=0; j<nsup; j++ ){
                int isup = sups[j];
                for( int i : nodes_sup[isup].branches ){
                    double   r2  = dist2limited( n, xs, nodes_sub[i].center, r2max );
                    if(r2>r2max) continue;
                    if(results )results[nfound]=i;
                    nsub++;
                    if(r2<rmin ){ rmin=r2; closest = i; }
                }
            }
            if (nsub==0) break;
            int* tmp = sups; sups = subs; sups = tmp; nsup=nsub;
            closest_ = closest;
        }
        level=ilevel;
        return nsub;
    }
    */

    int findNeighAtLevel( int n, double * xs, int& level, int& closest ){
        int ilevel=0;
        // ------- top level
        int nsup = 0; int imin  = -1; double rmin  = 1.0e+300d; double r2max = ROs[0];
        std::vector<SphereNodeND>& nodes_sup = nodes[0];
        for(int i=0; i<nodes_sup.size(); i++){
            nodes_sup[i].searchOp( i, n, xs, r2max, rmin, imin, nsup, sups );  // sups are level 1
        }
        closest = imin;
        if (nsup==0) return 0;
        // ------- lover levels
        for(ilevel=1;ilevel<level;ilevel++){
            int nsub  =  0; imin  = -1;  rmin  =  1.0e+300d; r2max = ROs[ilevel];
            std::vector<SphereNodeND>& nodes_sup = nodes[ilevel];
            std::vector<SphereNodeND>& nodes_sub = nodes[ilevel+1];
            for(int j=0; j<nsup; j++ ){
                for( int i : nodes_sup[sups[j]].branches ){
                    nodes_sub[i].searchOp( i, n, xs, r2max, rmin, imin, nsub, subs );
                }
            }
            if (nsub==0) break;
            int* tmp = sups; sups = subs; sups = tmp; nsup=nsub;
            closest = imin;
        }
        level=ilevel;
        return nsup;
    }

/*
    int findCloserThan( double r2max, int n, double * xs, double** results ){
        nsup=nodes.size(); for(int i=0; i<nsup; i++){ sups[i]=i; }
        int nsup = findNeighAtLevel( nLevel-1, n, xs );
        int nsub=0;
        for(int isup=0; isup<nsup; isup++ ){
            nsub += node[nLevel-1][sup]->findChildsCloserThan( ROs[ilevel], n, xs, &leafs[isup][0], sups+nnew );
        }
        for(int i=0; i<nsub; i++){ results[i] = leafs[sups[i]]; }
    }
*/

    int insert( int n, double * xs, int level, int isup ){
        SphereNodeND * sup;
        if(level==0){
            sup = new SphereNodeND( xs, isup );
            isup = nodes[0].size();
            nodes[0].push_back( *sup );
            level++;
        }else{
            sup = &nodes[level][isup];
        }
        for(int i=level; i<(nLevel-1); i++ ){
            SphereNodeND * node = new SphereNodeND( xs, isup );
            isup = nodes[level].size();
            sup->branches.push_back(isup);
            nodes[level] .push_back( *node );
            sup = node;
        }
        sup->branches.push_back( leafs.size() );
        leafs.push_back(xs);
    }

    int insert( int n, double * xs ){
        leafs.push_back( xs );
        int closest; int level = nLevel-1;
        int nsup = findNeighAtLevel( n, xs, level, closest );
        insert( n, xs, level, closest );
    }

    void init( int nLevel_, int nDim_, int nMaxFound_, double * RIs_ ){
        nDim  = nDim_; nLevel =  nLevel_;   RIs = RIs_; nMaxFound = nMaxFound_;
        nodes = new std::vector<SphereNodeND>[nLevel];
        //leafs = new std::vector<double*>     ;
        ROs   = new double[nLevel];
        subs  = new int[nMaxFound];
        sups  = new int[nMaxFound];
        double ro=0.0;
        for( int i=nLevel-1; i>=0; i--){
            ro    += RIs[i];
            ROs[i] = ro;
        }
    }

};


#endif

