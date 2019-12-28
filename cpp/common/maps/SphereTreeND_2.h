
#ifndef  SphereTreeND_2_h
#define  SphereTreeND_2_h

#include <vector>

//#include "VecN.h"

inline double dist2limited( int n, double * xs, double * x0s, double dist2max ){
    printf( "dist2limited (%g,%g) (%g,%g) \n", xs[0],xs[1],   x0s[0],x0s[1] );
    double d2sum = 0;
    for(int i=0; i<n; i++){ double d = xs[i]-x0s[i]; d*=d; d2sum+=d; }
    //for(int i=0; i<n; i++){ double d=center[i]-d; d*=d; d2sum+=d;  if(d2sum>dist2max) return 1e+300; }
    //printf( "dist2limited d2sum %g \n", d2sum );
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
        printf( "searchOp %i %g %g %g \n", i, r2, r2max, rmin  );
        if(r2>r2max) return;
        if(found)found[nfound]=i;
        nfound++;
        if(r2<rmin ){ rmin=r2; imin = i; }
        //printf( "--\n" );
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
    double * RIs  = NULL; // inner radius - child can be attached under this node ( but is not enclosed by it )
    double * ROs  = NULL; // outer radius - point can interact with some child of this node ( but not necessarily inside the node RI )
    double * R2Is = NULL;
    double * R2Os = NULL;
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
        printf( " ==== START findNeighAtLevel === \n" );
        int insert_level=  0;
        closest         = -1;
        // ------- top level
        printf(" --- level %i \n", 0 );
        int nsup = 0; int imin  = -1; double rmin  = 1.0e+300d; double r2max = R2Os[0];
        std::vector<SphereNodeND>& nodes_sup = nodes[0];
        printf(" nodes_sup.size() %i\n", nodes_sup.size() );
        for(int i=0; i<nodes_sup.size(); i++){
            nodes_sup[i].searchOp( i, n, xs, r2max, rmin, imin, nsup, sups );  // sups are level 0
        }
        printf(" findNeighAtLevel level 0 : nsup %i rmin %g \n", rmin, nsup );
        if( rmin < R2Is[0] ){ closest = imin; insert_level=1; }
        if (nsup>0){
            // ------- lover levels
            for(int ilevel=1;ilevel<level;ilevel++){
                printf(" --- level %i \n", ilevel );
                int nsub  =  0; imin  = -1;  rmin  =  1.0e+300d; r2max = R2Os[ilevel];
                std::vector<SphereNodeND>& nodes_sup = nodes[ilevel-1];
                std::vector<SphereNodeND>& nodes_sub = nodes[ilevel  ];
                for(int j=0; j<nsup; j++ ){
                    printf( "j %i nsup %i \n", j, nsup );
                    printf( "j %i nsup %i sups[j] %i \n", j, nsup, sups[j] );
                    printf( "j %i nsup %i sups[j] %i nbranch %i \n", j, nsup, sups[j], nodes_sup[sups[j]].branches.size() );
                    for( int i : nodes_sup[sups[j]].branches ){
                        nodes_sub[i].searchOp( i, n, xs, r2max, rmin, imin, nsub, subs );
                    }
                }
                printf(" nsub %i \n", nsub );
                if (nsub==0) break;
                int* tmp = sups; sups = subs; subs = tmp; nsup=nsub;
                if( rmin < R2Is[ilevel] ){ closest = imin; insert_level=ilevel+1; }
                printf(" --- END level %i \n", ilevel );
            }
        }
        printf( "findNeighAtLevel level %i closest %i\n", insert_level, closest );
        level=insert_level;
        printf( " ==== DONE findNeighAtLevel === \n" );
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

    void insert( int n, double * xs, int level, int isup ){
        //SphereNodeND * sup;
        printf( " insert level %i isup %i \n", level, isup );
        printf( "leafs.size() %i \n", leafs.size() );
        if(level==0){
            printf( "level %i isup %i \n", level, isup );
            SphereNodeND * node = new SphereNodeND( xs, isup );
            isup = nodes[0].size();
            nodes[0].push_back( *node );
            level++;
        }else{
            //sup = &nodes[level][isup];
        }
        for(level; level<nLevel; level++ ){
            printf( "--- level %i isup %i \n", level, isup );
            SphereNodeND * node = new SphereNodeND( xs, isup );
            int isup_ = nodes[level].size();
            printf( "isup_ %i \n", isup_  );
            nodes[level-1][isup].branches.push_back(isup_); //printf( "sup->branches.size() %i \n", sup->branches.size() );
            nodes[level  ].push_back( *node );
            isup = isup_;
        }
        printf( "--- last level %i \n", nLevel-1 );
        printf( "leafs.size() %i \n", leafs.size() );
        nodes[nLevel-1][isup].branches.push_back( leafs.size() );
        leafs.push_back(xs);
        printf( "leafs.size() %i \n", leafs.size() );
        printf( " ==== DONE insert === \n" );
        return ;
    }

    void insert( int n, double * xs ){
        int closest; int level = nLevel;
        int nsup = findNeighAtLevel( n, xs, level, closest );
        insert( n, xs, level, closest );
    }

    void init( int nLevel_, int nDim_, int nMaxFound_, double * RIs_ ){
        nDim  = nDim_; nLevel =  nLevel_;   RIs = RIs_; nMaxFound = nMaxFound_;
        nodes = new std::vector<SphereNodeND>[nLevel];
        //leafs = new std::vector<double*>     ;
        ROs   = new double[nLevel];
        R2Is  = new double[nLevel];
        R2Os  = new double[nLevel];
        subs  = new int[nMaxFound];
        sups  = new int[nMaxFound];
        double ro=0.0;
        for( int i=nLevel-1; i>=0; i--){
            double ri = RIs[i];
            R2Is[i]   = ri*ri;
            ro       += ri;
            ROs[i]    = ro;
            R2Os[i]   = ro*ro;
        }
    }

};


#endif

