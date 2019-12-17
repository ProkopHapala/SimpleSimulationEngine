
#ifndef  DistanceHierarchy_h
#define  DistanceHierarchy_h

#include <vector>

template <typename T>
class DHSearch{ public:
    double Rmax;
    int    nmax;
    T* p;
    T* ps;
    //int    nfound;
    //int*   found;
    inline double dist(int i){ return p.dist(ps[i]); };
};

template <typename T>
class DHNode{
    double Rcut       = 0;
    int pivot         =-1;
    DHNode*     parent= 0;
    std::vector<DHNode*> subs;
    std::vector<int>     leafs;

    DHNode* deepestContainer( const DHSearch<T>& srch ){ // not necesarily nearest
        double r = srch.dist( pivot );
        if( r>(Rcut+srch.Rmax) ) return 0;
        for( DHNode* s : subs ){
            DHNode* o = s.deepestContainer( srch );
            if( o ) return o;
        }
        return this;
    };

    void nearestContainer( const DHSearch<T>& srch, double& R, DHNode* best ){ // not necesarily nearest
        //double r = srch.dist( pivot );
        //if( r>(Rcut+srch.Rmax) ) return 0;
        for( DHNode* s : subs ){
            double r = srch.dist( s.pivot0 );
            if(r<R){ best=s; R=r; }
            if(r<s.Rcut){ nearestContainer( srch, R, best ); }
        }
    };

    //int findCloserThan( const T& p, double Rmax, const T* ps, int nmax, int* found ){
    int findCloserThan( const DHSearch<T>& srch, int& n, int* found ){
        double r = srch.dist( pivot );
        if( r<srch.Rmax ){ found[n]=pivot; n++; };
        if( r>Rcut ) return 0;
        for( int i : leafs ){
            r = srch.dist(i);
            if( r<srch.Rmax ){ found[n]=i; n++; };
        }
        for( DHNode* s : subs ){
            s.findCloserThan( srch, n, found );
        }
        return n;
    }

    //int findCloserThan( const T& p, double Rmax, const T* ps, int nmax, int* found ){
    int findFirstOverlap( const DHSearch<T>& srch ){
        double r = srch.dist( pivot );
        if     ( r<srch.Rmax ){ return pivot; }
        else if( r>Rcut      ){ return -1;    }
        for( int i : leafs ){
            r = srch.dist(i);
            if( r<srch.Rmax ){ return i; };
        }
        for( DHNode* s : subs ){
            int i = s.findFirstOverlap( srch );
            if ( i>0 ) return i;
        }
        return -1;
    }

};

class DHLevel{
    int    n;
    double R;
};

template<typename T>
class DistanceHierarchy : public DHNode<T>, public DHSearch<T> { public:

    std::vector<T> points;

    int nlevel;
    std::vector<DHLevel> levels;

    // === functions

    /*
    DHNode* findCloserThan( const T& p, double Rmax, int nmax, int* found ){
        double Rcut = ;
        for( int i : subs ){
            double r = p.dist( points[i] );
            if( r<Rmax ){

            };
        }
    }
    */

};

#endif

