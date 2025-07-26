
#ifndef  Selection_h
#define  Selection_h


#include <unordered_set>
#include <unordered_map>
#include <vector>

#include <functional>
#include <ranges>


#include "SDfuncs.h"



// template<typename Builder, typename DistFunc>
// int selectRectEdge( const Vec3d& p0, const Vec3d& p1, const Mat3d& rot, Builder& builder, DistFunc& distFunc ){ 
//     //printf( "Mesh::Builder2::selectRectEdge() p0(%g,%g,%g) p1(%g,%g,%g) \n", p0.x,p0.y,p0.z, p1.x,p1.y,p1.z );
//     Vec3d Tp0,Tp1;
//     //Mat3d rot = (Mat3d)cam.rot;
//     rot.dot_to(p0,Tp0);
//     rot.dot_to(p1,Tp1);
//     _order(Tp0.x,Tp1.x);
//     _order(Tp0.y,Tp1.y);
//     Tp0.z=-1e+300;
//     Tp1.z=+1e+300;
//     int nfound=0;
//     for(int i=0; i<edges.size(); i++ ){
//         Vec2i b = edges[i].lo;
//         Vec3d pa,pb;
//         rot.dot_to( verts[b.i].pos,pa);
//         rot.dot_to( verts[b.j].pos,pb);
//         if( pa.isBetween(Tp0,Tp1) && pb.isBetween(Tp0,Tp1) ){
//             curSelection->add(i);
//             nfound++;
//         }
//     }
//     return nfound;
// }

template<typename MESH, typename DistFunc>
int _selectRectVerts( MESH& mesh, const DistFunc& distFunc, double Rmax=1.0){ 
    double rmin   = 1e+300; 
    int    nfound = 0;
    int nverts = mesh._number_of_points();
    for(int i=0; i<nverts; i++ ){
        Vec3d p = mesh._get_point(i);
        double r = distFunc(p);
        if( r<rmin ) rmin = r;
        if( r<Rmax ){
            mesh._add_to_selection(i);
            nfound++;
        }
    }
    return nfound;
}


inline int make_unique(int n, int* ids){
    std::unordered_set<int> set; set.reserve(n);
    std::vector<int> ids_;       ids_.reserve(n);
    for(int i=0; i<n; i++){
        int id = ids[i];
        if( set.insert(id).second ){ ids_.push_back(id); } // if inserted successfully
    }
    for(int i=0;i<ids_.size();i++){ ids[i]=ids_[i]; }
    return ids_.size();
};

class Selection{public:
    int kind = -1;                   // kind of items selected 
    bool bUnique = true;
    std::vector<int>            vec; // this is ordered ( typically by order of insertion )
    std::unordered_map<int,int> map; // this is not ordered, but fast to check membership, ensure uniquness
    
    inline bool add     (int i)     { 
        if(bUnique){
            if(map.find(i)==map.end()){ map[i]=i; vec.push_back(i); return true; }
        }else{
            vec.push_back(i); return true; 
        }   
        return false; 
    };
    
    inline bool remove  (int i)     { 
        auto it = map.find(i); 
        if(it!=map.end()){
            int ind = it->second;
            //vec.erase(std::remove(vec.begin(), vec.end(), i), vec.end());  // this is costly
            vec[ind]=-1; // mark as removed
            map.erase(i); 
            return true;
        }
        return false;
    }

    inline int toggle( int i ){
        auto it = map.find(i); 
        if(it!=map.end()){ 
            int ind = it->second;
            vec[ind]=-1; // mark as removed
            map.erase(i); 
            return -1; 
        }else{
            map[i]=i; 
            vec.push_back(i); 
            return 1;
        }; 
        return 0;
    }

    inline bool contains(int i)const{ return map.count(i); }
    inline int  find    (int i)const{ auto it = map.find(i); if(it==map.end()) return -1; return it->second; }
    inline int  size()const{ return vec.size(); }
    inline int* data(){ return vec.data(); }

    inline void clear(){ map.clear(); vec.clear(); }
    //inline void unionWith(Selection& other        ){ for(int i : other.vec){ add(i); } }
    // inline void unionWith(std::co<int>& other ){ for(int i : other){ add(i); } }
    // inline void unionWith(std::m<int>& other ){ for(int i : other){ add(i); } }
    // inline void unionWith(std::vector<int>& other ){ for(int i : other){ add(i); } }

    template <typename Container> 
    inline void insert( const Container& other ){ for(int i : other){ add(i); } }
    inline void insert( int n, const int* is   ){ for(int i=0; i<n; i++){ add(is[i]); } }

    inline void invert( int nmax, int n=-1, const int* is=0 ){
        bool* mask = new bool[nmax];
        for(int i=0; i<nmax; i++){ mask[i]=true; }
        if(n<0){ n=vec.size(); is=vec.data(); }
        for(int i=0; i<n; i++){ mask[is[i]]=false; }
        clear();
        for(int i=0; i<nmax; i++){ if(mask[i]) add(i); }
        delete[] mask;
    }

    //inline void substract_remove    (Selection& other){ for(int i : vec){ if(other.contains(i)){ remove(i); } } }
    template <typename Container> 
    inline int substract(const Container& other, bool bInv=false){ 
        std::vector<int> new_vec;
        int nerased = 0;
        for(int i : vec){ 
            if( other.contains(i) != bInv ) { map.erase(i); nerased++; }
            else                            { new_vec.push_back(i);    } 
        }
        vec = new_vec;
        return nerased;
    }
    template <typename Container> 
    inline int intersectWith(const Container& other){ return substract(other, true); }


// General selection algorithm 
// see https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221OJl0tLGMMu_qhDAnOw4xvWJNMmgrFHmY%22%5D,%22action%22:%22open%22,%22userId%22:%22100958146796876347936%22,%22resourceKeys%22:%7B%7D%7D&usp=sharing    
template<typename Predicate>
int selectByPredicate( std::ranges::input_range auto&& range,  Predicate predicate ) {
    // Optional: A concept check to ensure the predicate is valid for our range.
    // This gives a much clearer error message if you pass the wrong kind of function.
    static_assert(std::is_invocable_r_v<bool, Predicate, const std::ranges::range_value_t<decltype(range)>&>);
    int index = 0;
    for (const auto& element : range) {
        if (predicate(element)) { 
            printf("Selected %d\n", index); 
            this->add(index); 
        }
        index++;
    }
    return vec.size(); // Return the number of items selected
}

};

class SelectionBanks{public:
    int icurSelection       = 0;
    Selection* curSelection = 0;
    std::vector<Selection> selections;
    SelectionBanks(int n){ selections.resize(n); pickSelection(0); }
    Selection* pickSelection(int i){ curSelection=&selections[i]; icurSelection=i; return curSelection; }
    Selection* nextSelection(){ icurSelection++; if(icurSelection>=selections.size()) icurSelection=0; return pickSelection(icurSelection); }
};

#endif
