
#ifndef  containers_h
#define  containers_h

#include <string>
#include <vector>
#include <unordered_map>

typedef std::unordered_map<int,int>         IntMap;
typedef std::unordered_map<std::string,int> Dictionary;

// === Convenience functions for std::unordered_map

template<typename K,typename V>
V get( const std::unordered_map<K,V>& map, const K& key, const V& def ){
    auto it = map.find(key);
    if( it != map.end() ){ return it->second; }else{ return def; }
}

template<typename K,typename V>
bool set( std::unordered_map<K,V>& map, const K& key, const V& val, bool bReplace=false ){
    auto it = map.find(key);
    if( it != map.end() ){ if(bReplace){it->second = val;} return true; }else{ map.insert({key,val}); return false; }
}

template<typename K,typename V>
V* setp( std::unordered_map<K,V*>& map, const K& key, const V* val, bool bReplace=false ){
    auto it = map.find(key);
    if( it != map.end() ){ V* ret=it->second; if(bReplace){it->second=val;} return ret; }else{ map.insert({key,val}); return 0; }
}

// ================== Ditionary
template<typename T>
class Dict{ public:
    std::unordered_map<std::string,int> map;
    std::vector<T*>                     vec;

    int getId( const char* name )const{
        auto it = map.find(name);
        if( it != map.end() ){ return it->second; }else{ return -1; }
    }

    T* get( const char* name )const{
        auto it = map.find(name);
        if( it != map.end() ){ return vec[it->second]; }else{ return 0; }
    }

    int add( T* mat, bool bDel=true ){
        auto it = map.find( mat->name );
        if( it == map.end() ){ int i=vec.size(); vec.push_back(mat); map.insert({mat->name,i}); return i; }else{ if(bDel)delete vec[it->second]; vec[it->second] = mat; return -1; }
    }
};


#endif
