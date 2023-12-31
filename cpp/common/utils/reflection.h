
#ifndef  reflection_h
#define  reflection_h

#include <string>
#include <vector>
#include <unordered_map>

std::vector<std::string>            ClassNames;
std::unordered_map<std::string,int> ClassDict;

int tryRegisterclass(const std::string& name){
    int ret = -1;
    auto it = ClassDict.find(name);
    if( it == ClassDict.end() ){
        ret = ClassNames.size();
        printf( "RegisterClass(%s) as id=%i \n", name.c_str(), ret );
        ClassNames.push_back(name);
        ClassDict.insert( {name,ret} );
    }else{
        ret = it->second;
    }
    return ret;
}

bool isThisClass( int instance_class_id_init, const std::string& name ){
    int id = instance_class_id_init;
    if( (id<0)||(id>=ClassNames.size()) ){ return false; }
    return (ClassNames[id] == name);
}

#define CLASS_REFLECTION_DATA    static int class_id=-1; int instance_class_id_init=-1;
#define REGISTER_CLASS           class_id = tryRegisterclass(typeid(*this).name() );

#endif
