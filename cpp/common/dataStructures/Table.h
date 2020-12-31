
#ifndef Table_h
#define Table_h

#include <stdio.h>
#include <vector>
#include <string>
#include <unordered_map>

enum class DataType{ Bool, Int, Float, Double, String };

#define _addColum(T,name,n,D)  { T.addColum( #name, &name, n, DataType::D); }
#define _addColum1d(T,name)    { T.addColum( #name, &name, 1, DataType::Double); }

/*
struct DataColumn{
    void*     data; // buffer where the items are stored
    DataType  type;
};
*/

/*
class Table{
    int n;
    std::unordered_map<std::string,int> name2column;
    std::vector<std::string> columnNames;
    std::vector<DataColumn>  columns;

    //???? getColumn(int i){ return dynamic_cast<columns[i].type>( columns[i].data ); }

    int* get(){ };
};
*/

struct Atribute{
    int      offset;
    int      nsub;    // number of sub units. e.g. Vec3
    //int      nbypes;  // number of bytes per item
    DataType type;    // type for conversion, necessarily exact - Float/Double can exchange

    Atribute() = default;
    Atribute(int offset_,int nsub_,DataType type_):offset(offset_),nsub(nsub_),type(type_){};
};

class Table{ public:
    int n;
    int   itemsize = 0; // number of bytes per item
    char* data     = 0; // pointer with data buffer of unknown type

    std::unordered_map<std::string,int> name2column;
    std::vector       <Atribute>        columns;

    void bind(void* data_, int itemsize_){
        data=(char*)data_;
        itemsize=itemsize_;
    }

    int addColum(void* ptr, int nsub, DataType type=DataType::Double){
        columns.push_back( Atribute( ((char*)ptr)-((char*)data), nsub, type ) );
        return columns.size()-1;
    }

    int addColum( std::string name, void* ptr, int nsub, DataType type=DataType::Double ){
        name2column.insert( {name,columns.size()} );
        return addColum(ptr, nsub, type);
    }

    /*
    void print(int i, int j){
        const Atribute& kind = columns[j];
        void* off = data+itemsize*i+kind.offset;
        int nb;
        switch(kind.type){
            case DataType::Bool   :{ bool*   arr=(bool  *)off; for(int i=0; i<kind.nsub; i++){ printf( "%i \n",  arr[i] ); }} break;
            case DataType::Int    :{ int*    arr=(int   *)off; for(int i=0; i<kind.nsub; i++){ printf( "%i \n",  arr[i] ); }} break;
            case DataType::Float  :{ float*  arr=(float *)off; for(int i=0; i<kind.nsub; i++){ printf( "%g \n",  arr[i] ); }} break;
            case DataType::Double :{ double* arr=(double*)off; for(int i=0; i<kind.nsub; i++){ printf( "%g \n",  arr[i] ); }} break;
            case DataType::String :{ char*   arr=(char  *)off; for(int i=0; i<kind.nsub; i++){ printf( "%s \n",  arr[i] ); }} break;
        }
    }
    */

    char* toStr(int i, int j, char* s){
        const Atribute& kind = columns[j];
        void* off = data+itemsize*i+kind.offset;
        switch(kind.type){
            case DataType::Bool   :{ bool*   arr=(bool  *)off; for(int i=0; i<kind.nsub; i++){ s+=sprintf(s,"%c ",  arr[i]?'T':'F' ); }} break;
            case DataType::Int    :{ int*    arr=(int   *)off; for(int i=0; i<kind.nsub; i++){ s+=sprintf(s,"%i ",  arr[i] ); }} break;
            case DataType::Float  :{ float*  arr=(float *)off; for(int i=0; i<kind.nsub; i++){ s+=sprintf(s,"%g ",  arr[i] ); }} break;
            case DataType::Double :{ double* arr=(double*)off; for(int i=0; i<kind.nsub; i++){ s+=sprintf(s,"%g ",  arr[i] ); }} break;
            case DataType::String :{ char*   arr=(char  *)off; for(int i=0; i<kind.nsub; i++){ s+=sprintf(s,"%s ",  arr[i] ); }} break;
        }
        return s;
    }

    void fromStr( char* s){
        //while(c){
        //
        //}
    }


    /*
    void loadCSV( char* fname, bool bRealloc=true ){
        if(bRealloc){
            delete [] data;

        }
        FILE * pFile;
        const int nbuff = 1024;
        char str[nbuff];
        pFile = fopen ( fname, "r");
        if (pFile == NULL){ printf("file not found: %s \n", fname ); return(-1); };
    }
    */

/*
    template<typename Func>
    void apply(int i, int j, Func func){
        const Atribute& kind = columns[j];
        void* off = data+itemsize*i+kind.offset;
        int nb;
        switch(kind.type){
            case DataType::Bool   : bool*   arr=(bool  *)off; for(int i=0; i<kind.nsub; i++){ func<bool  >( arr[i] ); } break;
            case DataType::Int    : int*    arr=(int   *)off; for(int i=0; i<kind.nsub; i++){ func<int   >( arr[i] ); } break;
            case DataType::Float  : float*  arr=(float *)off; for(int i=0; i<kind.nsub; i++){ func<float >( arr[i] ); } break;
            case DataType::Double : double* arr=(double*)off; for(int i=0; i<kind.nsub; i++){ func<double>( arr[i] ); } break;
            case DataType::String : char*   arr=(char  *)off; for(int i=0; i<kind.nsub; i++){ func<char* >( arr[i] ); } break;
        }
    }
*/

   /*
    template<typename Func>
    void apply(int i, int j, Func func){
        const Atribute& kind = columns[j];
        void* off = data+itemsize*i+kind.offset;
        int nb;
        for(int i=0; i<kind.nsub; i++){ func( arr[i] ); } break;
    }
    */

};



#endif
