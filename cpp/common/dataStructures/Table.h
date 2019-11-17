
#ifndef Table_h
#define Table_h

#include <vector>
#include <string>
#include <unordered_map>

enum class DataType{ Bool, Int, Float, Double, String };

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
};

class Table{
    int   itemsize = 0; // number of bytes per item
    void* data     = 0; // pointer with data buffer of unknown type

    std::unordered_map<std::string,int> name2column;
    std::vector       <Atribute>        columns;

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
