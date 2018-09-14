
#ifndef BatchBuff_h
#define BatchBuff_h

#include <vector>;

template<typename T>
class BatchBuffPow2{ public:
    int nbatches    = 0;
    int batch_pow   = 0;
    int batch_size  = 0;
    int batch_mask  = 0;
    std::vector<T*> buff;

    void setBatchSize(int batch_pow_){
        batch_pow   = batch_pow_;
        batch_size  = 1<<batch_pow;
        batch_mask  = batch_size-1;
        // TODO: maybe re-allocate buffer?
    }

    void set(int i, T& val){
        int ib = i>>batch_pow;
        if( ib<buff.size() ){;
            T* arr =  buff[ib];
            if(arr){ arr[i&batch_mask] = val; return; }
        }else{
            while( ib>buff.size() ){
                buff.push_back(nullptr);
            }
        }
        T* arr             = new T[batch_size];
        arr [i&batch_mask] = val;
        buff[ib]           = arr;
    }

    T* get(int i){
        T* arr =  buff[i>>batch_pow];
        if(arr){ return arr+(i&batch_mask); }
        else   { return 0; }
    }

    void reserve(int n){
        nbatches=(n>>batch_pow)+1;
        int nb = buff.size();
        for(int i=0;  i<nb;       i++ ){ if(buff[i]==nullptr)buff[i]=new T[batch_size]; }
        for(int i=nb; i<nbatches; i++ ){ buff.push_back( new T[batch_size] ); }
    }

    const T& operator[](int i)const{
        return buff[i>>batch_pow][i&batch_mask];
    }

    T&       operator[](int i){
        return buff[i>>batch_pow][i&batch_mask];
    }

};

template<typename T>
class BatchBuff{ public:
    int nbatches    = 0;
    int batch_size  = 0;
    std::vector<T*> buff;

    void setBatchSize(int batch_size_){
        batch_size  = batch_size_;
    }

    void set(int i, T& val){
        int ib = i/batch_size;
        int io = i%batch_size;
        if( ib<buff.size() ){;
            T* arr =  buff[ib];
            if(arr){ arr[io] = val; return; }
        }else{
            while( ib>buff.size() ){
                buff.push_back(nullptr);
            }
        }
        T* arr             = new T[batch_size];
        arr [io] = val;
        buff[ib]           = arr;
    }

    T* get(int i){
        int ib = i/batch_size;
        int io = i%batch_size;
        T* arr =  buff[ib];
        if(arr){ return arr+(io); }
        else   { return 0; }
    }

    void reserve(int n){
        nbatches=(n/batch_size)+1;
        int nb = buff.size();
        for(int i=0;  i<nb;       i++ ){ if(buff[i]==nullptr)buff[i]=new T[batch_size]; }
        for(int i=nb; i<nbatches; i++ ){ buff.push_back( new T[batch_size] ); }
    }

    const T& operator[](int i)const{
        return buff[i/batch_size][i%batch_size];
    }

    T&       operator[](int i){
        return buff[i/batch_size][i%batch_size];
    }

};

#endif









