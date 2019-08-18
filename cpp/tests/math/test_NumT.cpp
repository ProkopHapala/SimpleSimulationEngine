
#include "testUtils.h"

#include "Num.h"




using namespace NumT;

template<typename T>
T apply_func( int n, T* buff1, T* buff2, const std::function<void(T&,T&,T&)>& func ){
    T val;
    for(int i=0;i<n;i++){ func(buff2[i],buff1[i],val); };
    return val;
}

int main(){


    printf( "----\n" );
    (1<9_r).print();
    (1<9_s).print();

    //( 4_s | 1>>9/-2_s ).print();

    printf( "----\n" );
    (-5_s).print();
    (5/-1_s).print();
    (1<5/-1_s).print();
    //(1,9,-1_s).print();
    printf( "----\n" );
    //(3>>-6_r/2|Scan(4)|Scan(1,9,-1)).print();
    (3<-6_r/2 , 4_s , 1<9/-1_s).print();
    (3<-6_r/2 , 4_s , 1<9/-1_s, 3<19/3_s).print();


    //std::function<void(double&,double&,double&)> f{ _dot<double> };
    //std::function<void(double&,double&,double&)> f{ _dot };

    //Tuple_<int,3> t{1,2,3};
    Tuple<int,3> t{1,2,3};

    //for(int i=0;i<t.nDim; i++){ printf("%i \n",t[i]); }

    /*
    //double xs[] = {1.,2,3,4,};
    double* pxs   = new double[4]{1.,2,3,5};
    double  px[4] = reinterpret_cast<double[4]>(pxs);
    for( double x : (   (double[4])pxs) ){
        printf("x %g \n", x);
    }
    */

    //std::function<void(double&,double&,double&)> f = _dot<double>;
    VecNd v1{1.0d,1.3d,1.5d};
    VecNd v2{3.0d,2.3d,8.5d};
    //VecNd v{3.0d,2.3d,8.5d};
    VecNd v3(  VecNd{0.,1.2,3.66}  );
    VecNd v (5,0.);

    //VecNd v (5,new double[5]);

    //v1.apply_func(v2,f);

    //apply_func(1,v1.buff,v2.buff,f);
    //apply_func<double>(1,v1.buff,v2.buff,f);
    //apply_func<double>(1,v1.buff,v2.buff,_dot);
    //v1.apply_func(v1,v2, _dot<double>);

    println("v1 ",v1);
    println("v2 ",v2);
    println("v3 ",v3);
    println("v  ",v);

    //v1.print();
    //v2.print();
    //v1.scan( [](double& f,double _){ f=f*f; return 0; } );
    //std::function<void(double&,double&)> f = o_max<double>;
    //v1.apply_(f);
    //v1.scan(o_max<double>);
    //v1.apply_(o_max);
    //print( v.apply_( o_max ) );

    std::function<double(const double&,const double&)> f = std::plus<double>();
    //std::function<double(const double&,const double&)> f = std::plus<double>();
    //Func2<double> f = std::plus<double>();

    Func2<double> f2 = std::plus<double>();
    v.scan2(v1,f2 );
    //v.scan2(v1,v2, std::plus<double>() );
    //printf( "v.n %i \n", v.n);
    println("->v  ",v);

    //A* pa = new A{ A{1}+A{3}  };
    //pa->print();


    //(A{1},A{3}).print();

    //printf( "(%i,%i)\n", t.ts[0].i, t.ts[1].i );


    //(Scan(1)|Scan(1)).print();
    //VecN a{0.5,1,0,3,5,6,65848,8487,0.1};
    exit(0);

}
