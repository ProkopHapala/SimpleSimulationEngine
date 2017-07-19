

## Stack allocation and temp arrays

Today:
```
double * xs = { 1.0, 2.0, 6.0, 15.8 };
evalFunc( xs, [](double x)->double{ return sqrt(x); } );
```
Then:
```
evalFunc( double[]{ 1.0, 2.0, 6.0, 15.8 }, [](double x)->double{ return sqrt(x); } );
```
Notes:

- this works:
```
     double arr_func( int n, const double * xs ){ };
     arr_func( 3, (const double[]){1.0,2.0,3.0} );
```
- https://stackoverflow.com/questions/12907463/passing-static-array-to-function-in-c
- https://stackoverflow.com/questions/15458883/using-array-init-list-as-temporary-in-c11

## with struct operations

#### named bace assignemnt

Today:
```
body.pos={1.0,0.0,0.0};
body.vel.x=1.0;
body.mass=10.0;
```
Then:
```
body.{ pos={1.0,0.0,0.0}; vel.x=1.0; mass=10.0;}
```
Notes:
- https://stackoverflow.com/questions/6181715/convenient-c-struct-initialisation
- https://stackoverflow.com/questions/2279180/does-c-have-with-keyword-like-pascal

