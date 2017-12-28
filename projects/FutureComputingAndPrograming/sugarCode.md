
## Nicer function declaration and calling

It is quite common that you want to call one function by many different ways. In C or C++ you typically make overloaded functions, or you may use default or named parameters (although in C it is much more cumbersome than in python). Why not just to be able setup local variables within function body in clall?

```
void myAddFunc( ){
	double a, b; // locals, you can set them when calling
	return a+b;
};

myAddFunc( )(a=16.5, b=18.9);      // returns 16.5+18.9 = 35.4

```

#### Extend function at label

We often want slightly modified version of function which we already have; Templates and lambda function are some tool to solve this. However, they are not convenient, and forces modifie macrostructure of code.  

```
// function template
void myFunc{
	double c = sqrt(a*a + b*b); 
	@myLabel1 // here we can modify by adding some code
	double c=c + 10;
}

// make modified function
auto modifiedFunc = myFunc.{ 
	@myLabel1: 
		c *= 10.0;
	};

// run modified function
myFunc.{ 
	@myLabel1: 
		c *= 10.0;
	}();

```


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
- https://stackoverflow.com/questions/3016107/what-is-tagged-structure-initialization-syntax
- https://stackoverflow.com/questions/5790534/static-structure-initialization-with-tags-in-c

#### match variables by name

```
class RigidBody{ Vec3f pos; Mat3f rot;   Vec3f vel; Vec3f L; Mat3f invI; };
class BoundBox{ Vec3f pos; Mat3f rot;   Vec3f span; };

RigidBody rb; 
ClassBox  bb;

rb =.{pos,rot} bb;  // <=>      rb.pos=rb.pos; rb.rot=rb.rot;
rb =.{&*} bb; 

rb +=.{pos,rot} bb;  //    rb.pos+=rb.pos; rb.rot+=rb.rot;


```

#### Estension methods

- http://mariusbancila.ro/blog/2014/10/15/extension-methods-in-cpp/
