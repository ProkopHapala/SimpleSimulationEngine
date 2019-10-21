//module Vector_dlib;
module dlib.math.vector;

// Modified from dlib of Timur Gafarov.
// https://github.com/gecko0307/dlib
// https://gecko0307.github.io/dlib/



private{
    import std.conv;
    import std.math;
//    import std.random;
//    import std.range;
//    import std.format;
//    import std.traits;
//    import dlib.core.tuple;

    import dlib.math.utils;
    import dlib.math.matrix;
}

public:

struct Vector(T, int size){ public:
    //Vector constructor
    //Supports initializing from vector of arbitrary length and type
    this (T2, int size2)(Vector!(T2, size2) v){
        if (v.arrayof.length >= size){
            foreach(i; 0..size) arrayof[i] = cast(T)v.arrayof[i];
        }else{
            foreach(i; 0..v.arrayof.length) arrayof[i] = cast(T)v.arrayof[i];
        }
    }

    //Array constructor
    this (A)(A components) if (isDynamicArray!A && !isSomeString!A){
        if (components.length >= size){
            foreach(i; 0..size) arrayof[i] = cast(T)components[i];
        }else{
            foreach(i; 0..components.length) arrayof[i] = cast(T)components[i];
        }
    }

    // Static array constructor
    this (T2, size_t arrSize)(T2[arrSize] components){
        if (components.length >= size){
            foreach(i; 0..size)arrayof[i] = cast(T)components[i];
        } else {
            foreach(i; 0..components.length) arrayof[i] = cast(T)components[i];
        }
    }

   /*
    * Tuple constructor
    */
    this (F...)(F components)
    {
        foreach(i, v; components)
            static if (i < size)
                arrayof[i] = cast(T)v;
    }

   /*
    * String constructor
    */
    this (S)(S str) if (isSomeString!S)
    {
        arrayof = parse!(T[size])(str);
    }

   /*
    * Vector!(T,size) = Vector!(T,size2)
    */
    void opAssign(T2, int size2)(Vector!(T2,size2) v)
    {
        if (v.arrayof.length >= size)
        {
            foreach(i; 0..size)
                arrayof[i] = cast(T)v.arrayof[i];
        }
        else
        {
            foreach(i; 0..v.arrayof.length)
                arrayof[i] = cast(T)v.arrayof[i];
        }
    }

    // Same, but for enums
    void opAssign(T2)(T2 ev) if (is(T2 == enum))
    {
        opAssign(cast(OriginalType!T2)ev);
    }

   /*
    * -Vector!(T,size)
    */
    Vector!(T,size) opUnary(string s) () const if (s == "-")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = -arrayof[i];
        return res;
    }

   /*
    * +Vector!(T,size)
    */
    Vector!(T,size) opUnary(string s) () const if (s == "+")
    body
    {
        return Vector!(T,size)(this);
    }

   /*
    * Vector!(T,size) + Vector!(T,size)
    */
    const Vector!(T,size) opBinary(string op)(Vector!(T,size) v) if (op == "+")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] + v.arrayof[i]);
        return res;
    }

   /*
    * Vector!(T,size) - Vector!(T,size)
    */
    const Vector!(T,size) opBinary(string op)(Vector!(T,size) v) if (op == "-")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] - v.arrayof[i]);
        return res;
    }

   /*
    * Vector!(T,size) * Vector!(T,size)
    */
    const Vector!(T,size) opBinary(string op)(Vector!(T,size) v) if (op == "*")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] * v.arrayof[i]);
        return res;
    }

   /*
    * Vector!(T,size) / Vector!(T,size)
    */
    const Vector!(T,size) opBinary(string op)(Vector!(T,size) v) if (op == "/")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] / v.arrayof[i]);
        return res;
    }

   /*
    * Vector!(T,size) + T
    */
    const Vector!(T,size) opBinary(string op)(T t) if (op == "+")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] + t);
        return res;
    }

   /*
    * Vector!(T,size) - T
    */
    const Vector!(T,size) opBinary(string op)(T t) if (op == "-")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] - t);
        return res;
    }

   /*
    * Vector!(T,size) * T
    */
    const Vector!(T,size) opBinary(string op)(T t) if (op == "*")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] * t);
        return res;
    }

   /*
    * T * Vector!(T,size)
    */
    const Vector!(T,size) opBinaryRight(string op) (T t) if (op == "*" && isNumeric!T)
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] * t);
        return res;
    }

   /*
    * Vector!(T,size) / T
    */
    const Vector!(T,size) opBinary(string op)(T t) if (op == "/")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] / t);
        return res;
    }

   /*
    * Vector!(T,size) % T
    */
    Vector!(T,size) opBinary(string op, T2) (T2 t) const if (op == "%")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] % t);
        return res;
    }

   /*
    * Vector!(T,size) += Vector!(T,size)
    */
    Vector!(T,size) opOpAssign(string op)(Vector!(T,size) v) if (op == "+")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] += v.arrayof[i];
        return this;
    }

   /*
    * Vector!(T,size) -= Vector!(T,size)
    */
    Vector!(T,size) opOpAssign(string op)(Vector!(T,size) v) if (op == "-")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] -= v.arrayof[i];
        return this;
    }

   /*
    * Vector!(T,size) *= Vector!(T,size)
    */
    Vector!(T,size) opOpAssign(string op)(Vector!(T,size) v) if (op == "*")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] *= v.arrayof[i];
        return this;
    }

   /*
    * Vector!(T,size) /= Vector!(T,size)
    */
    Vector!(T,size) opOpAssign(string op)(Vector!(T,size) v) if (op == "/")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] /= v.arrayof[i];
        return this;
    }

   /*
    * Vector!(T,size) += T
    */
    Vector!(T,size) opOpAssign(string op)(T t) if (op == "+")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] += t;
        return this;
    }

   /*
    * Vector!(T,size) -= T
    */
    Vector!(T,size) opOpAssign(string op)(T t) if (op == "-")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] -= t;
        return this;
    }

   /*
    * Vector!(T,size) *= T
    */
    Vector!(T,size) opOpAssign(string op)(T t) if (op == "*")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] *= t;
        return this;
    }

   /*
    * Vector!(T,size) /= T
    */
    Vector!(T,size) opOpAssign(string op)(T t) if (op == "/")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] /= t;
        return this;
    }

   /*
    * Vector!(T,size) %= T
    */
    Vector!(T,size) opOpAssign(string op, T2)(T2 t) if (op == "%")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] %= t;
        return this;
    }

   /*
    * T = Vector!(T,size)[index]
    */
    auto ref T opIndex(this X)(size_t index)
    in
    {
        assert ((0 <= index) && (index < size),
            "Vector!(T,size).opIndex(int index): array index out of bounds");
    }
    body
    {
        return arrayof[index];
    }

   /*
    * Vector!(T,size)[index] = T
    */
    void opIndexAssign(T n, size_t index)
    in
    {
        assert ((0 <= index) && (index < size),
            "Vector!(T,size).opIndexAssign(int index): array index out of bounds");
    }
    body
    {
        arrayof[index] = n;
    }

   /*
    * T[] = Vector!(T,size)[index1..index2]
    */
    auto opSlice(this X)(size_t index1, size_t index2)
    in
    {
        assert ((0 <= index1) || (index1 < 3) || (0 <= index2) || (index2 < 3) || (index1 < index2),
            "Vector!(T,size).opSlice(int index1, int index2): array index out of bounds");
    }
    body
    {
        return arrayof[index1..index2];
    }

   /*
    * Vector!(T,size)[index1..index2] = T
    */
    T opSliceAssign(T t, size_t index1, size_t index2)
    in
    {
        assert ((0 <= index1) || (index1 < 3) || (0 <= index2) || (index2 < 3) || (index1 < index2),
            "Vector!(T,size).opSliceAssign(T t, int index1, int index2): array index out of bounds");
    }
    body
    {
        arrayof[index1..index2] = t;
        return t;
    }

   /*
    * T = Vector!(T,size)[]
    */
    auto opSlice(this X)()
    body
    {
        return arrayof[];
    }

   /*
    * Vector!(T,size)[] = T
    */
    T opSliceAssign(T t)
    body
    {
        //arrayof[] = t;
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] = t;
        return t;
    }

    static if (isNumeric!(T))
    {
       /*
        * Get vector length squared
        */
        @property T lengthsqr() const
        body
        {
            T res = 0;
            foreach (component; arrayof)
                res += component * component;
            return res;
        }

       /*
        * Get vector length
        */
        @property T length() const
        body
        {
            static if (isFloatingPoint!T)
            {
                T t = 0;
                foreach (component; arrayof)
                    t += component * component;
                return sqrt(t);
            }
            else
            {
                T t = 0;
                foreach (component; arrayof)
                    t += component * component;
                return cast(T)sqrt(cast(float)t);
            }
        }

       /*
        * Set vector length to 1
        */
        void normalize()
        body
        {
            static if (isFloatingPoint!T)
            {
                T lensqr = lengthsqr();
                if (lensqr > 0)
                {
                    T coef = 1.0 / sqrt(lensqr);
                    foreach (ref component; arrayof)
                        component *= coef;
                }
            }
            else
            {
                T lensqr = lengthsqr();
                if (lensqr > 0)
                {
                    float coef = 1.0 / sqrt(cast(float)lensqr);
                    foreach (ref component; arrayof)
                        component = cast(T)(component * coef);
                }
            }
        }

       /*
        * Return normalized copy
        */
        @property Vector!(T,size) normalized() const
        body
        {
            Vector!(T,size) res = this;
            res.normalize();
            return res;
        }

       /*
        * Return true if all components are zero
        */
        @property bool isZero() const
        body
        {
            foreach(i; RangeTuple!(0, size))
                if (arrayof[i] != 0)
                    return false;
            return true;
        }

       /*
        * Clamp components to min/max value
        */
        void clamp(T minv, T maxv)
        {
            foreach (ref v; arrayof)
                v = .clamp(v, minv, maxv);
        }
    }

   /*
    * Convert to string
    */
    @property string toString() const
    body
    {
        auto writer = appender!string();
        formattedWrite(writer, "%s", arrayof);
        return writer.data;
    }

   /*
    * Swizzling
    */
    template opDispatch(string s) if (valid(s))
    {
    /*
        static if (s.length == 1)
        {
            enum i = ["x":0, "y":1, "z":2, "w":3,
                      "r":0, "g":1, "b":2, "a":3,
                      "s":0, "t":1, "p":2, "q":3][s];

            @property auto ref opDispatch(this X)()
            {
                return arrayof[i];
            }

            @property auto ref opDispatch(this X, V)(auto ref V v)
            {
                return arrayof[i] = v;
            }
        }
        else*/
        static if (s.length <= 4)
        {
            @property auto ref opDispatch(this X)()
            {
                auto extend(string s)
                {
                    while (s.length < 4)
                        s ~= s[$-1];
                    return s;
                }

                enum p = extend(s);
                enum i = (char c) => ['x':0, 'y':1, 'z':2, 'w':3,
                                      'r':0, 'g':1, 'b':2, 'a':3,
                                      's':0, 't':1, 'p':2, 'q':3][c];
                enum i0 = i(p[0]),
                     i1 = i(p[1]),
                     i2 = i(p[2]),
                     i3 = i(p[3]);

                static if (s.length == 4)
                    return Vector!(T,4)(arrayof[i0], arrayof[i1], arrayof[i2], arrayof[i3]);
                else static if (s.length == 3)
                    return Vector!(T,3)(arrayof[i0], arrayof[i1], arrayof[i2]);
                else static if (s.length == 2)
                    return Vector!(T,2)(arrayof[i0], arrayof[i1]);
            }
        }
    }

    private static bool valid(string s)
    {
        if (s.length < 2)
            return false;

        foreach(c; s)
        {
            switch(c)
            {
                case 'w', 'a', 'q':
                    if (size < 4) return false;
                    else break;
                case 'z', 'b', 'p':
                    if (size < 3) return false;
                    else break;
                case 'y', 'g', 't':
                    if (size < 2) return false;
                    else break;
                case 'x', 'r', 's':
                    if (size < 1) return false;
                    else break;
                default:
                    return false;
            }
        }
        return true;
    }

   /*
    * Symbolic element access
    */
    private static string elements(string[4] letters) @property
    body
    {
        string res;
        foreach (i; 0..size)
        {
            res ~= "T " ~ letters[i] ~ "; ";
        }
        return res;
    }

   /*
    * Vector components
    */
    union
    {
        T[size] arrayof;
        static if (size < 5)
        {
            struct { mixin(elements(["x", "y", "z", "w"])); }
            struct { mixin(elements(["r", "g", "b", "a"])); }
            struct { mixin(elements(["s", "t", "p", "q"])); }
        }
    }
}

/*
 * Dot product
 */
T dot(T, int size) (Vector!(T,size) a, Vector!(T,size) b)
body
{
    static if (size == 1)
    {
        return a.x * b.x;
    }
    else
    static if (size == 2)
    {
        return ((a.x * b.x) + (a.y * b.y));
    }
    else
    static if (size == 3)
    {
        return ((a.x * b.x) + (a.y * b.y) + (a.z * b.z));
    }
    else
    {
        T d = 0;
        //foreach (i; 0..size)
        foreach(i; RangeTuple!(0, size))
            d += a[i] * b[i];
        return d;
    }
}

//
// Cross product
//
Vector!(T,size) cross(T, int size) (Vector!(T,size) a, Vector!(T,size) b) if (size == 3)body{
    return Vector!(T,size)(
        (a.y * b.z) - (a.z * b.y),
        (a.z * b.x) - (a.x * b.z),
        (a.x * b.y) - (a.y * b.x)
    );
}

Vector!(T,size) cross(T, int size) (Vector!(T,size) a, Vector!(T,size) b, Vector!(T,size) c) if (size == 4)body{
    return Vector!(T,size)(
        (a.y * b.z * c.w) - (a.y * b.w * c.z) + (a.z * b.w * c.y) - (a.z * b.y * c.w) + (a.w * b.y * c.z) - (a.w * b.z * c.y),
        (a.z * b.w * c.x) - (a.z * b.x * c.w) + (a.w * b.x * c.z) - (a.w * b.z * c.x) + (a.x * b.z * c.w) - (a.x * b.w * c.z),
        (a.w * b.x * c.y) - (a.w * b.y * c.x) + (a.x * b.y * c.w) - (a.x * b.w * c.y) + (a.y * b.w * c.x) - (a.y * b.x * c.w),
        (a.x * b.y * c.z) - (a.x * b.z * c.y) + (a.y * b.z * c.x) - (a.y * b.x * c.z) + (a.z * b.x * c.y) - (a.z * b.y * c.x)
    );
}

/*

TODO: Shou;ld go to matrix

//
// Tensor product
//
Matrix!(T,N) tensorProduct(T, size_t N) (Vector!(T,N) u, Vector!(T,N) v)body{
    Matrix!(T,N) res;
    foreach(i; 0..N) foreach(j; 0..N) {
        res[i, j] = u[i] * v[j];
    }
    return res;
}
*/

alias outerProduct = tensorProduct;

/*
 * Compute normal of a plane from three points
 */
Vector!(T,3) planeNormal(T) (Vector!(T,3) p1, Vector!(T,3) p2, Vector!(T,3) p3)body{
    Vector!(T,3) vec1 = Vector!(T,3)(p1 - p2);
    Vector!(T,3) vec2 = Vector!(T,3)(p1 - p3);

    Vector!(T,3) result = Vector!(T,3)(cross(vec1,vec2));
    result.normalize();

    return result;
}

void rotateAroundAxis(T) (ref Vector!(T,3) V, Vector!(T,3) P, Vector!(T,3) D, T angle){
    T axx,axy,axz,ax1;
    T ayx,ayy,ayz,ay1;
    T azx,azy,azz,az1;

    T u,v,w;
    T u2,v2,w2;
    T a,b,c;

    T sa,ca;

    sa = sin(angle);
    ca = cos(angle);

    u = D.x;
    v = D.y;
    w = D.z;

    u2 = u * u;
    v2 = v * v;
    w2 = w * w;

    a = P.x;
    b = P.y;
    c = P.z;

    axx = u2+(v2+w2)*ca;
    axy = u*v*(1-ca)-w*sa;
    axz = u*w*(1-ca)+v*sa;
    ax1 = a*(v2+w2)-u*(b*v+c*w)+(u*(b*v+c*w)-a*(v2+w2))*ca+(b*w-c*v)*sa;

    ayx = u*v*(1-ca)+w*sa;
    ayy = v2+(u2+w2)*ca;
    ayz = v*w*(1-ca)-u*sa;
    ay1 = b*(u2+w2)-v*(a*u+c*w)+(v*(a*u+c*w)-b*(u2+w2))*ca+(c*u-a*w)*sa;

    azx = u*w*(1-ca)-v*sa;
    azy = v*w*(1-ca)+u*sa;
    azz = w2+(u2+v2)*ca;

    az1 = c*(u2+v2)-w*(a*u+b*v)+(w*(a*u+b*v)-c*(u2+v2))*ca+(a*v-b*u)*sa;

    Vector!(T,3) W;
    W.x = axx * V.x + axy * V.y + axz * V.z + ax1;
    W.y = ayx * V.x + ayy * V.y + ayz * V.z + ay1;
    W.z = azx * V.x + azy * V.y + azz * V.z + az1;

    V = W;
}

/*
 * Compute distance between two points
 */
T distance(T) (Vector!(T,2) a, Vector!(T,2) b)
body
{
    T dx = a.x - b.x;
    T dy = a.y - b.y;
    return sqrt((dx * dx) + (dy * dy));
}

T distancesqr(T) (Vector!(T,2) a, Vector!(T,2) b) body{
    T dx = a.x - b.x;
    T dy = a.y - b.y;
    return ((dx * dx) + (dy * dy));
}

T distance(T) (Vector!(T,3) a, Vector!(T,3) b) body {
    T dx = a.x - b.x;
    T dy = a.y - b.y;
    T dz = a.z - b.z;
    return sqrt((dx * dx) + (dy * dy) + (dz * dz));
}

T distancesqr(T) (Vector!(T,3) a, Vector!(T,3) b) body {
    T dx = a.x - b.x;
    T dy = a.y - b.y;
    T dz = a.z - b.z;
    return ((dx * dx) + (dy * dy) + (dz * dz));
}

//
// Random unit length vectors
//
Vector!(T,2) randomUnitVector2(T)(){
    float azimuth = uniform(0.0, 1.0) * 2 * PI;
    return Vector!(T,2)(cos(azimuth), sin(azimuth));
}

Vector!(T,3) randomUnitVector3(T)(){
    float z = (2 * uniform(0.0, 1.0)) - 1;
    Vector!(T,2) planar = randomUnitVector2!(T)() * sqrt(1 - z * z);
    return Vector!(T,3)(planar.x, planar.y, z);
}

//
// Spherical linear interpolation
// (simple lerp is in dlib.math.interpolation)
//
Vector!(T,3) slerp(T) (Vector!(T,3) a, Vector!(T,3) b, T t){
    T dp = dot(a, b);
    dp = clamp(dp, -1.0, 1.0);
    T theta = acos(dp) * t;
    Vector!(T,3) relativeVec = b - a * dp;
    relativeVec.normalize();
    return ((a * cos(theta)) + (relativeVec * sin(theta)));
}

//
// Gradually decrease vector to zero length
//
Vector!(T,3) vectorDecreaseToZero(T) (Vector!(T,3) vector, T step){
    foreach (ref component; vector.arrayof){
        if (component > 0.0) component -= step;
        if (component < 0.0) component += step;
    }
    return vector;
}

//
// Almost zero
//
bool isAlmostZero(Vector3f v){ return (isConsiderZero(v.x) && isConsiderZero(v.y) && isConsiderZero(v.z)); }

Vector!(T,3) reflect    (T)( Vector!(T,3) I, Vector!(T,3) N ){ return I - N * dot(N, I) * 2.0; }
Vector!(T,3) faceforward(T)( Vector!(T,3) N, Vector!(T,3) I, Vector!(T,3) Nref ){    return dot(Nref, I) < 0.0 ? N : -N; }
Vector!(T,3) refract(T)    ( Vector!(T,3) I, Vector!(T,3) N, T r ){
    T d = 1.0 - r * r * (1.0 - dot(N, I) * dot(N, I));
    if (d < 0.0) return Vector!(T,3)(0.0, 0.0, 0.0);
    return I * r - N * (r * dot(N, I) + sqrt(d));
}



/*
 * Predefined vector types
 */
alias Vector2i = Vector!(int, 2);
alias Vector2u = Vector!(uint, 2);
alias Vector2f = Vector!(float, 2);
alias Vector2d = Vector!(double, 2);

alias Vector3i = Vector!(int, 3);
alias Vector3u = Vector!(uint, 3);
alias Vector3f = Vector!(float, 3);
alias Vector3d = Vector!(double, 3);

alias Vector4i = Vector!(int, 4);
alias Vector4u = Vector!(uint, 4);
alias Vector4f = Vector!(float, 4);
alias Vector4d = Vector!(double, 4);

/*
 * GLSL-like short aliases
 */
alias ivec2 = Vector2i;
alias uvec2 = Vector2u;
alias vec2  = Vector2f;
alias dvec2 = Vector2d;

alias ivec3 = Vector3i;
alias uvec3 = Vector3u;
alias vec3  = Vector3f;
alias dvec3 = Vector3d;

alias ivec4 = Vector4i;
alias uvec4 = Vector4u;
alias vec4  = Vector4f;
alias dvec4 = Vector4d;

//
// Axis vectors
//
static struct AxisVector
{
    Vector3f x = Vector3f(1.0f, 0.0f, 0.0f);
    Vector3f y = Vector3f(0.0f, 1.0f, 0.0f);
    Vector3f z = Vector3f(0.0f, 0.0f, 1.0f);
}

//
// Vector factory function
//
auto vectorf(T...)(T t) if (t.length > 0)
{
    return Vector!(float, t.length)(t);
}

//
// L-value pseudovector for assignment purposes.
//
// Usage example:
//
//  float a, b, c
//  lvector(a, b, c) = Vector3f(10, 4, 2);
//
auto lvector(T...)(ref T x)
{
    struct Result(T, uint size){
        T*[size] arrayof;

        void opAssign(int size2)(Vector!(T,size2) v){
            if (v.arrayof.length >= size)
                foreach(i; 0..size) *arrayof[i] = v.arrayof[i];
            else
                foreach(i; 0..v.arrayof.length) *arrayof[i] = v.arrayof[i];
        }
    }
    auto res = Result!(typeof(x[0]), x.length)();
    foreach(i, ref v; x) res.arrayof[i] = &v;
    return res;
}

unittest{
    {
        const vec3 a = vec3(10.5f, 20.0f, 33.12345f);
        const vec3 b = -a;
        const vec3 c = +a - b;
        const vec3 d = a * b / c;

        assert(isAlmostZero(to!vec3(c.toString()) - c));

        ivec2 ab = ivec2(5, 15);
        ab += ivec2(20, 30);
        ab *= 3;

        assert(ab[0] == 75 && ab[1] == 135);

        auto len = c.length();
        auto lensqr = c.lengthsqr();
        auto dist = distance(a, b);

        auto xy = a[0..1];
        auto n = a[];

        vec3 v1 = vec3(2.0f, 0.0f, 1.0f);
        ivec3 v2 = v1;
        assert(ivec3(v1) == ivec3(2, 0, 1));

        vec3 v3 = [0, 2, 3.5];
        assert(v3 == vec3(0.0f, 2.0f, 3.5f));

        ivec3 v4 = [7, 8, 3];
        v4 %= 2;
        assert(v4 == ivec3(1, 0, 1));
    }

    {
        Vector3f a = Vector3f(1, 2, 3);
        Vector2f b = Vector2f(a);
        assert(b == Vector2f(1, 2));
    }

    {
        Vector3f a = Vector3f([0, 1]);
        assert(isNaN(a.z));
    }

    {
        Vector3f a = Vector3f(0, 1, 2);
        a += 1;
        assert(a == Vector3f(1, 2, 3));
        a *= 2;
        assert(a == Vector3f(2, 4, 6));
        a -= 1;
        assert(a == Vector3f(1, 3, 5));
        a /= 3;
        assert(a.y == 1);
    }

    {
        Vector3f a;
        a[1] = 3;
        assert(a.y == 3);

        a[0..3] = 1;
        assert(a == Vector3f(1, 1, 1));

        a[] = 0;
        assert(a == Vector3f(0, 0, 0));
    }

    {
        Vector3i a = Vector3i(0, 0, 3);
        a = a.normalized;
        assert(a == Vector3i(0, 0, 1));
        assert(a.length == 1);
    }

    {
        Vector3f a = Vector3f(0, 0, 0);
        assert(a.isZero);
    }

    {
        Vector3f a = Vector3f(2, -3, 0);
        a.clamp(-1, 1);
        assert(a == Vector3f(1, -1, 0));
    }

    {
        Vector3f a = Vector3f(1, 2, 3, 4);
        assert(dot(a, a) == 14);
    }
}

