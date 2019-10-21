module Vector;

public:

struct Vec(T, int size){ public:

    // Vector components
    union{
        T[size] arrayof;
        static if (size < 5){
            struct { mixin(elements(["x", "y", "z", "w"])); }
            //struct { mixin(elements(["a", "b", "c", "d"])); }  // collision
            struct { mixin(elements(["i", "j", "k", "l"])); }
            struct { mixin(elements(["r", "g", "b", "a"])); }
            //struct { mixin(elements(["s", "t", "p", "q"])); }
        }
    }

    // =========== Swizzling

    template opDispatch(string s) if (valid(s)){
        static if (s.length <= 4){
            @property auto ref opDispatch(this X)(){
                auto extend(string s) {
                    while (s.length < 4) s ~= s[$-1];
                    return s;
                }
                enum p = extend(s);
                enum i = (char c) => ['x':0, 'y':1, 'z':2, 'w':3,
                 //                   'a':0, 'b':1, 'c':2, 'd':3,    // collision
                                      'i':0, 'j':1, 'k':2, 'l':3,
                                      'r':0, 'g':1, 'b':2, 'a':3,][c];
                enum i0 = i(p[0]), i1 = i(p[1]), i2 = i(p[2]), i3 = i(p[3]);
                static if (s.length == 4)      return Vector!(T,4)(arrayof[i0], arrayof[i1], arrayof[i2], arrayof[i3]);
                else static if (s.length == 3) return Vector!(T,3)(arrayof[i0], arrayof[i1], arrayof[i2]);
                else static if (s.length == 2) return Vector!(T,2)(arrayof[i0], arrayof[i1]);
            }
        }
    }

    private static bool valid(string s){
        if (s.length < 2)return false;
        foreach(c; s){
            switch(c){
                case 'w', 'a', 'l' : if (size < 4) return false; else break;
                case 'z', 'b', 'k' : if (size < 3) return false; else break;
                case 'y', 'g', 'j' : if (size < 2) return false; else break;
                case 'x', 'r', 'i' : if (size < 1) return false; else break;
                default: return false;
            }
        }
        return true;
    }

    // Symbolic element access
    private static string elements(string[4] letters) @property body {
        string res;
        foreach (i; 0..size) { res ~= "T " ~ letters[i] ~ "; "; }
        return res;
    }

    // ============= Constructors

    //Vector constructor
    //Supports initializing from vector of arbitrary length and type
    this (T2, int size2)(Vector!(T2, size2) v){
        if (v.arrayof.length >= size){ foreach(i; 0..size) arrayof[i] = cast(T)v.arrayof[i]; }
        else                         { foreach(i; 0..v.arrayof.length) arrayof[i] = cast(T)v.arrayof[i]; }
    }

    //Array constructor
    this (A)(A components) if (isDynamicArray!A && !isSomeString!A){
        if (components.length >= size){ foreach(i; 0..size) arrayof[i] = cast(T)components[i]; }
        else                          { foreach(i; 0..components.length) arrayof[i] = cast(T)components[i]; }
    }

    // Static array constructor
    this (T2, size_t arrSize)(T2[arrSize] components){
        if (components.length >= size){ foreach(i; 0..size)arrayof[i] = cast(T)components[i]; } 
        else                          { foreach(i; 0..components.length) arrayof[i] = cast(T)components[i]; }
    }

    this (F...)(F components){ foreach(i, v; components) static if (i < size) arrayof[i] = cast(T)v; } // Tuple constructor
    this (S)(S str) if (isSomeString!S){ arrayof = parse!(T[size])(str); }                             // String constructor


    // ============= Index

    // T = Vector!(T,size)[index]
    auto ref T opIndex(this X)(size_t index) in {
        assert ((0 <= index) && (index < size), "Vector!(T,size).opIndex(int index): array index out of bounds" );
    }body{
        return arrayof[index];
    }

    // Vector!(T,size)[index] = T
    void opIndexAssign(T n, size_t index) in { 
        assert ((0 <= index) && (index < size), "Vector!(T,size).opIndexAssign(int index): array index out of bounds"); 
    }body{
        arrayof[index] = n;
    }


    // T[] = Vector!(T,size)[index1..index2]
    auto opSlice(this X)(size_t index1, size_t index2) in {
        assert ((0 <= i1) || (i1 < 3) || (0 <= i2) || (i2 < 3) || (i1 < i2), "Vector!(T,size).opSlice(int index1, int index2): array index out of bounds");
    } body {
        return arrayof[i1..i2];
    }


    // Vector!(T,size)[index1..index2] = T
    T opSliceAssign(T t, size_t index1, size_t index2)in{
        assert ((0 <= index1) || (index1 < 3) || (0 <= index2) || (index2 < 3) || (index1 < index2), "Vector!(T,size).opSliceAssign(T t, int index1, int index2): array index out of bounds");
    } body {
        arrayof[index1..index2] = t;
        return t;
    }

    // T = Vector!(T,size)[]
    auto opSlice(this X)()body{
        return arrayof[];
    }

/*
    // Vector!(T,size)[] = T
    T opSliceAssign(T t)
    body{
        foreach(i; RangeTuple!(0, size)) arrayof[i] = t;
        return t;
    }
*/

    // ============= Assign

    // Vector!(T,size) = Vector!(T,size2)
    void opAssign(T2, int size2)(Vector!(T2,size2) v){
        if (v.arrayof.length >= size){ foreach(i; 0..size) arrayof[i] = cast(T)v.arrayof[i]; }
        else                         { foreach(i; 0..v.arrayof.length)arrayof[i] = cast(T)v.arrayof[i];  }
    }

    // Same, but for enums
    void opAssign(T2)(T2 ev) if (is(T2 == enum)) {
        opAssign(cast(OriginalType!T2)ev);

    }

    void opAssign(T2, int size2)(Vector!(T2,size2) v){
        if (v.arrayof.length >= size){ foreach(i; 0..size)             arrayof[i] = cast(T)v.arrayof[i]; }
        else                         { foreach(i; 0..v.arrayof.length) arrayof[i] = cast(T)v.arrayof[i]; }
    }

    void opAssign(T2)(T2 ev) if (is(T2 == enum)) { opAssign(cast(OriginalType!T2)ev); }

    void opAssign(string op)() if(op=="-"||op=="/"||op=="~"){
        foreach(i; 0..size) mixin( "arrayof[i]"+op+" = arrayof[i];" );
    }

    void opAssign(string op)(Vector!(T,size) v) if(op=="+" ||op=="-" ||op=="*" ||op=="/" ) {
        foreach(i; 0..size) mixin( "arrayof[i] "+op+"= cast(T)v.arrayof[i];" ); 
    }

    void opAssign(string op)(Vector!(T,size) v1, Vector!(T,size) v2) if(op=="+" ||op=="-" ||op=="*" ||op=="/" ) {
        foreach(i; 0..size) mixin( "arrayof[i] = cast(T)v1.arrayof[i] +"+op+" cast(T)v2.arrayof[i];" ); 
    }

    // ============= Operations

    Vector!(T,size) opUnary(string op) () const body{
        Vector!(T,size) res;
        foreach(i; 0..size) mixin( "res.arrayof[i] = "+op+" arrayof[i];" ); 
        return res;
    }

    Vector!(T,size) opUnary(string op)( Vector!(T,size) v ) const body{
        Vector!(T,size) res;
        res.opAssign!op( v );
        return res;
    }

    Vector!(T,size) opBinary(string op)( Vector!(T,size) v1, Vector!(T,size) v2 ) const body{
        Vector!(T,size) res;
        res.opAssign!op( v1, v2 );
        return res;
    }



/*

    // -Vector!(T,size)
    Vector!(T,size) opUnary(string s) () const if (s == "-") body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = -arrayof[i];
        return res;
    }

    // +Vector!(T,size)
    Vector!(T,size) opUnary(string s) () const if (s == "+")
    body
    {
        return Vector!(T,size)(this);
    }


    // Vector!(T,size) + Vector!(T,size)
    const Vector!(T,size) opBinary(string op)(Vector!(T,size) v) if (op == "+")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] + v.arrayof[i]);
        return res;
    }

    // Vector!(T,size) - Vector!(T,size)
    const Vector!(T,size) opBinary(string op)(Vector!(T,size) v) if (op == "-")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] - v.arrayof[i]);
        return res;
    }

    // Vector!(T,size) * Vector!(T,size)
    const Vector!(T,size) opBinary(string op)(Vector!(T,size) v) if (op == "*")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] * v.arrayof[i]);
        return res;
    }

    // Vector!(T,size) / Vector!(T,size)
    const Vector!(T,size) opBinary(string op)(Vector!(T,size) v) if (op == "/")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] / v.arrayof[i]);
        return res;
    }

    // Vector!(T,size) + T
    const Vector!(T,size) opBinary(string op)(T t) if (op == "+")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] + t);
        return res;
    }

    // Vector!(T,size) - T
    const Vector!(T,size) opBinary(string op)(T t) if (op == "-")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] - t);
        return res;
    }

    // Vector!(T,size) * T
    const Vector!(T,size) opBinary(string op)(T t) if (op == "*")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] * t);
        return res;
    }

    // T * Vector!(T,size)
    const Vector!(T,size) opBinaryRight(string op) (T t) if (op == "*" && isNumeric!T)
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] * t);
        return res;
    }

    // Vector!(T,size) / T
    const Vector!(T,size) opBinary(string op)(T t) if (op == "/")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] / t);
        return res;
    }

    // Vector!(T,size) % T
    Vector!(T,size) opBinary(string op, T2) (T2 t) const if (op == "%")
    body
    {
        Vector!(T,size) res;
        foreach(i; RangeTuple!(0, size))
            res.arrayof[i] = cast(T)(arrayof[i] % t);
        return res;
    }

    // Vector!(T,size) += Vector!(T,size)
    Vector!(T,size) opOpAssign(string op)(Vector!(T,size) v) if (op == "+")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] += v.arrayof[i];
        return this;
    }

    // Vector!(T,size) -= Vector!(T,size)
    Vector!(T,size) opOpAssign(string op)(Vector!(T,size) v) if (op == "-")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] -= v.arrayof[i];
        return this;
    }

    // Vector!(T,size) *= Vector!(T,size)
    Vector!(T,size) opOpAssign(string op)(Vector!(T,size) v) if (op == "*")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] *= v.arrayof[i];
        return this;
    }

    // Vector!(T,size) /= Vector!(T,size)
    Vector!(T,size) opOpAssign(string op)(Vector!(T,size) v) if (op == "/")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] /= v.arrayof[i];
        return this;
    }

    // Vector!(T,size) += T
    Vector!(T,size) opOpAssign(string op)(T t) if (op == "+")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] += t;
        return this;
    }

    // Vector!(T,size) -= T
    Vector!(T,size) opOpAssign(string op)(T t) if (op == "-")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] -= t;
        return this;
    }

    // Vector!(T,size) *= T
    Vector!(T,size) opOpAssign(string op)(T t) if (op == "*")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] *= t;
        return this;
    }

    // Vector!(T,size) /= T
    Vector!(T,size) opOpAssign(string op)(T t) if (op == "/")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] /= t;
        return this;
    }

    // Vector!(T,size) %= T
    Vector!(T,size) opOpAssign(string op, T2)(T2 t) if (op == "%")
    body
    {
        //foreach(i; 0..size)
        foreach(i; RangeTuple!(0, size))
            arrayof[i] %= t;
        return this;
    }

*/

   

/*
    static if (isNumeric!(T))
    {

        // Get vector length squared
        @property T lengthsqr() const
        body
        {
            T res = 0;
            foreach (component; arrayof)
                res += component * component;
            return res;
        }

        // Get vector length
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

        // Set vector length to 1
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

        // Return normalized copy
        @property Vector!(T,size) normalized() const
        body
        {
            Vector!(T,size) res = this;
            res.normalize();
            return res;
        }

       //Return true if all components are zero
        @property bool isZero() const
        body
        {
            foreach(i; RangeTuple!(0, size))
                if (arrayof[i] != 0)
                    return false;
            return true;
        }

       // Clamp components to min/max value
        void clamp(T minv, T maxv)
        {
            foreach (ref v; arrayof)
                v = .clamp(v, minv, maxv);
        }
    }
*/


/*
    // Convert to string
    @property string toString() const
    body
    {
        auto writer = appender!string();
        formattedWrite(writer, "%s", arrayof);
        return writer.data;
    }

*/





} // Vector


alias Vec2i = Vec!(int, 2);
alias Vec2u = Vec!(uint, 2);
alias Vec2f = Vec!(float, 2);
alias Vec2d = Vec!(double, 2);

alias Vec3i = Vec!(int, 3);
alias Vec3u = Vec!(uint, 3);
alias Vec3f = Vec!(float, 3);
alias Vec3d = Vec!(double, 3);

alias Vec4i = Vec!(int, 4);
alias Vec4u = Vec!(uint, 4);
alias Vec4f = Vec!(float, 4);
alias Vec4d = Vec!(double, 4);

/*
// GLSL-like short aliases
alias ivec2 = Vec2i;
alias uvec2 = Vec2u;
alias vec2  = Vec2f;
alias dvec2 = Vec2d;

alias ivec3 = Vec3i;
alias uvec3 = Vec3u;
alias vec3  = Vec3f;
alias dvec3 = Vec3d;

alias ivec4 = Vec4i;
alias uvec4 = Vec4u;
alias vec4  = Vec4f;
alias dvec4 = Vec4d;
*/


T dot(T, int size) (Vec!(T,size) a, Vec!(T,size) b)body{
    static      if (size == 1){ return  a.x * b.x; }
    else static if (size == 2){ return (a.x * b.x) + (a.y * b.y); }
    else static if (size == 3){ return (a.x * b.x) + (a.y * b.y) + (a.z * b.z); }
    else{
        T d = 0;
        foreach(i; 0..size) d += a[i] * b[i];
        return d;
    }
}
