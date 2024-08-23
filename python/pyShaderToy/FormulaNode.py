class Formula:
    def __init__(self, value):
        self.value = value
        self.left  = None
        self.right = None
        self.op    = None

    def __add__(self, other):
        return self._build_op_node(other, '+')

    def __sub__(self, other):
        return self._build_op_node(other, '-')

    def __mul__(self, other):
        return self._build_op_node(other, '*')

    def __truediv__(self, other):
        return self._build_op_node(other, '/')

    def _build_op_node(self, other, op):
        new_node       = Formula(None)
        new_node.left  = self
        new_node.right = other if isinstance(other, Formula) else Formula(other)
        new_node.op    = op
        return new_node

    def eval(self):
        if self.op is None:
            return str(self.value)
        else:
            left_str  = self.left.eval()
            right_str = self.right.eval()
            return f"({left_str} {self.op} {right_str})"
        
    def eval_CGS(self):
        if self.op is None:
            return str(self.value)
        else:
            left_str  = self.left.eval_CGS()
            right_str = self.right.eval_CGS()

            if( self.op == "+" ):
                return f" opU({left_str} , {right_str})"
            elif( self.op == "-" ):
                return f"opS({left_str} ,{right_str})"
            #elif( self.op == "*" ):
            #    return f"({left_str} * {right_str})")
            #elif( self.op == "/" ):
            #    return f"({left_str} / {right_str})")


def process_vec_arg( a ):
    if( isinstance( a, str ) ):
        return a
    elif( isinstance( a, float ) ):
        return str(a)
    elif ( isinstance( a, int ) ):
        return ("vec_data[%i]" %a)
    elif( isinstance( a, list ) ):
        return "vec3( " + ", ".join( process_vec_arg( x ) for x in a ) + " )"
    elif( isinstance( a, tuple ) ):
        return "vec3( " + ", ".join( process_vec_arg( x ) for x in a ) + " )"
    #elif( isinstance( a, np.ndarray ) ):
    #    return "vec3( " + ", ".join( process_vec_arg( x ) for x in a ) + " )"
    else:
        return str(a)


# float sdPlane    ( vec3 p, vec3 d )          { return dot(d,p); }
# float sdSphere   ( vec3 p, float s )         { return length(p)-s; }
# float sdBox      ( vec3 p, vec3 b )          { vec3 d = abs(p) - b; return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0)); }
# float sdEllipsoid( in vec3 p, in vec3 r )    { return (length( p/r ) - 1.0) * min(min(r.x,r.y),r.z); }
# float udRoundBox ( vec3 p, vec3 b, float r ) { return length(max(abs(p)-b,0.0))-r; }
# float sdTorus    ( vec3 p, vec2 t )          { return length( vec2(length(p.xz)-t.x,p.y) )-t.y; }
# float length2    ( vec2 p )                  { return sqrt( p.x*p.x + p.y*p.y ); }
# float length6    ( vec2 p )                  { p = p*p*p; p = p*p; return pow( p.x + p.y, 1.0/6.0 ); }
# float length8    ( vec2 p )                  { p = p*p; p = p*p; p = p*p; return pow( p.x + p.y, 1.0/8.0 ); }
# float sdTorus82  ( vec3 p, vec2 t )          { vec2 q = vec2(length2(p.xz)-t.x,p.y); return length8(q)-t.y; }
# float sdTorus88  ( vec3 p, vec2 t )          { vec2 q = vec2(length8(p.xz)-t.x,p.y); return length8(q)-t.y; }
# float sdCylinder6( vec3 p, vec2 h )          { return max( length6(p.xz)-h.x, abs(p.y)-h.y ); }

def Plane(p, d=(0.0,1.0,0.0) ):
    p = process_vec_arg(p,  )
    d = process_vec_arg(d)
    # sdPlane    ( vec3 p )
    return Formula( f"sdPlane2( pos-{p}, {d} )" )
def Sphere(p,r):
    p = process_vec_arg(p)
    # sdSphere   ( vec3 p, float r ) 
    return Formula( f"sdSphere( pos-{p}, {r} )" )       

def Box(p,b):
    p = process_vec_arg(p)
    b = process_vec_arg(b)
    # sdBox      ( vec3 p, vec3 b )
    return Formula( f"sdBox( pos-{p}, {b} )" )

def Ellipsoid(p,r):
    p = process_vec_arg(p)
    r = process_vec_arg(r)
    # sdEllipsoid( in vec3 p, in vec3 r )
    return Formula( f"sdEllipsoid( pos-{p}, {r} )" )

def RoundBox(p,b,r):
    p = process_vec_arg(p)
    b = process_vec_arg(b)
    # udRoundBox ( vec3 p, vec3 b, float r )
    return Formula( f"udRoundBox( pos-{p}, {b}, {r} )" )

def Torus(p,t):
    p = process_vec_arg(p)
    t = process_vec_arg(t)
    # sdTorus    ( vec3 p, vec2 t )
    return Formula( f"sdTorus( pos-{p}, {t} )" )



if __name__ == '__main__':
    # Example usage
    #a = Formula('a')
    #b = Formula('b')
    #c = Formula('c')

    #c = Sphere( [1.3,0.5,.2], 1.5 )
    #c = Sphere( 2, 1.5 )

    #print( c.eval_CGS() )

    S1 = Sphere( [1.3,0.5,.2], 1.5 )
    S2 = Sphere( [0.3,0.5,2.2], 0.5 )
    B1 = RoundBox( [-0.3,0.5,-2.2], [0.5,0.5,0.5], 1.5 )

    scene = (B1 - S1) + S2
    print( scene.eval_CGS() )

    # Build a formula: (a + b) * c
    #formula = (a + b) - c

    # Evaluate and print the formula
    #print(formula.eval())  # Outputs: ((a + b) * c)
    #print(formula.eval_CGS())