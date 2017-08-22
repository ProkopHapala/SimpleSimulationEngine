#!/usr/bin/python

import ExpressionParser as ep

def cleanComment( s, cCOMMNET="#" ):
    return "".join( [ l.split(cCOMMNET)[0].strip() for l in s.split("\n") ] ) 

def splitCommands( code, cCOMMNET="#", cNEXT=";" ):
    return [ cleanComment(s, cCOMMNET) for s in code.split( cNEXT ) ]

FUNTIONS = {
#  name       return type
"Sphere"    : ( "SDF",  ( "vec3", "float" ) ),
"Material"  : ( "Mat",  ( "vec3", "float" ) ),
"Object"    : ( "None", ( "SDF",  "Mat"   ) ),
}

VARS = []
# clear this by      del lst[:]     see https://stackoverflow.com/questions/1400608/how-to-empty-a-list-in-python

OPERATORS = { '+':"",'-':"",'*':""}

def checkFloatConst( s, cSEP="," ):
    print " ============= checkFloatConst"
    ss = s.split(cSEP)
    print ss
    vec = [ float(si) for si in ss ]
    print "vec = ", vec
    return  vec

def checkVar( op, st, level ):
    st = st.strip(" ,\t\n")
    if not ( st in VARS ):
        try:
            checkFloatConst(st)
        except Exception as error:
            print('caught this error: ' + repr(error))
            raise Exception('unknown variable: %s ' %st )

def checkFunc( op, st, level ):
    print "\t"*level,op,st
    st = st.strip()
    if st!="":
        if not ( st in FUNTIONS ):
            raise Exception('unknown function: %s ' %st )

def empty( op, st, level ):
    print "\t"*level,op,st
    return

def parse_expr( s ):
    iend,tree = ep.tokenTree( s, 0, len(s), 0 )
    print tree
    #print "::::: ", ep.tokenTree2string(s,tree, 0)
    print ep.tokenTree2string(s,tree, 0, indent=True)
    ep.walkTree_2( s, tree, 0, checkVar, empty, checkFunc )
    exit()

def parse_command( s, cASSIGN="=" ):
    print s
    sides = s.split( cASSIGN )
    if len(sides)>1:
        var_name = sides[0].strip
        if var_name in VARS:
            print "ERROR: redefinition of variable !"
            return None
        else:
            typ, expr = parse_expr( sides[1] )
            VARS[var_name] = typ
            return ( var_name, expr )
    else:
        typ, rhs = parse_expr( sides[0] )
        if typ is not None:
            print "ERROR: return value not stored !"
            return None
        return ( None, expr )

if __name__ == "__main__":
    s = '''
    s1 = Sphere( (1.0, 2.0, # some comment 
    3.0), 1.0 );         #  Sphere( vec3 pos, float R )
    s2 = Sphere( (3.0, 4.0, 5.0), 0.5 );
    s3 = Sphere( (3.0, 4.0, 5.0), 0.5 );
    m1 = Material( (0.5,1.0,0.5), 1.0, 1000.0 ); # Material( vec3 rgb, float cSpecular, float gloss )
    #m2 = Material( (0.5,0.5,0.5), 1.0, 15.0   );
    Object( s1 - ( s2 * s3 ), m1 );              # Object( sdf, Material m )
    '''
    
    print s
    print "=========="
    
    commands = splitCommands(s)
    print "\n".join(commands)
    print "=========="
    parse_command( commands[0] )

    
