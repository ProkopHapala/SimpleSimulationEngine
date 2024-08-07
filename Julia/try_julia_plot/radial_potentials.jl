
using Plots
include("plot_utils.jl")


Coulomb_const_eVA = 14.3996448915 

# ========== Functions

# ============ fast r^2 based polynomial potentials ============

function getRr3( r, R0, E0, Rc, Rf, K )
    """
        This function was developed to provide a smooth Energy similar to Lenard-Jones and Morse potentials but with finite cutoff other nice properties
         - has a minimum at R0
         - it is exactly parabolic around minimum R0 and for r<Rf with stiffness K
         - it has a cubic tail with finite cutoff at Rc
         - the function is continuous (and has continuous derivative) at Rf and Rc
         - curvature of the can be tuned with K, Rf=R0+dRf and cutoff Rc
         - the tail in the interval Rf<r<Rc is a cubic polynomial a3*x^3 + a2*x^2 matching the quadratic function around minimum
    """
    if r>Rc
        E = 0
        F = 0
    elseif r<Rf
        x = r - R0
        E   = K *    x^2 - E0
        F   = K * -2*x
    else
       
        # --- find coefs a3,a2 of cubic polynomial E(x) = a3*x^3 + a2*x^2
        Ef = K*(Rf-R0)^2 - E0
        Ff = K*2*(Rf-R0)
        xf = Rf - Rc
        a3 = -( 2*Ef - Ff*xf ) / xf^3
        a2 = +( 3*Ef - Ff*xf ) / xf^2

        #if x == 0 println( "a3 = ", a3, " a2 = ", a2 ) end
        #a2 = 0.025; a3 = 0.051

        x = r - Rc
        E   =    a3*x^3  +    a2*x^2
        F   = -3*a3*x^2  + -2*a2*x
    end
    return E,F
end


function getR_ir4r8( r, R0, E0, Rc, Rf, K )
    """
        This function is improved version of getRr3()
        In the region r>Rf the function does not need to evaluate r=sqrt( (pi-pj).norm2() ) because 
        it match the quadratic function around minimum using two basis functions (see getR4() and getR8() ) which depends only on r^2
        Otherwise it has similar properties as getRr3(), that is:
        - quadratic before the infex point Rf
        - finite cutoff at Rc
        - continuous Energy and Force at Rf and Rc
        - tunable stiffness around minimum using K
        - tunable curvature of the tail using K, Rf=R0+dRf and cutoff Rc
        - for speedup coefs a1,a2 can be precalculated and stored in a table
    """
    if r>Rc
        E = 0
        F = 0
    elseif r<Rf
        x = r - R0
        E   = K *    x^2 - E0
        F   = K * -2*x
    else
       
        # --- find coefs a3,a2 of cubic polynomial E(x) = a3*x^3 + a2*x^2
        Ef = K*(Rf-R0)^2 - E0
        Ff = K*2*(Rf-R0)
        a1 =        Rc^4 * ( Ff*( Rc^2 - Rf^2 ) + 8*Ef*Rf  ) / ( 4*Rf*( Rc^2 - Rf^2 )^2 )
        a2 = -0.25* Rc^8 * ( Ff*( Rc^2 - Rf^2 ) + 4*Ef*Rf  ) / (   Rf*( Rc^2 - Rf^2 )^4 )
        #if r == Rf println( "analyt a1 = ", a1, " a2 = ", a2 ) end

        #a1 = -1.8 * -0.5
        #a2 = -3.6 * (1.5)
        #a1 = -1.8 * 0.0
        #a2 = -3.6 * 1.0
        #if r == Rf println( "manual a1 = ", a1, " a2 = ", a2 ) end

        #a2 = 0.025; a3 = 0.051

        u  = 1-(r/Rc)^2
        u2 = u*u
        E   = ( a1  + a2*(   u2 ) ) * u2   
        F   = ( a1  + a2*( 2*u2 ) ) * u * 4*(r/(Rc^2))
        #F   = a1* 4*(r/(Rc^2))*( u )   +   a2* 8*(r/(Rc^2))*( u^3 )
    end
    return E,F
end


function getR4( r, Rc, E0 )
    """
        This simplest radial function with finite cutoff Rc which depends only on r^2, avoiding sqrt() for speedup of partile-particle interactions
    """
    if r>Rc
        E = 0
        F = 0
    else
        r2  = r*r
        u2  = 1-(r2/Rc^2)
        E   = E0 * ( 1-u2*u2 ) - E0
        F   = E0 * -4*(r/(Rc^2))*( u2*u2  ) 
    end
    return E,F
end

function getR4x2( r, Rf, E0, Rc, K )
    """
        This fast radial potential with minumum at around R0 and finite cutoff Rc 
        - for r>Rf it depends only on r^2, avoiding sqrt() for speedup of partile-particle interactions
        - The exact position of minimum is not exactly R0, but it can be fitted to be by tunig K and Rf
        - The minimum and repulsive part of the potential is created by adding a quadratic function with stiffness K around Rf for r<Rf
    """
    if r>Rc
        E = 0
        F = 0
    else
        r2  = r*r
        u2  = 1-(r2/Rc^2)
        E   = E0 * ( 1-u2*u2 ) - E0
        F   = E0 * -4*(r/(Rc^2))*( u2  ) 
        if r<Rf
            u =  Rf - r
            E += K * ( u*u )
            F += K * 2*u
        end
    end
    return E,F
end


function getR8( r, Rc, E0 )
    """
    This simple radial function with finite cutoff Rc which depends only on r^2, avoiding sqrt()
        - with respect to getR4() it has has better tail (more curved, simular to Morse) and approaches zero more smoothly at the cutoff Rc
    """
    if r>Rc
        E = 0
        F = 0
    else
        r2  = r*r
        u2  = 1-(r2/Rc^2)
        u4 = u2*u2
        u8 = u4*u4
        E   = E0 * ( 1-u8 ) - E0
        F   = E0 * -8*(r/(Rc^2))*( u4*u2  ) 
    end
    return E,F
end

function getR8x2( r, Rf, E0, Rc, K )
    """
    This fast radial potential with minumum at around R0 and finite cutoff Rc using getR8()
    - for r>Rf it depends only on r^2, avoiding sqrt() for speedup of partile-particle interactions
    - The exact position of minimum is not exactly R0, but it can be fitted to be by tunig K and Rf
    - The minimum and repulsive part of the potential is created by adding a quadratic function with stiffness K around Rf for r<Rf
    """
    if r>Rc
        E = 0
        F = 0
    else
        r2  = r*r
        u2  = 1-(r2/Rc^2)
        u4 = u2*u2
        u8 = u4*u4
        E   = E0 * ( 1-u8 ) - E0
        F   = E0 * -8*(r/(Rc^2))*( u4*u2  ) 
        if r<Rf
            u = Rf - r
            E += K * ( u*u )
            F += K * 2*u
        end
    end
    return E,F
end


# ============ exponetial based polynomial potentials ============

function getMorse( r, R0, E0, k )
    """
    Morse potential with minimum at R0, depth E0 and stiffness k
    - k~= 1.5-1.7 aproximates to Lennard-Jones potential
    """
    e  = exp( -k*(r-R0) )
    E  = E0 * ( e*e - 2*e )
    F  = E0 * 2*k* ( e*e -   e ) 
    return E,F
end

# ============  Lanard-Jones like 1/r^2 based reciprocal potnetials  ============

function getLJ( r, R0, E0 )
    """
    Lennard-Jones potential with minimum at R0 and depth E0
    """
    u  = R0/r
    u6 = u^6
    E  = E0 *    ( u6*u6 - 2*u6 )
    F  = E0 * 12*( u6*u6 -   u6 )/r 
    return E,F
end

function getLJr4( r, R0, E0 )
    """
    Alternative to Lenard-Jones potential lower degree of reciprocal power of r, therefore lower stiffness
     - ToDo: make poential which is linear combiation of LJ and LJr4 with smooth distance-dependent weigth function 
    """
    u  = R0/r
    u4 = u^4
    E  = E0 *   ( u4*u4 - 2*u4 )
    F  = E0 * 8*( u4*u4 -   u4 )/r 
    return E,F
end

function getLJr2( r, R0, E0 )
    u  = R0/r
    u2 = u^2
    E  = E0 *   ( u2*u2 - 2*u2 )
    F  = E0 * 8*( u2*u2 -   u2 )/r 
    return E,F
end

function getLJx2( r, R0, E0 )
    """
    Modified Lennard-Jones potential with quadratic function K*r^2 for r<R0
    - The stiffness of the quadratic function K is automatically adjusted to match curvature of the LJ part in the minimum (where derivative is equal to zero => also mathcing )
    - Notice that LJ potential does not need derivative for r>R0
    """
    if(r<R0)
        K = E0*72/R0^2
        E = K*(r-R0)^2 - E0
        F = -2*K*(r-R0)
    else
        u  = R0/r
        u6 = u^6
        E  = E0 *    ( u6*u6 - 2*u6 )
        F  = E0 * 12*( u6*u6 -   u6 )/r    
    end
    return E,F
end


function getLJs2( r, R0, E0 )
    """
    Modified Lennard-Jones potential with softer 1/r^4 for r<R0
    - The softer repulsive part is more useful for molecular simulations (avoiding hard collisions and high frequency)
    - The parametrs A,B of A*(u2*u2-2*u2)+B are found to match stiffness of the LJ part in the minimum (where derivative is equal to zero => also mathcing )
    - Notice that it does not need sqrt() only depends on powers of 1/r^2
    """
    u  = R0/r
    u2 = u^2
    if(r<R0)
        A = E0*(57.0/7.0)
        B = A - E0
        E = A * ( u2*u2 - 2*u2 ) + B
        F = A * 8*( u2*u2 -   u2 )/r 
    else
        u6 = u2*u2*u2
        E  = E0 *    ( u6*u6 - 2*u6 )
        F  = E0 * 12*( u6*u6 -   u6 )/r    
    end
    return E,F
end


# ============  Damped Coulomb potential  ============


function getCoulomb( r, Q )
    Q *= Coulomb_const_eVA
    E  = Q/r
    F  = Q/r^2
    return E,F
end

function getCoulomb_dampC2( r, Q, Rdamp )
    Q *= Coulomb_const_eVA
    r_ = sqrt( Rdamp^2 + r^2 )
    E  = Q/r_
    F  = Q*r/(r_^3)
    return E,F
end

function getCoulomb_dampR2( r, Q, Rdamp, ADamp )
    Q *= Coulomb_const_eVA
    if r<Rdamp
        D    = Rdamp - r
        dD   = -1  
        r_  = sqrt( r^2 + ADamp* D*D )
        dr_ =     ( r   + ADamp*dD*D )/r_
    else
        r_  = r
        dr_ = 1
    end
    E  = Q/r_
    F  = Q*dr_/(r_^2)
    return E,F
end

function getCoulomb_dampR3( r, Q, Rdamp, ADamp )
    Q *= Coulomb_const_eVA
    if r<Rdamp
        D    = Rdamp - r
        dD   = -1  
        r_  = sqrt( r^2 +     ADamp* D*D*D )
        dr_ =     ( r   + 1.5*ADamp*dD*D*D )/r_
    else
        r_  = r
        dr_ = 1
    end
    E  = Q/r_
    F  = Q*dr_/(r_^2)
    return E,F
end

function getCoulomb_dampR4( r, Q, Rdamp, ADamp )
    Q *= Coulomb_const_eVA
    if r<Rdamp
        D   = ( 1 - (r/Rdamp)^2 ) *Rdamp
        dD  = (   -2*r/Rdamp^2  ) *Rdamp

        #D   = ( Rdamp^2 - r^2  )*0.5
        #dD  = (        -2*r^2  )*0.5

        r_  = sqrt( r^2 + ADamp* D^2  )
        dr_ =     ( r   + ADamp*dD*D )/r_
    else
        r_  = r
        dr_ = 1
    end
    E  = Q/r_
    F  = Q*dr_/(r_^2)
    return E,F
end

function getCoulomb_dampSS( r, Q, Rdamp, ADamp, Rdam0 )
    Q *= Coulomb_const_eVA
    if r<Rdamp
        if r>Rdamp0
            x   = (r-Rdamp0)/(Rdamp-Rdam0)
            D   =  ADamp*( 1 - x * x * (3-2*x) )
            dD  = -ADamp*( 6*x*(1-x)/(Rdamp-Rdam0) )
        else
            D = ADamp
            dD  = 0
        end
        r_  = sqrt( D^2 + r^2 )
        dr_ = (r +dD*D)/r_ 
        #println( "r,x,D,dD = ", r," ",x," ",D," ",dD )
        #return D,dD
    else
        r_  = r
        dr_ = 1
    end
    E  = Q/r_
    F  = Q*dr_/(r_^2)
    return E,F
end

function getCoulomb_dampInv4( r, Q, ADamp )
    Q *= Coulomb_const_eVA

    ir2 = 1/r^2
    D   =    ir2*ir2
    dD  = -4*ir2*ir2*ir2

    #Ec =  1/r;
    #Fc = -1/(r^2)
    # e     = 1/( 1/f(x) + g(x) ) = f/( 1 + f*g )
    # de/dx = ( f' - f^2*g' )/(1+f*g) = ( f/(1+f*g) )^2  * ( f'/f^2 - g' ) = e^2 * ( f'/f^2 - g' ) 
    e  = 1/(  r +  D*ADamp );
    E = Q*e
    F = E*e*( 1 + r*dD*ADamp )
    return E,F
end    

function getCoulomb_dampInv2( r, Q, ADamp )
    Q *= Coulomb_const_eVA

    ir = 1/r
    D   =    ir*ir
    dD  = -2*ir*ir*ir

    #Ec =  1/r;
    #Fc = -1/(r^2)
    # e     = 1/( 1/f(x) + g(x) ) = f/( 1 + f*g )
    # de/dx = ( f' - f^2*g' )/(1+f*g) = ( f/(1+f*g) )^2  * ( f'/f^2 - g' ) = e^2 * ( f'/f^2 - g' ) 
    e  = 1/(  r +  D*ADamp );
    E = Q*e
    F = E*e*( 1 + dD*ADamp )
    return E,F
end    


function getCoulomb_x2( r, R0, E0, Q, Rf )
    #if r==Rf println("Q = ", Q, Q*Coulomb_const_eVA ) end
    Q *= Coulomb_const_eVA
    K   =  ( Q/Rf + E0 )/(Rf-R0)^2
    if r==Rf println("Q = ", Q, " K = ", K ) end
    if(r<Rf)
        # condition: Ef = Q/Rf = K*(Rf-R0)^2 - E0
        # => K = ( Q/Rf + E0 )/(Rf-R0)^2
        #K   =  ( Q/Rf + E0 )/(Rf-R0)^2
        E   =    K*(r-R0)^2 - E0
        F   = -2*K*(r-R0)
    else
        E = Q/r
        F = Q/r^2
    end
    return E,F
end    



function getCoulomb_x2smooth( r, R0, E0, Q, Rf, K )
    #if r==Rf println("Q = ", Q, " E0 ", E0 ," R0 ", R0, " Rf ", Rf ) end
    Q *= Coulomb_const_eVA

    #K   =  ( Q/Rf + E0 )/(Rf-R0)^2
    #if r==Rf println("Q = ", Q, " K = ", K ) end

    if(r<R0)
        E   =    K*(r-R0)^2 - E0
        F   = -2*K*(r-R0)
    elseif(r<Rf)
        x = (r-R0)/(Rf-R0)

        #smooth-step : https://en.wikipedia.org/wiki/Smoothstep
        s   =  1 - x*x*(3-2*x)
        ds  =  6*x*(1-x)/(Rf-R0)

        #smoother step : https://en.wikipedia.org/wiki/Smoothstep#Variations
        #s  =  1-x*x*x*( (6*x - 15)*x + 10 )
        #ds = 30*(x-1)*(x-1)*x*x /(Rf-R0)

        
        #s = (1-x)*(1-x)*(1-x)*(1-x)*x*x 
        #s = (1-x*x)*(1-x*x)/( 1 + x*x/(0.1) ) 
        s = (1-x*x)/( 1 + x*x/(0.1) ) 
        ds = 0
        
        return s, ds

        #s  = 1-x
        #ds = 1/(Rf-R0)

        #s  = 1-x*x
        #ds = 2*x/(Rf-R0)

        Eh   =    K*(r-R0)^2 - E0  
        Fh   = -2*K*(r-R0)         

        Ec = Q/r
        Fc = Q/r^2

        E   =    Eh*s  + Ec*(1-s)
        F   =  ( Fh*s  + Fc*(1-s)    + ds*Eh - ds*Ec  )
        #return s,ds
    else
        E = Q/r
        F = Q/r^2
    end
    return E,F
end 

function getCoulomb_x2lor( r, R0, E0, Q, Rf )
    #if r==Rf println("Q = ", Q, " E0 ", E0 ," R0 ", R0, " Rf ", Rf ) end
    Q *= Coulomb_const_eVA
    if(r<Rf)
        L  = (Rf-R0)   # lenght of the pseudo-Lorenz interval
        dx = 1/L       # dx/dr = 1/L
        Ef =  Q/Rf     # Energy at the end of the pseudo-Lorenz interval
        Ff = -Q/Rf^2   # Force at the end of the pseudo-Lorenz interval
        dE = E0 + Ef   # total energy difference covered by the pseudo-Lorenz interval
        w  =  2*dE*dx/Ff - 1   # width of the pseudo-Lorenz potential to match the derivative Ff ath the end of the pseudo-Lorenz interval
        if(r>R0)  # pseudo-Lorenz interval 
            x  = (r-R0)*dx
            x2 = x*x
            ix2w1 = 1/(x2*w + 1)
            w =  2*dE*dx/Ff - 1
            e   = -dE*ix2w1
            E   = e*(1-x2)  + Ef
            F   = e*2*x*dx*( (1-x2)*ix2w1*w + 1 )
            #E   = -dE*(1-x2)*ix2w1   + Ef
            #F   = -dE*(1-x2)*2*dx*w*x*ix2w1*ix2w1 - dE*ix2w1*2*x*dx
        else     # parabolic interval    
            # we fit stiffness of parabolic potential K*x^2 to stiffness of the pseudo-Lorenz potential
            K   =  dE*(w + 1) *dx*dx   #*L^2
            E   =    K*(r-R0)^2 - E0
            F   = -2*K*(r-R0)
        end
    else
        E = Q/r
        F = Q/r^2
    end
    return E,F
end 

#=
From Maxima:
==== L1x2

f         :  (1-x^2)/(w*x^2+1)
df  (x=1) : -2/(w+1)
ddf (x=0) : -2*w-2
df{S}     : -2* L^2 *  x      (w+1) 
ddf{S}    : -2* L^2 * (4*L-3)*(w+1)

==== L1x4   -This

f         : (1-x^4)/(w*x^2+1)
df  (x=1) : -4/(w+1)
ddf (x=0) : -2*w
df  {S}   : -2*L^2*x*(w*x^4+2*x^2+w)
ddf {S}   : -2*L^3*( w^2*x^6 + 3*w*x^4 - 3*w^2*x^2 + 6*x^2+w )

==== L0x4   -This

f         : (1-x^2)^2/(w*x^2+1)
df  (x=1) : 0
ddf (x=0) : (-2*w)-4
df  {S}   : 2*L^2 * (x-1) *x* (x+1) * (w*x^2+w+2)
ddf {S}   : 2*L^3 * ( w^2*x^6 + 3*w*x^4 + 3*w^2*x^2 + 6*w*x^2 + 6*x^2-w-2   )

=#

function getCoulomb_x2lor2( r, R0, E0, Q, Rf, K )
    #if r==Rf println("Q = ", Q, " E0 ", E0 ," R0 ", R0, " Rf ", Rf ) end
    Q *= Coulomb_const_eVA
    if(r<Rf)

        if(r>R0)  # pseudo-Lorenz interval 
            
            
            
            Ef =  Q/Rf     # Energy at the end of the pseudo-Lorenz interval
            Ff = -Q/Rf^2   # Force at the end of the pseudo-Lorenz interval
            dE = E0 + Ef   # total energy difference covered by the pseudo-Lorenz interval
            
            # transform variables from r@[R0,Rf] to x@[0,1]
            l  = (Rf-R0)
            dx = 1/l    
            x  = (r-R0)*dx
            Ff *= l   # resale force f=dE/dr to get derivative by dx 

            # System of equations to solve:
            #  derivative F(x=1) :  A1 * -4/(w+1)         = Ff         
            #  value      E(x=0) :  A1 + A2               = dE
            #  stiffness  K(x=0) :  A1*-2*w + A2*-2*(w+2) = K

            # solve the system of equations
            w = -( -4*dE - Ff - K ) / ( 2*dE+Ff )
            #w  = 4.55 
            A1 = Ff*(w+1)/4
            A2 = dE - A1

            #w  = 4.55 
            #A1 = dE
            #A2 = 0 
            #println("w,A1,A2 ", w, " ", A1/dE," ",A2/dE )

            x2 = x*x
            L  = 1/(x2*w+1)
            E1 = 1-x2*x2
            E2 = (1-x2)^2
            F1 = -( (w*x2 + 2)*x2 + w ) 
            F2 =  (x2-1) * (w*x2+w+2)
            E  = Ef - ( A1*E1 + A2*E2 )*L
            F  = (A1*F1 + A2*F2) * dx * 2*L^2 *x

            #E1 = L*(1-x^4)
            #E2 = L*(1-x^2)^2
            #F1 = -2*L^2 *x* (w*x^4 + 2*x^2 + w ) 
            #F2 =  2*L^2 *x* (x^2-1) * (w*x^2+w+2)
            #
            #F = (A1*F1 + A2*F2) * dx



            #E   = -dE*(1-x2)*ix2w1   + Ef
            #F   = -dE*(1-x2)*2*dx*w*x*ix2w1*ix2w1 - dE*ix2w1*2*x*dx
        else     # parabolic interval    
            # we fit stiffness of parabolic potential K*x^2 to stiffness of the pseudo-Lorenz potential
            #K   =  dE*(w + 1) *dx*dx   #*L^2
            E   =    K*(r-R0)^2 - E0
            F   = -2*K*(r-R0)
        end
    else
        E = Q/r
        F = Q/r^2
    end
    return E,F
end 


# ========== Body

# eval_forces = (position, velocity) -> eval_force_and_plot(position,velocity, plt, truss.bonds )

xs = xrange( 0.01, 0.01, 600 )

plt = plot( layout = (2, 1), size=(1000, 1000) )
mins = []

Rc    = 6.0
EvdW  = 0.01
RvdW  = 3.5
Ksr   = 1.7
Rf    = RvdW + 0.25

RHb = 2.0
EHb = 0.7

Q = -0.2*0.3   * 6.0
Rdamp0 = 0.5
Rdamp  = 2.0
Adamp  = 1.0

xlim = [1.0, 10.0]

push!( mins, plot_func( plt,  xs, (x)->getCoulomb(         x,Q     ),               clr=:black  , label="Coulomb"       ,  xlim=xlim ) )
#push!( mins, plot_func( plt, xs, (x)->getCoulomb_dampC2(   x,Q,Rdamp),              clr=:red   , label="Coulomb_C2"   ,  xlim=xlim, dnum=:true ) )
push!( mins, plot_func( plt, xs, (x)->getCoulomb_dampR2(   x,Q,Rdamp, Adamp),       clr=:blue   , label="Coulomb_R2"   ,  xlim=xlim, dnum=:true ) )
push!( mins, plot_func( plt, xs, (x)->getCoulomb_dampR3(   x,Q,Rdamp, Adamp),       clr=:green  , label="Coulomb_R3"   ,  xlim=xlim, dnum=:true ) )
push!( mins, plot_func( plt, xs, (x)->getCoulomb_dampR4(   x,Q,Rdamp, Adamp),       clr=:red    , label="Coulomb_R4"   ,  xlim=xlim, dnum=:true ) )
#push!( mins, plot_func( plt, xs, (x)->getCoulomb_dampSS(   x,Q,Rdamp, Adamp,Rdamp0), clr=:magenta , label="Coulomb_SS"  ,  xlim=xlim, dnum=:true ) )
#push!( mins, plot_func( plt, xs, (x)->getCoulomb_dampInv4( x,Q,Adamp),              clr=:green   , label="Coulomb_Inv4" ,  xlim=xlim, dnum=:true ) )
#push!( mins, plot_func( plt, xs, (x)->getCoulomb_dampInv2( x,Q,Adamp),              clr=:cyan   , label="Coulomb_Inv2" ,  xlim=xlim, dnum=:true ) )

#push!( mins, plot_func( plt, xs, (x)->getLJ(x,RvdW,EvdW),        clr=:black  , label="LJ"     ) )
#push!( mins, plot_func( plt, xs, (x)->getLJx2(x,RvdW,EvdW),      clr=:blue    , label="LJx2"   ) )
#push!( mins, plot_func( plt, xs, (x)->getLJs2(x,RvdW,EvdW),      clr=:red    , label="LJs2"   ) )

#push!( mins, plot_func( plt, xs, (x)->getLJr4(x,RvdW,EvdW),      clr=:green  , label="LJr4"   ) )
#push!( mins, plot_func( plt, xs, (x)->getLJr2(x,RvdW,EvdW),      clr=:red    , label="LJr4"   ) )
#push!( mins, plot_func( plt, xs, (x)->getMorse(x,RvdW,EvdW,1.7), clr=:black , label="Morse"  ) )

#push!( mins, plot_func( plt, xs, (x)->getRr3(     x, RvdW, EvdW, Rc, Rf, Ksr ), clr=:cyan , label="Rr3"  ) )

#push!( mins, plot_func( plt, xs, (x)->getR_ir4r8( x, RvdW, EvdW, Rc, Rf, Ksr ), clr=:magenta , label="R_ir4r8" ,  dnum=:true  ) )

#push!( mins, plot_func( plt, xs, (x)->getR4(x,Rc,E0),        clr=:black   , label="R4"     ) )
#push!( mins, plot_func( plt, xs, (x)->getR8(x,Rc,E0),        clr=:black  , label="R2"     ) )
#push!( mins, plot_func( plt, xs, (x)->getR4x2(x,R0+0.15,E0*1.85,Rc,Ksr),      clr=:red   , label="R4x2"  ,  dnum=:true  ) )
#push!( mins, plot_func( plt, xs, (x)->getR8x2(x,R0+0.17,E0*3.50,Rc,Ksr*2.0),  clr=:green  , label="R2x2" ,  dnum=:true) )


#push!( mins, plot_func( plt, xs, (x)->getCoulomb_dampR2(   x,Q,Rdamp, Adamp),                                                   clr=:blue    , label="Coulomb_R2" ,  xlim=xlim, dnum=:true ) )
#push!( mins, plot_func( plt, xs, (x)->getR_ir4r8( x, RHb, EHb, Rc, RHb+0,1, 10.0 ),                                             clr=:magenta , label="R_ir4r8"    ,  dnum=:true  ) )
#push!( mins, plot_func( plt, xs, (x)->(getR_ir4r8( x,RHb, EHb, Rc, RHb+0,1, 10.0 ) .+ getCoulomb_dampR2(   x,Q,Rdamp, Adamp)),  clr=:cyan    , label="R_ir4r8"    ,  dnum=:true  ) )


#push!( mins, plot_func( plt, xs, (x)->getCoulomb_x2( x, RHb, EHb, Q, RvdW ),                                                    clr=:blue    , label="Coulomb_x2" ,       xlim=xlim, dnum=:true ) )
#push!( mins, plot_func( plt, xs, (x)->getCoulomb_x2smooth( x, RHb, EHb, Q, RHb+0.8, 5.0 ),                                      clr=:blue    , label="Coulomb_x2smooth" , xlim=xlim, dnum=:true ) )
#push!( mins, plot_func( plt, xs, (x)->getCoulomb_x2lor( x, RHb, EHb, Q, RvdW-0.7 ),                                             clr=:blue    , label="Coulomb_x2lor" ,    xlim=xlim, dnum=:true ) )
#push!( mins, plot_func( plt, xs, (x)->getCoulomb_x2lor2( x, RHb, EHb, Q, RHb+0.8, 5.0 ),                                         clr=:blue    , label="Coulomb_x2lor2" ,  xlim=xlim, dnum=:true ) )

Emin = minimum( [min[1] for min in mins] ); ylims!( plt[1], Emin*1.5, -Emin )
Fmin = minimum( [min[2] for min in mins] ); ylims!( plt[2], Fmin*1.5, -Fmin )

#println( "Emin = ", Emin )

hline!( plt[1], [0.0], color=:black, label="", linestyle=:dash )
hline!( plt[2], [0.0], color=:black, label="",linestyle=:dash )

vline!( plt[1], [Rc],  color=:gray,  label="", linestyle=:dash )
vline!( plt[1], [RvdW],  color=:black, label="", linestyle=:dash )
vline!( plt[1], [RHb],  color=:black, label="", linestyle=:dash )
vline!( plt[1], [0.0], color=:black, label="", linestyle=:dash )

vline!( plt[2], [Rc],  color=:gray,  label="", linestyle=:dash )
vline!( plt[2], [RvdW],  color=:black, label="", linestyle=:dash )
vline!( plt[2], [RHb],  color=:black, label="", linestyle=:dash )
vline!( plt[2], [0.0], color=:black, label="", linestyle=:dash )

display(plt)

