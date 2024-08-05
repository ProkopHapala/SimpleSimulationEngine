
using Plots
include("plot_utils.jl")

# ========== Functions


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



# ========== Body

# eval_forces = (position, velocity) -> eval_force_and_plot(position,velocity, plt, truss.bonds )

xs = xrange( 2.5, 0.01, 600 )

plt = plot( layout = (2, 1), size=(1000, 1000) )
mins = []

E0  = 1.0
Rc  = 6.0
R0  = 3.0
Ksr = 1.7
Rf  = R0 + 0.25



push!( mins, plot_func( plt, xs, (x)->getLJ(x,R0,E0),        clr=:black  , label="LJ"     ) )
push!( mins, plot_func( plt, xs, (x)->getLJx2(x,R0,E0),      clr=:blue    , label="LJx2"   ) )

push!( mins, plot_func( plt, xs, (x)->getLJs2(x,R0,E0),      clr=:red    , label="LJs2"   ) )

#push!( mins, plot_func( plt, xs, (x)->getLJr4(x,R0,E0),      clr=:green  , label="LJr4"   ) )
#push!( mins, plot_func( plt, xs, (x)->getLJr2(x,R0,E0),      clr=:red    , label="LJr4"   ) )
#push!( mins, plot_func( plt, xs, (x)->getMorse(x,R0,E0,1.7), clr=:black , label="Morse"  ) )

#push!( mins, plot_func( plt, xs, (x)->getRr3(     x, R0, E0, Rc, Rf, Ksr ), clr=:cyan , label="Rr3"  ) )

#push!( mins, plot_func( plt, xs, (x)->getR_ir4r8( x, R0, E0, Rc, Rf, Ksr ), clr=:magenta , label="R_ir4r8" ,  dnum=:true  ) )

#push!( mins, plot_func( plt, xs, (x)->getR4(x,Rc,E0),        clr=:black   , label="R4"     ) )
#push!( mins, plot_func( plt, xs, (x)->getR8(x,Rc,E0),        clr=:black  , label="R2"     ) )
#push!( mins, plot_func( plt, xs, (x)->getR4x2(x,R0+0.15,E0*1.85,Rc,Ksr),      clr=:red   , label="R4x2"  ,  dnum=:true  ) )
#push!( mins, plot_func( plt, xs, (x)->getR8x2(x,R0+0.17,E0*3.50,Rc,Ksr*2.0),  clr=:green  , label="R2x2" ,  dnum=:true) )

Emin = minimum( [min[1] for min in mins] ); ylims!( plt[1], Emin*1.1, -Emin )
Fmin = minimum( [min[2] for min in mins] ); ylims!( plt[2], Fmin*1.1, -Fmin )

#println( "Emin = ", Emin )

hline!( plt[1], [0.0], color=:black, label="", linestyle=:dash )
hline!( plt[2], [0.0], color=:black, label="",linestyle=:dash )

vline!( plt[1], [Rc],  color=:gray,  label="", linestyle=:dash )
vline!( plt[1], [R0],  color=:black, label="", linestyle=:dash )
vline!( plt[1], [0.0], color=:black, label="", linestyle=:dash )

vline!( plt[2], [Rc],  color=:gray,  label="", linestyle=:dash )
vline!( plt[2], [R0],  color=:black, label="", linestyle=:dash )
vline!( plt[2], [0.0], color=:black, label="", linestyle=:dash )

display(plt)

