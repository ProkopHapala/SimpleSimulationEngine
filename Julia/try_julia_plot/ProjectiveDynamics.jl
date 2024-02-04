# Load the essential modules
#using LinearAlgebra
#using SparseArrays
#using Plots

"""
# NOTES_TO: https://doi.org/10.1145/2508363.2508406 Fast Simulation of Mass-Spring Systems, Liu, T., Bargteil, A. W., O’Brien, J. F., & Kavan, L. (2013). Fast simulation of mass-spring systems. ACM Transactions on Graphics, 32(6), 1–7.

to solve position of point q at discretized times ( h = dt is timestep ):
q_n+1 − 2*q_n + q_n−1 = h^2 * M^−1 f( q_n+1 )      eq.5
We substitute known positions as y and unkonwn future position as x:
x := qn+1
y := 2qn − qn−1
M(x−y) = h^2 * f(x)                                eq.7
this is solution to following minimization problem (  f = -Grad.E ):
g(x) = (1/2) (x−y)^T M (x−y) + h^2 * E(x)            eq.8
Grad g = 0 = M(x−y) + h^2 Grad.E(x)  =  M(x−y) - h^2 * f

Now for special case of linear springs
E = (1/2) Sum_i{ k_i*||x_i − x_j - l0_i ||^2 } - f_ext 
is solution of 
E = (1/2) Sum_i{ k_i*||x_i − x_j - d_i  ||^2 } - f_ext 
this follows from triangle inequality ( ||x_i − x_j - l0_i || <= ||x_i − x_j - d_i || )
E = min. (1/2) x^T*L*x - x^T*J*d + x^T f_ext         eq.13
L = Sum_i{ k_i * a_i * a_i^T } o I3x3         eq.12
J = Sum_i{ k_i * a_i * s_i^T } o I3x3         eq.12
where    
A_i is like ( 0,0,0,1,0,0,-1,0,0 ), selectes points belonging to the bond
S_i,j = delta_ij selects only this bond

By combination of eq.8 and eq.13 we get:
g(x) = (1/2) x^T *( M + h^2*L )* x  - x^T J d + x^T b
where b = f_ext - M * y     ( M*y is inertia resp. discretized momentum   m*v )

See also: Eq(14) https://doi.org/10.1145/3277644.3277779 Parallel iterative solvers for real-time elastic deformations. Fratarcangeli, M., Wang, H., & Yang, Y. (2018).  SIGGRAPH Asia 2018 Courses, 1–45.

eq.14   A*q=b
A = M/h^2 + Sum_c { w_c A_c^T A_c }
b = M/h^2 * s_t * q_n + Sum_c { w_c A_c^T A_c^T p_c }

s_t = q_t + h*v_t + h^2 M^−1 f_ext   :  is the expected position vector without internal forces.
p_c                                  :  is the projection of q into the energy-free space defined by c,

# NOTES_TO: As-Rigid-As-Possible Surface Modeling, Sorkine & Alexa, http://dx.doi.org/10.2312/SGP/SGP07/109-116

"""

# =========== Functions

# ======= Bulding system as dense matrix  

function make_PDMatrix( neighBs::Array{Vector{Float64},1}, bonds::Array{Tuple{Int,Int},1}, masses::Array{Float64,1}, dt::Float64, ks::Array{Float64,1} )
    # see:  1.3.3 A Simple Example,                 in https://doi.org/10.1145/3277644.3277779 Parallel iterative solvers for real-time elastic deformations.
    #       resp. eq.14 in the same paper, or eq.14 in https://doi.org/10.1145/2508363.2508406 Fast Simulation of Mass-Spring Systems
    np = length(masses)
    A = zeros( np, np)
    idt2 = 1. /dt^2
    for i = 1:np
        Aii = masses[i] * idt2
        for ib in neighBs[i]
            k = ks[ib]
            Aii += k
            (i_,j_) = bonds[ib]
            if     ( j_ > i )
                A[i,j_] = k
                A[j_,i] = k
            elseif ( i_ > i )
                A[i,i_] = k
                A[i_,i] = k
            end 
        end
        A[i,i] += Aii
    end
    return A
end 

function make_PD_rhs( neighBs::Array{Vector{Float64},1}, bonds::Array{Tuple{Int,Int},1}, masses::Array{Float64,1}, dt::Float64, ks::Array{Float64,1}, points::Array{Float64,2}, l0s::Array{Float64,1}, pnew::Array{Float64,2} )
    # see:  1.3.3 A Simple Example,                 in https://doi.org/10.1145/3277644.3277779 Parallel iterative solvers for real-time elastic deformations.
    #       resp. eq.14 in the same paper, or eq.14 in https://doi.org/10.1145/2508363.2508406 Fast Simulation of Mass-Spring Systems
    np = length(masses)
    b = zeros( np, 3)
    #print( "b : "); display(b)
    idt2 = 1. /dt^2
    for i = 1:np
        bi =  pnew[i,:] * (masses[i] * idt2)
        for ib in neighBs[i]
            k = ks[ib]
            (i_,j_) = bonds[ib]
            if i_ == i
                j = j_
            else
                j = i_
            end
            d = points[i,:] - points[j,:]
            d *= k * l0s[ib] / norm(d)
            bi += d        
        end
        b[i,:] = bi
    end
    #print( "b : "); display(b)
    return b
end


# ======= Bulding system as sparse solver