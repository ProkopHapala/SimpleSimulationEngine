# Load the essential modules
using LinearAlgebra
#using SparseArrays
#using Plots
#using Printf

include("Truss.jl")

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


function make_PD_Matrix( neighBs::Array{Vector{Int},1}, bonds::Array{Tuple{Int,Int},1}, masses::Array{Float64,1}, dt::Float64, ks::Array{Float64,1} )
    # see:  1.3.3 A Simple Example,                 in https://doi.org/10.1145/3277644.3277779 Parallel iterative solvers for real-time elastic deformations.
    #       resp. eq.14 in the same paper, or eq.14 in https://doi.org/10.1145/2508363.2508406 Fast Simulation of Mass-Spring Systems
    np = length(masses)
    A  = zeros( np, np)
    idt2 = 1. /dt^2
    for i = 1:np
        Aii = 0
        #println( "==i,i ", i," ", Aii," ", masses[i]," ", idt2 )
        for ib in neighBs[i]
            k = ks[ib]
            Aii += k
            (i_,j_) = bonds[ib]
            #println( "i,ib ",i," ",ib," ",  k  )
            # I think there is error in the paper https://doi.org/10.1145/3277644.3277779 Parallel iterative solvers for real-time elastic deformations. Fratarcangeli, M., Wang, H., & Yang, Y. (2018).  SIGGRAPH Asia 2018 Courses, 1–45.
            #  the off-diagonal elements should be -k, not +k
            if     ( j_ > i )
                A[i,j_] = -k
                A[j_,i] = -k
            elseif ( i_ > i )
                A[i,i_] = -k
                A[i_,i] = -k
            end 
        end
        A[i,i] += Aii
    end
    Mt  = zeros( np)
    for i = 1:np
        mti     = masses[i] * idt2
        Mt[i]   = mti
        A[i,i] += mti
    end
    #println( "make_PDMatrix() DONE !!!!!!!!!!!!!!!!!!!!! " )
    return A, Mt
end 

function make_PD_Matrix_bak( neighBs::Array{Vector{Int},1}, bonds::Array{Tuple{Int,Int},1}, masses::Array{Float64,1}, dt::Float64, ks::Array{Float64,1} )
    # see:  1.3.3 A Simple Example,                 in https://doi.org/10.1145/3277644.3277779 Parallel iterative solvers for real-time elastic deformations.
    #       resp. eq.14 in the same paper, or eq.14 in https://doi.org/10.1145/2508363.2508406 Fast Simulation of Mass-Spring Systems
    np = length(masses)
    A = zeros( np, np)
    idt2 = 1. /dt^2
    for i = 1:np
        Aii = masses[i] * idt2
        #println( "==i,i ", i," ", Aii," ", masses[i]," ", idt2 )
        for ib in neighBs[i]
            k = ks[ib]
            Aii += k
            (i_,j_) = bonds[ib]
            #println( "i,ib ",i," ",ib," ",  k  )
            # I think there is error in the paper https://doi.org/10.1145/3277644.3277779 Parallel iterative solvers for real-time elastic deformations. Fratarcangeli, M., Wang, H., & Yang, Y. (2018).  SIGGRAPH Asia 2018 Courses, 1–45.
            #  the off-diagonal elements should be -k, not +k
            if     ( j_ > i )
                A[i,j_] = -k
                A[j_,i] = -k
            elseif ( i_ > i )
                A[i,i_] = -k
                A[i_,i] = -k
            end 
        end
        A[i,i] += Aii
    end
    #println( "make_PDMatrix() DONE !!!!!!!!!!!!!!!!!!!!! " )
    return A
end 

function make_PD_rhs_decomp( neighBs::Array{Vector{Int},1}, bonds::Array{Tuple{Int,Int},1}, masses::Array{Float64,1}, dt::Float64, ks::Array{Float64,1}, points::Matrix{Float64}, l0s::Array{Float64,1}, pnew::Matrix{Float64} )
    # see:  1.3.3 A Simple Example,                 in https://doi.org/10.1145/3277644.3277779 Parallel iterative solvers for real-time elastic deformations.
    #       resp. eq.14 in the same paper, or eq.14 in https://doi.org/10.1145/2508363.2508406 Fast Simulation of Mass-Spring Systems
    np = length(masses)
    bK = zeros( np, 3)
    bM = zeros( np, 3)
    #print( "b : "); display(b)
    idt2 = 1. /dt^2
    for i = 1:np
        for ib in neighBs[i]
            (i_,j_) = bonds[ib]
            if i_ == i
                j = j_
            else
                j = i_
            end
            k = ks[ib]
            d = points[i,:] - points[j,:]
            d *= k * l0s[ib] / norm(d)
            bK[i,:] += d        
        end
    end
    for i = 1:np
        bM[i,:] =  pnew[i,:] * (masses[i] * idt2)
    end
    #print( "b : "); display(b)
    return bK, bM
end


function make_PD_rhs( neighBs::Array{Vector{Int},1}, bonds::Array{Tuple{Int,Int},1}, masses::Array{Float64,1}, dt::Float64, ks::Array{Float64,1}, points::Matrix{Float64}, l0s::Array{Float64,1}, pnew::Matrix{Float64} )
    # see:  1.3.3 A Simple Example,                 in https://doi.org/10.1145/3277644.3277779 Parallel iterative solvers for real-time elastic deformations.
    #       resp. eq.14 in the same paper, or eq.14 in https://doi.org/10.1145/2508363.2508406 Fast Simulation of Mass-Spring Systems
    np = length(masses)
    b = zeros( np, 3)
    #print( "b : "); display(b)
    idt2 = 1. /dt^2
    for i = 1:np
        bi =  pnew[i,:] * (masses[i] * idt2)
        ni = 0
        for ib in neighBs[i]
            k = ks[ib]
            (i_,j_) = bonds[ib]
            if i_ == i
                j = j_
            else
                j = i_
            end
            println( "rhs[",i-1,",ib=",ib-1,",j=",j-1,"] k=",k," l0=",l0s[ib] );
            #d = points[i,:] - points[j,:]   # using old positions seems to limit the propagation
            d = pnew[i,:] - pnew[j,:]        # using update positions seems to work much better !!!
            d *= k * l0s[ib] / norm(d)
            bi += d
            ni += 1        
        end
        println( "rhs[",i-1,"] ni=",ni," bi=",bi );
        b[i,:] = bi
    end
    #print( "b : "); display(b)
    return b
end

function make_PD_rhs_f( neighBs::Array{Vector{Int},1}, bonds::Array{Tuple{Int,Int},1}, masses::Array{Float32,1}, dt::Float32, ks::Array{Float32,1}, points::Matrix{Float32}, l0s::Array{Float32,1}, pnew::Matrix{Float32} )
    # see:  1.3.3 A Simple Example,                 in https://doi.org/10.1145/3277644.3277779 Parallel iterative solvers for real-time elastic deformations.
    #       resp. eq.14 in the same paper, or eq.14 in https://doi.org/10.1145/2508363.2508406 Fast Simulation of Mass-Spring Systems
    np = length(masses)
    b = zeros( Float32, np, 3 )
    #print( "b : "); display(b)
    idt2 = 1.f0 /dt^2
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
            #d = points[i,:] - points[j,:]   # using old positions seems to limit the propagation
            d = pnew[i,:] - pnew[j,:]        # using update positions seems to work much better !!!
            d *= k * l0s[ib] / norm(d)
            bi += d        
        end
        b[i,:] = bi
    end
    #print( "b : "); display(b)
    return b
end

struct Truss
    points::Matrix{Float64}
    bonds::Array{Tuple{Int, Int}, 1}
    masses::Array{Float64, 1}
    ks::Array{Float64, 1}
    l0s::Array{Float64, 1}
    fixed::Vector{Int}
    neighBs::Array{Vector{Int}, 1}
end

struct TrussSolution
    A::Matrix{Float64}
    L::Matrix{Float64}
    U::Matrix{Float64}
    neighsL::Vector{Vector{Int}}
    neighsU::Vector{Vector{Int}}
    LDLT_L::Matrix{Float64}
    LDLT_D::Vector{Float64}
    neighsLDLT::Vector{Vector{Int}}
end


struct LDLTsolution
    L::Matrix{Float64}
    D::Vector{Float64}
    neighs::Vector{Vector{Int}}
end

struct LDLTsolution_f
    L::Matrix{Float32}
    D::Vector{Float32}
    neighs::Vector{Vector{Int}}
end


struct Truss_f
    points::Matrix{Float32}
    bonds::Array{Tuple{Int, Int}, 1}
    masses::Array{Float32, 1}
    ks::Array{Float32, 1}
    l0s::Array{Float32, 1}
    fixed::Vector{Int}
    neighBs::Array{Vector{Int}, 1}
end

struct TrussSolution_f
    A::Matrix{Float32}
    L::Matrix{Float32}
    U::Matrix{Float32}
    neighsL::Vector{Vector{Int}}
    neighsU::Vector{Vector{Int}}
    LDLT_L::Matrix{Float32}
    LDLT_D::Vector{Float32}
    neighsLDLT::Vector{Vector{Int}}
end

function convert_to_truss_f(truss::Truss)
    return Truss_f(
        convert.(Float32, truss.points), 
        truss.bonds,  
        convert.(Float32, truss.masses), 
        convert.(Float32, truss.ks),  
        convert.(Float32, truss.l0s), 
        truss.fixed,  
        truss.neighBs 
    )
end

function convert_to_TrussSolution_f(sol::TrussSolution)
    return TrussSolution_f(
        convert.(Float32, sol.A), 
        convert.(Float32, sol.L),  
        convert.(Float32, sol.U), 
        sol.neighsL,  
        sol.neighsL,
        convert.(Float32, sol.LDLT_L),  
        convert.(Float32, sol.LDLT_D), 
        sol.neighsLDLT,
    )
end

function convert_to_LDLTsolution_f(sol::LDLTsolution)
    return LDLTsolution_f(
        convert.(Float32, sol.L),  
        convert.(Float32, sol.D), 
        sol.neighs,
    )
end

#===
function convert_Truss(T::Type, truss)
    return T(
        convert(Matrix{T}, truss.points),
        truss.bonds,
        convert(Vector{T}, truss.masses),
        convert(Vector{T}, truss.ks),
        convert(Vector{T}, truss.l0s),
        truss.fixed,
        truss.neighBs
    )
end

function convert_TrussSolution(T::Type, sol)
    return T(
        convert(Vector{T}, sol.A),
        convert(Vector{T}, sol.L),
        convert(Vector{T}, sol.U),
        sol.neighsU,
        sol.neighsU
    )
end
==#

function update_velocity_old( ps_cor::Matrix{Float64}, ps0::Matrix{Float64}, dt::Float64 )
    return (ps_cor .- ps0) / dt
end

function update_velocity( ps_cor::Matrix{Float64}, ps0::Matrix{Float64}, vel0::Matrix{Float64}, dt::Float64 )
    vel    = (ps_cor .- ps0) / dt
    #renorm = vecnorm(vel0, dims=2) / vecnorm(vel, dims=2)
    #vr  = sqrt.(  sum(abs2, vel,  dims=2) ); println( "vr  ", vr  )
    #vr0 = sqrt.(  sum(abs2, vel0, dims=2) ); println( "vr0 ", vr0 )
    #renorm = vr0 ./ vr
    renorm = sqrt.( sum(abs2, vel0, dims=2) ./ sum(abs2, vel, dims=2) )
    #println( "renorm ", renorm)
    vel   .*= renorm
    return vel
end

function run_solver_bak( truss::Truss, sol::LDLTsolution, velocity::Matrix{Float64}, eval_forces::Function; dt::Float64=0.1, niter::Int=100 ) 
    points = copy( truss.points )
    ps_cor = copy(points)
    for i=1:niter
        force                    = eval_forces( points, velocity )
        ps_pred                  = points .+ velocity*dt .+ force*(dt^2)
        ps_pred[truss.fixed,:]  .= truss.points[truss.fixed,:]
        b                        = make_PD_rhs( truss.neighBs, truss.bonds, truss.masses, dt, truss.ks, points, truss.l0s, ps_pred )  # ;print( "b : "); display(b)

        # ---- Method 2)
        #y      = forwardsub_sparse(sol.L,b,sol.neighsL)  #;print("y_ch : "); display(y_ch)
        #ps_cor = backsub_sparse(   sol.U,y,sol.neighsU)  #;print("x_ch : "); display(x_ch)

        ps_cor[:,1] = solve_LDLT( sol.L, sol.D, b[:,1] )
        ps_cor[:,2] = solve_LDLT( sol.L, sol.D, b[:,2] )
        ps_cor[:,3] = solve_LDLT( sol.L, sol.D, b[:,3] )

        #ps_cor[:,1] = solve_LDLT_sparse( sol.L, sol.D, sol.neighs, b[:,1] )
        #ps_cor[:,2] = solve_LDLT_sparse( sol.L, sol.D, sol.neighs, b[:,2] )
        #ps_cor[:,3] = solve_LDLT_sparse( sol.L, sol.D, sol.neighs, b[:,3] )

        #ps_cor = sol.A \ b

        
        # ---- Method 3)   using  forward-and-backsubstitution with Our Cholensky factorization
        #y    = forwardsub(L,b)
        #ps_cor = backsub(   U,y)
        # ---- Method 4)   using  forward-and-backsubstitution with Our Cholensky factorization from Julia's native algorithm
        #L = Matrix(AFact.L)
        #U = Matrix(AFact.U)
        #y    = forwardsub(L,b)
        #ps_cor = backsub(   U,y)
        # ---- Method 5b)   using factorization by Julia's native algorithm
        #y    = L \ b
        #y    = forwardsub(L,b)
        #ps_cor = U \ y
        #ps_cor    = forwardsub(U,y)
        #ps_cor = backsub(   U,y)

        # ---- Method 5)   using factorization by Julia's native algorithm
        #y    = AFact.L \ b
        #ps_cor = AFact.U \ y

        # ---- Method 6)   multiply by inverse matrix
        #ps_cor = AFact \ b

        # ---- Method 8)   multiply by inverse matrix
        #ps_cor = invA * b 

        # ---- Method 9)   direct sover
        #ps_cor = sol.A \ b

        #ps_cor = copy(ps_pred)

        # ---- residual
        res    = maximum( abs.(ps_cor-points) )   ;println("residual[$i] : ", res );
        # ---- update        
        velocity = update_velocity( ps_cor, points, velocity, dt )
        points[:,:] .= ps_cor[:,:]
        
        #   ToDo:   We need to update velocity and forces based on position update

        #points = points .+ x
        #plot!( plt, [points0[i,1],points[i,1]], [points0[i,2],points[i,2]], color=:blue, lw=0.5 )

        #plot_truss( plt, bonds, points, lw=1.0, c=:blue )
        #plot_truss( plt, bonds, points0, lw=1.0, c=:blue )
    end
    return points
end


#====== To improve numerical stability in Float32 use LDL^T ( root-free Cholesky decomposition) 

https://en.wikipedia.org/wiki/Cholesky_decomposition#LDL_decomposition

Alternative approach using two forward substitutions:
You're correct that there's a method to replace backward substitution with another forward substitution. This approach is sometimes called the "LDL^T factorization" or "root-free Cholesky decomposition". Here's how it works:
a. Instead of decomposing A into LL^T, decompose it into LDL^T, where L is lower triangular with 1's on the diagonal, and D is diagonal.
b. Solve the system in three steps:

Solve Lz = b (forward substitution)
Solve Dy = z (simple division)
Solve L^Tx = y (can be rewritten as a forward substitution)

===#

function run_solver( truss::Truss, sol::LDLTsolution, velocity::Matrix{Float64}, eval_forces::Function; dt::Float64=0.1, niter::Int=100, A_check::Matrix{Float64}, bRes::Bool=:true ) 
    n = size(truss.points,1)
    points  = copy( truss.points )
    ps_cor  = copy( points )
    ps_pred = copy( points )
    b       = copy( points )
    for iter=1:niter
        force                    = eval_forces( points, velocity )
        ps_pred[:,:]             = points .+ velocity*dt #.+ force*(dt^2)
        ps_pred[truss.fixed,:]  .= truss.points[truss.fixed,:]
        b[:,:]                   = make_PD_rhs( truss.neighBs, truss.bonds, truss.masses, dt, truss.ks, points, truss.l0s, ps_pred )  # ;print( "b : "); display(b)
        
        #for i=1:n
        #    #printf( "ps_pred[%i](%10.6f,%10.6f,%10.6f) v(%10.6f,%10.6f,%10.6f) p(%10.6f,%10.6f,%10.6f) dt=%g \n", ps_pred[i,1],ps_pred[i,2],ps_pred[i,3], velocity[i,1],velocity[i,2],velocity[i,3], points[i,1],points[i,2],points[i,3], dt );
        #    println( "ps_pred[",i,"]",ps_pred[i,:]," v", velocity[i,:],"  p", points[i,:]," dt= ", dt );
        #end

        #println( "run_solver().b=", b )

        ps_cor[:,1] = solve_LDLT( sol.L, sol.D, b[:,1] )
        ps_cor[:,2] = solve_LDLT( sol.L, sol.D, b[:,2] )
        ps_cor[:,3] = solve_LDLT( sol.L, sol.D, b[:,3] )
        #ps_cor[:,1] = solve_LDLT_sparse( sol.L, sol.D, sol.neighs, b[:,1] )
        #ps_cor[:,2] = solve_LDLT_sparse( sol.L, sol.D, sol.neighs, b[:,2] )
        #ps_cor[:,3] = solve_LDLT_sparse( sol.L, sol.D, sol.neighs, b[:,3] )
        # ---- residual
        if bRes
            #err = ps_cor - points
            #err = ps_cor - ps_pred
            #err = A_check*ps_cor - b     # check matrix solution
            err = process_bonds(truss.bonds,ps_cor)[2] - truss.l0s  # check bond lenghs 
            res = maximum( abs.(err) )
            println("residual[$iter] : ", res );
        end
        # ---- update        
        #velocity .+= update_velocity( ps_pred, ps_cor, dt )
        velocity = update_velocity( ps_cor, points, velocity, dt )
        points[:,:] .= ps_cor[:,:]
    end

    tmps = (ps_pred,ps_cor,b)
    return points, tmps
end

function run_solver_f( truss::Truss_f, sol::LDLTsolution_f, velocity::Matrix{Float32}, eval_forces::Function; dt::Float32=0.1f0, niter::Int=100 ) 
    points = copy( truss.points )
    ps_cor = copy( points )
    for i=1:niter
        force                    = eval_forces( points, velocity )
        ps_pred                  = points .+ velocity*dt .+ force*(dt^2)
        ps_pred[truss.fixed,:]  .= truss.points[truss.fixed,:]
        b                        = make_PD_rhs_f( truss.neighBs, truss.bonds, truss.masses, dt, truss.ks, points, truss.l0s, ps_pred )  # ;print( "b : "); display(b)
        
        println( "run_solver_f().b=", b )

        ps_cor[:,1] = solve_LDLT( sol.L, sol.D, b[:,1] )
        ps_cor[:,2] = solve_LDLT( sol.L, sol.D, b[:,2] )
        ps_cor[:,3] = solve_LDLT( sol.L, sol.D, b[:,3] )
        # ps_cor[:,1] = solve_LDLT_sparse( sol.L, sol.D, sol.neighs, b[:,1] )
        # ps_cor[:,2] = solve_LDLT_sparse( sol.L, sol.D, sol.neighs, b[:,2] )
        # ps_cor[:,3] = solve_LDLT_sparse( sol.L, sol.D, sol.neighs, b[:,3] )
        # ---- residual
        res    = maximum( abs.(ps_cor-points) )   ;println("residual[$i] : ", res );
        # ---- update        
        #velocity .+= update_velocity( ps_pred, ps_cor, dt )
        #velocity = update_velocity( ps_pred, points, velocity, dt )
        points[:,:] .= ps_cor[:,:]
    end
    return points
end






# ======= Bulding system as sparse solver