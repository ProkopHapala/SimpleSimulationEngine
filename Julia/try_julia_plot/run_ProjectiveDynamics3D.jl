using Random
using Revise # prevent problem with redefinition of types at runtipme (when using include() in REPL)

include("ProjectedGuassSeidel.jl")
include("plot_utils.jl")
include("IterativeSolvers.jl")
include("ProjectiveDynamics.jl")

# =========== Data

# struct Truss
#     points::Array{Float64, 2}
#     bonds::Array{Tuple{Int, Int}, 1}
#     masses::Array{Float64, 1}
#     ks::Array{Float64, 1}
#     l0s::Array{Float64, 1}
#     neighBs::Array{Vector{Int}, 1}
# end

# struct TrussSolution
#     A::Array{Float64,2}
#     U::Array{Float64,2}
#     L::Array{Float64,2}
#     neighsL::Vector{Vector{Int}}
#     neighsU::Vector{Vector{Int}}
# end

function eval_forces(position, velocity)
    force = -0.0 * position - 0.00 * velocity
    return force
end

###################
###    main()   ###
###################

function main( nx=5, ny=5, dt=5.0 )
    p0 =  Vector([0.0,0.0,0.0])
    p1 =  Vector([1.0,0.0,0.0])
    ax =  Vector([0.0,0.0,1.0])
    print("wheel()");         @time points, bonds, ks = wheel( p0, p1, ax, 8, 0.2 )

    #println("typeof(points) ", typeof(points) );

    print("process_bonds()"); @time   _, l0s          = process_bonds( bonds, points );
    #print("bonds :");  display(bonds)
    #print("points :"); display(points)
    #print("masses :"); display(masses)

    np = size(points,1)
    masses = ones(np)
    fixed  = []
    #masses[fixed] .= 100000.0

    print("point_neighs()");   @time neighs  = point_neighs(  bonds, np )     #;print("neighs : "); display(neighs)      # finds neighbors
    print("point_neighBs()");  @time neighBs = point_neighBs( bonds, np )     #;                                         # finds bonds to neighbors
    print("point_neighs2()");  @time neighs2 = point_neighs2( neighs    )     #; print("neighs2 : "); display(neighs2)   # finds neighbors of neighbors ( i.e. second order neighbors )
    print("make_PDMatrix()");  @time  A = make_PD_Matrix( neighBs, bonds, masses, dt, ks )  
    print("CholeskyDecomp_sparse()"); @time L,neighsCh  = CholeskyDecomp_sparse( A, neighs, neighs2 ) #;print("L: "); display(L)    # evaluate Cholensky decomposition for sparse matrix (efficient)                    # evaluate Cholensky decomposition using Cholesky–Crout algorithm
    U  = copy(L') 

    print("mapMatrixNeighs()"); @time neighsL = mapMatrixNeighs(L)    #;print("neighsL"); display(neighsL)   # find which elements of L are non-zero ( for each row    i map all non-zero Aij to columns j )
    print("invNeighs()");       @time neighsU = invNeighs( neighsL ) 

    # # ===========  Generate system matrix ( Projective Dynamics Matrix ) using neighborhood information

    # np = size(points,1)
    # print("point_neighs()");   @time neighs  = point_neighs(  bonds, np )     #;print("neighs : "); display(neighs)      # finds neighbors
    # print("point_neighBs()");  @time neighBs = point_neighBs( bonds, np )     #;                                         # finds bonds to neighbors
    # print("point_neighs2()");  @time neighs2 = point_neighs2( neighs    )     


    truss  = Truss(points, bonds, masses, ks, l0s, fixed, neighBs)
    sol    = TrussSolution(A,L,U,neighsL,neighsU)

    velocity = zeros(size(points))
    #eval_forces = (truss, sol, velocity) -> eval_forces( truss, sol, velocity, ks, l0s )

    points = run_solver( truss, sol, velocity, eval_forces; dt=0.1, niter=5 ) 

    #print("points : "); display(points)
    #print("bonds  : "); display(bonds)
    #print("x : "); display(x)
    plt = plot(legend=false, aspect_ratio = :equal )
    plot_truss( plt, truss.bonds, truss.points, lw=1.0, c=:blue, bLabel=true )
    plot_truss( plt, truss.bonds, points      , lw=1.0, c=:red )
    display( plt );
end # function main


main()