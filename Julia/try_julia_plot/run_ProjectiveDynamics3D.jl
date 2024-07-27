using Random
using Revise # prevent problem with redefinition of types at runtipme (when using include() in REPL)
using DelimitedFiles


include("Truss.jl")
#include("ProjectedGuassSeidel.jl")
include("plot_utils.jl")
#include("IterativeSolvers.jl")
include("Cholesky.jl")
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

function eval_force_and_plot( position, velocity,   plt, bonds )
    plot_truss( plt, bonds, position, lw=0.5, c=:black )
    force = -0.0 * position - 0.00 * velocity
    return force
end


function eval_force_and_plot_f( position::Matrix{Float32}, velocity::Matrix{Float32},   plt, bonds )
    plot_truss( plt, bonds, position, lw=0.5, c=:black )
    force = -0.0f0 * position - 0.0f0 * velocity
    return force
end

###################
###    main()   ###
###################

function main(;nseg=8, wheel_width=0.2, k=10000.0, dt=5.0, niter=5, omega=1.0 )

    # ===========  Generate Truss

    p0 =  Vector([0.0,0.0,0.0])
    p1 =  Vector([1.0,0.0,0.0])
    ax =  Vector([0.0,0.0,1.0])
    #print("wheel()");       @time points, bonds, ks = wheel( p0, p1, ax, wheel_width, n=nseg, kscale=k )
    print("ngonTruss()");    @time points, bonds, ks = ngonTruss( p0, p1, ax, n=3, k=k )
    plt = plot(legend=false, aspect_ratio = :equal )
    #println("typeof(points) ", typeof(points) );
    print("process_bonds()"); @time   _, l0s          = process_bonds( bonds, points );
    #print("bonds :");  display(bonds)
    #print("points :"); display(points)
    #print("masses :"); display(masses)
    np = size(points,1)
    masses = ones(np)
    fixed  = []
    #masses[fixed] .= 100000.0

    print("ks:       ", ks);  #display(ks     );
    #print("masses:  ", masses);  #display(masses );

    # ===========  Generate Neigbor Lists

    print("point_neighs()");   @time neighs  = point_neighs(  bonds, np )     #;print("neighs : "); display(neighs)      # finds neighbors
    print("point_neighBs()");  @time neighBs = point_neighBs( bonds, np )     #;                                         # finds bonds to neighbors
    print("point_neighs2()");  @time neighs2 = point_neighs2( neighs    )     #; print("neighs2 : "); display(neighs2)   # finds neighbors of neighbors ( i.e. second order neighbors )

    #print("neighBs: ");  display(neighBs);

    # ===========  Generate system matrix ( Projective Dynamics Matrix ) using neighborhood information

    print("make_PDMatrix()");  @time  A,Mt = make_PD_Matrix( neighBs, bonds, masses, dt, ks )  
    
    #print("CholeskyDecomp_sparse()"); @time L,neighsCh  = CholeskyDecomp_sparse( A, neighs, neighs2 ) #;print("L: "); display(L)    # evaluate Cholensky decomposition for sparse matrix (efficient)                    # evaluate Cholensky decomposition using Cholesky–Crout algorithm
    #U  = copy(L') 

    LDLT_L, LDLT_D              = CholeskyDecomp_LDLT(A)
    LDLT_L, LDLT_D, LDLT_neighs = CholeskyDecomp_LDLT_sparse(A,neighs)

    println("A  ", size(A    ) ); display(A    );
    println("Mt ", Mt );
    # println("LDLT_L ", size(LDLT_L ) ); display(LDLT_L);

    #writedlm("PDmat.csv",  A,      ',')
    #writedlm("LDLT_L.csv", LDLT_L, ',')
    #writedlm("PDmat.csv",  A,      ' ')
    #writedlm("LDLT_L.csv", LDLT_L, ' ')

    # ===========  Check C++ result vs Julia
    # A_cpp = readdlm("../../tests_bash/Orbital/PDmat.log",  Float64)
    # L_cpp = readdlm("../../tests_bash/Orbital/LDLT_L.log", Float64)
    # #println("A_cpp ", size(A_cpp) ); display(A_cpp);
    # dAcpp = A_cpp - A;   errAcpp = norm(dAcpp)
    # println("|A_cpp - A|=", errAcpp," size=", size(dAcpp),  ); #display(dAcpp);
    # #plot( dAcpp );
    # println("L_cpp  ", size(L_cpp  ) ); display(L_cpp );
    # println("LDLT_L ", size(LDLT_L ) ); display(LDLT_L);
    # dLcpp = L_cpp - LDLT_L; errLcpp = norm(dLcpp)
    # println("|L_cpp - LDLT_L|=", errLcpp," size=", size(dLcpp),  ); #display(dAcpp);

    reconstructed = LDLT_L * Diagonal(LDLT_D) * LDLT_L'
    diff = A - reconstructed
    println("Max.Error|A-LDL|=': ", maximum(abs.(diff)))

    #print("mapMatrixNeighs()"); @time neighsL = mapMatrixNeighs(L)    #;print("neighsL"); display(neighsL)   # find which elements of L are non-zero ( for each row    i map all non-zero Aij to columns j )
    #print("invNeighs()");       @time neighsU = invNeighs( neighsL ) 

    # np = size(points,1)
    # print("point_neighs()");   @time neighs  = point_neighs(  bonds, np )     #;print("neighs : "); display(neighs)      # finds neighbors
    # print("point_neighBs()");  @time neighBs = point_neighBs( bonds, np )     #;                                         # finds bonds to neighbors
    # print("point_neighs2()");  @time neighs2 = point_neighs2( neighs    )     

    truss   = Truss(points, bonds, masses, ks, l0s, fixed, neighBs)
    #sol    = TrussSolution(A,L,U,neighsL,neighsU, LDLT_L,LDLT_D,LDLT_neighs )
    sol     = LDLTsolution(LDLT_L,LDLT_D,LDLT_neighs)
    truss_f = convert_to_truss_f( truss )
    #sol_f  = convert_to_TrussSolution_f( sol  )
    sol_f   = convert_to_LDLTsolution_f( sol  )

    axis = [0.0, 0.0, 1.0]  # Example axis of rotation
    #velocity = Matrix( eachrow(points) .|> p -> cross(axis, p) )
    velocity = mapslices(p -> cross(axis, p)*omega, points, dims=2)
    #velocity = zeros(size(points))

    ps_pred  = points .+ velocity*dt
    bK, bMt  = make_PD_rhs_decomp( truss.neighBs, truss.bonds, truss.masses, dt, truss.ks, points, truss.l0s, ps_pred )  
    println("bK  ", bK  );
    println("bMt ", bMt );

    plot_truss( plt, truss.bonds, truss.points, lw=2.0, c=:blue, bLabel=true, legend="Initial" )
    plot_truss( plt, truss.bonds, ps_pred, lw=2.0, c=:red, legend="Predicted" )

    # --- Float64 solver 
    eval_forces = (position, velocity) -> eval_force_and_plot(position,velocity, plt, truss.bonds )
    points      = run_solver( truss, sol, velocity, eval_forces; dt=dt, niter=niter, A_check=A ) 
    plot_truss( plt, truss.bonds, points     , lw=2.0, c=:green, legend="Constrained"  )

    # --- Float32 solver
    #eval_forces_f = (position, velocity) -> eval_force_and_plot_f(position,velocity, plt, truss_f.bonds )
    #points_f      = run_solver_f( truss_f, sol_f, convert.(Float32, velocity), eval_forces_f; dt=Float32(dt), niter=niter ) 
    #plot_truss( plt, truss.bonds, points_f    , lw=2.0, c=:red )

    #print("points : "); display(points)
    #print("bonds  : "); display(bonds)
    #print("x : "); display(x)
    
    display( plt );
end # function main


#main()

main( dt=0.5, niter=5, k=1000000.0 )
#main( dt=0.1, niter=1 )

#main(dt=0.01 )
#main(dt=0.03)
#main(dt=0.05)
#main(dt=0.08)
#main(dt=0.05, nseg=6)
#main(dt=0.05, nseg=3)
#main(dt=0.1)
#main(dt=1.0)
#main(dt=5.0)