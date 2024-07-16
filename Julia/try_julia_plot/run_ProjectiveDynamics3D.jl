using Random

include("ProjectedGuassSeidel.jl")
include("plot_utils.jl")
include("IterativeSolvers.jl")
include("ProjectiveDynamics.jl")

# =========== Data


###################
###    main()   ###
###################

function main( nx=5, ny=5, dt=5.0 )
    p0 =  Vector([0.0,0.0,0.0])
    p1 =  Vector([1.0,0.0,0.0])
    ax =  Vector([0.0,0.0,1.0])
    points, bonds, ks = wheel( p0, p1, ax, 8, 0.2 )
    #print("points : "); display(points)
    #print("bonds  : "); display(bonds)
    #print("x : "); display(x)
    plt = plot(legend=false, aspect_ratio = :equal )
    plot_truss( plt, bonds, points, lw=1.0, c=:blue )
    display( plt );
end # function main


main()