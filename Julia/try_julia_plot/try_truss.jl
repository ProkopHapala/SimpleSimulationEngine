using Random

include("ProjectedGuassSeidel.jl")
include("plot_utils.jl")
include("IterativeSolvers.jl")
include("ProjectiveDynamics.jl")

# =========== Data

"""
# ==== Sytem 1 : triangle (3 points, 3 sticks)

masses = [10.0, 1.0, 1.5]

points = [
    0.0  0.0  0.0;
    1.0 -1.0  0.0;
   -1.0 -1.0  0.0
]

bonds = [
    (1, 2),
    (2, 3),
    (3, 1)
]
"""

# ==== Sytem 2 : rope with 10 segments (11 points, 10 sticks)

#bonds, points, masses = buildRope( 11 )
nx = 3; ny = 3
bonds, points0, masses, ks = buildGrid2D( nx,ny, kDiag=1.0 )
_, l0s = process_bonds( bonds, points0 );

#print("bonds :");  display(bonds)
#print("points :"); display(points)
#print("masses :"); display(masses)

plt = plot(legend=false, aspect_ratio = :equal )

#plot_truss( plt, bonds, points0, c=:blue )

# add random displacement to the points
Random.seed!(1234)  # Sets the seed to 1234
points = points0 + 0.1*randn(size(points0))
points[1   ,:] = points0[1   ,:]
points[nx+1,:] = points0[nx+1,:]

# =========== Main Body

neighBs = point_neighBs( bonds, points );
hs,ls   = process_bonds( bonds, points );

dl     = ls - l0s
strain = dl ./ l0s; #print("strain : "); display(strain)

#Ad      = build_PGSMatrix_direct( bonds, hs, l0s, masses, neighBs );

# linear solve for the displacement x = A_ \cdot b
#x = Ad \ dl
dt = 0.1
#print("x : "); display(x)
gravity = [0.0, 0.0, -9.81] 
pnew    = pnew = points .+ gravity' .* (dt^2)

dt = 0.1
A = make_PDMatrix( neighBs, bonds, masses, dt, ks )
b = make_PD_rhs(   neighBs, bonds, masses, dt, ks, points, l0s, pnew )
print( "b : "); display(b)
print( "A : "); display(A)
x0 = zeros(size(b)) 
x  = SolveIterative( A, b, x0, niter=10, tol=1.e-3, method=1 )

print("x : "); display(x)

# solve A*x=b using direct solver
x_ref = A \ b
print("x_ref : "); display(x_ref)


# A, J, M = build_PGSMatrix(bonds, points, masses)

# Output the final matrix J and A to the console
#print("J:"); display(J)
#print("A :"); display(A )
#print("A_:"); display(A_)

plot_truss( plt, bonds, points, strain=strain*10.0, lw=3.0 )
display( plt )
