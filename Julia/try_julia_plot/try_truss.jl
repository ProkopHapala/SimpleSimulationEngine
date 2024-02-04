using Random

include("ProjectedGuassSeidel.jl")
include("plot_utils.jl")
include("IterativeSolvers.jl")
include("ProjectiveDynamics.jl")

# =========== Data

function main()

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


# ==== Sytem 2 : rope with 10 segments (11 points, 10 sticks)

#bonds, points, masses = buildRope( 11 )
#nx = 3; ny = 3
nx = 3; ny = 3

#nx = 10; ny = 10

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
np = size(points,1)
neighs  = point_neighs( bonds,  np )     #;print("neighs : "); display(neighs)
neighBs = point_neighBs( bonds, np )      #;
hs,ls   = process_bonds( bonds, points )  #;

neighs2 = point_neighs2( neighs )    #; print("neighs2 : "); display(neighs2)

#return

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
A = make_PDMatrix( neighBs, bonds, masses, dt, ks )                     # ;print( "A : "); display(A)
b = make_PD_rhs(   neighBs, bonds, masses, dt, ks, points, l0s, pnew )  # ;print( "b : "); display(b)

#Ch = cholesky(A);             #;print("Cholesky.L : "); display(Ch.L)
#L  = CholeskyDecomp( A )       ;print("L : "); display(L)
#L_ = CholeskyDecomp_Crout( A ) ;print("L_ : "); display(L_)
#Ch = cholesky(A);              ;print("Cholesky.L : "); display(Ch.L)
@time begin
L_C = CholeskyDecomp_Crout( A )          #;print("L_C: "); display(L_C)
end
@time begin
L_   = CholeskyDecomp_sparse( A, neighs ) #;print("L_: "); display(L_)
end
@time begin
L  = CholeskyDecomp_sparse( A, neighs, neighs2 ) #;print("L: "); display(L)
end


plot_matrix_log(plt, L-L_C )
#plot_matrix_log(plt, L-L_ )

#display( plt ); return

U  = copy(L') 
println("CholeskyDecomp_sparse max_error=max(abs(L*L^T-A)) : ", maximum(abs.(L*U      - A)) );
println("CholeskyDecomp_Crout  max_error=max(abs(L*L^T-A)) : ", maximum(abs.(L_C*L_C' - A)) );

#plot_matrix_log(plt, Ch.L )
#plot_matrix_log(plt, L )
#plot_matrix_bwr( plt, err_Ch )
#plot_matrix_log(plt, Ch.U )

#y_ch = forwardsub_bak(L,b[:,1])  #;print("y_chb : "); display(y_ch)
#x_ch = backsub_bak(  U,y_ch)     ;print("x_chb : "); display(x_ch)


y_ch = forwardsub(L,b)  #;print("y_ch : "); display(y_ch)
x_ch = backsub(U,y_ch)  #;print("x_ch : "); display(x_ch)
x_ref = A \ b           #;print("x_ref : ");display(x_ref)

x_err = x_ch - x_ref    #;print("x_err : ");display(x_err)




#==


x0 = zeros(size(b)) 
x  = SolveIterative( A, b, x0, niter=10, tol=1.e-3, method=1 )   #; print("x : "); display(x)

# solve A*x=b using direct solver
x_ref = A \ b     #;print("x_ref : "); display(x_ref)


# invert matrix A_
invA = inv(A)      #;print("invA : "); display(invA)
x_inv = invA * b   #;print("x_inv : "); display(x_inv)

# if element of invA<trahold, set it to zero
threshold = 1.e-6
invA[ abs.(invA) .< threshold ] .= 0.0

# cholensky factorization of A

#plot_matrix_bwr(plt, A )
#plot_matrix_log(plt, A )
#plot_matrix_log(plt, invA )
plot_matrix_log(plt, L )

# A, J, M = build_PGSMatrix(bonds, points, masses)

# Output the final matrix J and A to the console
#print("J:"); display(J)
#print("A :"); display(A )
#print("A_:"); display(A_)

plot_truss( plt, bonds, points, strain=strain*10.0, lw=1.0 )
display( plt )

==#
display( plt )
end # function main


main()