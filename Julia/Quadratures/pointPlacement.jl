




#=
To Do : 

Use Static Arrays for better performance of 2D,3D points,vectors and matrices
https://github.com/JuliaArrays/StaticArrays.jl

Performance optimization tips:
https://docs.julialang.org/en/v1/manual/performance-tips/


* Optimization variables
  http://www.juliaopt.org/JuMP.jl/v0.19.0/variables/
   * Lazy broadcasting:
     http://www.cityinthesky.co.uk/2019/01/working-with-broadcasting-in-julia/

=#









#=

Julia Cheatsheet
https://juliadocs.github.io/Julia-Cheat-Sheet/
https://cheatsheets.quantecon.org/
https://learnxinyminutes.com/docs/julia/
https://juliabyexample.helpmanual.io/
https://julia.guide/broadcasting



b = [x for x in 1:0.5:5]'    # Row vector comprehension
a = [2, 3, 4, 5]             # Column vectro comprehension

Outer like operation
floor.(Int, sqrt.(a) .+ sqrt.(b))
@. floor(Int, sqrt(a) + sqrt(b))



ps = [ [ix*0.1 iy*0.1] for ix=1:10, iy=1:10 ]   # 2D comprenesion




in vscode : Open in Bach not Terminal Tools

// https://docs.julialang.org/en/v1/base/file/
pwd()
readdir()
include("pointPlacement.jl")

http://docs.juliaplots.org/latest/attributes/
scatter(y, marker = (:hexagon, 20, 0.6, :green, stroke(3, 0.2, :black, :dot)))
scatter(y, markershape = :hexagon,
    markersize = 20,
    markeralpha = 0.6,
    markercolor = :green,
    markerstrokewidth = 3,
    markerstrokealpha = 0.2,
    markerstrokecolor = :black,
    markerstrokestyle = :dot)

plot(y, xaxis = ("my label", (0,10), 0:0.5:10, :log, :flip, font(20, "Courier")))

plot(y, xlabel = "my label",
    xlims = (0,10),
    xticks = 0:0.5:10,
    xscale = :log,
    xflip = true,
    xtickfont = font(20, "Courier"))

=#

using LinearAlgebra

function hexGrid(n::Int,L::Float64)
    l::Float64  = L/n
    h::Float64  = l*sin(pi/3);
    #np = ( 2*n*n + numerator(n*(n-1)//2) )
    np::Int = ( (2*n+1)*(n+1) + n*n )
    #println( n," ",l," ",h," ",np )
    pts = Array{Float64}(undef,np,2)
    i::Int = 1
    x0::Float64 = -n*l*0.5
    y0::Float64 = n*h
    println( "n, x0, y0 ", n," ", x0," ", y0 )
    y = y0
    for iy::Int in 0:n
        x   = x0
        for ix::Int in 0:(n+iy)
            pts[i,1] = x
            pts[i,2] = y
            x += l
            i += 1
        end
        y  -= h
        x0 -= 0.5*l
    end
    x0 += 0.5*l
    for iy::Int in 0:(n-1)
        x0 += 0.5*l
        x   = x0
        for ix::Int in 0:(2*n-iy-1)
            pts[i,1] = x
            pts[i,2] = y
            x += l
            i += 1
        end
        y  -= h
    end
    #print( "iend = ", i )
    return pts
end

@inline function R2func(r2::Float64)
    r2 = 1-r2
    return r2*r2 
end

@inline function R2func_deriv(r2::Float64, invR2::Float64 )
    r2  *= invR2
    drf  = 1-r2
    return drf*drf, -4*drf*invR2
end

R2 = 0.04

function projectRadialToPoints( ps, centers, ys, params )
    np = size(ps,1)
    nc = size(centers,1)
    #println( "np, nc ", np," ", nc )
    for ip in 1:np
        pi = ps[ip,:]
        yi = 0.0::Float64
        for ic in 1:nc
            d  = centers[ic,:] - pi
            r2 = dot(d,d)            ::Float64
            #R2 = params[ic,1]
            #R2 = 0.125
            if r2<R2
                yi += R2func(r2/R2)
            end
        end
        ys[ip] = yi
    end
    return ys
end

function getForces( ps, centers, ys, params, fcs )
    np = size(ps,1)
    nc = size(centers,1)
    fc = zeros(Float64,2)
    dw = -1. /np
    for ic in 1:nc
        pci   = centers[ic,:] 
        #param = params[ic,:]
        fill!(fc,0.0)
        for ip in 1:np
            d  = ps[ip,:] - pci
            r2 = dot(d,d)            ::Float64
            #R2 = params[ic,1]
            #R2 = 0.125
            if r2<R2
                ytot  = ys[ip]
                y,dy  = R2func_deriv(r2, 1/R2 )
                fr    = (1-ytot)*dy
                fc   += d*fr
            end
        end
        fcs[ic,:] = fc*dw
        println( fcs ) 
    end
    return fcs
end

function move( ps, fs, vs; dt=0.1, damp=0.9 )
    #vs[:,:] *= damp ;
    #vs[:,:] += fs*dt   ;
    #ps[:,:] += vs*dt ;
    ps[:,:] += fs*dt ;
end

# ================= Main


nimg = 200
xmax = 1.0
dx   = 2*xmax/nimg

ps = [ [(-xmax+ix*dx) (-xmax+iy*dx)] for iy=1:nimg,ix=1:nimg ]
ps = vcat(ps...)
#xs = ps[:,1]
#ys = ps[:,2]

centers = hexGrid( 2, 0.5 )

#centers = [ [-0.1 0.]; [0.1 0.] ]

println( "centers : ", typeof(centers)," ", centers  )

np = size(ps,1)
nc = size(centers,1)
#;println("ncenters = ", nc)

centers .+= (rand( Float64, nc,2 ) .- 0.5)*0.1


params = rand( Float64, nc,2 )   #;println("params ", params)
ys     = Array{Float64}(undef,np)
fcs    = Array{Float64}(undef,nc,2)
vcs    = zeros(Float64,nc,2)


#println( fcs )


#using Plots

#=
x = 1:0.1:100
y = 1:0.1:100
f(x, y) = begin
        (3x + y ^ 2) * abs(sin(x) + cos(y))
    end
X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(y, 1, length(x))
Z = map(f, X, Y)
heatmap( Z, aspect_ratio=1)
#p1 = contour(x, y, f, fill=true)
#p2 = contour(x, y, Z)
#plot(p1, p2)
=#

Xs = reshape( ps[:,1], (nimg,nimg) )
Ys = reshape( ps[:,2], (nimg,nimg) )
Fs = reshape( ys,      (nimg,nimg) )
#println( Fs )
#p = contour( Xs, Ys, ys, fill=true)
#plot(p)
#@doc heatmap
#heatmap( Fs, aspect_ratio=1)

projectRadialToPoints( ps, centers, ys, params )
getForces( ps, centers, ys, params, fcs )

using PyPlot
imshow( Fs, origin="image", extent=(-1.,1.,-1.,1.) )
colorbar()
#x = range(0; stop=2*pi, length=1000); 
#y = sin.(3 * x + 4 * cos.(2 * x));
#plot(x, y, color="red", linewidth=2.0, linestyle="--")
#plot(x, y, ".-" )

#plot(ps[:,1], ps[:,2], ".",  markersize=1.0 )
plot(centers[:,1], centers[:,2], "." )
quiver( centers[:,1], centers[:,2], fcs[:,1], fcs[:,2], )


for iter in 0:100
    close()
    projectRadialToPoints( ps, centers, ys, params )
    getForces( ps, centers, ys, params, fcs )
    #println( fcs )
    #println( centers )
    println( maximum( fcs ) )
    
    imshow( Fs, origin="image", extent=(-1.,1.,-1.,1.) )
    colorbar()
    plot(centers[:,1], centers[:,2], "." )
    quiver( centers[:,1], centers[:,2], fcs[:,1], fcs[:,2], )
    #axis('equal')
    sleep(0.2)
    move( centers, fcs, vcs, dt=0.5, damp=0.9 )
end


#plot( centers[:,1], centers[:,2] )
#show()


#scatter!( ps[:,1], ps[:,2], markersize = 1, markerstrokewidth = 0, aspect_ratio=:equal )
#scatter!( centers[:,1], centers[:,2], markersize = 2, markerstrokewidth = 0, aspect_ratio=:equal )



#plot( pts[:,1], pts[:,2] )
#scatter( pts[:,1], pts[:,2] )
#scatter( pts[:,1], pts[:,2], markerstrokestyle = :dot )

#pts = hexGrid( 2, 1.0 )
#scatter( pts[:,1], pts[:,2], markersize = 2, markerstrokewidth = 0, aspect_ratio=:equal )
#plot( pts[:,2] )

#show()
#close()