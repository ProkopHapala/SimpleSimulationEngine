
using Plots


function getLJ( r, R, e )
    u  = R/r
    u6 = u^6
    E  = e *    ( u6*u6 - 2*u6 )
    F  = e * 12*( u6*u6 - 2*u6 )/r 
    return E,F
end



function getLJx2( r, R, e )
    if(r<R)
        E  = (r-R)^2 - e
        F  = 0
    else
        u  = R/r
        u6 = u^6
        E  = e *    ( u6*u6 - 2*u6 )
        F  = e * 12*( u6*u6 - 2*u6 )/r    
    end
    return E,F
end

function getLJx2( r, R, e )
    if(r<R)
        E  = (r-R)^2 - e
        F  = 0
    else
        u  = R/r
        u6 = u^6
        E  = e *    ( u6*u6 - 2*u6 )
        F  = e * 12*( u6*u6 - 2*u6 )/r    
    end
    return E,F
end


# eval_forces = (position, velocity) -> eval_force_and_plot(position,velocity, plt, truss.bonds )

function xrange( x0, dx, n )
    xs = Vector{Float64}(undef,n)
    for i in 1:n 
        xs[i] = x0 + dx*(i-1)    
    end
    return xs
end

function plot_func( plt, xs, func::Function; clr=nothing )
    #p = plot()
    # Add each edge as a line segment to the plot
    n = size(xs,1)
    E = Vector{Float64}(undef,n)
    F = Vector{Float64}(undef,n)
    for i in 1:n 
        E[i],F[i] = func( xs[i] )    
    end
    plot!(plt, xs, E, seriestype=:path, line=clr, legend=false)
end

xs = xrange( 2.5, 0.1, 300 )


plt = plot()
plot_func( plt, xs, (x)->getLJ(x,3.0,1.0),   clr=:black  )
plot_func( plt, xs, (x)->getLJx2(x,3.0,1.0), clr=:blue   )

# Display the plot
display(plt)

