using Plots
# inport get_cmap


#function Plot_Structure( bonds::Array{Int64,2}, points::Array{Float64,2} )
function plot_truss( plt, bonds, points; axes=(1,2), c=:black, strain=nothing, lw=1.0, bPoints=false, bLabel=false )
    # println("points: "); display(points)
    #plt = plot(legend=false, aspect_ratio = :equal )
    #plt = plot(legend=false, aspect_ratio = :equal )
    nb    = length(bonds)
    ix,iy = axes
    #cmap = get_cmap("viridis")
    #cmap = get_color_gradient(:viridis)
    #cmap = cgrad(:bwr)
    if strain !== nothing
        cmap = cgrad([:blue,:black,:red], [-1.0, 0.0, 1.0])
    end
    for ib=1:nb
        (i,j) = bonds[ib]
        pi = points[i,:]
        pj = points[j,:]
        if strain !== nothing
            #c = cmap( strain[ib] )
            c = get(cmap, strain[ib])
            #println("c : ", c)
        end
        plot!( plt, [pi[ix],pj[ix]], [pi[iy],pj[iy]], color=c, legend=false, lw=lw )
    end 
    # Plot each point
    if bPoints
        pxs, pys = points[:, ix], points[:, iy]
        scatter!(plt, pxs, pys, color=c, marker=:circle)
    end
    # Add labels to points
    if bLabel
        for i in 1:size(points, 1)
            pi = points[i, :]
            annotate!(plt, pi[ix], pi[iy], text(string(i), :black, :center))
        end
    end
    #display( plt )
end

#plot matrix using heatmap
function plot_matrix_bwr( plt, A )
    # select color
    cmap = cgrad(:bwr)
    vmax = max( abs(minimum(A)), abs(maximum(A)) )
    return heatmap!( plt, A, aspect_ratio = :equal , color=cmap, colorbar=true, clim=(-vmax,vmax) )
end

#plot log of absolute value of matrix using heatmap
function plot_matrix_log( plt, A, vrange=6.0, cmap=:default )
    logA = log10.( abs.(A) .+ 1.e-300 )
    vmax = maximum(logA)
    return heatmap!(plt, logA, aspect_ratio = :equal, colorbar=true, clim=(vmax-vrange,vmax) )
end