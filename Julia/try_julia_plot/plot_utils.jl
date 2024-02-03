using Plots
# inport get_cmap


#function Plot_Structure( bonds::Array{Int64,2}, points::Array{Float64,2} )
function plot_truss( plt, bonds, points; axes=(1,2), c=:black, strain=nothing, lw=1.0 )
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
    pxs, pys = points[:, ix], points[:, iy]
    scatter!(plt, pxs, pys, color=c, marker=:circle)
    #display( plt )
end