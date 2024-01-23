# Docs:
# https://docs.juliaplots.org/latest/generated/unitfulext_examples/#Heatmaps
# https://docs.juliaplots.org/latest/generated/unitfulext_plots/
# array slicing in Julia:  https://www.matecdev.com/posts/julia-array-indexing.html
using Plots
#using Plots.Heatmaps  # Add this line to import the Heatmaps module

function plot_bonds( ps, edges )
    # Initialize an empty plot
    p = plot()
    # Add each edge as a line segment to the plot
    for (i, j) in edges
        pi = ps[i]
        pj = ps[j]
        plot!(p, [pi[1], pj[1]], [pi[2], pj[2]], seriestype=:path, line=:black, legend=false)
    end
end

# Example data
ps = [ [1,2], [2,4], [3,1], [4,3]]
# Define edges (pairs of vertex indices)
#edges = [(1, 2), (2, 3), (3, 4), (4, 1), (4, 1), (4, 1) ]
edges = [(1, 2), (1, 3), (2, 4), (3, 4), (1,4) ]

plot_bonds( ps, edges )

# Display the plot
display(p)