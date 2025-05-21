function quiver(map::PerlinNoiseMap{N, T}) where {N,T}
    plt = plot()
    quiver!(plt, map)
    return plt
end

"""Generate a quiver (vector field) plot of the noise map."""
function quiver!(plt, map::PerlinNoiseMap{N, T}) where {N,T}
    dimensions = length(dims(map))
    @assert 2 <= dimensions <= 3 "Can only plot 2D or 3D noise maps."

    tails::Vector{Vector{T}} = [ [] for _ in 1:dimensions ]
    deltas::Vector{Vector{T}} = [ [] for _ in 1:dimensions ]
    gridsizes = dims(map)

    # Extract the matrix index for each gridpoint.
    gridpoints = Iterators.product([
        1:gridsizes[dim] for dim in 1:dimensions
    ]...)

    # Generate an arrow for each (gridpoint, gradient) pair.
    for gridpoint in gridpoints
        gradient = map.gradients[gridpoint...]
        for dim in 1:dimensions
            # Julia is 1-indexed, but cartesian coordinates start from 0.
            push!(tails[dim], gridpoint[dim] - 1)
            # Scale the gradients to avoid cluttering the plot.
            gradscalar::Float64 = 1 / 5
            push!(deltas[dim], gradient[dim] * gradscalar)

        end
    end

    Plots.quiver!(plt, tails..., quiver=Tuple(deltas), arrow=(:closed,:head,2.0,2.0))
end

"""Generate a heatmap plot of the noise map."""
function image2d(map::PerlinNoiseMap{N, T}, width::Integer, height::Integer;
        pscatter::Bool=false, pquiver::Bool=false) where {N,T}
    dimensions = length(dims(map))
    @assert dimensions == 2 "Can only plot 2D noise maps."

    xmax, ymax = gridmax(map)
    image = Array{Float64, 2}(undef, width, height)
    samples::Vector{Tuple{Float64, Float64}} = []
    for x in 1:width
        for y in 1:height
            sample = (
                (x - 1) / (width - 1) * xmax,
                (y - 1) / (height - 1) * ymax
            )
            push!(samples, sample)
            image[x,y] = noise(map, sample)
        end
    end
    xticks = (0:width-1)/(width - 1) * gridmax(map)[1]
    yticks = (0:height-1)/(height - 1) * gridmax(map)[2]
    p = heatmap(xticks, yticks, image, color = :greys)
    if pscatter begin scatter!(p, samples, markersize=1) end end
    if pquiver begin quiver!(p, map) end end

    return p
end
