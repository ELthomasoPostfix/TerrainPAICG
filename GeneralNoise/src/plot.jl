function quiver(map::PerlinNoiseMap{N, T}) where {N,T}
    plt = plot()
    quiver!(plt, map)
    return plt
end

"""Generate a quiver (vector field) plot of the noise map."""
function quiver!(plt, map::PerlinNoiseMap{N, T}) where {N,T}
    dimensions = dim(map)
    @assert 1 <= dimensions <= 2 "Can only plot 1D or 2D noise maps."

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

    # For a 1D map, fix the y-coordinates of all gradient tails to zero.
    if dimensions == 1
        push!(tails, zero(tails[1]))
        pushfirst!(deltas, zero(deltas[1]))
    end

    Plots.quiver!(plt, tails..., quiver=Tuple(deltas), arrow=(:closed,:head,2.0,2.0))
end

"""Generate a noise curve plot of the noise map."""
function image1d(map::PerlinNoiseMap{N, T}, width::Integer;
        pscatter::Bool=false, pquiver::Bool=false) where {N,T}
    dimensions = dim(map)
    @assert dimensions == 1 "Can only plot 1D noise maps."
    @assert width >= 1 "The width must be at least 1."

    xmax, = gridmax(map)
    samples::Vector{Tuple{Float64, Float64}} = []
    for x in 1:width
        sample = (x - 1) / width * xmax
        pnoise = noise(map, (sample,))
        push!(samples, (sample, pnoise))
    end
    p = plot(samples, color = :greys)
    if pscatter begin scatter!(p, samples, markersize=1) end end
    if pquiver begin quiver!(p, map) end end
    plot!(p, yticks=[-1, 0, 1], ylims=(-1,1))

    return p
end

"""Generate a heatmap plot of the noise map."""
function image2d(map::PerlinNoiseMap{N, T}, width::Integer, height::Integer;
        pscatter::Bool=false, pquiver::Bool=false) where {N,T}
    dimensions = dim(map)
    @assert dimensions == 2 "Can only plot 2D noise maps."
    @assert width >= 1 && height >= 1 "The dimensions must be at least 1x1."

    xmax, ymax = gridmax(map)
    image = Array{Float64, 2}(undef, width, height)
    samples::Vector{Tuple{Float64, Float64}} = []
    for x in 1:width
        for y in 1:height
            sample = (
                (x - 1) / width * xmax,
                (y - 1) / height * ymax
            )
            push!(samples, sample)
            image[x,y] = noise(map, sample)
        end
    end
    xticks = (0:width-1) / width * xmax
    yticks = (0:height-1) / height * ymax
    p = heatmap(xticks, yticks, image, color = :greys)
    if pscatter begin scatter!(p, samples, markersize=1) end end
    if pquiver begin quiver!(p, map) end end

    return p
end
