using LinearAlgebra
using Plots



"""An N-dimensional perlin noise map.

A 2D noise map of dimensions nxm contains all grid points in the cartesian product
`(0:1:n) x (0:1:m)`.

    (0, 0) (1, 0) ... (n, 0)
    (0, 1) (1, 1) ... (n, 1)
    ...    ...    ... ...
    (0, m) (1, m) ... (n, m)
"""
struct PerlinNoiseMap{N, T<:Number}
    """The N-dimension gradient vector map."""
    gradients :: Array{NTuple{N, T}, N}

    """PerlinNoiseMap inner constructor.

    @param[in] dims The number of gridpoints of each dimension.
    """
    function PerlinNoiseMap{N, T}(
            dims :: NTuple{N, I}) where {N, T<:Number, I<:Integer}

        @assert N > 1 "The noise map must be at least 2-dimensional."
        @assert all( dims .> 1 ) "Require 2 or more gridpoints per dimension."

        # Each corner of a grid cell needs an associated gradient vector.
        gradients = Array{NTuple{N, T}, N}(undef, dims)
        # Generate unit gradient vectors for each gridpoint.
        for idx in eachindex(gradients)
            # Normalization is strictly required for vectors of > 1 dimensions.
            gradients[idx] = Tuple(normalize(rand(-1:0.01:1, N)))
        end

        return new{N, T}(gradients)
    end

end

"""Determine the number of gridcells in each dimension."""
function dims(map::PerlinNoiseMap{N, T}) where {N, T}
    return size(map.gradients)
end

"""Compute the maximum coordinate inbound the grid for each dimension."""
function gridmax(map::PerlinNoiseMap{N, T}) where {N, T}
    # Assume each gridcell is a cube with sides of length 1.
    # Gradients correspond to gridcell corners, so there is one less
    # gradient in each dimension than there are cells.
    return size(map.gradients) .- 1
end

"""Check that the gridcell is in bounds of the noise map."""
function inbounds(map::PerlinNoiseMap{N, T}, cell::NTuple{N, Integer}) where {N, T}
    # If a dimension spans x gridpoints (gradients), then it spans x-1 gridcells.
    # Then, the valid gridcell indices are 0, 1, ..., x-2 for that dimension.
    return all(0 .<= cell .&& cell .<= gridmax(map))
end

"""Find the gridcell that the point is contained in.

@post The corresponding gridcell must be inbounds of the noise map.
"""
function gridcell(map::PerlinNoiseMap{N, T}, point::NTuple{N, T}) where {N, T}
    # Assume each gridcell is a cube with sides of length 1.
    cell::NTuple{N, Integer} = Tuple( @. Integer(floor(point)) )
    # Given a gridcell with left bound lo and right bound hi.
    # Then consider a point x in bounds of the gridcell if lo <= x < hi.
    # But, for the rightmost gridcell in the grid, instead lo <= x <= hi holds.
    # The use of `floor` means points on the rightmost boundary would naturally
    # fall outside the grid points. Manually take this edge case into account.
    cell = cell .- (cell .== gridmax(map))
    @assert inbounds(map, cell) "The gridcell falls outside the noise map."
    return cell
end

"""The smoothstep sigmoid-like smoothing function.

@param[in] x The value to smooth.
@return A smoothed value in [0, 1].
"""
function smoothstep(x::T) where {T<:Number}
    @assert 0.0 <= x <= 1.0 "$x in [0, 1] does not hold."
    # Perform sigmoid-like smoothing: 3(x^2) - 2(x^3) = (3 - 2x) * x^2
    return (3 - 2*x) * x^2
end

"""Check if x = 2^n for any n = 0, 1, ... holds."""
function ispow2(x::Integer)
    return (x & (x - 1)) == 0
end

"""Linearly interpolate between a and b."""
function lerp(frac::T, a::T, b::T) where {T}
    @assert 0.0 <= frac <= 1.0 "The interpolation fraction must be in [0, 1]."
    return a + frac*(b - a)
end

"""Obtain the noise value for the given point in space."""
function noise(map::PerlinNoiseMap{N, T}, samplepoint::NTuple{N, T}) where {N,T}
    # Determine which grid cell the space point falls in.
    cell::NTuple{N, Integer} = gridcell(map, samplepoint)

    # Given the cell indexes, the corresponding gradients are found at
    # index and index+1 of the gradient map of each dimension.
    indexes = Iterators.product([
        [cellidx, cellidx+1] for cellidx in cell
    ]...)

    # FIXME: Julia is 1-indexed, so increment EVERY index.
    Jindexes = [ cell .+ 1 for cell in reduce(vcat, indexes) ]

    # Obtain one noise value for each gridpoint.
    # FIXME: The noise values are not in the range [-1, 1]???
    #        Determine the proper noise values, so that lerped it is in [-1, 1].
    noise::Vector{T} = [
        dot(map.gradients[Jcorner...], samplepoint .- corner)
        # FIXME: Julia has 1-indexed matrices.
        for (corner, Jcorner) in zip(indexes, Jindexes)
    ]

    # Linearly interpolate the noise value sequence.
    # The i-th iteration lerps the i-th dimension.
    # Each iteration halves the length of the sequence.
    for dim in 1:N
        # The lerp sequence must be of length 2^n.
        @assert ispow2(length(noise))

        # The projection of the sample point onto the x-axis is fraction t
        # along the line segment going from corner c_1 to corner c_2.
        t = samplepoint[dim] - cell[dim]
        # Make the linear interpolation sigmoid-like.
        t = smoothstep(t)

        noise = [
            lerp(t, noise[i], noise[i+1])
            for i in 1:2:length(noise)
        ]
    end

    # Finally, extract the lerped noise value.
    @assert length(noise) == 1
    return noise[1]
end

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
function image2d(map::PerlinNoiseMap{N, T}, width::Integer, height::Integer) where {N,T}
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
    # scatter!(p, samples, markersize=1) # FIXME: uncomment?
    quiver!(p, map)

    return p
end

map = PerlinNoiseMap{2,Float64}((5,5))
println(typeof(map.gradients)); println(map.gradients); println()



# # TODO: If a 2D noise map represents a terrain topology, then
# #          a 3D noise map represents the evolution of a topoly
# #       throughout time. I.e. if you take the third
# #       noise map dimension to represent the time axis, and assign
# #       it a suitable small time step, then you can animate the
# #       evolution of topoly through time (make a gif of this).


p = image2d(map, 400, 400)

display(p)
println("Press ENTER to continue.")
readline()
