"""An N-dimensional perlin noise map.

A 2D noise map of dimensions nxm contains all grid points in the cartesian product
`(0:1:n) x (0:1:m)`.

    (0, 0)  (1, 0)  ...  (n, 0)
    (0, 1)  (1, 1)  ...  (n, 1)
    ...     ...     ...  ...
    (0, m)  (1, m)  ...  (n, m)
"""
struct PerlinNoiseMap{N, T<:Number}
    """The N-dimensional gradient vector map."""
    gradients :: Array{NTuple{N, T}, N}

    """PerlinNoiseMap inner constructor.

    @param[in] cells The number of grid cells in each dimension.
    """
    function PerlinNoiseMap{N, T}(cells :: NTuple{N, I}) where {N, T<:Number, I<:Integer}

        @assert N > 1 "The noise map must be at least 2-dimensional."
        @assert all( cells .> 0 ) "Require 1 or more gridcells in each dimension."

        # Each corner of a grid cell needs an associated gradient vector.
        # There is one more corner/gradiant in each dimension than gridcells.
        gradients = Array{NTuple{N, T}, N}(undef, cells .+ 1)
        # Generate unit gradient vectors for each gridpoint.
        for idx in eachindex(gradients)
            # Normalization is strictly required for vectors of > 1 dimensions.
            gradients[idx] = Tuple(normalize(rand(-1:0.01:1, N)))
        end

        return new{N, T}(gradients)
    end

end

"""PerlinNoiseMap outer constructor."""
function PerlinNoiseMap{T}(cells :: NTuple{N, I}) where {N, T<:Number, I}
    return PerlinNoiseMap{N, T}(cells)
end

"""Determine the number of gradients in each dimension."""
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
    return all(0 .<= cell .&& cell .< gridmax(map))
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
    cell = cell .- (point .== gridmax(map))
    @assert inbounds(map, cell) "The gridcell falls outside the noise map."
    return cell
end

"""Compute the indices of the corners of the gridcell."""
function corners(cell::NTuple{N, Integer}) where {N, T}
    # Given the cell indexes, the corresponding gradients are found at
    # index and index+1 of the gradient map of each dimension.
    return reduce(vcat, Iterators.product([
        [cellidx, cellidx+1] for cellidx in cell
    ]...))
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

    # Obtain one noise value for each gridpoint.
    # FIXME: The noise values are not in the range [-1, 1]???
    #        Determine the proper noise values, so that lerped it is in [-1, 1].
    # Note: Julia is 1-indexed, so increment EVERY gradient corner index.
    noise::Vector{T} = [
        dot(map.gradients[(corner .+ 1)...], samplepoint .- corner)
        # FIXME: Julia has 1-indexed matrices.
        for corner in corners(cell)
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
